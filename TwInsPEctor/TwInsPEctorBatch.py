import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import shlex
import shutil
import subprocess
import sys
import pandas as pd


BOOL_ARGS = {
    "recoding_mode",
    "plot_full_reads",
    "no_collapse_displayed_alleles",
    "no_alignment_adjustments",
    "ignore_extraspacer_deletions",
    "no_figures",
    "no_summary_figures",
    "no_per_base_figures",
    "no_mutation_figures",
    "plot_extended_tables",
    "save_pdf",
    "no_save_png",
    "no_rerun",
    "keep_crispresso_outputs",
    "verbose",
}

INT_ARGS = {
    "min_num_base_edits",
    "default_min_aln_score",
    "max_n_rows",
    "max_n_alleles_to_write",
    "cleavage_offset_a",
    "cleavage_offset_b",
    "n_processes",
}

FLOAT_ARGS = {
    "min_frequency_alleles",
}

EQUALS_STYLE_ARGS = {
    "crispresso_args",
    "cleavage_offset_a",
    "cleavage_offset_b",
}

BATCH_ONLY_ARGS = {
    "batch_file",
    "sample_name",
    "dry_run",
}

REQUIRED_ARGS = ["fastq_r1", "wt_seq", "tpe_seq", "peg_spacers"]

COLUMN_ALIASES = {
    "twin_seq": "tpe_seq",
}


def _is_missing(value):
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):
        pass
    if isinstance(value, str) and value.strip() == "":
        return True
    return False


def _parse_bool(value, col_name):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in ("true", "1", "yes"):
            return True
        if lowered in ("false", "0", "no"):
            return False
    raise ValueError(f"Invalid boolean value '{value}' for column '{col_name}'.")


def _fastq_stem(path):
    pattern = r'([^/]+?)(?=(?:\.fastq|\.fq)?(?:\.gzip|\.gz|\.bz2|\.bz|\.xz|\.lzma)?$)'
    match = re.search(pattern, str(path))
    return match.group(1) if match else os.path.basename(str(path))


def _resolve_twinspector_executable():
    exe = shutil.which("TwInsPEctor")
    if exe:
        return [exe]
    return [sys.executable, "-m", "TwInsPEctor.TwInsPEctor"]


def build_sample_commands(args, batch_df, per_sample_processes):
    batch_df = batch_df.copy()

    # Map legacy column aliases (e.g. twin_seq -> tpe_seq) if the canonical column is absent
    for alias, canonical in COLUMN_ALIASES.items():
        if alias in batch_df.columns:
            if canonical not in batch_df.columns:
                print(f"Note: Renaming legacy batch column '{alias}' to '{canonical}'.")
                batch_df = batch_df.rename(columns={alias: canonical})
            else:
                print(f"Warning: Both '{alias}' and '{canonical}' found in batch file; using '{canonical}'.")

    valid_columns = set(vars(args).keys()) | {"sample_name"}
    for col in batch_df.columns:
        if col not in valid_columns:
            print(f"Warning: Column '{col}' in batch file is not a recognized argument and will be ignored.")

    twinspector_cmd = _resolve_twinspector_executable()
    commands = []

    for row_idx, row in batch_df.iterrows():
        sample_label = f"row {row_idx + 1}"
        sample_name = None
        if "sample_name" in batch_df.columns and not _is_missing(row["sample_name"]):
            sample_name = str(row["sample_name"]).strip()
            sample_label = sample_name

        resolved = {}
        for arg_name, cli_val in vars(args).items():
            if arg_name in BATCH_ONLY_ARGS:
                continue

            row_has_val = arg_name in batch_df.columns and not _is_missing(row[arg_name])
            val = row[arg_name] if row_has_val else cli_val

            if arg_name in BOOL_ARGS:
                val = _parse_bool(val, arg_name) if not _is_missing(val) else False
            elif not _is_missing(val):
                if arg_name in INT_ARGS:
                    val = int(val)
                elif arg_name in FLOAT_ARGS:
                    val = float(val)
                else:
                    val = str(val).strip()
            else:
                val = None

            resolved[arg_name] = val

        # Validate required arguments for this sample
        missing_required = [req for req in REQUIRED_ARGS if _is_missing(resolved.get(req))]
        if missing_required:
            raise ValueError(
                f"Sample '{sample_label}' is missing required argument(s): {', '.join(missing_required)}."
            )

        # Resolve per-sample output_root so batch samples do not overwrite each other
        row_has_output_root = "output_root" in batch_df.columns and not _is_missing(row["output_root"])
        if not row_has_output_root:
            if sample_name:
                if args.output_root:
                    resolved["output_root"] = os.path.join(args.output_root, sample_name)
                else:
                    resolved["output_root"] = sample_name
            elif args.output_root:
                r1_stem = _fastq_stem(resolved["fastq_r1"])
                if resolved.get("fastq_r2"):
                    r2_stem = _fastq_stem(resolved["fastq_r2"])
                    subdir = f"TwInsPEctor_on_{r1_stem}_{r2_stem}"
                else:
                    subdir = f"TwInsPEctor_on_{r1_stem}"
                resolved["output_root"] = os.path.join(args.output_root, subdir)

        # Apply computed per-sample process budget unless explicitly overridden per row
        row_has_n_processes = "n_processes" in batch_df.columns and not _is_missing(row["n_processes"])
        if not row_has_n_processes:
            resolved["n_processes"] = per_sample_processes

        cmd = list(twinspector_cmd)
        for arg_name, val in resolved.items():
            if arg_name in BOOL_ARGS:
                if val:
                    cmd.append(f"--{arg_name}")
            elif val is not None:
                if arg_name in EQUALS_STYLE_ARGS:
                    cmd.append(f"--{arg_name}={val}")
                else:
                    cmd.extend([f"--{arg_name}", str(val)])

        commands.append((sample_label, cmd))

    return commands


def _format_cmd(cmd):
    return " ".join(shlex.quote(tok) for tok in cmd)


def _run_sample(sample_label, cmd):
    print(f"[{sample_label}] Running: {_format_cmd(cmd)}", flush=True)
    result = subprocess.run(cmd, check=False)
    return sample_label, result.returncode


def main(args=None):
    if args is None:
        args = parse_args()

    input_batch_file = args.batch_file
    if not os.path.isfile(input_batch_file):
        print(f"Error: Batch file '{input_batch_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    batch_df = pd.read_csv(input_batch_file, sep="\t")
    if batch_df.empty:
        print(f"Error: Batch file '{input_batch_file}' contains no samples.", file=sys.stderr)
        sys.exit(1)

    batch_workers = min(len(batch_df), args.n_processes)
    per_sample_processes = max(1, args.n_processes // batch_workers)

    try:
        commands = build_sample_commands(args, batch_df, per_sample_processes)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        for sample_label, cmd in commands:
            print(f"[{sample_label}] {_format_cmd(cmd)}")
        return

    failures = []
    if batch_workers == 1:
        for sample_label, cmd in commands:
            label, returncode = _run_sample(sample_label, cmd)
            if returncode != 0:
                failures.append((label, returncode))
    else:
        with ThreadPoolExecutor(max_workers=batch_workers) as executor:
            futures = [
                executor.submit(_run_sample, sample_label, cmd)
                for sample_label, cmd in commands
            ]
            for future in as_completed(futures):
                label, returncode = future.result()
                if returncode != 0:
                    failures.append((label, returncode))

    if failures:
        failed_summary = ", ".join(f"{label} (exit {code})" for label, code in failures)
        print(
            f"Finished with {len(failures)}/{len(commands)} failed sample(s): {failed_summary}",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Finished running TwInsPEctor on {len(commands)} sample(s).")


def parse_args(argv=None):
    prog = "TwInsPEctorBatch"
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Analyzes Twin Prime Editing outcomes in batches by running TwInsPEctor on each sample.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-f", "--batch_file",
        required=True,
        help="Path to a tab-delimited batch file (columns can include sample_name, fastq_r1, fastq_r2, wt_seq, tpe_seq, peg_spacers, rt_templates, and any other TwInsPEctor CLI argument names).",
    )
    parser.add_argument("-r1", "--fastq_r1", type=str, default=None, help="Path to FASTQ R1 file (if shared across batch).")
    parser.add_argument("-r2", "--fastq_r2", type=str, default=None, help="Path to FASTQ R2 file for paired-end data (if shared across batch).")
    parser.add_argument("-w", "--wt_seq", type=str, default=None, help="Full wild-type reference amplicon sequence including spacers.")
    parser.add_argument("-t", "--tpe_seq", type=str, default=None, help="Full Twin prime edited reference amplicon sequence with 5' & 3' ends identical to wildtype reference amplicon.")
    parser.add_argument("-g", "--peg_spacers", type=str, default=None, help="Comma-separated pegRNA spacer sequences: <spacer A>,<spacer B>.")
    parser.add_argument("-rt", "--rt_templates", type=str, default=None, help="Comma-separated pegRNA reverse transcriptase templates: <RT template A>,<RT template B>.")
    parser.add_argument("-o", "--output_root", type=str, default=None, help="Parent output folder for batch results. Each sample is placed in its own subfolder.")
    parser.add_argument("-rcm", "--recoding_mode", action="store_true", help="Run in recoding mode when wild-type and twin prime edited sequences have the same length.")
    parser.add_argument("-ne", "--min_num_base_edits", type=int, default=None, help="Minimum number of base changes required for a read to be considered edited.")
    parser.add_argument("-dmas", "--default_min_aln_score", type=int, default=None, help="Default minimum homology score for a read to align to the compound reference amplicon.")
    parser.add_argument("-pfr", "--plot_full_reads", action="store_true", help="Display full read sequences in allele tables.")
    parser.add_argument("-ncda", "--no_collapse_displayed_alleles", action="store_true", help="Do not combine alleles that become identical in the displayed allele table window.")
    parser.add_argument("-naa", "--no_alignment_adjustments", action="store_true", help="Do not visually adjust homologies in the allele tables.")
    parser.add_argument("-ied", "--ignore_extraspacer_deletions", action="store_true", help="Classification ignores deletions occurring beyond the spacers (outside edit window).")
    parser.add_argument("-nf", "--no_figures", action="store_true", help="Skip all figures if only text outputs are desired.")
    parser.add_argument("-nsf", "--no_summary_figures", action="store_true", help="Skip summary barplots.")
    parser.add_argument("-nbf", "--no_per_base_figures", action="store_true", help="Skip per-base barplots.")
    parser.add_argument("-nmf", "--no_mutation_figures", action="store_true", help="Skip mutation barplots.")
    parser.add_argument("-pet", "--plot_extended_tables", action="store_true", help="Generate separate allele tables for each category.")
    parser.add_argument("-pdf", "--save_pdf", action="store_true", help="Save PDF versions of all plots.")
    parser.add_argument("-npng", "--no_save_png", action="store_true", help="Do not save PNG versions of all plots.")
    parser.add_argument("-mfa", "--min_frequency_alleles", type=float, default=None, help="Minimum percent read frequency required to report an allele in the allele tables.")
    parser.add_argument("-mnr", "--max_n_rows", type=int, default=None, help="Maximum number of allele rows to display in the allele tables by category.")
    parser.add_argument("-mna", "--max_n_alleles_to_write", type=int, default=None, help="Maximum number of alleles per category to write to the f7 text file.")
    parser.add_argument("-nrr", "--no_rerun", action="store_true", help="Don't rerun CRISPResso2 if a run using the same parameters has already been finished.")
    parser.add_argument("-kco", "--keep_crispresso_outputs", action="store_true", help="Don't delete CRISPResso2 output folders after analysis.")
    parser.add_argument("--crispresso_args", type=str, default=None, help='Additional arguments to pass to CRISPResso2 (wrapped in quotes, use equal sign); do not include --n_processes.')
    parser.add_argument("-coa", "--cleavage_offset_a", type=int, default=None, help="Cleavage offset for pegRNA spacer A.")
    parser.add_argument("-cob", "--cleavage_offset_b", type=int, default=None, help="Cleavage offset for pegRNA spacer B.")
    parser.add_argument("-p", "--n_processes", "--n_threads", dest="n_processes", type=int, default=1, metavar="N", help="Total process budget across batch samples and TwInsPEctor runs (default: 1).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose CRISPResso2 output.")
    parser.add_argument("--dry_run", action="store_true", help="Print the TwInsPEctor commands that would be executed without running them.")
    parser.add_argument("-V", "--version", action="version", version="%(prog)s 0.2.0")

    args = parser.parse_args(argv)

    if args.n_processes < 1:
        parser.error("--n_processes/--n_threads must be at least 1")

    return args


if __name__ == "__main__":
    main()
