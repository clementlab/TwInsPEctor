import os
import re
import sys
import html
import shutil
import zipfile
import argparse
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

# Recent CRISPresso2 update causes global style changes that break plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib import colors as colors_mpl
import matplotlib.ticker
# snapshot current rcParams before importing CRISPResso2
_RC_BEFORE_CRISP = matplotlib.rcParams.copy()
from CRISPResso2 import CRISPRessoShared
# restore rcParams to the snapshot to undo CRISPResso2 global style changes
matplotlib.rcParams.update(_RC_BEFORE_CRISP)


CATEGORY_COLORS = {
    "WT": "#b2182b",
    "WT Indel": "#d6604d",
    "Imperfect WT": "#f4a582",
    "Right Flap": "#d9d9d9",
    "Left Flap": "#9e9e9e",
    "Imperfect TPE": "#92c5de",
    "TPE Indel": "#4393c3",
    "Perfect TPE": "#2166ac",
    "Uncategorized": "#93de92"
}

CATEGORY_ORDER = [
    "WT", 
    "WT Indel",
    "Imperfect WT",
    "Right Flap",
    "Left Flap",
    "Imperfect TPE", 
    "TPE Indel",
    "Perfect TPE"
]

ALPHA = 0.4

BASE_COLORS = {
    "A": (127/255, 201/255, 127/255, ALPHA),
    "T": (190/255, 174/255, 212/255, ALPHA),
    "C": (253/255, 192/255, 134/255, ALPHA),
    "G": (255/255, 255/255, 153/255, ALPHA),
    "N": (200/255, 200/255, 200/255, ALPHA)
}


def get_folder_names(args):
    r1 = args.fastq_r1
    r2 = args.fastq_r2 if args.fastq_r2 else None
    pattern = r'([^/]+?)(?=(?:\.fastq|\.fq)?(?:\.gzip|\.gz|\.bz2|\.bz|\.xz|\.lzma)?$)'
    r1m = re.search(pattern, r1)
    r2m = re.search(pattern, r2) if r2 else None
    parent_folder = None
    crispresso_wt = None
    crispresso_tpe = None
    crispresso_composite_a = None
    crispresso_composite_b = None

    # If output_root provided use as parent for CRISPResso and TwInsPEctor results directories
    if args.output_root:
        parent_folder = os.path.join(os.getcwd(), args.output_root.rstrip("/"))
        # Mimic CRISPResso output folder naming conventions to get correct path
        if r1m and r2m:
                crispresso_wt = os.path.join(parent_folder, "CRISPResso_wt", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
                crispresso_tpe = os.path.join(parent_folder, "CRISPResso_tpe", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
                crispresso_composite_a = os.path.join(parent_folder, "CRISPResso_composite_a", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
                crispresso_composite_b = os.path.join(parent_folder, "CRISPResso_composite_b", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
        elif r1m and not r2m:
                crispresso_wt = os.path.join(parent_folder, "CRISPResso_wt", f"CRISPResso_on_{r1m.group(1)}")
                crispresso_tpe = os.path.join(parent_folder, "CRISPResso_tpe", f"CRISPResso_on_{r1m.group(1)}")
                crispresso_composite_a = os.path.join(parent_folder, "CRISPResso_composite_a", f"CRISPResso_on_{r1m.group(1)}")
                crispresso_composite_b = os.path.join(parent_folder, "CRISPResso_composite_b", f"CRISPResso_on_{r1m.group(1)}")
        else:   
            raise ValueError("Could not find fastq file(s).")
    # If output_root not provided, create one based on fastq file names
    else:
        if r1m and r2m:
            parent_folder = os.path.join(os.getcwd(), f"TwInsPEctor_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_wt = os.path.join(parent_folder, "CRISPResso_wt", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_tpe = os.path.join(parent_folder, "CRISPResso_tpe", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_composite_a = os.path.join(parent_folder, "CRISPResso_composite_a", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_composite_b = os.path.join(parent_folder, "CRISPResso_composite_b", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
        elif r1m and not r2m:
            parent_folder = os.path.join(os.getcwd(), f"TwInsPEctor_on_{r1m.group(1)}")
            crispresso_wt = os.path.join(parent_folder, "CRISPResso_wt", f"CRISPResso_on_{r1m.group(1)}")
            crispresso_tpe = os.path.join(parent_folder, "CRISPResso_tpe", f"CRISPResso_on_{r1m.group(1)}")
            crispresso_composite_a = os.path.join(parent_folder, "CRISPResso_composite_a", f"CRISPResso_on_{r1m.group(1)}")
            crispresso_composite_b = os.path.join(parent_folder, "CRISPResso_composite_b", f"CRISPResso_on_{r1m.group(1)}")
        else:
            raise ValueError("Could not find fastq file(s).")
        
    twinspector_results_folder = os.path.join(parent_folder, "TwInsPEctor_results")

    return parent_folder, crispresso_wt, crispresso_tpe, crispresso_composite_a, crispresso_composite_b, twinspector_results_folder


def get_spacer_seqs(peg_spacers_arg):
    spacers = [s.strip() for s in peg_spacers_arg.split(",")]
    if len(spacers) != 2:
        raise ValueError("pegRNA spacers must be entered as two comma-separated sequences")
    spacer_a, spacer_b = spacers
    spacer_a = spacer_a.upper()
    spacer_b = spacer_b.upper()

    return spacer_a, spacer_b


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    return ''.join(complement.get(base, base) for base in reversed(seq.upper()))


def analyze_references(wt_seq, tpe_seq, spacer_a, spacer_b, cleavage_offset_a, cleavage_offset_b, output_root, recoding_mode=False):
    """
    Validate sequence inputs, build composite references, return detailed info 
    Requires wt_seq and tpe_seq inputs share identical 5' and 3' anchor sequences.
    """
    # Find spacer A in WT sequence
    if wt_seq.find(spacer_a) != -1:
        spacer_a_start_wt = wt_seq.find(spacer_a)
        is_spacer_a_rc = False
    elif wt_seq.find(reverse_complement(spacer_a)) != -1:
        spacer_a_start_wt = wt_seq.find(reverse_complement(spacer_a))
        is_spacer_a_rc = True
    else:
        raise ValueError("Could not find pegRNA spacer A or its reverse complement in WT sequence")
    
    # Find spacer B in WT sequence
    if wt_seq.find(spacer_b) != -1:
        spacer_b_start_wt = wt_seq.find(spacer_b)
        is_spacer_b_rc = False
    elif wt_seq.find(reverse_complement(spacer_b)) != -1:
        spacer_b_start_wt = wt_seq.find(reverse_complement(spacer_b))
        is_spacer_b_rc = True
    else:
        raise ValueError("Could not find pegRNA spacer B or its reverse complement in WT sequence")
    
    # Check that the spacers are in the correct order
    if spacer_a_start_wt > spacer_b_start_wt:
        raise ValueError("pegRNA spacer B should be located downstream of pegRNA spacer A in WT sequence")
    
    # Find spacer A in TPE sequence
    if not is_spacer_a_rc:
        spacer_a_tpe = spacer_a[:cleavage_offset_a if cleavage_offset_a < 0 else len(spacer_a)]
        spacer_a_start_tpe = tpe_seq.find(spacer_a_tpe)
    else:
        spacer_a_tpe = spacer_a[-cleavage_offset_a if cleavage_offset_a < 0 else 0:]
        spacer_a_start_tpe = tpe_seq.find(reverse_complement(spacer_a)[:cleavage_offset_a if cleavage_offset_a < 0 else len(spacer_a)])
    if spacer_a_start_tpe == -1:
        raise ValueError("Could not find pegRNA spacer A in TPE sequence")

    # Find spacer B in TPE sequence
    if not is_spacer_b_rc:
        spacer_b_tpe = spacer_b[-cleavage_offset_b if cleavage_offset_b < 0 else 0:]
        spacer_b_start_tpe = tpe_seq.find(spacer_b_tpe)
    else:
        spacer_b_tpe = spacer_b[:cleavage_offset_b if cleavage_offset_b < 0 else len(spacer_b)]
        spacer_b_start_tpe = tpe_seq.find(reverse_complement(spacer_b)[-cleavage_offset_b if cleavage_offset_b < 0 else 0:])
    if spacer_b_start_tpe == -1:
        raise ValueError("Could not find pegRNA spacer B in TPE sequence")

    # Get positions of nick sites in WT and TPE sequences
    spacer_a_nick_site_wt = (spacer_a_start_wt + len(spacer_a) + cleavage_offset_a - 1)
    spacer_b_nick_site_wt = (spacer_b_start_wt - cleavage_offset_b - 1)
    spacer_a_nick_site_tpe = (spacer_a_start_tpe + len(spacer_a) + cleavage_offset_a - 1)
    spacer_b_nick_site_tpe = spacer_b_start_tpe - 1

    # Extract sequence components
    prefix_seq = wt_seq[:spacer_a_nick_site_wt + 1]
    suffix_seq = wt_seq[spacer_b_nick_site_wt + 1:]
    wt_deleted_seq = wt_seq[spacer_a_nick_site_wt + 1:spacer_b_nick_site_wt + 1]
    tpe_inserted_seq = tpe_seq[spacer_a_nick_site_tpe + 1:spacer_b_nick_site_tpe + 1]

    # Build composite reference sequences
    composite_a_ref_seq = prefix_seq + wt_deleted_seq + tpe_inserted_seq + suffix_seq
    composite_b_ref_seq = prefix_seq + tpe_inserted_seq + wt_deleted_seq + suffix_seq
    
    # Align standard reference sequences
    wt_aln_seq_a = prefix_seq + wt_deleted_seq + len(tpe_inserted_seq) * '-' + suffix_seq
    wt_aln_seq_b = prefix_seq + len(tpe_inserted_seq) * '-' + wt_deleted_seq + suffix_seq
    tpe_aln_seq_a = prefix_seq + len(wt_deleted_seq) * '-' + tpe_inserted_seq + suffix_seq
    tpe_aln_seq_b = prefix_seq + tpe_inserted_seq + len(wt_deleted_seq) * '-' + suffix_seq
    
    # Check lengths
    if not (len(composite_a_ref_seq) == len(wt_aln_seq_a) == len(tpe_aln_seq_a) == len(composite_b_ref_seq) == len(wt_aln_seq_b) == len(tpe_aln_seq_b)):
        raise ValueError("Composite references, WT alignments, and Twin alignments are not the same length")
    
    # Find spacer A in Composite A reference sequence
    if not is_spacer_a_rc:
        spacer_a_start_composite_a = composite_a_ref_seq.find(spacer_a)
    else:
        spacer_a_start_composite_a = composite_a_ref_seq.find(reverse_complement(spacer_a))

    # Find spacer B in Composite A reference sequence
    if not is_spacer_b_rc:
        spacer_b_start_composite_a = composite_a_ref_seq.find(spacer_b_tpe)
    else:
        spacer_b_start_composite_a = composite_a_ref_seq.find(reverse_complement(spacer_b_tpe))

    # Find spacer A in Composite B reference sequence
    if not is_spacer_a_rc:
        spacer_a_start_composite_b = composite_b_ref_seq.find(spacer_a_tpe)
    else:
        spacer_a_start_composite_b = composite_b_ref_seq.find(reverse_complement(spacer_a_tpe))

    # Find spacer B in Composite B reference sequence
    if not is_spacer_b_rc:
        spacer_b_start_composite_b = composite_b_ref_seq.find(spacer_b)
    else:
        spacer_b_start_composite_b = composite_b_ref_seq.find(reverse_complement(spacer_b))

    # pegRNA positions in reference sequences
    pegRNA_intervals_wt = [[spacer_a_start_wt, spacer_a_start_wt + len(spacer_a) - 1], [spacer_b_start_wt, spacer_b_start_wt + len(spacer_b) - 1]]
    pegRNA_intervals_tpe = [[spacer_a_start_tpe, spacer_a_start_tpe + len(spacer_a_tpe) - 1], [spacer_b_start_tpe, spacer_b_start_tpe + len(spacer_b_tpe) - 1]]
    pegRNA_intervals_composite_a = [[spacer_a_start_composite_a, spacer_a_start_composite_a + len(spacer_a) - 1], [spacer_b_start_composite_a, spacer_b_start_composite_a + len(spacer_b_tpe) - 1]]
    pegRNA_intervals_composite_b = [[spacer_a_start_composite_b, spacer_a_start_composite_b + len(spacer_a_tpe) - 1], [spacer_b_start_composite_b, spacer_b_start_composite_b + len(spacer_b) - 1]]

    # Get positions of nick sites in Composite reference sequences
    spacer_b_nick_site_comp_a = spacer_b_start_composite_a - 1
    spacer_b_nick_site_comp_b = spacer_b_start_composite_b - cleavage_offset_b - 1

    # Get deletion and insertion info for Composite A reference sequence
    composite_a_del_start = spacer_a_nick_site_wt + 1
    composite_a_del_end = composite_a_del_start + len(wt_deleted_seq) - 1
    composite_a_ins_start = composite_a_del_end + 1
    composite_a_ins_end = composite_a_ins_start + len(tpe_inserted_seq) - 1

    # Get deletion and insertion info for Composite B reference sequence
    composite_b_ins_start = spacer_a_nick_site_tpe + 1
    composite_b_ins_end = composite_b_ins_start + len(tpe_inserted_seq) - 1
    composite_b_del_start = composite_b_ins_end + 1
    composite_b_del_end = composite_b_del_start + len(wt_deleted_seq) - 1

    # Special recoding mode composite references for comparing base substitutions
    composite_wt = None
    composite_tpe = None
    if recoding_mode:
        composite_wt = prefix_seq + wt_deleted_seq + wt_deleted_seq + suffix_seq
        composite_tpe = prefix_seq + tpe_inserted_seq + tpe_inserted_seq + suffix_seq

    # Special replacement mode standard references for comparing bases
    tpe_seq_replacement_bp_changes = None
    wt_seq_replacement_bp_changes = None
    if not recoding_mode:
        tpe_seq_replacement_bp_changes = prefix_seq + len(wt_deleted_seq) * '-' + suffix_seq
        wt_seq_replacement_bp_changes = prefix_seq + len(tpe_inserted_seq) * '-' + suffix_seq

    with open(
        os.path.join(output_root, "c4.reference_sequences.txt"), "w"
    ) as fout:
        fout.write(f"@Sequence inputs\n")
        fout.write(f">Wildtype reference sequence\n{wt_seq}\n")
        fout.write(f">TwinPE reference sequence\n{tpe_seq}\n")
        fout.write(f">pegRNA spacer a sequence\n{spacer_a}\n")
        fout.write(f">pegRNA spacer b sequence\n{spacer_b}\n\n")
        fout.write(f"@Composite A alignments\n")
        fout.write(f">Wildtype reference sequence alignment\n{wt_aln_seq_a}\n")
        fout.write(f">Composite A reference sequence\n{composite_a_ref_seq}\n")
        fout.write(f">TwinPE reference sequence alignment\n{tpe_aln_seq_a}\n\n")
        fout.write(f"@Composite B alignments\n")
        fout.write(f">Wildtype reference sequence alignment\n{wt_aln_seq_b}\n")
        fout.write(f">Composite B reference sequence\n{composite_b_ref_seq}\n")
        fout.write(f">TwinPE reference sequence alignment\n{tpe_aln_seq_b}\n\n")
        if recoding_mode:
            fout.write(f"@Recoding mode base change reference sequences\n")
            fout.write(f">Composite WT reference sequence\n{composite_wt}\n")
            fout.write(f">Composite TPE reference sequence\n{composite_tpe}\n\n")
        fout.write(f"@Modified pegRNA spacer sequences\n")
        fout.write(f">pegRNA spacer a for twinPE and composite b reference sequences\n{spacer_a_tpe}\n")
        fout.write(f">pegRNA spacer b for twinPE and composite a reference sequences\n{spacer_b_tpe}\n")

    reference_info = {
        "composite_a_ref_seq": composite_a_ref_seq, 
        "wt_aln_seq_a": wt_aln_seq_a,
        "tpe_aln_seq_a": tpe_aln_seq_a,
        "composite_b_ref_seq": composite_b_ref_seq,
        "wt_aln_seq_b": wt_aln_seq_b,
        "tpe_aln_seq_b": tpe_aln_seq_b,
        "composite_wt": composite_wt,
        "composite_tpe": composite_tpe,
        "spacer_a_wt": spacer_a,
        "spacer_b_wt": spacer_b,
        "spacer_a_tpe": spacer_a_tpe,
        "spacer_b_tpe": spacer_b_tpe,
        "spacer_a_composite_a": spacer_a,
        "spacer_b_composite_a": spacer_b_tpe,
        "spacer_a_composite_b": spacer_a_tpe,
        "spacer_b_composite_b": spacer_b,
        "spacer_a_start_wt": spacer_a_start_wt,
        "spacer_b_start_wt": spacer_b_start_wt,
        "spacer_a_start_tpe": spacer_a_start_tpe,
        "spacer_b_start_tpe": spacer_b_start_tpe, 
        "spacer_a_start_composite_a": spacer_a_start_composite_a,
        "spacer_b_start_composite_a": spacer_b_start_composite_a,
        "spacer_a_start_composite_b": spacer_a_start_composite_b,
        "spacer_b_start_composite_b": spacer_b_start_composite_b,
        "pegRNA_intervals_wt": pegRNA_intervals_wt,
        "pegRNA_intervals_tpe": pegRNA_intervals_tpe,
        "pegRNA_intervals_composite_a": pegRNA_intervals_composite_a,
        "pegRNA_intervals_composite_b": pegRNA_intervals_composite_b,
        "is_spacer_a_rc": is_spacer_a_rc,
        "is_spacer_b_rc": is_spacer_b_rc,
        "cleavage_offset_a": cleavage_offset_a,
        "cleavage_offset_b": cleavage_offset_b, 
        "cut_points_wt": [spacer_a_nick_site_wt, spacer_b_nick_site_wt],
        "cut_points_tpe": [spacer_a_nick_site_tpe, spacer_b_nick_site_tpe],
        "cut_points_composite_a": [spacer_a_nick_site_wt, spacer_b_nick_site_comp_a],
        "cut_points_composite_b": [spacer_a_nick_site_tpe, spacer_b_nick_site_comp_b], 
        "composite_a_del_start": composite_a_del_start,
        "composite_a_del_end": composite_a_del_end,
        "composite_a_ins_start": composite_a_ins_start,
        "composite_a_ins_end": composite_a_ins_end,
        "composite_b_del_start": composite_b_del_start,
        "composite_b_del_end": composite_b_del_end,
        "composite_b_ins_start": composite_b_ins_start,
        "composite_b_ins_end": composite_b_ins_end,
        "inserted_seq": tpe_inserted_seq, 
        "deleted_seq": wt_deleted_seq,
        "ins_region_len": len(tpe_inserted_seq),
        "del_region_len": len(wt_deleted_seq), 
        "tpe_seq_replacement_bp_changes": tpe_seq_replacement_bp_changes, 
        "wt_seq_replacement_bp_changes": wt_seq_replacement_bp_changes
    }

    return reference_info


def get_crispresso_command(args=None, ref_seq=None, ref_name=None, spacer_a=None, spacer_b=None, crispresso_output_folder=None, twinspector_results_folder=None, append=True):
    # May need to wrap in more CRISPResso options
    cmd = [
        "CRISPResso",
        "--fastq_r1", args.fastq_r1,
        "--amplicon_seq", f"{ref_seq}",
        "--amplicon_name", f"{ref_name}",
        "--guide_seq", f"{spacer_a},{spacer_b}", 
        # "--guide_name", f"pegRNA a,pegRNA b",
        "--default_min_aln_score", "0", 
        # "--min_frequency_alleles_around_cut_to_plot", str(args.min_frequency_alleles), 
        # "--max_rows_alleles_around_cut_to_plot", str(args.max_n_rows), 
        "--output_folder", os.path.dirname(crispresso_output_folder),
        "--write_detailed_allele_table",
        "--n_processes", "1"
    ]

    if args.fastq_r2:
        cmd.extend(["--fastq_r2", args.fastq_r2])
    if args.no_rerun:
        cmd.append("--no_rerun")
    if args.trim_string:
        cmd.extend(["--trim_sequences"])
        cmd.extend(["--fastp_options_string", f"{args.trim_string}"])
    if args.fastp_command:
        cmd.extend(["--fastp_command", args.fastp_command])
    if args.debug:
        cmd.append("--debug")

    with open(
        os.path.join(twinspector_results_folder, "c5.crispresso2_commands.txt"), 
        "a" if append else "w" 
    ) as fout:
        fout.write(" ".join(cmd) + "\n")

    return cmd


def run_crispresso_command(cmd, verbose=False):
    if verbose:
        print("Running CRISPResso2 with command:\n", " ".join(cmd), "\n")
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_crispresso_commands_parallel(crispresso_tasks, verbose=False):
    # Using CRISPRessoBatch with multiple processes instead of this wrapper is likely best
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(run_crispresso_command, cmd, verbose) for cmd in crispresso_tasks]
        for future in as_completed(futures):
            future.result()


def get_allele_df_keys(df):
    # Create sequence_key and take the lexicographically smaller of the forward and reverse complement
    df = df.copy()
    df['sequence_key_fw'] = df['Aligned_Sequence'].str.replace('-', '', regex=False)
    df['sequence_key_rc'] = df['sequence_key_fw'].apply(reverse_complement)
    df['sequence_key'] = df[['sequence_key_fw', 'sequence_key_rc']].min(axis=1)

    return df


def load_allele_table(folder, crispresso_info, suffix):
    zip_path = os.path.join(
        folder,
        crispresso_info["running_info"]["allele_frequency_table_zip_filename"],
    )
    table_name = crispresso_info["running_info"]["allele_frequency_table_filename"]

    with zipfile.ZipFile(zip_path) as z:
        with z.open(table_name) as zf:
            df = pd.read_csv(zf, sep="\t")
            ref_name = crispresso_info["results"]["ref_names"][0]
            ref = crispresso_info["results"]["refs"][ref_name]
            pegrna_info = {
                "pegRNA_cut_points": ref["sgRNA_cut_points"],
                # "pegRNA_plot_cut_points": ref["sgRNA_plot_cut_points"],
                "pegRNA_intervals": ref["sgRNA_intervals"],
                "pegRNA_mismatches": ref["sgRNA_mismatches"],
                # "pegRNA_names": ref["sgRNA_names"], 
            }

    df = get_allele_df_keys(df)

    df_merged = (
        df[
            [
                "sequence_key",
                "#Reads",
                "Aligned_Sequence",
                "Reference_Sequence",
                "Aligned_Reference_Scores", 
                "%Reads"
            ]
        ]
        .rename(columns={
            "#Reads": f"#Reads_{suffix}",
            "Aligned_Sequence": f"Aligned_Sequence_{suffix}",
            "Reference_Sequence": f"Reference_Sequence_{suffix}",
            "Aligned_Reference_Scores": f"Aligned_Reference_Score_{suffix}", 
            "%Reads": f"%Reads_{suffix}",
        })
    )

    return df_merged, pegrna_info

def merge_crispresso_allele_tables(crispresso_wt=None, crispresso_tpe=None, crispresso_composite_a=None, crispresso_composite_b=None):
    crispresso2_wt_info = CRISPRessoShared.load_crispresso_info(crispresso_wt)
    crispresso2_tpe_info = CRISPRessoShared.load_crispresso_info(crispresso_tpe)
    crispresso2_composite_a_info = CRISPRessoShared.load_crispresso_info(crispresso_composite_a)
    crispresso2_composite_b_info = CRISPRessoShared.load_crispresso_info(crispresso_composite_b)
    
    df_alleles_wt, pegrna_info_wt = load_allele_table(crispresso_wt, crispresso2_wt_info, "wt")
    df_alleles_tpe, pegrna_info_tpe = load_allele_table(crispresso_tpe, crispresso2_tpe_info, "tpe")
    df_alleles_comp_a, pegrna_info_comp_a = load_allele_table(crispresso_composite_a, crispresso2_composite_a_info, "comp_a")
    df_alleles_comp_b, pegrna_info_comp_b = load_allele_table(crispresso_composite_b, crispresso2_composite_b_info, "comp_b")

    # pegrna_info = {
    #     "wt": pegrna_info_wt,
    #     "tpe": pegrna_info_tpe,
    #     "comp_a": pegrna_info_comp_a,
    #     "comp_b": pegrna_info_comp_b,
    # }

    # Check for duplicate sequence keys
    for name, df in {
        "wt": df_alleles_wt,
        "tpe": df_alleles_tpe,
        "comp_a": df_alleles_comp_a,
        "comp_b": df_alleles_comp_b,
    }.items():
        if df["sequence_key"].duplicated().any():
            raise ValueError(f"Duplicate sequence_key in {name}")
    
    df_merged = (
        df_alleles_wt
        .merge(df_alleles_tpe, on='sequence_key', how='outer', validate='one_to_one')
        .merge(df_alleles_comp_a, on='sequence_key', how='outer', validate='one_to_one')
        .merge(df_alleles_comp_b, on='sequence_key', how='outer', validate='one_to_one')
    )

    # Ensure no missing values in dataframe
    if df_merged.isnull().values.any():
        raise ValueError("Merged allele dataframe contains missing values.")


    # Ensure #Reads columns agree across references for each sequence_key
    reads_cols = [
        "#Reads_wt",
        "#Reads_tpe",
        "#Reads_comp_a",
        "#Reads_comp_b",
    ]
    mismatch = df_merged[reads_cols].nunique(axis=1) > 1
    if mismatch.any():
        raise ValueError(
            f"Read counts do not agree for sequence_keys:\n"
            f"{df_merged.loc[mismatch, ['sequence_key'] + reads_cols]}"
        )

    # Drop duplicate #Reads columns, keep only one and rename to "#Reads"
    df_merged = df_merged.rename(columns={"#Reads_wt": "#Reads"})
    df_merged = df_merged.rename(columns={"%Reads_wt": "%Reads"})
    df_merged = df_merged.drop(columns=["#Reads_tpe", "#Reads_comp_a", "#Reads_comp_b"])
    df_merged = df_merged.drop(columns=["%Reads_tpe", "%Reads_comp_a", "%Reads_comp_b"])

    
    return df_merged  # , pegrna_info


def get_refpos_values(ref_aln_seq, read_aln_seq):
    """
    Given a reference alignment this returns a dictionary such that refpos_dict[ind] is the value of the read at the position corresponding to the ind'th base in the reference
    Any additional bases in the read (gaps in the ref) are assigned to the first position of the ref (i.e. refpos_dict[0])
    For other additional bases in the ref (gaps in the read), the value is appended to the last position of the ref that had a non-gap base (to the left)
    For example:
    ref_seq =  '--A-TGC-'
    read_seq = 'GGAGTCGA'
    get_refpos_values(ref_seq, read_seq)
    {0: 'GGAG', 1: 'T', 2: 'C', 3: 'GA'}
    Args:
    - ref_aln_seq: str, reference alignment sequence
    - read_aln_seq: str, read alignment sequence
    Returns:
    - refpos_dict: dict, dictionary such that refpos_dict[ind] is the value of the read at the position corresponding to the ind'th base in the reference
    """
    refpos_dict = defaultdict(str)

    # First, if there are insertions in read, add those to the first position in ref
    if ref_aln_seq[0] == '-':
        aln_index = 0
        read_start_bases = ""
        while aln_index < len(ref_aln_seq) and ref_aln_seq[aln_index] == '-':
            read_start_bases += read_aln_seq[aln_index]
            aln_index += 1
        refpos_dict[0] = read_start_bases
        ref_aln_seq = ref_aln_seq[aln_index:]
        read_aln_seq = read_aln_seq[aln_index:]
        
    ref_pos = 0
    last_nongap_ref_pos = 0
    for ind in range(len(ref_aln_seq)):
        ref_base = ref_aln_seq[ind]
        read_base = read_aln_seq[ind]
        if ref_base == '-':
            refpos_dict[last_nongap_ref_pos] += read_base
        else:
            refpos_dict[ref_pos] += read_base
            last_nongap_ref_pos = ref_pos
            ref_pos += 1
    return refpos_dict


def get_mutations(allele_map=None, ref_seq=None, cut_points=None):
    """
    Determines sub, del, ins information from allele map.
    """
    # Get substitution positions and check if within edit window
    all_sub_pos = [pos for pos, base in allele_map.items() if len(base) == 1 and base != '-' and base != ref_seq[pos]]
    sub_between_cuts = False
    for sub_pos in all_sub_pos:
        if sub_pos >= cut_points[0]+1 and sub_pos <= cut_points[1]:  # +1 to cut_points[0] since cut is after the base position
            sub_between_cuts = True
            break

    # Get deletion positions and check if within edit window
    all_del_pos = [pos for pos, base in allele_map.items() if base == '-']
    del_between_cuts = False
    for del_pos in all_del_pos:
        if del_pos >= cut_points[0]+1 and del_pos <= cut_points[1]:
            del_between_cuts = True
            break

    # Get insertion positions
    all_ins_pos = [pos for pos, base in allele_map.items() if len(base) > 1]

    has_substitutions = len(all_sub_pos) > 0
    has_deletions = len(all_del_pos) > 0
    has_insertions = len(all_ins_pos) > 0

    return all_sub_pos, sub_between_cuts, all_del_pos, del_between_cuts, all_ins_pos, has_substitutions, has_deletions, has_insertions


def get_replacement_base_changes(comp_ref_seq=None, wt_aln_seq=None, twin_aln_seq=None):
    bp_changes_arr = []
    for idx in range(len(comp_ref_seq)):
        wt_base = wt_aln_seq[idx]
        twin_base = twin_aln_seq[idx]
        if wt_base != '-' and twin_base == '-':
            bp_changes_arr.append((idx, wt_base, twin_base))
        elif wt_base == '-' and twin_base != '-':
            bp_changes_arr.append((idx, wt_base, twin_base))
        elif wt_base != twin_base:
            raise ValueError('Substitution detected in Replacement mode.')

    return bp_changes_arr


def get_recoding_base_changes(wt_seq=None, twin_seq=None, composite_wt=None, composite_tpe=None):
    std_bp_changes_arr = []
    for idx, (wt_base, twin_base) in enumerate(zip(wt_seq, twin_seq)):
        if wt_base != twin_base:
            std_bp_changes_arr.append((idx, wt_base, twin_base))

    comp_bp_changes_arr = []
    for idx, (wt_base, twin_base) in enumerate(zip(composite_wt, composite_tpe)):
        if wt_base != twin_base:
            comp_bp_changes_arr.append((idx, wt_base, twin_base))
    
    return std_bp_changes_arr, comp_bp_changes_arr


def get_allele_match_array(bp_changes_arr, allele_map, del_start, del_end, ins_start, ins_end):
    match_arr = ["0"] * len(bp_changes_arr)
    # full_ins_arr = ["0"] * len(bp_changes_arr)
    # full_del_arr = ["0"] * len(bp_changes_arr)
    # full_sub_arr = ["0"] * len(bp_changes_arr)
    for ind, (comp_ind, wt_base, tpe_base) in enumerate(bp_changes_arr):
        allele_base = allele_map.get(comp_ind, "")
        if allele_base == wt_base:
            match_arr[ind] = "W" # matches WT base
        elif allele_base == tpe_base:
            match_arr[ind] = "T" # matches TwinPE base
        elif len(allele_base) > 1:
            # full_ins_arr[ind] = "I" # non-programmed insertion
            if allele_base[0] == wt_base:
                match_arr[ind] = "WI" # matches WT base with insertion after
            elif allele_base[0] == tpe_base:
                match_arr[ind] = "TI" # matches TwinPE base with insertion after
            elif allele_base[0] in {"A", "C", "G", "T"}:
                match_arr[ind] = "SI" # non-programmed substitution with insertion after
                # full_sub_arr[ind] = "S" # non-programmed substitution with insertion after
            else:
                match_arr[ind] = "NI" # ambiguous base with insertion after
        elif allele_base in {"A", "C", "G", "T"}:
            match_arr[ind] = "S" # non-programmed substitution
            # full_sub_arr[ind] = "S" # non-programmed substitution
        elif allele_base == "-":
            # full_del_arr[ind] = "D" # deletion relative to wt/tpe references aligned to composite references (replacement - not possible) or relative to composite wt/tpe references (recoding - possible)
            match_arr[ind] = "D" # "-"
        else:
            match_arr[ind] = "N" # ambiguous base

    # Split match_arr by del and ins regions for plotting
    del_match_arr = []
    ins_match_arr = []
    for match, (comp_ind, wt_base, tpe_base) in zip(match_arr, bp_changes_arr):
        if comp_ind >= del_start and comp_ind <= del_end:
            del_match_arr.append(match)
        # elif comp_ind >= ins_start and comp_ind <= ins_end:
        else:
            ins_match_arr.append(match)

    return match_arr, ins_match_arr, del_match_arr #  , full_ins_arr, full_sub_arr, full_del_arr


def check_indel_positions(all_insertion_left_positions, all_deletion_positions, del_start, del_end, ins_start, ins_end, ignore_extraspacer_deletions, pegRNA_intervals):
    has_any_ins_byproduct = False
    has_del_in_spacer_window = False
    has_any_del_byproduct = False

    # Check for insertions anywhere in read
    if all_insertion_left_positions != []:
        has_any_ins_byproduct = True
    else:
        if del_start < ins_start:
            edit_range = range(del_start, ins_end + 1)
        else:
            edit_range = range(ins_start, del_end + 1)
        # Ignore deletions beyond spacers if flagged
        if ignore_extraspacer_deletions:
            for del_ind in all_deletion_positions:
                if del_ind >= pegRNA_intervals[0][0] and del_ind <= pegRNA_intervals[1][1] and del_ind not in edit_range:
                    has_del_in_spacer_window = True
                    break
        # Check for deletions anywhere outside of the edit region if not flagged
        else:
            for del_ind in all_deletion_positions:
                if del_ind not in edit_range:
                    has_any_del_byproduct = True
                    break

    # Set has_indel based on flag
    if ignore_extraspacer_deletions:
        has_indel = has_any_ins_byproduct or has_del_in_spacer_window
    else:
        has_indel = has_any_ins_byproduct or has_any_del_byproduct

    return has_indel


def resolve_composite_categories(
    category_a,
    score_a,
    category_b,
    score_b,
    flap_score_delta=99,
):
    """
    Resolve conflicting Composite A and Composite B classifications.

    Returns
    -------
    tuple
        (winning_category, classified_by)
    """

    CATEGORY_PRIORITY = {
        "Perfect_TPE": 1,
        "TPE_Indel": 2,
        "WT": 3,
        "WT_Indel": 4,
        "Imperfect_TPE": 5,
        "Left_Flap": 6,
        "Right_Flap": 6,
        "Imperfect_WT": 7,
    }

    DEFAULT_PRIORITY = 99

    # Left flap vs Right flap disagreement
    flap_disagreement = (
        (category_a == "Left_Flap" and category_b == "Right_Flap")
        or
        (category_a == "Right_Flap" and category_b == "Left_Flap")
    )

    if flap_disagreement and abs(score_a - score_b) <= flap_score_delta:
        return "Imperfect_TPE", "Composite_A&B"

    # Missing categories
    if not category_a:
        return category_b, "Composite_B"

    if not category_b:
        return category_a, "Composite_A"

    # Rank categories
    rank_a = CATEGORY_PRIORITY.get(category_a, DEFAULT_PRIORITY)
    rank_b = CATEGORY_PRIORITY.get(category_b, DEFAULT_PRIORITY)

    # Lower rank wins
    if rank_b < rank_a:
        return category_b, "Composite_B"

    if rank_a < rank_b:
        return category_a, "Composite_A"

    # Tie -> higher alignment score wins
    if score_b > score_a:
        return category_b, "Composite_B"

    return category_a, "Composite_A"


def categorize_alleles(
        df_merged=None, 
        wt_seq=None, 
        tpe_seq=None, 
        reference_info=None, 
        # pegrna_info=None, 
        num_changes_to_check=2, 
        ignore_extraspacer_deletions=False, 
        default_min_aln_score=30, 
        recoding_mode=False
    ):
    comp_a_ref_seq=reference_info["composite_a_ref_seq"]
    comp_b_ref_seq=reference_info["composite_b_ref_seq"]
    wt_aln_seq_comp_a=reference_info["wt_aln_seq_a"]
    tpe_aln_seq_comp_a=reference_info["tpe_aln_seq_a"]
    wt_aln_seq_comp_b=reference_info["wt_aln_seq_b"]
    tpe_aln_seq_comp_b=reference_info["tpe_aln_seq_b"]
    composite_wt=reference_info["composite_wt"]
    composite_tpe=reference_info["composite_tpe"]

    # Check for agreement between reference_info and pegrna_info
    # if reference_info["cut_points_wt"] != pegrna_info["wt"]["pegRNA_cut_points"]:
    #     raise ValueError("Inconsistent cut points between reference_info and pegrna_info for WT")
    # if reference_info["cut_points_tpe"] != pegrna_info["tpe"]["pegRNA_cut_points"]:
    #     raise ValueError("Inconsistent cut points between reference_info and pegrna_info for TPE")
    # if reference_info["cut_points_composite_a"] != pegrna_info["comp_a"]["pegRNA_cut_points"]:
    #     raise ValueError("Inconsistent cut points between reference_info and pegrna_info for Composite A")
    # if reference_info["cut_points_composite_b"] != pegrna_info["comp_b"]["pegRNA_cut_points"]:
    #     raise ValueError("Inconsistent cut points between reference_info and pegrna_info for Composite B")

    # Drop alleles with insufficient alignment scores
    aln_score_cols = [
        "Aligned_Reference_Score_wt",
        "Aligned_Reference_Score_tpe",
        "Aligned_Reference_Score_comp_a",
        "Aligned_Reference_Score_comp_b",
    ]
    keep_mask = df_merged[aln_score_cols].apply(
        lambda row: any(pd.notna(val) and val >= default_min_aln_score for val in row),
        axis=1,
    )
    df_merged = df_merged.loc[keep_mask].copy()

    df_merged['Category_wt'] = ""
    df_merged['Category_tpe'] = ""
    df_merged['Category_comp_a'] = ""
    df_merged['Category_comp_b'] = ""
    df_merged['tpe_ins_match_arr'] = ""
    df_merged['comp_a_ins_match_arr'] = ""
    df_merged['comp_b_ins_match_arr'] = ""
    df_merged['wt_del_match_arr'] = ""
    df_merged['comp_a_del_match_arr'] = ""
    df_merged['comp_b_del_match_arr'] = ""

    # Base changes
    if recoding_mode:
        # Uses composite wt/tpe references (and standard wt/tpe references for plotting only)
        std_bp_changes_arr, comp_bp_changes_arr = get_recoding_base_changes(wt_seq, tpe_seq, composite_wt, composite_tpe)
        bp_changes_arrs = {"std_bp_changes_arr": std_bp_changes_arr, "comp_bp_changes_arr": comp_bp_changes_arr}
    else:
        # Uses composite a/b references (and altered standard wt/tpe references that do not directly compare bases for plotting only)
        comp_a_bp_changes_arr = get_replacement_base_changes(comp_a_ref_seq, wt_aln_seq_comp_a, tpe_aln_seq_comp_a)
        comp_b_bp_changes_arr = get_replacement_base_changes(comp_b_ref_seq, wt_aln_seq_comp_b, tpe_aln_seq_comp_b)
        wt_bp_changes_arr = get_replacement_base_changes(comp_ref_seq=wt_seq, wt_aln_seq=wt_seq, twin_aln_seq=reference_info["tpe_seq_replacement_bp_changes"])
        tpe_bp_changes_arr = get_replacement_base_changes(comp_ref_seq=tpe_seq, wt_aln_seq=reference_info["wt_seq_replacement_bp_changes"], twin_aln_seq=tpe_seq)
        bp_changes_arrs = {"comp_a_bp_changes_arr": comp_a_bp_changes_arr, "comp_b_bp_changes_arr": comp_b_bp_changes_arr, "wt_bp_changes_arr": wt_bp_changes_arr, "tpe_bp_changes_arr": tpe_bp_changes_arr}

    for idx, allele in df_merged.iterrows():
        # Classify perfect WT alleles and some Imperfect WT alleles using WT reference alignment
        wt_seq_aln_allele = allele.Reference_Sequence_wt
        allele_seq_aln_wt = allele.Aligned_Sequence_wt

        wt_map = get_refpos_values(wt_seq_aln_allele, allele_seq_aln_wt)

        if recoding_mode:
            wt_del_match_arr, _, _ = get_allele_match_array(std_bp_changes_arr, wt_map, del_start=0, del_end=0, ins_start=0, ins_end=0)           
        else:
            wt_del_match_arr, _, _ = get_allele_match_array(wt_bp_changes_arr, wt_map, del_start=0, del_end=0, ins_start=0, ins_end=0)           

        wt_all_sub_pos, wt_sub_between_cuts, wt_all_del_pos, wt_del_between_cuts, wt_all_ins_pos, wt_has_substitutions, wt_has_deletions, wt_has_insertions = get_mutations(wt_map, wt_seq, reference_info["cut_points_wt"])

        if ignore_extraspacer_deletions:
            wt_has_all_extraspacer_deletions = all(x < reference_info["spacer_a_start_wt"] or x > (reference_info["spacer_b_start_wt"]+len(reference_info["spacer_b_wt"])) for x in wt_all_del_pos)

        if allele['Aligned_Reference_Score_wt'] == 100.0:
            df_merged.at[idx, 'Category_wt'] = 'WT'
        elif not wt_has_insertions and not wt_has_deletions and not wt_sub_between_cuts:  # WTs with substitutions outside of cut sites, can add check for all wt bases in edit window if needed
            df_merged.at[idx, 'Category_wt'] = 'WT'
        elif ignore_extraspacer_deletions and not wt_has_insertions and len(wt_all_sub_pos) <= 1 and wt_has_deletions and wt_has_all_extraspacer_deletions:
            if not wt_sub_between_cuts:
                df_merged.at[idx, 'Category_wt'] = 'WT'
            else:
                df_merged.at[idx, 'Category_wt'] = 'Imperfect_WT'
        elif not wt_has_insertions and not wt_has_deletions and len(wt_all_sub_pos) <= 1:
            if not wt_sub_between_cuts:
                df_merged.at[idx, 'Category_wt'] = 'WT'
            else:
                df_merged.at[idx, 'Category_wt'] = 'Imperfect_WT'

        # Classify perfect TPE alleles and some Imperfect TPE alleles using TPE reference alignment
        tpe_seq_aln_allele = allele.Reference_Sequence_tpe
        allele_seq_aln_tpe = allele.Aligned_Sequence_tpe

        tpe_map = get_refpos_values(tpe_seq_aln_allele, allele_seq_aln_tpe)

        if recoding_mode:
            tpe_ins_match_arr, _, _ = get_allele_match_array(std_bp_changes_arr, tpe_map, del_start=0, del_end=0, ins_start=0, ins_end=0)           
        else:
            tpe_ins_match_arr, _, _ = get_allele_match_array(tpe_bp_changes_arr, tpe_map, del_start=0, del_end=0, ins_start=0, ins_end=0)

        tpe_all_sub_pos, tpe_sub_between_cuts, tpe_all_del_pos, tpe_del_between_cuts, tpe_all_ins_pos, tpe_has_substitutions, tpe_has_deletions, tpe_has_insertions = get_mutations(tpe_map, tpe_seq, reference_info["cut_points_tpe"])

        if ignore_extraspacer_deletions:
            tpe_has_all_extraspacer_deletions = all(x < reference_info["spacer_a_start_tpe"] or x > (reference_info["spacer_b_start_tpe"]+len(reference_info["spacer_b_tpe"])) for x in tpe_all_del_pos)

        if allele['Aligned_Reference_Score_tpe'] == 100.0:
            df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'
        elif not tpe_has_insertions and not tpe_has_deletions and not tpe_sub_between_cuts:  # TPEs with substitutions outside of cut sites, can add check for all tpe bases in edit window if needed
            df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'
        elif ignore_extraspacer_deletions and not tpe_has_insertions and len(tpe_all_sub_pos) <= 1 and tpe_has_deletions and tpe_has_all_extraspacer_deletions:
            if not tpe_sub_between_cuts:
                df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'
            else:
                df_merged.at[idx, 'Category_tpe'] = 'Imperfect_TPE'
        elif not tpe_has_insertions and not tpe_has_deletions and len(tpe_all_sub_pos) <= 1:
            if not tpe_sub_between_cuts:
                df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'
            else:
                df_merged.at[idx, 'Category_tpe'] = 'Imperfect_TPE'

        # Classify all alleles using Composite A reference alignment
        comp_a_seq_aln_allele = allele.Reference_Sequence_comp_a
        allele_seq_aln_comp_a = allele.Aligned_Sequence_comp_a

        comp_a_map = get_refpos_values(comp_a_seq_aln_allele, allele_seq_aln_comp_a)

        if recoding_mode:
            comp_a_match_arr, comp_a_ins_match_arr, comp_a_del_match_arr = get_allele_match_array(comp_bp_changes_arr, comp_a_map, reference_info["composite_a_del_start"], reference_info["composite_a_del_end"], reference_info["composite_a_ins_start"], reference_info["composite_a_ins_end"])
        else:
            comp_a_match_arr, comp_a_ins_match_arr, comp_a_del_match_arr = get_allele_match_array(comp_a_bp_changes_arr, comp_a_map, reference_info["composite_a_del_start"], reference_info["composite_a_del_end"], reference_info["composite_a_ins_start"], reference_info["composite_a_ins_end"])

        comp_a_all_sub_pos, comp_a_sub_between_cuts, comp_a_all_del_pos, comp_a_del_between_cuts, comp_a_all_ins_pos, comp_a_has_substitutions, comp_a_has_deletions, comp_a_has_insertions = get_mutations(comp_a_map, comp_a_ref_seq, reference_info["cut_points_composite_a"])

        comp_a_has_indel = check_indel_positions(comp_a_all_ins_pos, comp_a_all_del_pos, reference_info["composite_a_del_start"], reference_info["composite_a_del_end"], reference_info["composite_a_ins_start"], reference_info["composite_a_ins_end"], ignore_extraspacer_deletions, reference_info["pegRNA_intervals_composite_a"])

        comp_a_total_TPE_count = comp_a_match_arr.count("T") + comp_a_match_arr.count("TI")
        comp_a_has_all_TPE = (comp_a_total_TPE_count == len(comp_a_match_arr))
        comp_a_has_any_TPE = (comp_a_total_TPE_count >= num_changes_to_check)
        comp_a_has_any_TPE_in_insertion = (comp_a_ins_match_arr.count("T") + comp_a_ins_match_arr.count("TI") >= num_changes_to_check)

        comp_a_total_WT_count = comp_a_match_arr.count("W") + comp_a_match_arr.count("WI")
        comp_a_has_all_WT = (comp_a_total_WT_count == len(comp_a_match_arr))
        comp_a_has_any_WT = (comp_a_total_WT_count > 0)

        comp_a_has_left_flap = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[:num_changes_to_check])
        comp_a_has_right_flap = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[-num_changes_to_check:])

        if comp_a_has_all_TPE and comp_a_has_indel:
            df_merged.at[idx,'Category_comp_a'] = "TPE_Indel"
        elif comp_a_has_all_TPE:
            df_merged.at[idx,'Category_comp_a'] = "Perfect_TPE"
        elif comp_a_has_all_WT and comp_a_has_indel:
            df_merged.at[idx,'Category_comp_a'] = "WT_Indel"
        elif comp_a_has_all_WT:
            df_merged.at[idx,'Category_comp_a'] = "WT"
        elif comp_a_has_left_flap and not comp_a_has_right_flap:
            df_merged.at[idx,'Category_comp_a'] = "Left_Flap"
        elif comp_a_has_right_flap and not comp_a_has_left_flap:
            df_merged.at[idx,'Category_comp_a'] = "Right_Flap"
        elif comp_a_has_any_WT and not comp_a_has_any_TPE_in_insertion:
            df_merged.at[idx,'Category_comp_a'] = "Imperfect_WT"
        elif comp_a_has_any_TPE:
            df_merged.at[idx,'Category_comp_a'] = "Imperfect_TPE"
        elif not comp_a_has_left_flap and not comp_a_has_right_flap and not comp_a_has_any_TPE:
            df_merged.at[idx,'Category_comp_a'] = "Imperfect_WT"
        else:
            df_merged.at[idx,'Category_comp_a'] = "Uncategorized"

        # Classify all alleles using Composite B reference alignment
        comp_b_seq_aln_allele = allele.Reference_Sequence_comp_b
        allele_seq_aln_comp_b = allele.Aligned_Sequence_comp_b

        comp_b_map = get_refpos_values(comp_b_seq_aln_allele, allele_seq_aln_comp_b)

        if recoding_mode:
            comp_b_match_arr, comp_b_ins_match_arr, comp_b_del_match_arr = get_allele_match_array(comp_bp_changes_arr, comp_b_map, reference_info["composite_b_del_start"], reference_info["composite_b_del_end"], reference_info["composite_b_ins_start"], reference_info["composite_b_ins_end"])
        else:
            comp_b_match_arr, comp_b_ins_match_arr, comp_b_del_match_arr = get_allele_match_array(comp_b_bp_changes_arr, comp_b_map, reference_info["composite_b_del_start"], reference_info["composite_b_del_end"], reference_info["composite_b_ins_start"], reference_info["composite_b_ins_end"])

        comp_b_all_sub_pos, comp_b_sub_between_cuts, comp_b_all_del_pos, comp_b_del_between_cuts, comp_b_all_ins_pos, comp_b_has_substitutions, comp_b_has_deletions, comp_b_has_insertions = get_mutations(comp_b_map, comp_b_ref_seq, reference_info["cut_points_composite_b"])

        comp_b_has_indel = check_indel_positions(comp_b_all_ins_pos, comp_b_all_del_pos, reference_info["composite_b_del_start"], reference_info["composite_b_del_end"], reference_info["composite_b_ins_start"], reference_info["composite_b_ins_end"], ignore_extraspacer_deletions, reference_info["pegRNA_intervals_composite_b"])

        comp_b_total_TPE_count = comp_b_match_arr.count("T") + comp_b_match_arr.count("TI")
        comp_b_has_all_TPE = (comp_b_total_TPE_count == len(comp_b_match_arr))
        comp_b_has_any_TPE = (comp_b_total_TPE_count >= num_changes_to_check)
        comp_b_has_any_TPE_in_insertion = (comp_b_ins_match_arr.count("T") + comp_b_ins_match_arr.count("TI") >= num_changes_to_check)

        comp_b_total_WT_count = comp_b_match_arr.count("W") + comp_b_match_arr.count("WI")
        comp_b_has_all_WT = (comp_b_total_WT_count == len(comp_b_match_arr))
        comp_b_has_any_WT = (comp_b_total_WT_count > 0)

        comp_b_has_left_flap = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[:num_changes_to_check])
        comp_b_has_right_flap = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[-num_changes_to_check:])

        if comp_b_has_all_TPE and comp_b_has_indel:
            df_merged.at[idx,'Category_comp_b'] = "TPE_Indel"
        elif comp_b_has_all_TPE:
            df_merged.at[idx,'Category_comp_b'] = "Perfect_TPE"
        elif comp_b_has_all_WT and comp_b_has_indel:
            df_merged.at[idx,'Category_comp_b'] = "WT_Indel"
        elif comp_b_has_all_WT:
            df_merged.at[idx,'Category_comp_b'] = "WT"
        elif comp_b_has_left_flap and not comp_b_has_right_flap:
            df_merged.at[idx,'Category_comp_b'] = "Left_Flap"
        elif comp_b_has_right_flap and not comp_b_has_left_flap:
            df_merged.at[idx,'Category_comp_b'] = "Right_Flap"
        elif comp_b_has_any_WT and not comp_b_has_any_TPE_in_insertion:
            df_merged.at[idx,'Category_comp_b'] = "Imperfect_WT"
        elif comp_b_has_any_TPE:
            df_merged.at[idx,'Category_comp_b'] = "Imperfect_TPE"
        elif not comp_b_has_left_flap and not comp_b_has_right_flap and not comp_b_has_any_TPE:
            df_merged.at[idx,'Category_comp_b'] = "Imperfect_WT"
        else:
            df_merged.at[idx,'Category_comp_b'] = "Uncategorized"

        # Write classification results to df_merged
        df_merged.at[idx,'tpe_ins_match_arr'] = tpe_ins_match_arr
        df_merged.at[idx,'comp_a_ins_match_arr'] = comp_a_ins_match_arr
        df_merged.at[idx,'comp_b_ins_match_arr'] = comp_b_ins_match_arr
        df_merged.at[idx,'wt_del_match_arr'] = wt_del_match_arr
        df_merged.at[idx,'comp_a_del_match_arr'] = comp_a_del_match_arr
        df_merged.at[idx,'comp_b_del_match_arr'] = comp_b_del_match_arr

    

    # Resolve category conflicts
    df_merged["Category_final"] = ""
    df_merged["Classified_by"] = ""

    for allele in df_merged.itertuples():

        # WT/TPE classifications always take precedence
        if allele.Category_tpe:
            df_merged.at[allele.Index, "Category_final"] = allele.Category_tpe
            df_merged.at[allele.Index, "Classified_by"] = "TPE"

        elif allele.Category_wt:
            df_merged.at[allele.Index, "Category_final"] = allele.Category_wt
            df_merged.at[allele.Index, "Classified_by"] = "WT"
        else:
            category, source = resolve_composite_categories(
                allele.Category_comp_a,
                allele.Aligned_Reference_Score_comp_a,
                allele.Category_comp_b,
                allele.Aligned_Reference_Score_comp_b,
            )

            df_merged.at[allele.Index, "Category_final"] = category
            df_merged.at[allele.Index, "Classified_by"] = source

    # print("Allele distribution:")
    # print(df_merged["Category_final"].value_counts())

    # print("df_merged columns:")
    # print(df_merged.columns)

    # df_merged.to_csv("/uufs/chpc.utah.edu/common/home/u0493285/clement/projects/20240216_rowley_hdr/analysis/01_20240718_jules_twinedit/03_nate_practicum/14_4ref_alignments/df_merged.csv", index=False)

    return df_merged, bp_changes_arrs


#### Plotting Functions ####
def get_plotting_stats(df=None, reference_info=None, bp_changes_arrs=None,twinspector_results_folder=None, recoding_mode=False):

    if recoding_mode:
        inserted_seq = "".join(tpe for _, _, tpe in bp_changes_arrs["std_bp_changes_arr"])
        deleted_seq = wt_seq = "".join(wt for _, wt, _ in bp_changes_arrs["std_bp_changes_arr"])
        ins_region_len = len(inserted_seq)
        del_region_len = len(deleted_seq)
    else:
        inserted_seq = reference_info["inserted_seq"]
        deleted_seq = reference_info["deleted_seq"]
        ins_region_len = reference_info["ins_region_len"]
        del_region_len = reference_info["del_region_len"]

    # Stats for summary barplots
    wt_count = 0
    wt_indel_count = 0
    imperfect_wt_count = 0
    right_flap_count = 0
    left_flap_count = 0
    imperfect_tpe_count = 0
    tpe_indels_count = 0
    perfect_tpe_count = 0
    uncategorized_count = 0
    # Stats for per-base barplots
    total_reads = 0
    edit_counts_arr = [0] * ins_region_len
    from_right_all_edit_counts_arr = [0] * ins_region_len
    from_left_all_edit_counts_arr = [0] * ins_region_len
    perfect_edit_counts = 0
    edit_counts_del_region_arr = [0] * del_region_len
    from_right_all_edit_counts_del_region_arr = [0] * del_region_len
    from_left_all_edit_counts_del_region_arr = [0] * del_region_len

    for idx, allele in df.iterrows():

        # Stats for summary barplots
        if allele.Category_final == "WT":
            wt_count += allele["#Reads"]
        elif allele.Category_final == "WT_Indel":
            wt_indel_count += allele["#Reads"]
        elif allele.Category_final == "Imperfect_WT":
            imperfect_wt_count += allele["#Reads"]
        elif allele.Category_final == "Right_Flap":
            right_flap_count += allele["#Reads"]
        elif allele.Category_final == "Left_Flap":
            left_flap_count += allele["#Reads"]
        elif allele.Category_final == "Imperfect_TPE":
            imperfect_tpe_count += allele["#Reads"]
        elif allele.Category_final == "TPE_Indel":
            tpe_indels_count += allele["#Reads"]
        elif allele.Category_final == "Perfect_TPE":
            perfect_tpe_count += allele["#Reads"]
        else:
            uncategorized_count += allele["#Reads"]

        # Stats for per-base barplots
        del_match_arr = None
        ins_match_arr = None
        if allele["Classified_by"] == "WT":
            del_match_arr = allele["wt_del_match_arr"]
        if allele["Classified_by"] == "TPE":
            ins_match_arr = allele["tpe_ins_match_arr"]
        if allele["Classified_by"] == "Composite_A":
            ins_match_arr = allele["comp_a_ins_match_arr"]
            del_match_arr = allele["comp_a_del_match_arr"]
        if allele["Classified_by"] == "Composite_B":
            ins_match_arr = allele["comp_b_ins_match_arr"]
            del_match_arr = allele["comp_b_del_match_arr"]

        if ins_match_arr:
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match in {"T", "TI"}:
                    edit_counts_arr[pos_idx] += allele['#Reads']

            for pos_idx, match in zip(reversed(range(len(ins_match_arr))), reversed(ins_match_arr)):
                if match in {"T", "TI"}:
                    from_right_all_edit_counts_arr[pos_idx] += allele['#Reads']
                else:
                    break

            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match in {"T", "TI"}:
                    from_left_all_edit_counts_arr[pos_idx] += allele['#Reads']
                else:
                    break

            if all(base in {"T", "TI"} for base in ins_match_arr):
                perfect_edit_counts += allele['#Reads']

        if del_match_arr:    
            for pos_idx, match in zip(range(len(del_match_arr)), del_match_arr):
                if match in {"T", "TI"}:
                    edit_counts_del_region_arr[pos_idx] += allele['#Reads']

            for pos_idx, match in zip(range(len(del_match_arr)), del_match_arr):
                if match in {"T", "TI"}:
                    from_left_all_edit_counts_del_region_arr[pos_idx] += allele['#Reads']
                else:
                    break

            for pos_idx, match in zip(reversed(range(len(del_match_arr))), reversed(del_match_arr)):
                if match in {"T", "TI"}:
                    from_right_all_edit_counts_del_region_arr[pos_idx] += allele['#Reads']
                else:
                    break

        total_reads += allele["#Reads"]

    total_reads_arr = [total_reads] * ins_region_len
    perfect_edit_counts_arr = [perfect_edit_counts] * ins_region_len
    total_reads_del_region_arr = [total_reads] * del_region_len
    perfect_edit_counts_del_region_arr = [perfect_edit_counts] * del_region_len
    edit_counts_del_region_arr = [x + perfect_edit_counts for x in edit_counts_del_region_arr]
    from_left_all_edit_counts_del_region_arr = [x + perfect_edit_counts for x in from_left_all_edit_counts_del_region_arr]
    from_right_all_edit_counts_del_region_arr = [x + perfect_edit_counts for x in from_right_all_edit_counts_del_region_arr]

    # Stats for summary barplots
    category_counts = {
        "Perfect TPE": perfect_tpe_count,
        "TPE Indel": tpe_indels_count,
        "Imperfect TPE": imperfect_tpe_count,
        "Left Flap": left_flap_count,
        "Right Flap": right_flap_count,
        "Imperfect WT": imperfect_wt_count,
        "WT Indel": wt_indel_count,
        "WT": wt_count,
        "Uncategorized": uncategorized_count,
    }

    # All stats
    stats = {
        "category_counts": category_counts,
        "total_reads_arr": total_reads_arr, 
        "edit_counts_arr": edit_counts_arr, 
        "from_right_all_edit_counts_arr": from_right_all_edit_counts_arr,
        "from_left_all_edit_counts_arr": from_left_all_edit_counts_arr,
        "perfect_edit_counts_arr": perfect_edit_counts_arr, 
        "inserted_seq": inserted_seq, 
        "total_reads_del_region_arr": total_reads_del_region_arr,
        "perfect_edit_counts_del_region_arr": perfect_edit_counts_del_region_arr,
        "edit_counts_del_region_arr": edit_counts_del_region_arr,
        "from_right_all_edit_counts_del_region_arr": from_right_all_edit_counts_del_region_arr,
        "from_left_all_edit_counts_del_region_arr": from_left_all_edit_counts_del_region_arr, 
        "deleted_seq": deleted_seq
    }

    # Write to files
    with open(twinspector_results_folder + "/c1.category_counts.txt", "w") as fout:
        fout.write("\t".join(["Perfect TPE", "TPE Indel", "Imperfect TPE", "Left Flap", "Right Flap", "Imperfect WT", "WT Indel", "WT", "Uncategorized"]) + "\n")
        fout.write("\t".join([str(category_counts[cat]) for cat in ["Perfect TPE", "TPE Indel", "Imperfect TPE", "Left Flap", "Right Flap", "Imperfect WT", "WT Indel", "WT", "Uncategorized"]]) + "\n")

    with open(twinspector_results_folder + "/c2.top_alleles_by_category.txt", "w") as fout:
        for c in CATEGORY_ORDER + ["Uncategorized"]:
            fc = c.replace(" ", "_")
            if fc not in df['Category_final'].values:
                continue
            fout.write(f"Category: {c}, Total Reads: {category_counts[c]}\n")
            for idx, row in df[df['Category_final'] == fc].sort_values(by='#Reads', ascending=False).head(50).iterrows():
                fout.write(f"Read: {idx}  count: {row['#Reads']}  Classified by: {row['Classified_by']}  Alignment Scores: {row['Aligned_Reference_Score_wt']} (WT), {row['Aligned_Reference_Score_tpe']} (TPE), {row['Aligned_Reference_Score_comp_a']} (Composite A), {row['Aligned_Reference_Score_comp_b']} (Composite B)\n")
                if row['Classified_by'] == 'WT':
                    fout.write(f"{row['Aligned_Sequence_wt']}\n")
                    fout.write(f"{row['Reference_Sequence_wt']}\n")
                elif row['Classified_by'] == 'TPE':
                    fout.write(f"{row['Aligned_Sequence_tpe']}\n")
                    fout.write(f"{row['Reference_Sequence_tpe']}\n")
                elif row['Classified_by'] == 'Composite_A':
                    fout.write(f"{row['Aligned_Sequence_comp_a']}\n")
                    fout.write(f"{row['Reference_Sequence_comp_a']}\n")
                elif row['Classified_by'] == 'Composite_B':
                    fout.write(f"{row['Aligned_Sequence_comp_b']}\n")
                    fout.write(f"{row['Reference_Sequence_comp_b']}\n")
                elif row['Classified_by'] == 'Composite_A&B':
                    fout.write(f"{row['Aligned_Sequence_comp_a']}\n")
                    fout.write(f"{row['Reference_Sequence_comp_a']}\n")
                    fout.write(f"{row['Aligned_Sequence_comp_b']}\n")
                    fout.write(f"{row['Reference_Sequence_comp_b']}\n")
            fout.write("\n")

    with open(twinspector_results_folder + "/c3.base_counts.txt", "w") as fout:
        fout.write("inserted_tpe_bases\t" + "\t".join([str(x) for x in inserted_seq]) + "\n")
        fout.write("perfectly_edited_base_counts\t" + '\t'.join([str(x) for x in perfect_edit_counts_arr]) + "\n")
        fout.write("imperfectly_edited_base_counts\t" + '\t'.join([str(x) for x in edit_counts_arr]) + "\n")
        fout.write("continuously_edited_base_counts_from_left\t" + '\t'.join([str(x) for x in from_left_all_edit_counts_arr]) + "\n")
        fout.write("continuously_edited_base_counts_from_right\t" + '\t'.join([str(x) for x in from_right_all_edit_counts_arr]) + "\n")
        fout.write("total_read_counts\t" + '\t'.join([str(x) for x in total_reads_arr]) + "\n\n")
        fout.write("removed_wt_bases\t" + '\t'.join([str(x) for x in deleted_seq]) + "\n")
        fout.write("perfectly_edited_base_counts\t" + '\t'.join([str(x) for x in perfect_edit_counts_del_region_arr]) + "\n")
        fout.write("imperfectly_edited_base_counts\t" + '\t'.join([str(x) for x in edit_counts_del_region_arr]) + "\n")
        fout.write("continuously_edited_base_counts_from_left\t" + '\t'.join([str(x) for x in from_left_all_edit_counts_del_region_arr]) + "\n")
        fout.write("continuously_edited_base_counts_from_right\t" + '\t'.join([str(x) for x in from_right_all_edit_counts_del_region_arr]) + "\n")
        fout.write("total_read_counts\t" + '\t'.join([str(x) for x in total_reads_del_region_arr]) + "\n")

    return stats


def setBarMatplotlibDefaults():
    matplotlib.rcParams["font.sans-serif"] = [
        "Arial",
        "Liberation Sans",
        "Bitstream Vera Sans",
    ]
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["axes.facecolor"] = "white"
    plt.ioff()


def setAlleleMatplotlibDefaults():
    font = {"size": 22}
    matplotlib.rc("font", **font)
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    sns.set(style="white", font_scale=2.2)


#### Summary barplots ####
def plot_reads_input_summary_barplot(crispresso_output_folder, counts_dict, fig_root=None, produce_pdf=False):

    crispresso_mapping_statistics_file = os.path.join(crispresso_output_folder, 'CRISPResso_mapping_statistics.txt')
    read_stats = pd.read_csv(crispresso_mapping_statistics_file, sep="\t")
    # Update Reads Aligned count for post-CRISPResso homology filtering
    counts = [read_stats['READS IN INPUTS'][0], read_stats['READS AFTER PREPROCESSING'][0], sum(counts_dict.values())]  # read_stats['READS ALIGNED'][0]]
    labels = ["Input", "After Preprocessing", "Analyzed"]
    total = read_stats['READS IN INPUTS'][0]

    sorted_pairs = sorted(zip(labels, counts), key=lambda x: x[1], reverse=True)
    sorted_labels = [p[0] for p in sorted_pairs]
    sorted_values = [p[1] for p in sorted_pairs]

    percent_labels = [f"{lab}\n({val/total*100:.1f}%)" for lab, val in sorted_pairs]

    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)

    bars = ax.bar(sorted_labels, sorted_values, 
                color="lightgrey", edgecolor="white", linewidth=0.6)

    for bar in bars:
        height = bar.get_height()
            
        ax.text(bar.get_x() + bar.get_width()/2.,
                height,
                f"{height:,}", 
                ha='center',
                va='bottom',
                fontsize=10,
                color='black'
            )
        
    ax.set_ylabel("Read Counts", fontsize=10)

    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(percent_labels, rotation=0, ha="center")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ax.minorticks_on()
    # ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    # ax.tick_params(axis="x", which="minor", bottom=False)

    plt.tight_layout()

    plt.savefig(fig_root + "/a1.Reads_input_summary.png", bbox_inches='tight', dpi=300)
    if produce_pdf:
        plt.savefig(fig_root + "/a1.Reads_input_summary.pdf", bbox_inches='tight')


def plot_category_stacked_summary_barplot(crispresso_output_folder, counts_dict, fig_root=None, produce_pdf=False, category_colors=None):

    aligned_labels = list(counts_dict.keys())
    aligned_counts = list(counts_dict.values())
    crispresso_mapping_statistics_file = os.path.join(crispresso_output_folder, 'CRISPResso_mapping_statistics.txt')
    read_stats = pd.read_csv(crispresso_mapping_statistics_file, sep="\t")
    unaligned_count = read_stats['READS AFTER PREPROCESSING'][0] - sum(counts_dict.values())  # read_stats['READS ALIGNED'][0]

    # fixed ordering
    custom_order = [
        "WT",
        "WT Indel",
        "Imperfect WT",
        "Right Flap",
        "Left Flap",
        "Imperfect TPE",
        "TPE Indel",
        "Perfect TPE"
    ]
    count_dict = dict(zip(aligned_labels, aligned_counts))
    plot_order = custom_order[::-1]
    sorted_counts_labels = [
        (label, count_dict.get(label, 0))
        for label in plot_order
    ]

    # Order by counts
    # sorted_counts_labels = sorted(zip(aligned_labels, aligned_counts), key=lambda x: x[1], reverse=True)

    total = sum(aligned_counts)
    legend_labels = [f"{lab} ({val/total*100:.1f}%)" for lab, val in sorted_counts_labels]

    fig, axes = plt.subplots(1, 2, figsize=(3, 6), dpi=300, sharey=True, gridspec_kw={'wspace': 0}) # dpi=300

    x = [0]
    bottom = 0
    for (lab, val), legend_label in zip(sorted_counts_labels, legend_labels):
        color = category_colors.get(lab, "black")
        axes[0].bar(
            x, val, 
            bottom=bottom, 
            label=legend_label, 
            color=color, 
            edgecolor='white',
            linewidth=.3
        )
        bottom += val

    axes[0].text(x[0], total, f"{total:,}", ha='center', va='bottom', fontsize=8)

    bar2 = axes[1].bar([0], unaligned_count, color='black')

    axes[1].text([0][0], unaligned_count, f"{unaligned_count:,}", ha='center', va='bottom', fontsize=8)

    handles, labels = axes[0].get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    handles.append(bar2[0])
    axes[0].legend(
        handles, labels,
        bbox_to_anchor=(2.03, 0.5),
        loc="center left",
        borderaxespad=0.25,
        fontsize=8
    )

    axes[0].set_xticks([0], labels=["Analyzed"], fontsize=8)
    # axes[0].set_xlabel("Analyzed", fontsize=8)
    axes[0].set_ylabel("Read Counts", fontsize=8)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].minorticks_on()
    axes[0].tick_params(axis="y", which="minor", length=3, width=0.8)
    axes[0].tick_params(axis="y", which="major", length=5, width=1)
    axes[0].tick_params(axis="y", labelsize=8)

    axes[1].set_xticks([0], labels=["Discarded"], fontsize=8)
    # axes[1].set_xlabel("Discarded", fontsize=8)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[1].minorticks_on()
    axes[1].tick_params(axis="y", which="minor", length=3, width=0.8, left=False)
    axes[1].tick_params(axis="y", which="major", length=5, width=1, left=False)

    # plt.tight_layout()
    plt.savefig(fig_root + "/a2.Category_summary_stacked.png", bbox_inches='tight', dpi=300)
    if produce_pdf:
        plt.savefig(fig_root + "/a2.Category_summary_stacked.pdf", bbox_inches='tight')


def plot_category_summary_barplot(counts_dict, fig_root=None, produce_pdf=False, category_colors=None):
    labels = list(counts_dict.keys())
    counts = list(counts_dict.values())

    # Fixed ordering
    custom_order = [
        "WT",
        "WT Indel",
        "Imperfect WT",
        "Right Flap",
        "Left Flap",
        "Imperfect TPE",
        "TPE Indel",
        "Perfect TPE"
    ]
    count_dict = dict(zip(labels, counts))

    sorted_pairs = [
        (label, count_dict[label])
        for label in reversed(custom_order)
        if label in count_dict
    ]

    sorted_labels = [p[0] for p in sorted_pairs]
    sorted_values = [p[1] for p in sorted_pairs]

    # Ordering by counts    
    # sorted_pairs = sorted(zip(labels, counts), key=lambda x: x[1], reverse=True)
    # sorted_labels = [p[0] for p in sorted_pairs]
    # sorted_values = [p[1] for p in sorted_pairs]

    total = sum(sorted_values)

    percent_labels = [f"{lab}\n({val/total*100:.1f}%)" for lab, val in sorted_pairs]

    colors = [category_colors.get(lab, "black") for lab in sorted_labels]

    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)

    bars = ax.bar(sorted_labels, sorted_values, 
                color=colors, edgecolor="white", linewidth=0.6)
    
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.,
                height,
                f"{height:,}",
                ha='center',
                va='bottom',
                fontsize=8,
            )

    ax.set_ylabel("Read Counts", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)

    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(percent_labels, rotation=0, ha="center", fontsize=7)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ax.minorticks_on()
    # ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    ax.tick_params(axis="x", which="minor", bottom=False)

    plt.tight_layout()

    plt.savefig(fig_root + "/a3.Category_summary.png", bbox_inches='tight', dpi=300)
    if produce_pdf:
        plt.savefig(fig_root + "/a3.Category_summary.pdf", bbox_inches='tight')


def plot_summary_barplots(category_counts, crispresso_output_folder, twinspector_results_folder, produce_pdf):

    setBarMatplotlibDefaults()
    
    plot_reads_input_summary_barplot(
        crispresso_output_folder,
        category_counts,
        fig_root=twinspector_results_folder,
        produce_pdf=produce_pdf
    )

    plot_category_stacked_summary_barplot(
        crispresso_output_folder,
        category_counts, 
        fig_root=twinspector_results_folder, 
        produce_pdf=produce_pdf, 
        category_colors=CATEGORY_COLORS
    )

    plot_category_summary_barplot(
        category_counts,
        fig_root=twinspector_results_folder, 
        produce_pdf=produce_pdf, 
        category_colors=CATEGORY_COLORS
    )


#### Per-base barplots ####
def plot_flap_incorporation(
    total_counts,
    edit_counts,
    from_right_all_edit_counts,
    from_left_all_edit_counts,
    perfect_edit_counts, 
    insert_sequence,
    show_total_reads=True, 
    recoding_mode=False, 
    fig_root=None,
    produce_pdf=False, 
    category_colors=None
):
    n = len(insert_sequence)
    indices = np.arange(n)

    physical_bar_width = 0.25 
    physical_gap_width = 0.02 

    width_per_base = physical_bar_width + physical_gap_width
    bar_width = physical_bar_width / width_per_base

    spine_gap = bar_width / 2 + 0.15
    data_range = (n - 1 + spine_gap) - (-spine_gap)
    axes_width = data_range * width_per_base

    min_fig_width = 16
    fig_height = 7
    y_label_space = 1.0  
    right_margin = 0.5   

    natural_fig_width = axes_width + y_label_space + right_margin
    fig_width = max(min_fig_width, natural_fig_width)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

    extra_padding = fig_width - natural_fig_width
    axes_left = y_label_space + (extra_padding / 2)

    ax = fig.add_axes([axes_left / fig_width, 0.20, axes_width / fig_width, 0.72])

    total_reads = max(total_counts)

    total_counts_pct = np.array(total_counts) / total_reads * 100
    edit_counts_pct = np.array(edit_counts) / total_reads * 100
    from_right_all_edit_counts_pct = np.array(from_right_all_edit_counts) / total_reads * 100
    from_left_all_edit_counts_pct = np.array(from_left_all_edit_counts) / total_reads * 100
    perfect_edit_counts_pct = np.array(perfect_edit_counts) / total_reads * 100

    if show_total_reads:
        ax.bar(indices, total_counts_pct, width=bar_width, label="Total Reads", color=category_colors["WT"], alpha=1.0)
        plot_max = 100
    else:
        plot_max = min(100, max(edit_counts_pct) * 1.05)

    ax.bar(indices, edit_counts_pct, width=bar_width, label="Total TPEs", color=category_colors["Imperfect TPE"], alpha=1.0)
    ax.bar(indices, from_left_all_edit_counts_pct, width=bar_width, color=category_colors["Left Flap"], label="Continuous 3'-Flap A Integration", alpha=1.0)
    ax.bar(indices, from_right_all_edit_counts_pct, width=bar_width, color=category_colors["Right Flap"], label="Continuous 3'-Flap B Integration", alpha=0.75)
    ax.bar(indices, perfect_edit_counts_pct, width=bar_width, label="Perfect TPE Reads", color=category_colors["Perfect TPE"], alpha=1.0)

    ax.set_ylabel("% of analyzed reads", fontsize=16)
    ax.set_xlim(-spine_gap, n - 1 + spine_gap)
    ax.set_ylim(0, plot_max)

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.tick_params(axis='x', length=0)

    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', length=3, width=0.8)
    ax.tick_params(axis='y', which='major', labelsize=16, length=6, width=1)

    max_height = plot_max
    gap_inches = 0.08
    rect_height_inches = 0.25
    fig_height_in = fig.get_size_inches()[1]
    ax_pos = ax.get_position()
    ax_height_in = fig_height_in * ax_pos.height

    gap_data = gap_inches / ax_height_in * max_height
    rect_height = rect_height_inches / ax_height_in * max_height
    y_base = -(gap_data + rect_height)

    rect_y_offset = 0.0085 * max_height
    text_y_offset = 0.0055 * max_height
    region_text_offset = 0.014 * max_height

    for i, base in enumerate(insert_sequence):
        rect = patches.Rectangle(
            (i - bar_width/2, y_base + rect_y_offset),
            bar_width,
            rect_height,
            facecolor=BASE_COLORS.get(base, "#ffffff"),
            edgecolor="none",
            clip_on=False
        )
        ax.add_patch(rect)
        ax.text(
            i,
            y_base + rect_height/2 + text_y_offset,
            base,
            ha="center",
            va="center",
            fontsize=12,
            clip_on=False
        )

    region_gap_inches = 0.05
    region_height_inches = 0.12

    region_gap_data = region_gap_inches / ax_height_in * max_height
    region_height = region_height_inches / ax_height_in * max_height
    y_region = y_base - region_gap_data - region_height

    label_text = "Edited Bases" if recoding_mode else "Edited Sequence"
    ax.text(
        (n - 1)/2,
        y_region + region_text_offset,
        label_text,
        ha="center",
        va="top",
        fontsize=16,
        color="black",
        clip_on=False
    )

    legend_fontsize = 16
    points_per_inch = 72
    perfect_handle_length = (physical_bar_width * points_per_inch) / legend_fontsize

    fig.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 0.1),
        ncol=5 if n > 33 else 3,
        frameon=False,
        fontsize=legend_fontsize,
        handlelength=perfect_handle_length,
    )

    file_name = "a4.3'_flap_integration_all_reads" if show_total_reads else "a5.3'_flap_integration_edited_reads"
    
    plt.savefig(f"{fig_root}/{file_name}.png", dpi=300)
    if produce_pdf:
        plt.savefig(f"{fig_root}/{file_name}.pdf")


def plot_flap_removal(
    total_counts,
    edit_counts,
    from_right_all_edit_counts_del_region,
    from_left_all_edit_counts_del_region,
    perfect_edit_counts_del_region, 
    deleted_sequence,
    show_total_reads=True, 
    recoding_mode=False, 
    fig_root=None,
    produce_pdf=False, 
    category_colors=None
):
    n = len(deleted_sequence)
    indices = np.arange(n)

    physical_bar_width = 0.25 
    physical_gap_width = 0.02 

    width_per_base = physical_bar_width + physical_gap_width
    bar_width = physical_bar_width / width_per_base

    spine_gap = bar_width / 2 + 0.15
    data_range = (n - 1 + spine_gap) - (-spine_gap)
    axes_width = data_range * width_per_base

    min_fig_width = 16
    fig_height = 7
    y_label_space = 1.0  
    right_margin = 0.5   

    natural_fig_width = axes_width + y_label_space + right_margin
    fig_width = max(min_fig_width, natural_fig_width)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

    extra_padding = fig_width - natural_fig_width
    axes_left = y_label_space + (extra_padding / 2)

    ax = fig.add_axes([axes_left / fig_width, 0.20, axes_width / fig_width, 0.72])

    total_reads = max(total_counts)

    total_counts_pct = np.array(total_counts) / total_reads * 100
    edit_counts_pct = np.array(edit_counts) / total_reads * 100
    from_right_all_edit_counts_pct = np.array(from_right_all_edit_counts_del_region) / total_reads * 100
    from_left_all_edit_counts_pct = np.array(from_left_all_edit_counts_del_region) / total_reads * 100
    perfect_edit_counts_pct = np.array(perfect_edit_counts_del_region) / total_reads * 100

    if show_total_reads:
        ax.bar(indices, total_counts_pct, width=bar_width, label="Total Reads", color=category_colors["WT"], alpha=1.0)
        plot_max = 100
    else:
        plot_max = max(edit_counts_pct) * 1.05
        
    ax.bar(indices, edit_counts_pct, width=bar_width, label="Total TPEs", color=category_colors["Imperfect TPE"], alpha=1.0)
    ax.bar(indices, from_left_all_edit_counts_pct, width=bar_width, color=category_colors["Left Flap"], label="Continuous 5'-Flap A Removal", alpha=1.0)
    ax.bar(indices, from_right_all_edit_counts_pct, width=bar_width, color=category_colors["Right Flap"], label="Continuous 5'-Flap B Removal", alpha=0.75)
    ax.bar(indices, perfect_edit_counts_pct, width=bar_width, label="Perfect TPE Reads", color=category_colors["Perfect TPE"], alpha=1.0)

    ax.set_ylabel("% of analyzed reads", fontsize=16)
    ax.set_xlim(-spine_gap, n - 1 + spine_gap)
    ax.set_ylim(0, plot_max)

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.tick_params(axis='x', length=0)

    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', length=3, width=0.8)
    ax.tick_params(axis='y', which='major', labelsize=16, length=6, width=1)

    max_height = plot_max
    gap_inches = 0.08
    rect_height_inches = 0.25
    fig_height_in = fig.get_size_inches()[1]
    ax_pos = ax.get_position()
    ax_height_in = fig_height_in * ax_pos.height

    gap_data = gap_inches / ax_height_in * max_height
    rect_height = rect_height_inches / ax_height_in * max_height
    y_base = -(gap_data + rect_height)

    rect_y_offset = 0.0085 * max_height
    text_y_offset = 0.0055 * max_height
    region_text_offset = 0.014 * max_height

    for i, base in enumerate(deleted_sequence):
        rect = patches.Rectangle(
            (i - bar_width/2, y_base + rect_y_offset),
            bar_width,
            rect_height,
            facecolor=BASE_COLORS.get(base, "#ffffff"),
            edgecolor="none",
            clip_on=False
        )
        ax.add_patch(rect)
        ax.text(
            i,
            y_base + rect_height/2 + text_y_offset,
            base,
            ha="center",
            va="center",
            fontsize=12,
            clip_on=False
        )

    region_gap_inches = 0.05
    region_height_inches = 0.12

    region_gap_data = region_gap_inches / ax_height_in * max_height
    region_height = region_height_inches / ax_height_in * max_height
    y_region = y_base - region_gap_data - region_height

    label_text = "Wild-type Bases" if recoding_mode else "Wild-type Sequence"
    ax.text(
        (n - 1)/2,
        y_region + region_text_offset,
        label_text,
        ha="center",
        va="top",
        fontsize=16,
        color="black",
        clip_on=False
    )

    legend_fontsize = 16
    points_per_inch = 72
    perfect_handle_length = (physical_bar_width * points_per_inch) / legend_fontsize

    fig.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 0.1),
        ncol=5 if n > 33 else 3,
        frameon=False,
        fontsize=legend_fontsize,
        handlelength=perfect_handle_length,
    )

    file_name = "a6.5'_flap_removal_all_reads" if show_total_reads else "a7.5'_flap_removal_edited_reads"
    
    plt.savefig(f"{fig_root}/{file_name}.png", dpi=300)
    if produce_pdf:
        plt.savefig(f"{fig_root}/{file_name}.pdf")


def plot_per_base_pos_barplots(plotting_info, twinspector_results_folder=None, produce_pdf=False, recoding_mode=False):

    setBarMatplotlibDefaults()

    plot_flap_incorporation(
        total_counts=plotting_info["total_reads_arr"],  
        edit_counts=plotting_info["edit_counts_arr"], 
        from_right_all_edit_counts=plotting_info["from_right_all_edit_counts_arr"], 
        from_left_all_edit_counts=plotting_info["from_left_all_edit_counts_arr"], 
        perfect_edit_counts=plotting_info["perfect_edit_counts_arr"], 
        insert_sequence=plotting_info["inserted_seq"], 
        show_total_reads=True,
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder,
        produce_pdf=produce_pdf, 
        category_colors=CATEGORY_COLORS 
    )

    plot_flap_incorporation(
        total_counts=plotting_info["total_reads_arr"], 
        edit_counts=plotting_info["edit_counts_arr"],
        from_right_all_edit_counts=plotting_info["from_right_all_edit_counts_arr"],
        from_left_all_edit_counts=plotting_info["from_left_all_edit_counts_arr"],
        perfect_edit_counts=plotting_info["perfect_edit_counts_arr"],
        insert_sequence=plotting_info["inserted_seq"], 
        show_total_reads=False, 
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder,
        produce_pdf=produce_pdf, 
        category_colors=CATEGORY_COLORS
    )

    plot_flap_removal(
        total_counts=plotting_info["total_reads_del_region_arr"], 
        edit_counts=plotting_info["edit_counts_del_region_arr"],
        from_right_all_edit_counts_del_region=plotting_info["from_right_all_edit_counts_del_region_arr"],
        from_left_all_edit_counts_del_region=plotting_info["from_left_all_edit_counts_del_region_arr"],
        perfect_edit_counts_del_region=plotting_info["perfect_edit_counts_del_region_arr"],
        deleted_sequence=plotting_info["deleted_seq"], 
        show_total_reads=True, 
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder,
        produce_pdf=produce_pdf, 
        category_colors=CATEGORY_COLORS
    )
    
    plot_flap_removal(
        total_counts=plotting_info["total_reads_del_region_arr"], 
        edit_counts=plotting_info["edit_counts_del_region_arr"],
        from_right_all_edit_counts_del_region=plotting_info["from_right_all_edit_counts_del_region_arr"],
        from_left_all_edit_counts_del_region=plotting_info["from_left_all_edit_counts_del_region_arr"],
        perfect_edit_counts_del_region=plotting_info["perfect_edit_counts_del_region_arr"],
        deleted_sequence=plotting_info["deleted_seq"], 
        show_total_reads=False, 
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder,
        produce_pdf=produce_pdf, 
        category_colors=CATEGORY_COLORS
    )


#### Allele tables ####
def get_dataframe_allele_region(
    df_alleles,
    pegRNA_intervals,
    ref_seq,
    ref_aln_seq,
    twin_aln_seq,
    cut_points,
    window_by_intervals=False,
    left_pad=20,
    right_pad=20
):
    """
    Return aligned sequences trimmed so all arrays match.
    If window_by_intervals is True, first trim to a window spanning pegRNA intervals
    with left_pad bases before the first interval start and right_pad bases after
    the last interval end; then harmonize to the shortest length within that window.
    """
    if df_alleles.shape[0] == 0:
        return df_alleles, ref_seq, ref_aln_seq, twin_aln_seq, cut_points, pegRNA_intervals

    ordered_intervals = sorted(list(pegRNA_intervals or []), key=lambda x: x[0])

    # Compute initial window bounds (on aligned string indices)
    start_idx = 0
    stop_idx_excl = len(ref_aln_seq)
    if window_by_intervals and len(ordered_intervals) >= 1:
        first_start = ordered_intervals[0][0]
        last_end_incl = ordered_intervals[-1][1]  # inclusive end
        # Apply padding and clamp
        left_bound = max(0, int(first_start) - int(left_pad))
        right_bound_incl = min(len(ref_aln_seq) - 1, int(last_end_incl) + int(right_pad))
        start_idx = left_bound
        stop_idx_excl = right_bound_incl + 1

    # Slice all aligned strings to the window first
    ref_aln_seq_win  = ref_aln_seq[start_idx:stop_idx_excl] if isinstance(ref_aln_seq, str) else ref_aln_seq
    twin_aln_seq_win = twin_aln_seq[start_idx:stop_idx_excl] if isinstance(twin_aln_seq, str) else twin_aln_seq

    # Aggregate counts but DO NOT slice by arbitrary lengths yet
    df = df_alleles.copy()
    # Renames columns to match
    df = df.rename(
        columns={
            col: "Reference_Sequence" if "Reference_Sequence" in col 
            else "Aligned_Sequence" if "Aligned_Sequence" in col 
            else col 
            for col in df.columns
        }
    )
    # df = (
    #     df.groupby(
    #         ['Aligned_Sequence', 'Reference_Sequence'],  # , 'Read_Status', 'n_deleted', 'n_inserted', 'n_mutated'],
    #         dropna=False
    #     ).sum().reset_index()
    # )

    # Slice allele strings to the window (if present)
    if 'Aligned_Sequence' in df.columns:
        df['Aligned_Sequence'] = df['Aligned_Sequence'].astype(str).str.slice(start_idx, stop_idx_excl)
    if 'Reference_Sequence' in df.columns:
        df['Reference_Sequence'] = df['Reference_Sequence'].astype(str).str.slice(start_idx, stop_idx_excl)

    # Determine smallest length across aligned strings in the window
    min_len_candidates = []
    if isinstance(ref_aln_seq_win, str):
        min_len_candidates.append(len(ref_aln_seq_win))
    if isinstance(twin_aln_seq_win, str):
        min_len_candidates.append(len(twin_aln_seq_win))
    for col in ('Aligned_Sequence', 'Reference_Sequence'):
        if col in df.columns:
            col_series = df[col].dropna()
            if len(col_series) > 0:
                lens = col_series[col_series.map(lambda x: isinstance(x, str))].map(len)
                if len(lens) > 0:
                    min_len_candidates.append(int(lens.min()))
    min_len = max(0, min(min_len_candidates)) if min_len_candidates else len(ref_aln_seq_win)

    # Final trim to shortest length
    ref_aln_seq_trim  = ref_aln_seq_win[:min_len] if isinstance(ref_aln_seq_win, str) else ref_aln_seq_win
    twin_aln_seq_trim = twin_aln_seq_win[:min_len] if isinstance(twin_aln_seq_win, str) else twin_aln_seq_win
    if 'Aligned_Sequence' in df.columns:
        df['Aligned_Sequence'] = df['Aligned_Sequence'].astype(str).str.slice(0, min_len)
    if 'Reference_Sequence' in df.columns:
        df['Reference_Sequence'] = df['Reference_Sequence'].astype(str).str.slice(0, min_len)

    # Adjust cut points from global to window coords, then clamp to min_len
    cut_points_window = None
    if cut_points is not None:
        try:
            cps = list(cut_points)
        except TypeError:
            cps = [cut_points]
        # shift into window coordinates
        cps = [cp - start_idx for cp in cps]
        cut_points_window = [max(0, min(int(cp), max(0, min_len - 1))) for cp in cps]

    # Adjust intervals from global to window coords, then clamp to [0, min_len]
    pegRNA_intervals_window = []
    for (s, e) in ordered_intervals:
        s_win = int(s) - start_idx
        e_win = int(e) - start_idx
        s_clamped = max(0, min(s_win, max(0, min_len - 1)))
        e_clamped = max(s_clamped, min(e_win, max(0, min_len)))
        pegRNA_intervals_window.append((s_clamped, e_clamped))

    df = df.set_index('Aligned_Sequence')
    df.sort_values(
        by=["#Reads", "sequence_key"],
        ascending=[False, True],
        inplace=True,
    )
    # df["Unedited"] = df["Read_Status"].eq("UNMODIFIED")

    return (
        df,
        ref_seq,                 # keep raw ref_seq unchanged
        ref_aln_seq_trim,        # windowed+trimmed aligned reference
        twin_aln_seq_trim,       # windowed+trimmed aligned twin
        cut_points_window if cut_points_window is not None else cut_points,
        pegRNA_intervals_window
    )


def prep_alleles_table(
    df_alleles,
    reference_seq,
    ref_aln_seq_region,
    tpe_aln_seq_region,
    MAX_N_ROWS,
    MIN_FREQUENCY,
    pegRNA_intervals, 
):
    """
    Prepare matrices and metadata required to render an allele heatmap.

    This function converts a subset of rows from a CRISPResso allele table into:
      - X: numeric matrix encoding per-position allele bases for each displayed allele
      - annot: parallel matrix of string characters for on-cell annotation
      - y_labels: formatted labels per allele row (percent and read count)
      - insertion_dict: mapping of row index to a list of (start,end) x-spans that denote
        insertion events (identified by runs of '-' in the reference alignment)
      - per_element_annot_kws: per-cell Text style dictionaries (e.g., bold substitutions)
      - is_reference: boolean flags indicating whether the allele is identical to the
        provided reference sequence without indels

    Selection is limited to the top MAX_N_ROWS rows that meet MIN_FREQUENCY based on
    df_alleles['%Reads'].

    Args:
        df_alleles (pd.DataFrame): Allele table indexed by Aligned_Sequence (or similar),
            containing at least columns: 'Reference_Sequence', '#Reads', '%Reads'.
        reference_seq (str): Ungapped reference sequence used to determine reference rows.
        MAX_N_ROWS (int): Maximum number of allele rows to include.
        MIN_FREQUENCY (float): Minimum percentage (df_alleles['%Reads']) to include a row.

    Returns:
        tuple:
            X (list[list[int]]): Numeric-encoded alleles using dna_to_numbers mapping.
            annot (list[list[str]]): Same shape as X with literal characters per cell.
            y_labels (list[str]): Display labels per row (e.g., "12.34% (123 reads)").
            insertion_dict (defaultdict(list)): row_index -> list of (start, end) insertion spans.
            per_element_annot_kws (np.ndarray): Per-cell dicts of Text properties used for styling.
            is_reference (list[bool]): True if row exactly matches reference_seq and has no indels.
    """
    dna_to_numbers = {"-": 0, "A": 1, "T": 2, "C": 3, "G": 4, "N": 5, " ": 6}
    seq_to_numbers = lambda seq: [dna_to_numbers[x] for x in seq]
    X = []
    annot = []
    y_labels = []
    insertion_dict = defaultdict(list)
    per_element_annot_kws = []
    is_reference = []
    num_blanks = 2

    # Regex to find contiguous gap runs ('-') in the reference alignment; these correspond to
    # insertions in the read relative to the reference.
    re_find_indels = re.compile("(-*-)")
    idx_row = 0
    for idx, row in df_alleles[df_alleles["%Reads"] >= MIN_FREQUENCY][
        :MAX_N_ROWS
    ].iterrows():
        # Encode the allele (index) sequence
        X.append(seq_to_numbers(idx.upper()))
        annot.append(list(idx))

        # Track insertion spans based on gaps in the reference sequence
        has_indels = False
        for p in re_find_indels.finditer(row["Reference_Sequence"]):
            has_indels = True
            insertion_dict[idx_row].append((p.start(), p.end()))

        # Build y-axis labels with percentage and read count
        y_labels.append("%.2f%% (%d reads)" % (row["%Reads"], row["#Reads"]))
        if idx == reference_seq and not has_indels:
            is_reference.append(True)
        else:
            is_reference.append(False)

        idx_row += 1

        # Detect substitutions (non-gap mismatches) to style them in bold/black
        idxs_sub = [
            i_sub
            for i_sub in range(len(idx))
            if (row["Reference_Sequence"][i_sub] != idx[i_sub])
            and (row["Reference_Sequence"][i_sub] != "-")
            and (idx[i_sub] != "-")
            # If only the bases that differ from both displayed refs (not just the aligned ref) should be bolded then
            # check WT and TPE Seqs as well (for standard ref alingments, for composite refs also need to check alternate composite - tpe_aln_seq_region_a/b and wt_aln_seq_region_a/b)
            # and (idx[i_sub] != ref_aln_seq_region[i_sub])
            # and (idx[i_sub] != tpe_aln_seq_region[i_sub])
        ]
        to_append = np.array([{}] * len(idx), dtype=object)
        to_append[idxs_sub] = {"weight": "bold", "color": "black", "size": 16}
        per_element_annot_kws.append(to_append)

    if tpe_aln_seq_region is not None and len(tpe_aln_seq_region) > 0:
        # Append two blank rows to the end of outputs
        for _ in range(num_blanks):
            blank_row = " " * len(tpe_aln_seq_region)
            X.append(seq_to_numbers(blank_row))
            annot.append(list(blank_row))
            y_labels.append("")
            is_reference.append(False)
            to_append = np.array([{"color": "white"}] * len(blank_row), dtype=object)
            per_element_annot_kws.append(to_append)

        # Append tpe_aln_seq_region to the end of outputs
        X.append(seq_to_numbers(tpe_aln_seq_region.upper()))
        annot.append(list(tpe_aln_seq_region))
        y_labels.append("TwinPE Reference")
        # if len(pegRNA_intervals) > 1:
        #     y_labels.append("TwinPE Reference")
        # else:
        #     y_labels.append("PE Reference")
        is_reference.append(False)
        # # detect substitutions (non-gap mismatches) between tpe_aln_seq_region and ref_aln_seq_region to style them in bold/black
        # if recoding_mode:
        #     idxs_sub = [
        #         i_sub
        #         for i_sub in range(len(tpe_aln_seq_region))
        #         if (ref_aln_seq_region[i_sub] != tpe_aln_seq_region[i_sub])
        #         and (ref_aln_seq_region[i_sub] != "-")
        #         and (tpe_aln_seq_region[i_sub] != "-")
        #     ]
        to_append = np.array([{}] * len(tpe_aln_seq_region), dtype=object)
        # if recoding_mode:
        #     to_append[idxs_sub] = {"weight": "bold", "color": "black", "size": 16}
        per_element_annot_kws.append(to_append)
        # # track insertions in tpe_aln_seq_region relative to ref_aln_seq_region
        # for p in re_find_indels.finditer(ref_aln_seq_region):
        #     insertion_dict[idx_row + num_blanks].append((p.start(), p.end()))

    return X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference


def get_nuc_color(nuc, alpha):
    """
    Return a consistent RGBA color tuple for a nucleotide or special token.

    Args:
        nuc (str): One of {'A','T','C','G','N','INS','DEL','-'} or any other string.
            'N' denotes ambiguous; '-' denotes gap. 'INS'/'DEL' share the same color to
            visually group indels. Any unknown token results in a deterministic color
            derived from its character codes.
        alpha (float): Alpha channel in [0.0, 1.0] controlling transparency.

    Returns:
        tuple: (r, g, b, a) with values in [0.0, 1.0].
    """
    get_color = lambda x, y, z: (x / 255.0, y / 255.0, z / 255.0, alpha)
    if nuc == "A":
        return get_color(127, 201, 127)
    elif nuc == "T":
        return get_color(190, 174, 212)
    elif nuc == "C":
        return get_color(253, 192, 134)
    elif nuc == "G":
        return get_color(255, 255, 153)
    elif nuc == "N":
        return get_color(200, 200, 200)
    elif nuc == "INS":
        #        return get_color(185,219,253)
        #        return get_color(177,125,76)
        return get_color(193, 129, 114)
    elif nuc == "DEL":
        # return get_color(177,125,76)
        #        return get_color(202,109,87)
        return get_color(193, 129, 114)
    elif nuc == "-":
        # return get_color(177,125,76)
        #        return get_color(202,109,87)
        return get_color(30, 30, 30)
    elif nuc == " ":
        # white space for padding
        return get_color(255, 255, 255)
    else:  # return a random color (that is based on the nucleotide given)
        charSum = 0
        for char in nuc.upper():
            thisval = ord(char) - 65  #'A' is 65
            if thisval < 0 or thisval > 90:
                thisval = 0
            charSum += thisval
        charSum = (charSum / len(nuc)) / 90.0

        return (charSum, (1 - charSum), (2 * charSum * (1 - charSum)))


def get_rows_for_sgRNA_annotation(sgRNA_intervals, amp_len):
    """
    Assign a vertical "row" for each sgRNA interval so that overlapping intervals
    are staggered and do not visually collide when drawn.

    The algorithm greedily places each interval on the top-most row that does not
    already contain any of its covered x positions. Occupancy is tracked per integer
    x position between the interval's (clipped) start and end.

    Args:
        sgRNA_intervals (list[tuple[int,int]]): List of (start, end) sgRNA spans in reference
            coordinates. Intervals are clipped to [0, amp_len-1] for overlap calculations.
        amp_len (int): Amplicon/reference length for clipping.

    Returns:
        np.ndarray: Row indices (int) per sgRNA, where 0 is the top-most row. The rows are
        inverted (highest row index becomes 0) so that earlier rows appear visually higher
        when drawn relative to negative y offsets.
    """
    # figure out how many rows are needed to show all sgRNAs
    sgRNA_plot_rows = [0] * len(
        sgRNA_intervals
    )  # which row each sgRNA should be plotted on
    sgRNA_plot_occupancy = []  # which idxs are already filled on each row
    sgRNA_plot_occupancy.append([])
    for idx, sgRNA_int in enumerate(sgRNA_intervals):
        this_sgRNA_start = max(0, sgRNA_int[0])
        this_sgRNA_end = min(sgRNA_int[1], amp_len - 1)
        curr_row = 0
        if this_sgRNA_start > amp_len or this_sgRNA_end < 0:
            # Interval entirely outside; place on row 0 and continue
            sgRNA_plot_rows[idx] = curr_row
            continue
        # Bump the row until there is no position overlap with already-placed intervals
        while (
            len(
                np.intersect1d(
                    sgRNA_plot_occupancy[curr_row],
                    range(this_sgRNA_start, this_sgRNA_end),
                )
            )
            > 0
        ):
            next_row = curr_row + 1
            if not next_row in sgRNA_plot_occupancy:
                sgRNA_plot_occupancy.append([])
            curr_row = next_row
        sgRNA_plot_rows[idx] = curr_row
        # Mark occupancy for the chosen row
        sgRNA_plot_occupancy[curr_row].extend(range(this_sgRNA_start, this_sgRNA_end))
    # Invert rows so that the last-created (lowest) row is drawn lowest when using negative offsets
    return np.subtract(max(sgRNA_plot_rows), sgRNA_plot_rows)


class Custom_HeatMapper(sns.matrix._HeatMapper):
    """
    Extension of seaborn's private _HeatMapper to support per-element annotation style
    (per-element text properties) and to suppress the colorbar.

    This utility mirrors seaborn.heatmap internals while allowing an "annot" matrix to be
    styled cell-by-cell via a parallel matrix of dictionaries (per_element_annot_kws), where
    each dictionary can specify matplotlib.text.Text properties for the corresponding cell.

    Caution: sns.matrix._HeatMapper is a private API and may change across seaborn versions.
    """

    def __init__(
        self,
        data,
        vmin,
        vmax,
        cmap,
        center,
        robust,
        annot,
        fmt,
        annot_kws,
        per_element_annot_kws,
        cbar,
        cbar_kws,
        xticklabels=True,
        yticklabels=True,
        mask=None,
    ):
        """
        Initialize the heatmap plotter and capture optional per-element annotation styles.

        Args mirror seaborn.heatmap/_HeatMapper with the following addition:
            per_element_annot_kws (np.ndarray | list | None): Same shape as `annot` where each
                element is a dict of Text properties applied to that cell's annotation.
                If None, an empty dict is used for every cell.
        """
        super(Custom_HeatMapper, self).__init__(
            data,
            vmin,
            vmax,
            cmap,
            center,
            robust,
            annot,
            fmt,
            annot_kws,
            cbar,
            cbar_kws,
            xticklabels,
            yticklabels,
            mask,
        )

        # Prepare a mirror structure for per-element annotation keyword arguments
        if annot is not None:
            if per_element_annot_kws is None:
                self.per_element_annot_kws = np.empty_like(annot, dtype=object)
                self.per_element_annot_kws[:] = dict()
            else:
                self.per_element_annot_kws = per_element_annot_kws

    # add per element dict to style the annotation
    def _annotate_heatmap(self, ax, mesh):
        """Add textual labels with the value in each cell.

        This override allows passing a per-cell dictionary of Text properties to fine-tune
        the appearance (e.g., bold substitutions) while preserving seaborn's luminance-based
        foreground color choice.
        """
        mesh.update_scalarmappable()
        xpos, ypos = np.meshgrid(ax.get_xticks(), ax.get_yticks())

        # Iterate the mesh values, facecolors, annotations, and per-cell styles in lock-step
        for x, y, m, color, val, per_element_dict in zip(
            xpos.flat,
            ypos.flat,
            mesh.get_array().flat,
            mesh.get_facecolors(),
            self.annot_data.flat,
            self.per_element_annot_kws.flat,
        ):
            # print per_element_dict
            if m is not np.ma.masked:
                l = sns.utils.relative_luminance(color)
                text_color = ".15" if l > 0.408 else "w"
                annotation = ("{:" + self.fmt + "}").format(str(val))
                text_kwargs = dict(color=text_color, ha="center", va="center")
                text_kwargs.update(self.annot_kws)
                text_kwargs.update(per_element_dict)

                ax.text(x, y, annotation, **text_kwargs)

    # removed the colorbar
    def plot(self, ax, cax, kws):
        """Draw the heatmap body and tick labels on the provided Axes.

        This version deliberately avoids attaching a colorbar and leaves any colorbar
        management to the caller.
        """
        # Remove all the Axes spines for a cleaner matrix look
        sns.utils.despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap as a pcolormesh for performance on large matrices
        # If a Normalize is supplied, do NOT also pass vmin/vmax.
        if "norm" in kws and kws["norm"] is not None:
            mesh = ax.pcolormesh(self.plot_data, cmap=self.cmap, **kws)
        else:
            mesh = ax.pcolormesh(
                self.plot_data, vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, **kws
            )

        # Set axis limits to span the matrix exactly
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Add row and column labels
        ax.set(xticks=self.xticks, yticks=self.yticks)
        xtl = ax.set_xticklabels(self.xticklabels)
        ytl = ax.set_yticklabels(self.yticklabels, rotation="vertical", va="center")

        # Possibly rotate them if they overlap after layout
        plt.draw()
        if sns.utils.axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if sns.utils.axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, mesh)


def custom_heatmap(
    data,
    vmin=None,
    vmax=None,
    cmap=None,
    center=None,
    robust=False,
    annot=None,
    fmt=".2g",
    annot_kws=None,
    per_element_annot_kws=None,
    linewidths=0,
    linecolor="white",
    cbar=True,
    cbar_kws=None,
    cbar_ax=None,
    square=False,
    ax=None,
    xticklabels=True,
    yticklabels=True,
    mask=None,
    **kwargs,
):
    """
    Convenience wrapper around Custom_HeatMapper to draw a heatmap matrix with optional
    per-element annotation styling and without a colorbar by default.

    Args:
        data (np.ndarray | list): 2D array of numeric values to visualize.
        vmin, vmax (float | None): Colormap scaling bounds.
        cmap (matplotlib.colors.Colormap | str | None): Colormap to use.
        center (float | None): If set, shift the colormap center to this value.
        robust (bool): If True, use robust quantiles rather than min/max for colormap scaling.
        annot (np.ndarray | list | None): 2D array of values/strings to annotate each cell.
        fmt (str): Format string applied to annotations.
        annot_kws (dict | None): Global matplotlib.text.Text properties applied to all annotations.
        per_element_annot_kws (np.ndarray | list | None): Same shape as annot; each element is a
            dict of Text properties applied to that cell, allowing per-cell styles.
        linewidths (float): Line width between cells (pcolormesh edge widths).
        linecolor (str): Line color between cells.
        cbar (bool): Present for API parity; colorbar is not created by this function.
        cbar_kws (dict | None): Ignored here; reserved for compatibility.
        cbar_ax (matplotlib.axes.Axes | None): Ignored here; reserved for compatibility.
        square (bool): If True, set aspect to equal so each cell is square.
        ax (matplotlib.axes.Axes | None): Axes to draw into; if None, uses current axes.
        xticklabels, yticklabels: Tick label configuration as in seaborn.heatmap.
        mask (np.ndarray | None): Boolean mask specifying cells not to plot.
        **kwargs: Additional arguments passed to Axes.pcolormesh (e.g., shading, antialiased).

    Returns:
        matplotlib.axes.Axes: The Axes containing the heatmap.
    """

    # Initialize the plotter object
    plotter = Custom_HeatMapper(
        data,
        vmin,
        vmax,
        cmap,
        center,
        robust,
        annot,
        fmt,
        annot_kws,
        per_element_annot_kws,
        cbar,
        cbar_kws,
        xticklabels,
        yticklabels,
        mask,
    )

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax


def add_sgRNA_to_ax(ax,
                    sgRNA_intervals,
                    sgRNA_y_start,
                    sgRNA_y_height,
                    amp_len,
                    x_offset=0,
                    sgRNA_mismatches=None,
                    sgRNA_names=None,
                    sgRNA_rows=None,
                    font_size=None,
                    clip_on=True,
                    label_at_zero=False,
                    sgRNA_label_sides=None,
                    ref_row_seq=None,                 # WT ref (with '-')
                    extend_left_non_gap=None,           # per-sgRNA left extension in non-gap bases
                    extend_right_non_gap=None):       # per-sgRNA right extension in non-gap bases
    """
    Draw one or more sgRNA annotations on a Matplotlib Axes.

    Each sgRNA is a semi-transparent rectangle (start..end). Optional mismatches
    are red sub-rectangles. Labels can be placed per sgRNA either to the left
    or right of the rectangle via sgRNA_label_sides.

    Notes:
        - Right-side labels are anchored just beyond the original rectangle end (end+1).
        - If extend_left_non_gap is set, the visual rectangle start is shifted
          left by that many non-gap bases (skipping '-' in ref_row_seq).
        - If ref_row_seq is provided, the rectangle is split into segments and drawn
          only over non-gap columns (skipping any '-' within the span).
        - Mismatch blocks remain aligned to the original (unextended) sgRNA start.
    """
    if font_size is None:
        font_size = matplotlib.rcParams['font.size']

    if sgRNA_rows is None:
        sgRNA_rows = [0]*len(sgRNA_intervals)
    max_sgRNA_row = max(sgRNA_rows)+1
    this_sgRNA_y_height = sgRNA_y_height/float(max_sgRNA_row)

    # Normalize label sides
    if sgRNA_label_sides is None:
        sgRNA_label_sides = ['left'] * len(sgRNA_intervals)
    else:
        if len(sgRNA_label_sides) < len(sgRNA_intervals):
            sgRNA_label_sides = sgRNA_label_sides + ['left']*(len(sgRNA_intervals)-len(sgRNA_label_sides))
        else:
            sgRNA_label_sides = sgRNA_label_sides[:len(sgRNA_intervals)]
        sgRNA_label_sides = [('right' if s.lower()=='right' else 'left') for s in sgRNA_label_sides]

    # Normalize extension list
    if extend_left_non_gap is None:
        extend_left_non_gap = [0]*len(sgRNA_intervals)
    else:
        if len(extend_left_non_gap) < len(sgRNA_intervals):
            extend_left_non_gap = extend_left_non_gap + [0]*(len(sgRNA_intervals)-len(extend_left_non_gap))
        else:
            extend_left_non_gap = extend_left_non_gap[:len(sgRNA_intervals)]

    if extend_right_non_gap is None:
        extend_right_non_gap = [0]*len(sgRNA_intervals)
    else:
        if len(extend_right_non_gap) < len(sgRNA_intervals):
            extend_right_non_gap = extend_right_non_gap + [0]*(len(sgRNA_intervals)-len(extend_right_non_gap))
        else:
            extend_right_non_gap = extend_right_non_gap[:len(sgRNA_intervals)]

    def _left_shift_by_non_gaps(row_seq, start_idx, n_non_gaps):
        """Return how many columns to shift left to include n_non_gaps non-'-' bases."""
        if row_seq is None or n_non_gaps <= 0 or start_idx <= 0:
            return 0
        count = 0
        steps = 0
        i = int(start_idx) - 1
        while i >= 0 and count < n_non_gaps:
            if row_seq[i] != '-':
                count += 1
            steps += 1
            i -= 1
        return steps
    
    def _right_shift_by_non_gaps(row_seq, end_idx, n_non_gaps):
        if row_seq is None or n_non_gaps <= 0 or end_idx >= len(row_seq)-1:
            return 0
        count = 0
        steps = 0
        i = int(end_idx) + 1
        while i < len(row_seq) and count < n_non_gaps:
            if row_seq[i] != '-':
                count += 1
            steps += 1
            i += 1
        return steps

    def _non_gap_runs(row_seq, start_idx, end_idx):
        """Return list of (run_start, run_end) contiguous non-gap segments in [start_idx, end_idx]."""
        if row_seq is None:
            return [(int(start_idx), int(end_idx))]
        if len(row_seq) == 0:
            return []
        s = max(0, int(start_idx))
        e = min(int(end_idx), len(row_seq) - 1)
        runs = []
        i = s
        while i <= e:
            while i <= e and row_seq[i] == '-':
                i += 1
            if i > e:
                break
            run_start = i
            while i <= e and row_seq[i] != '-':
                i += 1
            run_end = i - 1
            runs.append((run_start, run_end))
        return runs

    min_sgRNA_x = None
    label_left_sgRNA = True

    for idx, sgRNA_int in enumerate(sgRNA_intervals):
        # Original clipped interval (for mismatch alignment and caps)
        this_sgRNA_start = max(0, sgRNA_int[0])
        this_sgRNA_end   = min(sgRNA_int[1], amp_len - 1)
        if this_sgRNA_start > amp_len or this_sgRNA_end < 0:
            continue

        this_sgRNA_y_row_start = sgRNA_y_start + this_sgRNA_y_height*sgRNA_rows[idx]

        # Visual start extended left by N non-gap bases (skip '-')
        left_extra_cols = _left_shift_by_non_gaps(ref_row_seq, this_sgRNA_start, extend_left_non_gap[idx])
        right_extra_cols = _right_shift_by_non_gaps(ref_row_seq, this_sgRNA_end, extend_right_non_gap[idx])
        display_start = max(0, this_sgRNA_start - left_extra_cols)
        display_end = min(amp_len-1, this_sgRNA_end + right_extra_cols)

        # Draw as multiple rectangles over non-gap runs only
        # runs = _non_gap_runs(ref_row_seq, display_start, this_sgRNA_end)
        runs = _non_gap_runs(ref_row_seq, display_start, display_end)
        for seg_start, seg_end in runs:
            if seg_start > seg_end:
                continue
            ax.add_patch(
                patches.Rectangle(
                    (x_offset + seg_start, this_sgRNA_y_row_start),
                    1 + seg_end - seg_start,
                    this_sgRNA_y_height,
                    facecolor=(0, 0, 0, 0.15),
                    clip_on=clip_on
                )
            )

        # Clip caps (based on original interval vs clipping)
        if this_sgRNA_start != sgRNA_int[0]:
            ax.add_patch(
                patches.Rectangle(
                    (x_offset + 0.1 + this_sgRNA_start, this_sgRNA_y_row_start),
                    0.1,
                    this_sgRNA_y_height,
                    facecolor='w',
                    clip_on=clip_on
                )
            )
        if this_sgRNA_end != sgRNA_int[1]:
            ax.add_patch(
                patches.Rectangle(
                    (x_offset + 0.8 + this_sgRNA_end, this_sgRNA_y_row_start),
                    0.1,
                    this_sgRNA_y_height,
                    facecolor='w',
                    clip_on=clip_on
                )
            )

        # Mismatches (relative to original sgRNA start)
        if sgRNA_mismatches is not None and idx < len(sgRNA_mismatches):
            for mismatch in sgRNA_mismatches[idx]:
                mismatch_plot_pos = sgRNA_int[0] + mismatch
                if 0 <= mismatch_plot_pos < amp_len:
                    ax.add_patch(
                        patches.Rectangle(
                            (x_offset + mismatch_plot_pos, this_sgRNA_y_row_start),
                            1,
                            this_sgRNA_y_height,
                            facecolor='r',
                            clip_on=clip_on
                        )
                    )

        # For left-anchored label and min_x heuristic, use first visible segment start
        leftmost_visible = runs[0][0] if runs else display_start
        if min_sgRNA_x is None or leftmost_visible < min_sgRNA_x:
            min_sgRNA_x = leftmost_visible

        # Label
        if sgRNA_names is not None and idx < len(sgRNA_names) and sgRNA_names[idx] != "":
            side = sgRNA_label_sides[idx]
            if side == 'left':
                anchor_x = x_offset + leftmost_visible
                if (label_at_zero and anchor_x < len(sgRNA_names[idx])*0.66):
                    ax.text(
                        0,
                        this_sgRNA_y_row_start + this_sgRNA_y_height/2,
                        sgRNA_names[idx] + " ",
                        horizontalalignment='left',
                        verticalalignment='center',
                        fontsize=font_size
                    )
                else:
                    ax.text(
                        anchor_x,
                        this_sgRNA_y_row_start + this_sgRNA_y_height/2,
                        sgRNA_names[idx] + " ",
                        horizontalalignment='right',
                        verticalalignment='center',
                        fontsize=font_size
                    )
            else:  # right
                label_x = x_offset + this_sgRNA_end + 1.0
                ax.text(
                    label_x,
                    this_sgRNA_y_row_start + this_sgRNA_y_height/2,
                    " " + sgRNA_names[idx],
                    horizontalalignment='left',
                    verticalalignment='center',
                    fontsize=font_size
                )
            label_left_sgRNA = False

    if min_sgRNA_x is not None and label_left_sgRNA:
        if (label_at_zero and x_offset + min_sgRNA_x < 5):
            ax.text(0,
                    this_sgRNA_y_row_start + this_sgRNA_y_height/2,
                    'sgRNA ',
                    horizontalalignment='left',
                    verticalalignment='center',
                    fontsize=font_size)
        else:
            ax.text(x_offset+min_sgRNA_x,
                    this_sgRNA_y_row_start + this_sgRNA_y_height/2,
                    'sgRNA ',
                    horizontalalignment='right',
                    verticalalignment='center',
                    fontsize=font_size)


def plot_alleles_heatmap(
        reference_aln_seq,
        X,
        annot,
        y_labels,
        insertion_dict,
        per_element_annot_kws,
        fig_filename_root=None,
        custom_colors=None,
        SAVE_ALSO_PDF=False,
        plot_cut_point=True,
        cut_point_ind=None,
        sgRNA_intervals=None,
        sgRNA_names=None,
        sgRNA_mismatches=None,
        extend_left_non_gap=None, 
        extend_right_non_gap=None, 
        category=None,
        **kwargs):
    """
    Plots alleles in a heatmap (nucleotides color-coded for easy visualization)
    input:
    -reference_seq: sequence of reference allele to plot
    -X: list of numbers representing nucleotides of the allele
    -annot: list of nucleotides (letters) of the allele
    -y_labels: list of labels for each row/allele
    -insertion_dict: locations of insertions -- red squares will be drawn around these
    -per_element_annot_kws: annotations for each cell (e.g. bold for substitutions, etc.)
    -fig_filename_root: figure filename to plot (not including '.pdf' or '.png'). If None, plots are shown interactively.
    -custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    -SAVE_ALSO_PDF: whether to write pdf file as well
    -plot_cut_point: if false, won't draw 'predicted cleavage' line
    -cut_point_ind: index of cut point (if None, will be plot in the middle calculated as len(reference_seq)/2)
    -sgRNA_intervals: locations where sgRNA is located
    -sgRNA_mismatches: array (for each sgRNA_interval) of locations in sgRNA where there are mismatches
    -sgRNA_names: array (for each sgRNA_interval) of names of sgRNAs (otherwise empty)
    -custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    """
    plot_nuc_len=len(reference_aln_seq)

    # make a color map of fixed colors
    alpha=0.4
    A_color=get_nuc_color('A', alpha)
    T_color=get_nuc_color('T', alpha)
    C_color=get_nuc_color('C', alpha)
    G_color=get_nuc_color('G', alpha)
    INDEL_color = get_nuc_color('N', alpha)
    blank_color = get_nuc_color(' ', alpha)

    if custom_colors is not None:
        hex_alpha = '66'  # this is equivalent to 40% in hexadecimal
        if 'A' in custom_colors:
            A_color = custom_colors['A'] + hex_alpha
        if 'T' in custom_colors:
            T_color = custom_colors['T'] + hex_alpha
        if 'C' in custom_colors:
            C_color = custom_colors['C'] + hex_alpha
        if 'G' in custom_colors:
            G_color = custom_colors['G'] + hex_alpha
        if 'N' in custom_colors:
            INDEL_color = custom_colors['N'] + hex_alpha
        if ' ' in custom_colors:
            blank_color = custom_colors[' '] + hex_alpha

    dna_to_numbers={'-':0,'A':1,'T':2,'C':3,'G':4,'N':5, ' ':6}
    seq_to_numbers= lambda seq: [dna_to_numbers[x] for x in seq]

    # cmap = colors_mpl.ListedColormap([INDEL_color, A_color, T_color, C_color, G_color, INDEL_color, blank_color])
    # cmap.set_over(blank_color)
    # norm = colors_mpl.Normalize(vmin=0, vmax=5)

    # New: 7-color colormap + discrete norm for 0..6
    cmap = colors_mpl.ListedColormap([INDEL_color, A_color, T_color, C_color, G_color, INDEL_color, blank_color])
    bnorm = colors_mpl.BoundaryNorm(np.arange(-0.5, 7.5, 1), cmap.N, clip=False)

    #ref_seq_around_cut=reference_seq[max(0,cut_point-plot_nuc_len/2+1):min(len(reference_seq),cut_point+plot_nuc_len/2+1)]

#    print('per element anoot kws: ' + per_element_annot_kws)
    if len(per_element_annot_kws) > 1:
        per_element_annot_kws=np.vstack(per_element_annot_kws[::-1])
    else:
        per_element_annot_kws=np.array(per_element_annot_kws)
    ref_seq_hm=np.expand_dims(seq_to_numbers(reference_aln_seq), 1).T
    ref_seq_annot_hm=np.expand_dims(list(reference_aln_seq), 1).T

    annot=annot[::-1]
    X=X[::-1]

    N_ROWS=len(X)
    N_COLUMNS=plot_nuc_len

    # if N_ROWS < 1:
    #     fig, ax = plt.subplots()
    #     fig.text(0.5, 0.5, 'No Alleles', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    #     ax.set_clip_on(False)

    #     # if fig_filename_root is None:
    #     #     plt.show()
    #     # else:
    #     fig.savefig(fig_filename_root+'.png', bbox_inches='tight', dpi=300)
    #     if SAVE_ALSO_PDF:
    #         fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight')
    #     plt.close(fig)
    #     return

    sgRNA_rows = []
    num_sgRNA_rows = 0

    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        sgRNA_rows = get_rows_for_sgRNA_annotation(sgRNA_intervals, plot_nuc_len)
        num_sgRNA_rows = max(sgRNA_rows) + 1
        fig=plt.figure(figsize=(plot_nuc_len*0.3, (N_ROWS+1 + num_sgRNA_rows)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+2, N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+2, N_COLUMNS)
        #ax_hm_ref heatmap for the reference
        ax_hm_ref=fig.add_subplot(gs1[0:1,:])
        ax_hm=fig.add_subplot(gs2[2:,:])
    else:
        fig=plt.figure(figsize=(plot_nuc_len*0.3, (N_ROWS+1)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+1, N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+1, N_COLUMNS)
        #ax_hm_ref heatmap for the reference
        ax_hm_ref=fig.add_subplot(gs1[0,:])
        ax_hm=fig.add_subplot(gs2[1:,:])


    custom_heatmap(ref_seq_hm, annot=ref_seq_annot_hm, annot_kws={'size':16}, cmap=cmap, fmt='s', ax=ax_hm_ref, norm=bnorm, vmin=None, vmax=None, square=True)
    custom_heatmap(X, annot=np.array(annot), annot_kws={'size':16}, cmap=cmap, fmt='s', ax=ax_hm, norm=bnorm, vmin=None, vmax=None, square=True, per_element_annot_kws=per_element_annot_kws)
    

    # place ticks at row centers
    ax_hm.yaxis.tick_right()
    ax_hm.set_yticks(np.arange(N_ROWS) + 0.5)
    # apply labels in reverse order
    ax_hm.set_yticklabels(y_labels[::-1], rotation=True, va='center')
    # turn off only the 2nd and 3rd tick marks (Python is 0-based)
    for idx in [1, 2]:
        tick = ax_hm.yaxis.get_major_ticks()[idx]
        tick.tick1line.set_visible(False)  # left tick
        tick.tick2line.set_visible(False)  # right tick
    ax_hm.xaxis.set_ticks([])

    # Newer method - includes right-side labels and extended grey pegRNA spacerB span by caller
    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        this_sgRNA_y_start = -1*num_sgRNA_rows
        this_sgRNA_y_height = num_sgRNA_rows - 0.3

        # Normalize left extension to a list the length of sgRNA_intervals
        if extend_left_non_gap is None:
            left_extend = [0]*len(sgRNA_intervals)
        elif isinstance(extend_left_non_gap, int):
            left_extend = [extend_left_non_gap]*len(sgRNA_intervals)
        elif isinstance(extend_left_non_gap, dict):
            left_extend = [0]*len(sgRNA_intervals)
            for i in range(len(sgRNA_intervals)):
                left_extend[i] = extend_left_non_gap.get(
                    sgRNA_names[i] if sgRNA_names and i < len(sgRNA_names) else i,
                    extend_left_non_gap.get(i, 0)
                )
        else:
            left_extend = list(extend_left_non_gap)
            if len(left_extend) < len(sgRNA_intervals):
                left_extend += [0]*(len(sgRNA_intervals)-len(left_extend))
            else:
                left_extend = left_extend[:len(sgRNA_intervals)]

        # Normalize right extension to a list the length of sgRNA_intervals
        if extend_right_non_gap is None:
            right_extend = [0]*len(sgRNA_intervals)
        elif isinstance(extend_right_non_gap, int):
            right_extend = [extend_right_non_gap]*len(sgRNA_intervals)
        elif isinstance(extend_right_non_gap, dict):
            right_extend = [0]*len(sgRNA_intervals)
            for i in range(len(sgRNA_intervals)):
                right_extend[i] = extend_right_non_gap.get(
                    sgRNA_names[i] if sgRNA_names and i < len(sgRNA_names) else i,
                    extend_right_non_gap.get(i, 0)
                )
        else:
            right_extend = list(extend_right_non_gap)
            if len(right_extend) < len(sgRNA_intervals):
                right_extend += [0]*(len(sgRNA_intervals)-len(right_extend))
            else:
                right_extend = right_extend[:len(sgRNA_intervals)]      

        add_sgRNA_to_ax(
            ax_hm_ref,
            sgRNA_intervals,
            sgRNA_y_start=this_sgRNA_y_start,
            sgRNA_y_height=this_sgRNA_y_height,
            amp_len=plot_nuc_len,
            font_size='small',
            clip_on=False,
            sgRNA_names=sgRNA_names,
            sgRNA_mismatches=sgRNA_mismatches,
            x_offset=0,
            label_at_zero=True,
            sgRNA_rows=sgRNA_rows,
            sgRNA_label_sides=['left','right'],
            # ref_row_seq=ref_aln_seq_region,
            ref_row_seq=reference_aln_seq,
            extend_left_non_gap=left_extend, 
            extend_right_non_gap=right_extend
        )
    # End of Newer method

    #create boxes for ins
    for idx, lss in insertion_dict.items():
        for ls in lss:
            ax_hm.add_patch(patches.Rectangle((ls[0], N_ROWS-idx-1), ls[1]-ls[0], 1, linewidth=3, edgecolor='r', fill=False))

    # Newer method to draw cut point lines at correct position and for WT and TwinPE refs rows
    # cut point vertical line
    if plot_cut_point:
        if cut_point_ind is None:
            cut_point_ind = [plot_nuc_len / 2]

        def _is_blank_row(chars):
            return all(c == " " for c in chars)

        N_ROWS = len(annot)

        n_bottom_blank = 0
        for i in range(1, min(3, N_ROWS)):
            if _is_blank_row(annot[i]):
                n_bottom_blank += 1
            else:
                break

        # convert base indices to boundary positions (draw line AFTER the specified base).
        # if cut_point_ind values are 0-based indices of the cut base, we add +1 so the
        # vertical line appears between base i and i+1.
        raw_points = cut_point_ind if isinstance(cut_point_ind, (list, tuple, np.ndarray)) else [cut_point_ind]
        xs = [cp + 1 for cp in raw_points if cp is not None]

        # on allele table (ax_hm): draw over TwinPE row
        if N_ROWS >= 1:
            ax_hm.vlines(xs, 0, 1, linestyles="dashed")

        # skip the blank spacer rows and draw through remaining allele rows
        y_start = 1 + n_bottom_blank
        if N_ROWS > y_start:
            ax_hm.vlines(xs, y_start, N_ROWS, linestyles="dashed")

        # on WT Reference (ax_hm_ref): single row
        ax_hm_ref.vlines(xs, 0, 1, linestyles="dashed")
    # End of New method
    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(['WT Reference'], rotation=True, va='center')

    gs2.update(left=0, right=1, hspace=0.05, wspace=0, top=1*(((N_ROWS)*1.13))/(N_ROWS))
    gs1.update(left=0, right=1, hspace=0.05, wspace=0,)

    sns.set_context(rc={'axes.facecolor':'white','lines.markeredgewidth': 1,'mathtext.fontset' : 'stix','text.usetex':True,'text.latex.unicode':True} )

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                    mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'), ms=18),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='r', marker='s', ms=8, markeredgewidth=2.5),
              matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='black', marker='_', ms=2,)]
    descriptions=['Substitutions', 'Insertions', 'Deletions']

    if plot_cut_point:
        proxies.append(
              matplotlib.lines.Line2D([0], [1], linestyle='--', c='black', ms=6))
        descriptions.append('Predicted cleavage position')

    if category:
        category = category.replace("_", r"\ ")
        proxies.append(patches.Patch(color='none'))
        descriptions.append(rf"$\bf{{{category}\ Alleles}}$")

    #ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)
    lgd = ax_hm.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='upper center', bbox_to_anchor=(0.5, 0), ncol=1, fancybox=True, shadow=False)

    # if fig_filename_root is None:
    #     plt.show()
    # else:
    fig.savefig(fig_filename_root+'.png', bbox_inches='tight', bbox_extra_artists=(lgd,), dpi=300)
    if SAVE_ALSO_PDF:
        fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
    plt.close(fig)


def plot_categorical_ref_allele_tables(
    df_alleles, 
    reference_info, 
    ref_seq,
    ref_aln_seq,
    tpe_aln_seq,
    min_frequency,
    max_n_rows, 
    plot_full_reads=False, 
    fig_root=None,
    produce_pdf=False, 
    ref_type=None
):
    if ref_type == "wt":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_wt", None)
        cut_points = reference_info.get("cut_points_wt", None)
        extend_left_non_gap = [0, 0]
        extend_right_non_gap = [0, 0]
    elif ref_type == "tpe":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_tpe", None)
        cut_points = reference_info.get("cut_points_tpe", None)
        extend_left_non_gap = [0, 0]
        extend_right_non_gap = [0, 0]
    elif ref_type == "comp_a":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_composite_a", None)
        cut_points = reference_info.get("cut_points_composite_a", None)
        extend_left_non_gap = [0, abs(reference_info['cleavage_offset_b'])]
        extend_right_non_gap = [0, 0]
    elif ref_type == "comp_b":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_composite_b", None)
        cut_points = reference_info.get("cut_points_composite_b", None)
        extend_left_non_gap = [0, 0]
        extend_right_non_gap = [abs(reference_info['cleavage_offset_a']), 0]

    for cat, name_comp_a, name_comp_b, name_tpe, name_wt in [
        ("Perfect_TPE", "b1.Perfect_TPE.aligned_to_composite_a", "b1.Perfect_TPE.aligned_to_composite_b", "b1.Perfect_TPE.aligned_to_tpe", "b1.Perfect_TPE.aligned_to_wt"), 
        ("TPE_Indel", "b2.TPE_Indel.aligned_to_composite_a", "b2.TPE_Indel.aligned_to_composite_b", "b2.TPE_Indel.aligned_to_tpe", "b2.TPE_Indel.aligned_to_wt"), 
        ("Imperfect_TPE", "b3.Imperfect_TPE.aligned_to_composite_a", "b3.Imperfect_TPE.aligned_to_composite_b", "b3.Imperfect_TPE.aligned_to_tpe", "b3.Imperfect_TPE.aligned_to_wt"),
        ("Left_Flap", "b4.Left_Flap.aligned_to_composite_a", "b4.Left_Flap.aligned_to_composite_b", "b4.Left_Flap.aligned_to_tpe", "b4.Left_Flap.aligned_to_wt"),
        ("Right_Flap", "b5.Right_Flap.aligned_to_composite_a", "b5.Right_Flap.aligned_to_composite_b", "b5.Right_Flap.aligned_to_tpe", "b5.Right_Flap.aligned_to_wt"),
        ("Imperfect_WT", "b6.Imperfect_WT.aligned_to_composite_a", "b6.Imperfect_WT.aligned_to_composite_b", "b6.Imperfect_WT.aligned_to_tpe", "b6.Imperfect_WT.aligned_to_wt"),
        ("WT_Indel", "b7.WT_Indel.aligned_to_composite_a", "b7.WT_Indel.aligned_to_composite_b", "b7.WT_Indel.aligned_to_tpe", "b7.WT_Indel.aligned_to_wt"),
        ("WT", "b8.WT.aligned_to_composite_a", "b8.WT.aligned_to_composite_b", "b8.WT.aligned_to_tpe", "b8.WT.aligned_to_wt"),
        ("Uncategorized", "b9.Uncategorized.aligned_to_composite_a", "b9.Uncategorized.aligned_to_composite_b", "b9.Uncategorized.aligned_to_tpe", "b9.Uncategorized.aligned_to_wt"),
    ]:
        
        df_alleles_cat = df_alleles[df_alleles["Category_final"] == cat]

        if df_alleles_cat.empty or not (df_alleles_cat["%Reads"] >= min_frequency).any():
            continue
        elif ref_type == "comp_a":
            name = name_comp_a
        elif ref_type == "comp_b":
            name = name_comp_b
        else:
            if ref_type == "wt":
                name = name_wt
            elif ref_type == "tpe":
                name = name_tpe

        df_alleles_around_region, ref_seq_region, ref_aln_seq_region, tpe_aln_seq_region, cut_points_window, pegRNA_intervals_region = get_dataframe_allele_region(
            df_alleles_cat,
            pegRNA_intervals,
            ref_seq,
            ref_aln_seq,
            tpe_aln_seq,
            cut_points, 
            window_by_intervals=not plot_full_reads, 
            # left_pad=6, 
            # right_pad=6
        )

        X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = prep_alleles_table(
            df_alleles_around_region,
            ref_seq_region,
            ref_aln_seq_region,
            tpe_aln_seq_region,
            max_n_rows,
            min_frequency,
            pegRNA_intervals
        )

        plot_alleles_heatmap(
            reference_aln_seq=ref_aln_seq_region,
            X=X,
            annot=annot,
            y_labels=y_labels,
            insertion_dict=insertion_dict,
            per_element_annot_kws=per_element_annot_kws,
            fig_filename_root=fig_root + f"/{name}",
            SAVE_ALSO_PDF=produce_pdf,
            plot_cut_point=[True, True],
            cut_point_ind=cut_points_window,
            sgRNA_intervals=pegRNA_intervals_region,
            sgRNA_names=["pegRNA a", "pegRNA b"],
            sgRNA_mismatches=[[], []], 
            category=cat,
            extend_left_non_gap=extend_left_non_gap, 
            extend_right_non_gap=extend_right_non_gap, 
        )


def plot_one_allele_table(task):
    df_alleles, reference_info, args, twinspector_results_folder, ref_type, ref_aln_seq, tpe_aln_seq = task

    setAlleleMatplotlibDefaults()

    plot_categorical_ref_allele_tables(
        df_alleles=df_alleles,
        reference_info=reference_info,
        ref_seq=args.wt_seq,
        ref_aln_seq=ref_aln_seq,
        tpe_aln_seq=tpe_aln_seq,
        min_frequency=args.min_frequency_alleles,
        max_n_rows=args.max_n_rows,
        plot_full_reads=args.plot_full_reads,
        fig_root=twinspector_results_folder,
        produce_pdf=args.produce_pdf,
        ref_type=ref_type 
    )
    return ref_type


#### Parallelized version of plot_ref_allele_tables using ThreadPoolExecutor ####
def plot_ref_allele_tables_parallel_tpe(args, df_categorized, reference_info, twinspector_results_folder):
    ref_specs = [
        ("comp_a", ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_comp_a', 'Reference_Sequence_comp_a', 'Category_final', 'Classified_by'], reference_info['wt_aln_seq_a'], reference_info['tpe_aln_seq_a']),
        ("comp_b", ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_comp_b', 'Reference_Sequence_comp_b', 'Category_final', 'Classified_by'], reference_info['wt_aln_seq_b'], reference_info['tpe_aln_seq_b']),
        ("tpe",    ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_tpe', 'Reference_Sequence_tpe', 'Category_final', 'Classified_by'], args.tpe_seq, args.tpe_seq),
        ("wt",     ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_wt', 'Reference_Sequence_wt', 'Category_final', 'Classified_by'], args.wt_seq, args.tpe_seq),
    ]

    tasks = [
        (df_categorized[cols], reference_info, args, twinspector_results_folder, ref_type, ref_aln_seq, tpe_aln_seq)
        for ref_type, cols, ref_aln_seq, tpe_aln_seq in ref_specs
    ]

    with ThreadPoolExecutor(max_workers=min(4, len(tasks))) as executor:
        futures = [executor.submit(plot_one_allele_table, task) for task in tasks]
        for future in as_completed(futures):
            future.result()


#### Parallelized version of plot_ref_allele_tables using ProcessPoolExecutor ####
def plot_ref_allele_tables_parallel_ppe(args, df_categorized, reference_info, twinspector_results_folder):
    ref_specs = [
        (
            "comp_a",
            ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_comp_a', 'Reference_Sequence_comp_a', 'Category_final', 'Classified_by'],
            reference_info['wt_aln_seq_a'],
            reference_info['tpe_aln_seq_a'],
        ),
        (
            "comp_b",
            ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_comp_b', 'Reference_Sequence_comp_b', 'Category_final', 'Classified_by'],
            reference_info['wt_aln_seq_b'],
            reference_info['tpe_aln_seq_b'],
        ),
        (
            "tpe",
            ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_tpe', 'Reference_Sequence_tpe', 'Category_final', 'Classified_by'],
            args.wt_seq,
            args.tpe_seq,
        ),
        (
            "wt",
            ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_wt', 'Reference_Sequence_wt', 'Category_final', 'Classified_by'],
            args.wt_seq,
            args.tpe_seq,
        ),
    ]

    tasks = [
        (df_categorized[cols], reference_info, args, twinspector_results_folder, ref_type, ref_aln_seq, tpe_aln_seq)
        for ref_type, cols, ref_aln_seq, tpe_aln_seq in ref_specs
    ]

    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(plot_one_allele_table, task) for task in tasks]
        for future in as_completed(futures):
            future.result()


#### Sequential version of plot_ref_allele_tables ####
def plot_ref_allele_tables_sequential(args, df_categorized, reference_info, twinspector_results_folder):

    wt_cols =     ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_wt','Reference_Sequence_wt', 'Category_final', 'Classified_by']
    tpe_cols =    ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_tpe', 'Reference_Sequence_tpe', 'Category_final', 'Classified_by']
    comp_a_cols = ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_comp_a','Reference_Sequence_comp_a', 'Category_final', 'Classified_by']
    comp_b_cols = ['sequence_key', '#Reads', '%Reads', 'Aligned_Sequence_comp_b', 'Reference_Sequence_comp_b','Category_final', 'Classified_by']

    setAlleleMatplotlibDefaults()

    plot_categorical_ref_allele_tables(
        df_alleles=df_categorized[comp_a_cols], 
        reference_info=reference_info, 
        ref_seq=args.wt_seq, 
        ref_aln_seq=reference_info['wt_aln_seq_a'], 
        tpe_aln_seq=reference_info['tpe_aln_seq_a'], 
        min_frequency=args.min_frequency_alleles,
        max_n_rows=args.max_n_rows,
        plot_full_reads=args.plot_full_reads,
        fig_root=twinspector_results_folder,
        produce_pdf=args.produce_pdf, 
        ref_type="comp_a" 
    )

    plot_categorical_ref_allele_tables(
        df_alleles=df_categorized[comp_b_cols], 
        reference_info=reference_info, 
        ref_seq=args.wt_seq, 
        ref_aln_seq=reference_info['wt_aln_seq_b'], 
        tpe_aln_seq=reference_info['tpe_aln_seq_b'], 
        min_frequency=args.min_frequency_alleles,
        max_n_rows=args.max_n_rows,
        plot_full_reads=args.plot_full_reads,
        fig_root=twinspector_results_folder,
        produce_pdf=args.produce_pdf, 
        ref_type="comp_b" 
    )

    plot_categorical_ref_allele_tables(
        df_alleles=df_categorized[tpe_cols], 
        reference_info=reference_info, 
        ref_seq=args.wt_seq, 
        ref_aln_seq=args.wt_seq, 
        tpe_aln_seq=args.tpe_seq, 
        min_frequency=args.min_frequency_alleles,
        max_n_rows=args.max_n_rows,
        plot_full_reads=args.plot_full_reads,
        fig_root=twinspector_results_folder,
        produce_pdf=args.produce_pdf, 
        ref_type="tpe" 
    )

    plot_categorical_ref_allele_tables(
        df_alleles=df_categorized[wt_cols], 
        reference_info=reference_info, 
        ref_seq=args.wt_seq, 
        ref_aln_seq=args.wt_seq, 
        tpe_aln_seq=args.tpe_seq, 
        min_frequency=args.min_frequency_alleles,
        max_n_rows=args.max_n_rows,
        plot_full_reads=args.plot_full_reads,
        fig_root=twinspector_results_folder,
        produce_pdf=args.produce_pdf, 
        ref_type="wt" 
    )


def generate_report(twinspector_results_folder, parent_folder, html_filename="TwInsPEctor_report.html"):
    """
    Bundle every PNG inside `twinspector_results_folder` into a single HTML report.
    """

    if not os.path.isdir(twinspector_results_folder):
        raise ValueError(f"{twinspector_results_folder} is not a valid directory")

    images = sorted(
        f for f in os.listdir(twinspector_results_folder)
        if f.lower().endswith(".png")
    )

    if not images:
        raise ValueError("No PNG files found in folder.")

    # Group allele-table alignment variants
    groups = {}
    alignment_groups = {}

    toggle_regex = re.compile(
        r"^(?P<stem>.+?)\.aligned_to_(?P<target>wt|tpe|composite_a|composite_b)\.png$"
    )

    for img in images:
        m = toggle_regex.match(img)

        if m:
            key = m.group("stem")              # e.g. b1.Perfect_TPE
            target = m.group("target")         # wt/tpe/composite_a/composite_b
            alignment_groups.setdefault(key, {})[target] = img
        else:
            groups[img] = {
                "title": img,
                "single": img,
            }

    # Helper to build any multi-image card
    def _make_special_group(group_key, title, files, button_texts):
        for f in files:
            groups.pop(f, None)
        variants = {}
        for i, fname in enumerate(files):
            v_key = chr(ord("a") + i)  # 'a', 'b', 'c',...
            variants[v_key] = {
                "img": fname,
                "text": button_texts[i] if i < len(button_texts) else v_key.upper(),
            }
        groups[group_key] = {
            "title": title,
            "variants": variants,
        }

    # Custom titles for allele table cards
    alignment_titles = {
        "b1.Perfect_TPE": "Perfect TPE Alleles",
        "b2.TPE_Indel": "TPE Alleles with Indels",
        "b3.Left_Flap": "Left Flap TPE Alleles",
        "b4.Right_Flap": "Right Flap TPE Alleles",
        "b5.Imperfect_TPE": "Imperfect TPE Alleles",
        "b6.Imperfect_WT": "Imperfect WT Alleles",
        "b7.WT_Indel": "WT Alleles with Indels",
        "b8.WT": "WT Alleles",
    }

    alignment_order = [
        ("composite_a", "Composite A"),
        ("composite_b", "Composite B"), 
        ("tpe", "TPE"), 
        ("wt", "WT")
    ]

    for key, target_map in alignment_groups.items():

        present = [t for t, _ in alignment_order if t in target_map]

        # If only one image exists, treat it like a normal card
        if len(present) == 1:
            img = target_map[present[0]]
            groups[img] = {
                "title": alignment_titles.get(key, key),
                "single": img,
            }
            continue

        files = [target_map[t] for t in present]
        button_texts = [
            label
            for t, label in alignment_order
            if t in target_map
        ]

        _make_special_group(
            key,
            alignment_titles.get(key, key),
            files,
            button_texts,
        )

    # First special card: a1, a2, a3
    group1_files = [
        "a3.Category_summary.png", 
        "a2.Category_summary_stacked.png", 
        "a1.Reads_input_summary.png"
    ]
    existing_group1 = [f for f in group1_files if f in images]
    if existing_group1:
        _make_special_group(
            "a1-3.summary",
            "Summary",
            existing_group1,
            ["Categories", "Categories stacked", "Reads"][:len(existing_group1)],
        )

    # Second special card: a5, a6, a10.a
    group2_files = [
        "a4.3'_flap_integration_all_reads.png",
        "a5.3'_flap_integration_edited_reads.png",
        "a6.5'_flap_removal_all_reads.png",
        "a7.5'_flap_removal_edited_reads.png",
        "a8.a.Edit_type_summary.png",
    ]

    existing_group2 = [f for f in group2_files if f in images]
    if existing_group2:
        _make_special_group(
            "a4-8.TPE_profiles",
            "Editing Per Base",
            existing_group2,
            ["3' Flaps - All reads", "3' Flaps - Edited reads", "5' Flaps - All reads", "5' Flaps - Edited reads", "Edit type - Full edit"][:len(existing_group2)],
        )

    group3_files = [
        "a9.Total_TPEs_vs_indels.png",
        "a10.Left_TPEs_vs_indels.png",
        "a11.Right_TPEs_vs_indels.png",
    ]
    existing_group3 = [f for f in group3_files if f in images]
    if existing_group3:
        _make_special_group(
            "a7-9.TPEs_vs_indels",
            "Indels Per Base",
            existing_group3,
            ["Total", "From Left", "From Right"][:len(existing_group3)],
        )

    cards_html = []
    uid = 0

    for key in sorted(groups):
        entry = groups[key]
        raw_title = entry["title"]
        if raw_title.lower().endswith(".png"):
            raw_title = os.path.splitext(raw_title)[0]
        title = html.escape(raw_title)

        if "single" in entry:
            # point HTML to Outputs/<filename>
            img_rel = os.path.join("Outputs", entry["single"])
            img_rel = html.escape(img_rel)
            cards_html.append(f"""
            <div class="card mb-4 shadow-sm">
              <div class="card-header fw-semibold">{title}</div>
              <div class="card-body text-center">
                <div class="img-panel">
                  <img src="{img_rel}" class="report-img"/>
                </div>
              </div>
            </div>
            """)
        else:
            uid += 1
            buttons = []
            panels = []

            variant_items = sorted(entry["variants"].items(), key=lambda kv: kv[0])
            for idx, (v_key, v_info) in enumerate(variant_items):
                # point HTML to Outputs/<filename>
                img_rel = os.path.join("TwInsPEctor_results", v_info["img"])
                img_rel = html.escape(img_rel)
                panel_id = f"panel-{uid}-{v_key}"
                is_first = (idx == 0)
                active = "active" if is_first else ""
                hidden = "" if is_first else "d-none"
                label_text = html.escape(v_info.get("text", v_key.upper()))

                buttons.append(
                    f'<button type="button" '
                    f'class="btn btn-outline-primary {active}" '
                    f'data-target="{panel_id}">{label_text}</button>'
                )

                panels.append(f"""
                    <div id="{panel_id}" class="img-panel {hidden}">
                      <img src="{img_rel}" class="report-img"/>
                    </div>
                """)

            cards_html.append(f"""
            <div class="card mb-4 shadow-sm">
              <div class="card-header fw-semibold">{title}</div>
              <div class="card-body text-center">
                <div class="btn-group mb-3" role="group">
                  {''.join(buttons)}
                </div>
                {''.join(panels)}
              </div>
            </div>
            """)

    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>TwInsPEctor</title>
  <link rel="stylesheet"
        href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.3/css/bootstrap.min.css">
  <style>
    body {{
      font-family: 'Montserrat', sans-serif;
      background-color: #f7f9fc;
    }}

    .report-wrapper {{
      width: 100%;
      margin: 0 auto;
    }}

    .card {{
        width: 100%;
    }}

    .img-panel {{
      background-color: #ffffff;
      border-radius: 6px;
      padding: 0rem;
      height: 68vh;
      display: flex;
      justify-content: center;
      align-items: center;
      overflow: hidden;
    }}

    .img-panel.zoomed-panel {{
      overflow: auto;
    #   cursor: grab;
    }}

    .report-img {{
      display: block;
      max-width: 100%;
      max-height: 100%;
      width: auto;
      height: auto;
      border: 1px solid #e0e0e0;
      border-radius: 4px;
      object-fit: contain;
      cursor: zoom-in;
      transition: transform 0.2s ease-out;
    }}

    .report-img.zoomed {{
      cursor: zoom-out;
    }}

    .report-path {{
      font-size: 0.75rem;
      word-break: break-all;
      margin-top: 0.25rem;
      margin-bottom: 0rem;
    }}

    /* Logo styles */
    .wrapper {{
      text-align: center;
      margin-top: -1.5rem;
      margin-bottom: 0.25rem;
    }}

    .logo {{
      font-size: 42px;
      font-weight: 350;
      letter-spacing: 1px;
      display: inline-flex;
      align-items: baseline;
    }}

    .twinpe {{
      color: black;
    }}

    .spector {{
      color: gray;
      font-weight: 300;
    }}

    .magnifier {{
      width: 50px;
      height: 60px;
      margin: 0 -8px;
      overflow: visible;
    }}

    .subtitle {{
      font-size: 8px;
      letter-spacing: 2px;
      color: gray;
      margin-top: -8px;
    }}
  </style>
</head>
<body class="p-4">
  <div class="report-wrapper">
    <div class="wrapper">
      <div class="logo">
        <span class="twinpe">TwIn</span>
        <span class="spector">s</span>

        <svg class="magnifier" viewBox="0 0 94 120">
          <line x1="19" y1="59" x2="19" y2="125"
                stroke="black"
                stroke-width="5"
                stroke-linecap="round"/>
          <circle cx="50" cy="60" r="27"
                  stroke="black"
                  stroke-width="5"
                  fill="none"/>
          <circle cx="55" cy="55" r="16" fill="gray" opacity="0.72"/>
          <circle cx="55" cy="56" r="6" fill="black"/>
          <path d="
            M24 27
            C37 16, 47 13, 57 19
            S75 30, 80 28
            C83 28, 85 32, 81 34
            C71 41, 62 31, 53 25
            S40 19, 26 26
            Z
          " fill="black"/>
        </svg>

        <span class="twinpe">E</span>
        <span class="spector">ctor</span>
      </div>
      <div class="subtitle">TWIN PRIME EDITING ANALYSIS</div>
      <p class="text-muted report-path">{html.escape(parent_folder)}</p>
    </div>

    {''.join(cards_html)}
  </div>

  <script>
    document.querySelectorAll('.btn-group button').forEach(function(btn) {{
      btn.addEventListener('click', function() {{
        const group = btn.closest('.card-body');
        group.querySelectorAll('.btn-group button')
             .forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        group.querySelectorAll('.img-panel')
             .forEach(p => p.classList.add('d-none'));
        const target = document.getElementById(btn.dataset.target);
        if (target) target.classList.remove('d-none');
      }});
    }});

    const ZOOM_LEVELS = [1, 1.7, 3];

    document.querySelectorAll('.report-img').forEach(function(img) {{
      img.dataset.zoom = '1';
      img.addEventListener('click', function() {{
        const panel = img.closest('.img-panel');
        let current = parseFloat(img.dataset.zoom || '1');
        let idx = ZOOM_LEVELS.indexOf(current);
        if (idx === -1) idx = 0;
        idx = (idx + 1) % ZOOM_LEVELS.length;
        const next = ZOOM_LEVELS[idx];
        img.dataset.zoom = String(next);

        if (next === 1) {{
          img.classList.remove('zoomed');
          img.style.transform = '';
          img.style.transformOrigin = '';
          if (panel) {{
            panel.classList.remove('zoomed-panel');
            panel.scrollTop = 0;
            panel.scrollLeft = 0;
          }}
        }} else {{
          img.classList.add('zoomed');
          img.style.transformOrigin = 'top left';
          img.style.transform = 'scale(' + next + ')';
          if (panel) {{
            panel.classList.add('zoomed-panel');
          }}
        }}
      }});
    }});
  </script>
</body>
</html>
"""

    out_path = os.path.join(parent_folder, html_filename)
    with open(out_path, "w", encoding="utf-8") as fout:
        fout.write(html_out)
    return out_path


def parse_args():
    prog = 'TwInsPEctor'
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Analyzes amplicon sequencing data from twin prime editing experiments by aligning reads to standard and composite references using CRISPResso2. Allele outcomes are classified into 8 categories with detailed visualizations.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            f'Example: {prog}'
            "-r1 <fastq R1 file> -r2 <fastq R2 file> "
            "-w <full wild-type amplicon sequence> -t <full twin prime edited amplicon sequence> "
            "-g <pegRNA-a spacer sequence>,<pegRNA-b spacer sequence>\n\n"
            "----Category Definitions----\n"
            "Perfect TPE: complete programmed edit without indels.\n"
            "TPE Indel: complete programmed edit with indels.\n"
            "Left Flap: at least N consecutive programmed bases starting from the left (pegRNA-a) but not from the right.\n"
            "Right Flap: at least N consecutive programmed bases starting from the right (pegRNA-b) but not from the left.\n"
            "Imperfect TPE: incomplete programmed edit (contains TPE sequence but meets neither or both flap criteria).\n"
            "Imperfect WT: incomplete wildtype sequence and none of the programmed edit.\n"
            "WT Indel: complete wildtype sequence with indels and none of the programmed edit.\n"
            "WT: complete wildtype sequence without indels and none of the programmed edit.\n"
            "Uncategorized: does not fit into any category.\n"
        )
    )

    parser.add_argument("-r1", "--fastq_r1", type=str, required=True, help="Path to FASTQ R1 file.")
    parser.add_argument("-r2", "--fastq_r2", type=str, required=False, help="Path to FASTQ R2 file for paired-end data.")
    parser.add_argument("-w", "--wt_seq", type=str, required=True, help="Full wild-type reference amplicon sequence including spacers.")
    parser.add_argument("-t", "--tpe_seq", type=str, required=True, help="Full Twin prime edited reference amplicon sequence with 5' & 3' ends identical to wildtype reference amplicon.")
    parser.add_argument("-g", "--peg_spacers", type=str, required=True, help="Comma-separated pegRNA spacer sequences: <spacer A>,<spacer B>. Should include bases immediately adjacent to but not including the PAM sequence (usually 20nt 5' of NGG).")
    parser.add_argument("-o", "--output_root", type=str, default=None, help="Root output folder for CRISPResso2 and TwInsPEctor results. If not provided, a folder will be created in the current working directory based on the input fastq file names.")
    parser.add_argument("-ne", "--num_changes_to_check", type=int, default=2, help="Minimum number of programmed bases that must be edited for read to be classified.")
    parser.add_argument("-rcm", "--recoding_mode", action="store_true", help="Run in recoding mode if the wild-type and twin prime edited sequences are the same length and should be evaluated as having only base substitutions.")
    parser.add_argument("-dmas", "--default_min_aln_score", type=int, default=30, help="Default minimum homology score for a read to align to the compound reference amplicon")
    parser.add_argument("-pfr", "--plot_full_reads", action="store_true", help="Allele tables will display full read sequences.")
    parser.add_argument("-ied", "--ignore_extraspacer_deletions", action="store_true", help="Classification ignores deletions occurring beyond the spacers (outside edit window).")
    parser.add_argument("-np", "--no_plots", action="store_true", help="Skip all figures if only text outputs are desired.")
    parser.add_argument("-nsp", "--no_summary_plots", action="store_true", help="Skip summary barplots if they are not desired.")
    parser.add_argument("-npp", "--no_per_base_plots", action="store_true", help="Skip per-base barplots if they are not desired.")
    parser.add_argument("-nat", "--no_allele_tables", action="store_true", help="Skip allele tables if they are not desired.")
    parser.add_argument("-pdf", "--produce_pdf", action="store_true", help="Produce PDF versions of all plots in addition to PNG versions.")
    parser.add_argument("-mfa", "--min_frequency_alleles", type=float, default=0.2, help="Minimum percent read frequency required to report an allele in the alleles tables.")
    parser.add_argument("-mnr", "--max_n_rows", type=int, default=30, help="Maximum number of allele rows to display in the allele tables.")
    parser.add_argument("-nrr", "--no_rerun", action="store_true", help="Don't rerun CRISPResso2 if a run using the same parameters has already been finished.")
    parser.add_argument("-kco", "--keep_crispresso_outputs", action="store_true", help="Don't delete CRISPResso2 output folders after analysis.")
    parser.add_argument("-ts", "--trim_string", type=str, default=None, help="String to trim reads using fastp with override options within CRISPResso2 before analysis.")
    parser.add_argument("-fp", "--fastp_command", type=str, default=None, help="Command to run fastp for read trimming within CRISPResso2 before analysis.")
    parser.add_argument("-coa", "--cleavage_offset_a", type=int, default=-3, help="Cleavage offset for pegRNA spacer A (default: -3).")
    parser.add_argument("-cob", "--cleavage_offset_b", type=int, default=-3, help="Cleavage offset for pegRNA spacer B (default: -3).")
    parser.add_argument("-nmt", "--no_multithreading", action="store_true", help="Run CRISPResso2 and allele plotting sequentially instead of in parallel (significantly increases runtime).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose CRISPResso2 output.")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-V", "--version", action="version", version="%(prog)s 1.0.0")

    args = parser.parse_args()

    args.wt_seq = args.wt_seq.upper()
    args.tpe_seq = args.tpe_seq.upper()

    return args


def main():
    print("\nStarting TwInsPEction...")
    args = parse_args()

    parent_folder, crispresso_wt, crispresso_tpe, crispresso_composite_a, crispresso_composite_b, twinspector_results_folder = get_folder_names(args)
    os.makedirs(twinspector_results_folder, exist_ok=True)
    os.makedirs(crispresso_wt, exist_ok=True)
    os.makedirs(crispresso_tpe, exist_ok=True)
    os.makedirs(crispresso_composite_a, exist_ok=True)
    os.makedirs(crispresso_composite_b, exist_ok=True)

    print("Analyzing reference inputs...")
    spacer_a, spacer_b = get_spacer_seqs(args.peg_spacers)
    reference_info = analyze_references(args.wt_seq, args.tpe_seq, spacer_a, spacer_b, cleavage_offset_a=args.cleavage_offset_a, cleavage_offset_b=args.cleavage_offset_b, output_root=twinspector_results_folder, recoding_mode=args.recoding_mode)

    crispresso_cmd_wt = get_crispresso_command(args=args, ref_seq=args.wt_seq, ref_name="WT", spacer_a=reference_info["spacer_a_wt"], spacer_b=reference_info["spacer_b_wt"], crispresso_output_folder=crispresso_wt, twinspector_results_folder=twinspector_results_folder, append=False)
    crispresso_cmd_tpe = get_crispresso_command(args=args, ref_seq=args.tpe_seq, ref_name="TPE", spacer_a=reference_info["spacer_a_tpe"], spacer_b=reference_info["spacer_b_tpe"], crispresso_output_folder=crispresso_tpe, twinspector_results_folder=twinspector_results_folder)
    crispresso_cmd_composite_a = get_crispresso_command(args=args, ref_seq=reference_info["composite_a_ref_seq"], ref_name="Composite_A", spacer_a=reference_info["spacer_a_composite_a"], spacer_b=reference_info["spacer_b_composite_a"], crispresso_output_folder=crispresso_composite_a, twinspector_results_folder=twinspector_results_folder)
    crispresso_cmd_composite_b = get_crispresso_command(args=args, ref_seq=reference_info["composite_b_ref_seq"], ref_name="Composite_B", spacer_a=reference_info["spacer_a_composite_b"], spacer_b=reference_info["spacer_b_composite_b"], crispresso_output_folder=crispresso_composite_b, twinspector_results_folder=twinspector_results_folder)

    if args.no_multithreading:
        print("Running CRISPResso2...")
        run_crispresso_command(crispresso_cmd_wt, verbose=args.verbose)
        run_crispresso_command(crispresso_cmd_tpe, verbose=args.verbose)
        run_crispresso_command(crispresso_cmd_composite_a, verbose=args.verbose)
        run_crispresso_command(crispresso_cmd_composite_b, verbose=args.verbose)
    else:
        print("Running CRISPResso2 in parallel...")
        crispresso_tasks = [crispresso_cmd_wt, crispresso_cmd_tpe, crispresso_cmd_composite_a, crispresso_cmd_composite_b]
        run_crispresso_commands_parallel(crispresso_tasks, verbose=args.verbose)

    print("Analyzing CRISPResso2 outputs...")
    df_merged = merge_crispresso_allele_tables(crispresso_wt=crispresso_wt, crispresso_tpe=crispresso_tpe, crispresso_composite_a=crispresso_composite_a, crispresso_composite_b=crispresso_composite_b)
    df_categorized, bp_changes_arrs = categorize_alleles(df_merged=df_merged, wt_seq=args.wt_seq, tpe_seq=args.tpe_seq, reference_info=reference_info, num_changes_to_check=args.num_changes_to_check, ignore_extraspacer_deletions=args.ignore_extraspacer_deletions, default_min_aln_score=args.default_min_aln_score, recoding_mode=args.recoding_mode)
    plotting_info = get_plotting_stats(df=df_categorized, reference_info=reference_info, bp_changes_arrs=bp_changes_arrs, twinspector_results_folder=twinspector_results_folder, recoding_mode=args.recoding_mode)

    if not args.no_plots:
        if not args.no_summary_plots:
            print("Generating summary plots...")
            plot_summary_barplots(category_counts=plotting_info["category_counts"], crispresso_output_folder=crispresso_wt, twinspector_results_folder=twinspector_results_folder, produce_pdf=args.produce_pdf)

        if not args.no_per_base_plots:
            print("Generating per-base plots...")
            plot_per_base_pos_barplots(plotting_info=plotting_info, twinspector_results_folder=twinspector_results_folder, produce_pdf=args.produce_pdf, recoding_mode=args.recoding_mode)

        if not args.no_allele_tables:
            if args.no_multithreading:
                print("Generating allele tables...")
                plot_ref_allele_tables_sequential(args=args, df_categorized=df_categorized, reference_info=reference_info, twinspector_results_folder=twinspector_results_folder)
            else:
                print("Generating allele tables in parallel...")
                plot_ref_allele_tables_parallel_ppe(args=args, df_categorized=df_categorized, reference_info=reference_info, twinspector_results_folder=twinspector_results_folder)
                # plot_ref_allele_tables_parallel_tpe(args=args, df_categorized=df_categorized, reference_info=reference_info, twinspector_results_folder=twinspector_results_folder)

        print("Generating report...")
        generate_report(twinspector_results_folder, parent_folder)

    # Safe deletion of CRISPResso outputs
    if not args.keep_crispresso_outputs:
        parent_folder = os.path.abspath(parent_folder)
        for name in ["CRISPResso_wt", "CRISPResso_tpe", "CRISPResso_composite_a", "CRISPResso_composite_b"]:
            folder = os.path.join(parent_folder, name)
            if os.path.isdir(folder):
                shutil.rmtree(folder)

    print("TwInsPEction complete!")
    sys.exit(0)


if __name__ == "__main__":
    main()
