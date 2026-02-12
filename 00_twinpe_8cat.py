import os
import re
import sys
import zipfile
import argparse
import matplotlib
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from CRISPResso2 import CRISPRessoShared
# from CRISPResso2 import CRISPResso2Align
from matplotlib import colors as colors_mpl


# TODO:
# Add recoding mode.
# Add mismatched prefix/suffix support. 
# Add multi-sample support.


CATEGORY_COLORS = {
    "Perfect PE":  "#1f77b4", 
    "PE Indel": "#9467bd", 
    "Left Flap": "#ff7f0e", 
    "Right Flap": "#2ca02c", 
    "Imperfect PE": "#d62728", 
    "Imperfect WT": "#e377c2", 
    "WT Indel":   "#bcbd22", 
    "WT": "#8c564b", 
    "Uncategorized": "#7f7f7f"
}


def main():
    args = parse_args()

    parent_folder, crispresso_output_folder_a, crispresso_output_folder_b, twinpe_8cat_results_folder = resolve_output_folders(args)
    crispresso_run_a_folder = os.path.join(parent_folder, "Crispresso_outputs", "Run_a")
    crispresso_run_b_folder = os.path.join(parent_folder, "Crispresso_outputs", "Run_b")
    # twinpe_8cat_a_folder = os.path.join(twinpe_8cat_results_folder, "Run_a")
    # twinpe_8cat_b_folder = os.path.join(twinpe_8cat_results_folder, "Run_b")

    # os.makedirs(twinpe_8cat_a_folder, exist_ok=True)
    # os.makedirs(twinpe_8cat_b_folder, exist_ok=True)
    os.makedirs(twinpe_8cat_results_folder, exist_ok=True)
    os.makedirs(crispresso_run_a_folder, exist_ok=True)
    os.makedirs(crispresso_run_b_folder, exist_ok=True)

    spacer_info = find_spacers_in_references(args.wt_seq, args.twin_seq, args.peg_spacers[0], args.peg_spacers[1])
    comp_ref_seq_a, wt_aln_seq_a, twin_aln_seq_a, comp_ref_seq_b, wt_aln_seq_b, twin_aln_seq_b = build_compound_reference_alignments(args.wt_seq, args.twin_seq, args.peg_spacers[0], args.peg_spacers[1], spacer_info, twinpe_8cat_results_folder)
    # compound_ref_seq = build_compound_reference(args.wt_seq, args.twin_seq)

    crispresso_cmd_a = build_crispresso_command(args, comp_ref_seq_a, crispresso_run_a_folder, spacer_info, twinpe_8cat_results_folder, run_label="a")
    crispresso_cmd_b = build_crispresso_command(args, comp_ref_seq_b, crispresso_run_b_folder, spacer_info, twinpe_8cat_results_folder, run_label="b")

    print("Running CRISPResso with command:\n", " ".join(crispresso_cmd_a), "\n")
    subprocess.run(crispresso_cmd_a, check=True)
    print("Running CRISPResso with command:\n", " ".join(crispresso_cmd_b), "\n")
    subprocess.run(crispresso_cmd_b, check=True)

    print("Analyzing CRISPResso output...")
    analyze_visualize_sample(twinpe_8cat_results_folder, crispresso_output_folder_a, crispresso_output_folder_b, args, comp_ref_seq_a, wt_aln_seq_a, twin_aln_seq_a, comp_ref_seq_b, wt_aln_seq_b, twin_aln_seq_b, spacer_info) # twinpe_8cat_a_folder, twinpe_8cat_b_folder, 
    print("Finished TwinPE analysis!")

    sys.exit(0)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="04_twinpe_8cat_compound_crispresso.py",
        description="Analyzes Twin Prime Editing outcomes by running CRISPResso2 on raw fastq sequencing files with alignment to a WT-TwinPE compound reference. Classifies reads into 8 categories and provides detailed visualizations.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Example: python 04_twinpe_8cat_compound_crispresso.py "
            "-r1 <fastq R1 file> -r2 <fastq R2 file> "
            "-w <full wildtype sequence> -t <full twinpe sequence> "
            "-g <peg spacer A sequence> <peg spacer B sequence>\n\n"
            "----Category Definitions----\n"
            "Perfect PE: complete programmed edit without indels anywhere.\n"
            "PE Indel: complete programmed edit with indels anywhere.\n"
            "Left Flap: at least N consecutive programmed bases starting from the left and not from the right.\n"
            "Right Flap: at least N consecutive programmed bases starting from the right and not from the left.\n"
            "Imperfect PE: incomplete or incorrect programmed edit.\n"
            "Imperfect WT: some but not all of the wildtype sequence and none of the programmed edit.\n"
            "WT Indel: full wildtype sequence with indels anywhere and none of the programmed edit.\n"
            "WT: full wildtype sequence without indels anywhere and none of the programmed edit.\n"
            "Uncategorized: does not fit into any category.\n"
        )
    )

    parser.add_argument("-r1", "--fastq_r1", type=str, required=True, help="Path to FASTQ R1 file.")
    parser.add_argument("-r2", "--fastq_r2", type=str, required=False, help="Path to FASTQ R2 file (optional).")
    parser.add_argument("-w", "--wt_seq", type=str, required=True, help="Full wildtype reference amplicon sequence.")
    parser.add_argument("-t", "--twin_seq", type=str, required=True, help="Full TwinPE reference amplicon sequence with 5' & 3' ends identical to wildtype reference amplicon.")
    parser.add_argument("-g", "--peg_spacers", type=str, nargs=2, required=True, help="Space-separated pegRNA spacer sequences: <spacer A> <spacer B>")
    parser.add_argument("-o", "--output_root", type=str, default=None, help="Root output folder for CRISPResso2 and TwinPE 8cat results. If not provided, a folder will be created in the current working directory based on the input fastq file names.")
    parser.add_argument("-n", "--num_changes_to_check", type=int, default=2, help="Minimum number of programmed base edits required for read to be classified as left/right flap allele.")
    parser.add_argument("-dmas", "--default_min_aln_score", type=int, default=50, help="Default minimum homology score for a read to align to the compound reference amplicon")
    parser.add_argument("-f", "--plot_full_reads", action="store_true", help="Allele tables will display full read sequences.")
    parser.add_argument("--ignore_extraspacer_deletions", action="store_true", help="Classification ignores deletions occurring beyond the spacers (outside edit window).")
    parser.add_argument("--produce_png", action="store_true", help="Produce PNG versions of all plots in addition to PDF versions.")
    parser.add_argument("--min_frequency_alleles", type=float, default=0.0, help="Minimum percent read frequency required to report an allele in the alleles table plot.")
    parser.add_argument("--max_n_rows", type=int, default=50, help="Maximum number of allele rows to display in the allele table plot.")
    parser.add_argument("--no_rerun", action="store_true", help="Don't rerun CRISPResso2 if a run using the same parameters has already been finished.")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-V", "--version", action="version", version="%(prog)s 1.0")

    args = parser.parse_args()

    args.wt_seq = args.wt_seq.upper()
    args.twin_seq = args.twin_seq.upper()
    args.peg_spacers = [s.upper() for s in args.peg_spacers]

    return args


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    return ''.join(complement.get(base, base) for base in reversed(seq.upper()))


def find_spacers_in_references(wt_seq, twin_seq, spacer_a, spacer_b):
    """
    Searches reference sequences for spacers or their reverse complements.
    Returns indices of spacers in both sequences, whether they were found as reverse complements,
    and number of bases removed cleaved by nicking.
    """
    # May need to make sure spacers are found at different indices for the cases where they could be the same sequence
    # Find spacers in wt_seq
    spacer_a_rc = False
    spacer_a_index_wt_f = wt_seq.find(spacer_a)
    spacer_a_index_wt_rc = -1
    if spacer_a_index_wt_f == -1:
        spacer_a_index_wt_rc = wt_seq.find(reverse_complement(spacer_a))
        if spacer_a_index_wt_rc != -1:
            spacer_a_rc = True
    if spacer_a_index_wt_f == -1 and spacer_a_index_wt_rc == -1:
        raise ValueError("Could not find peg spacer A in WT Sequence")
    spacer_b_rc = False
    spacer_b_index_wt_f = wt_seq.find(spacer_b)
    spacer_b_index_wt_rc = -1
    if spacer_b_index_wt_f == -1:
        spacer_b_index_wt_rc = wt_seq.find(reverse_complement(spacer_b))
        if spacer_b_index_wt_rc != -1:
            spacer_b_rc = True
    if spacer_b_index_wt_f == -1 and spacer_b_index_wt_rc == -1:
        raise ValueError("Could not find peg spacer B in WT Sequence")
    spacer_a_index_wt = spacer_a_index_wt_f if spacer_a_index_wt_f != -1 else spacer_a_index_wt_rc
    spacer_b_index_wt = spacer_b_index_wt_rc if spacer_b_index_wt_rc != -1 else spacer_b_index_wt_f

    # Find spacers in twin_seq
    spacer_a_num_bases_removed = 0
    for idx in range(len(spacer_a)):
        if not spacer_a_rc:
            spacer_a_subseq = spacer_a[:len(spacer_a)-idx]
        else:
            spacer_a_subseq = reverse_complement(spacer_a)[idx:len(spacer_a)]
        spacer_a_index_twin = twin_seq.find(spacer_a_subseq)
        if spacer_a_index_twin != -1:
            break
        spacer_a_num_bases_removed += 1
        if len(spacer_a)-spacer_a_num_bases_removed < 6: # better length cutoff?
            raise ValueError("Could not find peg spacer A in Twin Sequence")
    spacer_b_num_bases_removed = 0
    for idx in range(len(spacer_b)):
        if not spacer_b_rc:
            spacer_b_subseq = spacer_b[:len(spacer_b)-idx]
        else:
            spacer_b_subseq = reverse_complement(spacer_b)[idx:len(spacer_b)]
        spacer_b_index_twin = twin_seq.find(spacer_b_subseq)
        if spacer_b_index_twin != -1:
            break
        spacer_b_num_bases_removed += 1
        if len(spacer_b)-spacer_b_num_bases_removed < 6: # better length cutoff?
            raise ValueError("Could not find peg spacer B in Twin Sequence")

    results = {
        "spacer_a_index_wt": spacer_a_index_wt,
        "spacer_b_index_wt": spacer_b_index_wt,
        "spacer_a_index_twin": spacer_a_index_twin,
        "spacer_b_index_twin": spacer_b_index_twin,
        "spacer_a_rc": spacer_a_rc,
        "spacer_b_rc": spacer_b_rc,
        "spacer_a_num_bases_removed": spacer_a_num_bases_removed,
        "spacer_b_num_bases_removed": spacer_b_num_bases_removed
    }

    return results


def build_compound_reference_alignments(wt_seq, twin_seq, spacer_a, spacer_b, spacer_info, output_root):
    """
    Requires that wt_seq and twin_seq share identical 5' and 3' anchors.
    """
    spacer_a_coverage_wt = range(spacer_info['spacer_a_index_wt'], spacer_info['spacer_a_index_wt'] + len(spacer_a))
    spacer_b_coverage_wt = range(spacer_info['spacer_b_index_wt'], spacer_info['spacer_b_index_wt'] + len(spacer_b))
    prefix = ''
    for wt_base, twin_base in zip(wt_seq, twin_seq):
        if wt_base == twin_base:
            prefix += wt_base
        else:
            break
    # validate that prefix ends within spacer_a region
    # if len(prefix) not in spacer_a_coverage_wt and len(prefix) not in spacer_b_coverage_wt:
    if len(prefix) not in spacer_a_coverage_wt:
        raise ValueError("5' anchor does not terminate within Spacer A region")
    suffix = ''
    for wt_base, twin_base in zip(wt_seq[::-1], twin_seq[::-1]):
        if wt_base == twin_base:
            suffix += wt_base
        else:
            break
    suffix = suffix[::-1]
    # validate that suffix starts within spacer_b region
    # if wt_seq.find(suffix) not in spacer_a_coverage_wt or wt_seq.find(suffix) not in spacer_b_coverage_wt:
    if wt_seq.find(suffix) not in spacer_b_coverage_wt:
        raise ValueError("3' anchor does not terminate within Spacer B region")

    prefix_and_deletion = wt_seq[:wt_seq.find(suffix)]
    insertion_and_suffix = twin_seq[len(prefix):]
    comp_ref_seq_a = prefix_and_deletion + insertion_and_suffix

    prefix_and_insertion = twin_seq[:twin_seq.find(suffix)]
    deletion_and_suffix = wt_seq[len(prefix):]
    comp_ref_seq_b = prefix_and_insertion + deletion_and_suffix

    # build WT and TwinPE reference alignments to Compound Reference
    wt_aln_seq_a = prefix_and_deletion + (len(insertion_and_suffix)-len(suffix)) * '-' + suffix
    twin_aln_seq_a = prefix + (len(prefix_and_deletion)-len(prefix)) * '-' + insertion_and_suffix
    wt_aln_seq_b = prefix + (len(prefix_and_insertion)-len(prefix)) * '-' + deletion_and_suffix
    twin_aln_seq_b = prefix_and_insertion + (len(deletion_and_suffix)-len(suffix)) * '-' + suffix

    # validate lengths
    if not (len(comp_ref_seq_a) == len(wt_aln_seq_a) == len(twin_aln_seq_a)):
        raise ValueError("Compound reference, WT alignment, and Twin alignment sequences are not the same length")
    
    with open(
        os.path.join(output_root, "c5.aligned_reference_sequences.txt"), "w"
    ) as fout:
        fout.write(f"@Alignment A\n")
        fout.write(f">Wildtype_reference_sequence\n{wt_aln_seq_a}\n")
        fout.write(f">Compound_reference_sequence\n{comp_ref_seq_a}\n")
        fout.write(f">TwinPE_reference_sequence\n{twin_aln_seq_a}\n\n")
        fout.write(f"@Alignment B\n")
        fout.write(f">Wildtype_reference_sequence\n{wt_aln_seq_b}\n")
        fout.write(f">Compound_reference_sequence\n{comp_ref_seq_b}\n")
        fout.write(f">TwinPE_reference_sequence\n{twin_aln_seq_b}\n\n")
        fout.write(f"@pegRNA Spacers\n")
        fout.write(f">SpacerA_sequence\n{spacer_a}\n")
        fout.write(f">SpacerB_sequence\n{spacer_b}\n")

    return comp_ref_seq_a, wt_aln_seq_a, twin_aln_seq_a, comp_ref_seq_b, wt_aln_seq_b, twin_aln_seq_b

# def build_compound_reference(wt_seq, twin_seq):
    # """
    # Requires that wt_seq and twin_seq share identical 5' and 3' anchor sequences.
    # """
#     prefix = ''
#     for wt_base, twin_base in zip(wt_seq, twin_seq):
#         if wt_base == twin_base:
#             prefix += wt_base
#         else:
#             break
#     if prefix == '':
#         raise ValueError("Could not find matching prefix anchor between provided WT and TwinPE sequences")
#     suffix = ''
#     for wt_base, twin_base in zip(wt_seq[::-1], twin_seq[::-1]):
#         if wt_base == twin_base:
#             suffix += wt_base
#         else:
#             break
#     if suffix == '':
#         raise ValueError("Could not find matching suffix anchor between provided WT and TwinPE sequences")
#     suffix = suffix[::-1]

#     prefix_and_deletion = wt_seq[:wt_seq.find(suffix)]
#     insertion_and_suffix = twin_seq[len(prefix):]
#     compound_ref_seq = prefix_and_deletion + insertion_and_suffix

#     return compound_ref_seq


def resolve_output_folders(args):
    r1 = args.fastq_r1
    r2 = args.fastq_r2 if args.fastq_r2 else None
    # any other extension types?
    pattern = r'([^/]+?)(?=(?:\.fastq|\.fq)?(?:\.gzip|\.gz|\.bz2|\.bz|\.xz|\.lzma)?$)'
    r1m = re.search(pattern, r1)
    r2m = re.search(pattern, r2) if r2 else None
    if args.output_root:
        # If output_folder is provided in command args
        # Use that as the parent folder for CRISPResso output and TwinPE 8cat results
        parent_folder = os.path.join(os.getcwd(), args.output_root.rstrip("/"))
        # Mimic CRISPResso output folder naming conventions to get correct path
        if r1m and r2m:
            crispresso_output_folder_a = os.path.join(parent_folder, "Crispresso_outputs","Run_a", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_output_folder_b = os.path.join(parent_folder, "Crispresso_outputs","Run_b", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
        elif r1m and not r2m:
            crispresso_output_folder_a = os.path.join(parent_folder, "Crispresso_outputs","Run_a", f"CRISPResso_on_{r1m.group(1)}")
            crispresso_output_folder_b = os.path.join(parent_folder, "Crispresso_outputs","Run_b", f"CRISPResso_on_{r1m.group(1)}")
        else:   
            raise ValueError("Could not parse fastq file names for output folder naming.")
    else:
        # If output_folder not provided, create own
        if r1m and r2m:
            parent_folder = os.path.join(os.getcwd(), f"TwinPE_8cat_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_output_folder_a = os.path.join(parent_folder, "Crispresso_outputs","Run_a", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
            crispresso_output_folder_b = os.path.join(parent_folder, "Crispresso_outputs","Run_b", f"CRISPResso_on_{r1m.group(1)}_{r2m.group(1)}")
        elif r1m and not r2m:
            parent_folder = os.path.join(os.getcwd(), f"TwinPE_8cat_on_{r1m.group(1)}")
            crispresso_output_folder_a = os.path.join(parent_folder, "Crispresso_outputs","Run_a", f"CRISPResso_on_{r1m.group(1)}")
            crispresso_output_folder_b = os.path.join(parent_folder, "Crispresso_outputs","Run_b", f"CRISPResso_on_{r1m.group(1)}")
        else:
            raise ValueError("Could not parse fastq file names for output folder naming.")
        
    twinpe_8cat_results_folder = os.path.join(parent_folder, "TwinPE_8cat_results")

    return parent_folder, crispresso_output_folder_a, crispresso_output_folder_b, twinpe_8cat_results_folder


def build_crispresso_command(args, compound_ref_seq, crispresso_output_folder, spacer_info, twinpe_8cat_output_folder, run_label):
    if spacer_info['spacer_a_index_wt'] < spacer_info['spacer_b_index_wt']:
        if run_label == "a":
            first_spacer = args.peg_spacers[0]
            second_spacer = args.peg_spacers[1][:len(args.peg_spacers[1])-spacer_info['spacer_b_num_bases_removed']]
        else:
            first_spacer = args.peg_spacers[0][:len(args.peg_spacers[0])-spacer_info['spacer_a_num_bases_removed']]
            second_spacer = args.peg_spacers[1]
    else:
        if run_label == "a":
            first_spacer = args.peg_spacers[1]
            second_spacer = args.peg_spacers[0][:len(args.peg_spacers[0])-spacer_info['spacer_a_num_bases_removed']]
        else:
            first_spacer = args.peg_spacers[1][:len(args.peg_spacers[1])-spacer_info['spacer_b_num_bases_removed']]
            second_spacer = args.peg_spacers[0]
    # May want to wrap in all CRISPResso parameters somehow or remove some of the hardcoded ones below
    cmd = [
        "CRISPResso",
        "--fastq_r1", args.fastq_r1,
        "--amplicon_seq", compound_ref_seq,
        "--amplicon_name", "Compound",
        "--guide_seq", f"{first_spacer},{second_spacer}", 
        "--default_min_aln_score", str(args.default_min_aln_score), 
        "--min_frequency_alleles_around_cut_to_plot", str(args.min_frequency_alleles), 
        "--max_rows_alleles_around_cut_to_plot", str(args.max_n_rows), 
        "--write_detailed_allele_table",
    ]

    cmd.extend(["--output_folder", crispresso_output_folder])

    if args.fastq_r2:
        cmd.extend(["--fastq_r2", args.fastq_r2])
    if args.no_rerun:
        cmd.append("--no_rerun")
    if args.debug:
        cmd.append("--debug")

    with open(
        os.path.join(twinpe_8cat_output_folder, "c6.crispresso2_command_a.txt" if run_label == "a" else "c7.crispresso2_command_b.txt"), "w"
    ) as fout:
        fout.write(" ".join(cmd) + "\n")

    return cmd


def analyze_visualize_sample(
        twinpe_8cat_results_folder, 
        # twinpe_8cat_a_folder, 
        # twinpe_8cat_b_folder, 
        crispresso_output_folder_a, 
        crispresso_output_folder_b, 
        args, 
        comp_ref_seq_a, 
        wt_aln_seq_a, 
        twin_aln_seq_a, 
        comp_ref_seq_b, 
        wt_aln_seq_b, 
        twin_aln_seq_b, 
        spacer_info,
    ):
    """
    Runs classification and plotting functions for a single sample.
    """
    df_alleles_a_results = categorize_alleles(crispresso_output_folder_a, comp_ref_seq_a, twin_aln_seq_a, wt_aln_seq_a, args.num_changes_to_check, args.ignore_extraspacer_deletions)

    df_alleles_b_results = categorize_alleles(crispresso_output_folder_b, comp_ref_seq_b, twin_aln_seq_b, wt_aln_seq_b, args.num_changes_to_check, args.ignore_extraspacer_deletions)

    df_alleles_final = collapse_allele_categories(df_alleles_a_results["df_alleles"], df_alleles_b_results["df_alleles"])

    results_final = analyze_collapsed_categorized_alleles(
        df_alleles_final, 
        # crispresso_output_folder=crispresso_output_folder_a,
        twinpe_8cat_results_folder=twinpe_8cat_results_folder, 
        # wt_seq=args.wt_seq,
        # twin_seq=args.twin_seq,
        # wt_aln_seq_a=wt_aln_seq_a, 
        # twin_aln_seq_a=twin_aln_seq_a,
        # comp_ref_seq_a=comp_ref_seq_a, 
        # wt_aln_seq_b=wt_aln_seq_b, 
        # twin_aln_seq_b=twin_aln_seq_b,
        # comp_ref_seq_b=comp_ref_seq_b, 
        bp_changes_arr_a=df_alleles_a_results["bp_changes_arr"], 
        del_start_a=df_alleles_a_results["del_start"], 
        del_end_a=df_alleles_a_results["del_end"], 
        ins_start_a=df_alleles_a_results["ins_start"], 
        ins_end_a=df_alleles_a_results["ins_end"], 
        ins_region_len_a=df_alleles_a_results["ins_region_len"], 
        pegRNA_intervals_a=df_alleles_a_results["pegRNA_intervals"], 
        bp_changes_arr_b=df_alleles_b_results["bp_changes_arr"], 
        del_start_b=df_alleles_b_results["del_start"], 
        del_end_b=df_alleles_b_results["del_end"], 
        ins_start_b=df_alleles_b_results["ins_start"], 
        ins_end_b=df_alleles_b_results["ins_end"], 
        ins_region_len_b=df_alleles_b_results["ins_region_len"], 
        pegRNA_intervals_b=df_alleles_b_results["pegRNA_intervals"], 
        # num_changes_to_check=args.num_changes_to_check,
        ignore_extraspacer_deletions=args.ignore_extraspacer_deletions, 
        # produce_png=args.produce_png
    )

    setBarMatplotlibDefaults()
    
    # for results, crispresso_output_folder, twinpe_8cat_results_folder in [(results_a["folder_category_counts"], crispresso_output_folder_a, twinpe_8cat_a_folder), (results_b["folder_category_counts"], crispresso_output_folder_b, twinpe_8cat_b_folder)]:
    plot_summary_barplots(
        results_final["folder_category_counts"],  
        crispresso_output_folder_a, 
        twinpe_8cat_results_folder, 
        args.produce_png
    )

    plot_per_base_pos_barplots(
        results_final, 
        df_alleles_a_results["ins_start"], 
        df_alleles_a_results["ins_end"], 
        df_alleles_a_results["del_start"], 
        twinpe_8cat_results_folder, 
        args.produce_png
    )

    setAlleleMatplotlibDefaults()

    # for results, twinpe_8cat_results_folder in [(results_a, twinpe_8cat_a_folder)]: # , (results_b, twinpe_8cat_b_folder)]: # Need to make many changes for this to work
    plot_categorical_allele_tables(
        args.min_frequency_alleles,
        args.max_n_rows, 
        df_alleles_final,
        args.wt_seq,
        wt_aln_seq_a,
        twin_aln_seq_a,
        df_alleles_a_results["pegRNA_cut_points"],
        df_alleles_a_results["pegRNA_plot_cut_points"],
        df_alleles_a_results["pegRNA_intervals"],
        df_alleles_a_results["pegRNA_mismatches"],
        df_alleles_a_results["pegRNA_names"],
        spacer_info=spacer_info,
        fig_root=twinpe_8cat_results_folder,
        produce_png=args.produce_png,
        plot_full_reads=args.plot_full_reads
    )


def plot_summary_barplots(folder_category_counts, crispresso_output_folder_a, twinpe_8cat_results_folder, produce_png):
    
    plot_reads_input_summary_barplot(
        crispresso_output_folder_a,
        fig_root=twinpe_8cat_results_folder,
        produce_png=produce_png
    )

    plot_category_stacked_summary_barplot(
        crispresso_output_folder_a,
        folder_category_counts, 
        fig_root=twinpe_8cat_results_folder, 
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS
    )

    # plot_stacked_summary_barplot(
    #     results["folder_category_counts"], 
    #     fig_root=twinpe_8cat_results_folder, 
    #     produce_png=args.produce_png
    # )

    plot_category_summary_barplot(
        folder_category_counts, 
        fig_root=twinpe_8cat_results_folder, 
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS
    )
    

def plot_per_base_pos_barplots(results_final, ins_start, ins_end, del_start, twinpe_8cat_results_folder, produce_png):
    
    plot_successful_twin_edit_counts_by_category(
        # results["bp_changes_arr"],
        results_final["edit_counts"],
        results_final["cat_perfect_pe_count_arr"],
        results_final["cat_left_flap_count_arr"],
        results_final["cat_right_flap_count_arr"],
        results_final["cat_imperfect_pe_count_arr"],
        # results["cat_imperfect_wt_count_arr"],
        results_final["cat_pe_indels_count_arr"],
        # results["cat_wt_indel_count_arr"],
        # results["cat_wt_count_arr"],
        # results["cat_uncategorized_count_arr"],
        ins_start,
        ins_end,
        fig_root=twinpe_8cat_results_folder,
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS
    )

    plot_total_read_counts(
        # results["bp_changes_arr"],
        results_final["total_counts"],
        results_final["edit_counts"],
        results_final["from_right_all_edit_counts"],
        results_final["from_left_all_edit_counts"],
        results_final["perfect_edit_counts"],
        # results["deletion_counts"],
        # results["insertion_counts"],
        # results["substitution_counts"],
        ins_start,
        ins_end,
        fig_root=twinpe_8cat_results_folder,
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS 
    )

    plot_edit_read_counts(
        # results["bp_changes_arr"],
        results_final["edit_counts"],
        results_final["from_right_all_edit_counts"],
        results_final["from_left_all_edit_counts"],
        results_final["perfect_edit_counts"],
        # results["deletion_counts"],
        # results["insertion_counts"],
        # results["substitution_counts"],
        ins_start,
        ins_end,
        fig_root=twinpe_8cat_results_folder,
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS
    )

    plot_edit_read_counts_with_indels(
        results_final["edit_counts"],
        # results["bp_changes_arr"],
        results_final["edit_counts_with_indels"],
        results_final["from_right_all_edit_counts"],
        results_final["from_right_all_edit_counts_with_indels"],
        results_final["from_left_all_edit_counts"],
        results_final["from_left_all_edit_counts_with_indels"],
        ins_start,
        ins_end,
        fig_root=twinpe_8cat_results_folder,
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS
    )

    plot_editing_summary(
        results_final["full_deletion_counts"], 
        results_final["full_insertion_counts"], 
        results_final["full_substitution_counts"], 
        results_final["full_edit_counts"], 
        results_final["full_total_counts"], 
        # results["ins_start"], 
        ins_end,
        del_start,
        # results["del_end"],
        fig_root=twinpe_8cat_results_folder, 
        produce_png=produce_png, 
        category_colors=CATEGORY_COLORS
    )

    # plot_nonprogrammed_edit_counts(
    #     results["full_deletion_counts"],
    #     results["full_insertion_counts"],
    #     results["full_substitution_counts"],
    #     ins_start=results["ins_start"],
    #     ins_end=results["ins_end"], 
    #     del_start=results["del_start"],
    #     del_end=results["del_end"],
    #     fig_root=twinpe_8cat_results_folder,
    #     produce_png=args.produce_png
    # )


def get_ref_base_changes(comp_ref_seq, wt_aln_seq, twin_aln_seq):
    # works for both compound variants
    bp_changes_arr = []
    del_start = None
    del_end = None
    ins_start = None
    ins_end = None
    for idx in range(len(comp_ref_seq)):
        wt_base = wt_aln_seq[idx]
        twin_base = twin_aln_seq[idx]
        if wt_base != '-' and twin_base == '-':
            bp_changes_arr.append((idx, wt_base, twin_base))
            if del_start is None:
                del_start = idx
            del_end = idx
        elif wt_base == '-' and twin_base != '-':
            bp_changes_arr.append((idx, wt_base, twin_base))
            if ins_start is None:
                ins_start = idx
            ins_end = idx

    # Debug check
    del_region_len = del_end - del_start + 1
    ins_region_len = ins_end - ins_start + 1
    if len(bp_changes_arr) != del_region_len + ins_region_len:
        raise Exception('Number of total base changes does not equal sum of deletions and insertions.')
    
    return bp_changes_arr, del_start, del_end, ins_start, ins_end, ins_region_len


def get_read_match_array(bp_changes_arr, read_map, del_start, del_end, ins_start, ins_end):
    match_arr = []
    # for ind, (comp_ind, wt_base, twin_base) in enumerate(bp_changes_arr):
    for (comp_ind, wt_base, twin_base) in bp_changes_arr:
        read_base = read_map.get(comp_ind, "")
        if read_base == wt_base:
            match_arr.append("W") # matches WT base
        elif read_base == twin_base:
            match_arr.append("T") # matches TwinPE base
        # this will never trigger bc every position in bp_changes_arr has a '-' in either wt or twin
        elif read_base == "-":
            match_arr.append("D") # non-programmed deletion
        elif len(read_base) > 1:
            match_arr.append("I") # non-programmed insertion
        elif read_base in ["A", "C", "G", "T"]:
            match_arr.append("S") # non-programmed substitution
        else:
            match_arr.append("N") # ambiguous base
        # match_arr_to_comp_base_index[ind] = comp_ind

    if del_start < ins_start:
        del_match_arr = match_arr[:del_end-del_start+1]
        ins_match_arr = match_arr[del_end-del_start+1:]
    else:
        del_match_arr = match_arr[ins_end-ins_start+1:]
        ins_match_arr = match_arr[:ins_end-ins_start+1]

    return match_arr, ins_match_arr, del_match_arr


def indel_checking(all_insertion_left_positions, all_deletion_positions, del_start, del_end, ins_start, ins_end, ignore_extraspacer_deletions, pegRNA_intervals):
    # Updated Indel checking to include flag 
    has_any_ins_byproduct = False
    has_del_in_spacer_window = False
    has_any_del_byproduct = False

    # Check for insertions anywhere in read
    if all_insertion_left_positions != "[]":
        has_any_ins_byproduct = True
    else:
        if del_start < ins_start:
            edit_range = range(del_start, ins_end + 1)
        else:
            edit_range = range(ins_start, del_end + 1)
        all_del_pos = [int(x) for x in all_deletion_positions.strip("[]").split(",") if x.strip().isdigit()]
        # Ignore deletions beyond spacers if flagged (spacer_info['spacer_a_num_bases_removed'])
        if ignore_extraspacer_deletions:
            for del_ind in all_del_pos:
                if del_ind >= pegRNA_intervals[0][0] and del_ind <= pegRNA_intervals[1][1] and del_ind not in edit_range:  # allele table stores exact position of '-' for deletions. Uses ungapped ref coords which is same coord system as pegRNA_intervals.
                    has_del_in_spacer_window = True
                    break
        # Check for deletions anywhere outside of the edit region if not flagged (del_ind to comp_ref index mapping not necessary)
        else:
            for del_ind in all_del_pos:
                if del_ind not in edit_range:
                    has_any_del_byproduct = True
                    break

    # Set has_indel based on flag
    if ignore_extraspacer_deletions:
        has_indel = has_any_ins_byproduct or has_del_in_spacer_window
    else:
        has_indel = has_any_ins_byproduct or has_any_del_byproduct

    return has_indel


def categorize_alleles(crispresso_output_folder=None, comp_ref_seq=None, twin_aln_seq=None, wt_aln_seq=None,  num_changes_to_check=2, ignore_extraspacer_deletions=False):
    # Load and validate CRISPResso2 outputs
    crispresso_info_file = os.path.join(crispresso_output_folder, "CRISPResso2_info.json")
    if not os.path.exists(crispresso_info_file):
        sys.exit(f"CRISPResso2 output missing: {crispresso_info_file}")
    try:
        crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_output_folder)
    except Exception as e:
        sys.exit(f"Could not open CRISPResso2 info file: {e}")
        
    ref_names = crispresso2_info["results"].get("ref_names", [])
    if len(ref_names) != 1 or ref_names[0] != "Compound":
        sys.exit(f"CRISPResso2 was not run against the 'Compound' reference only - did CRISPResso2 run complete successfully?\nFound reference: {ref_names}")

    # Load all necessary info once
    ref = crispresso2_info["results"]["refs"]["Compound"]
    pegRNA_intervals = ref["sgRNA_intervals"]
    # pegRNA_cut_points = ref["sgRNA_cut_points"]
    # pegRNA_plot_cut_points = ref["sgRNA_plot_cut_points"]
    # pegRNA_mismatches = ref["sgRNA_mismatches"]
    # pegRNA_names = ref["sgRNA_names"]

    # TODO: need to add supoort to get_ref_base_changes() for determining ins del region in A vs B
    # Get insertion/deletion regions and base changes
    bp_changes_arr, del_start, del_end, ins_start, ins_end, ins_region_len = get_ref_base_changes(comp_ref_seq, wt_aln_seq, twin_aln_seq)

    # del_region_coords = (del_start, del_end)
    # ins_region_coords = (ins_start, ins_end)

    # Map Compound ref bases to WT and TwinPE ref bases
    # wt_map = get_refpos_values(comp_ref_seq, wt_aln_seq)
    # pe_map = get_refpos_values(comp_ref_seq, twin_aln_seq)

    z = zipfile.ZipFile(
        os.path.join(
            crispresso_output_folder,
            crispresso2_info["running_info"]["allele_frequency_table_zip_filename"],
        )
    )
    zf = z.open(crispresso2_info["running_info"]["allele_frequency_table_filename"])
    df_alleles = pd.read_csv(zf, sep="\t")

    for idx, allele in df_alleles.iterrows():
        comp_aln_seq_read = allele.Reference_Sequence
        read_aln_seq = allele.Aligned_Sequence

        # Map this reads aligned bases to Compound ref aligned bases 
        read_map = get_refpos_values(comp_aln_seq_read, read_aln_seq)

        # Map this alleles match_arr index to ref_seq index
        # match_arr_to_comp_base_index = {}

        # Determine reads base matching for insertion/deletion region
        match_arr, ins_match_arr, del_match_arr = get_read_match_array(bp_changes_arr, read_map, del_start, del_end, ins_start, ins_end)

        # check match_arr for indels - repetitive if checking allele table
        # has_indel_in_match_arr = any(m in ["I", "D"] for m in match_arr)
        # if has_indel_in_match_arr:
        #     has_indel_in_match_arr_count += 1

        has_indel = indel_checking(allele.all_insertion_left_positions, allele.all_deletion_positions, del_start, del_end, ins_start, ins_end, ignore_extraspacer_deletions, pegRNA_intervals)

        total_TE_count = match_arr.count("T")
        has_all_TE = (total_TE_count == len(match_arr))

        # Changed to require at least num_changes_to_check
        # has_any_TE = (total_TE_count > 0)
        has_any_TE = (total_TE_count >= num_changes_to_check)
        
        # changed to require at least num_changes_to_check
        # without has_any_non_WT_in_insertion its allowing reads with non-TE matches in insertion region to be categorized as Imperfect_WT - with it causes more issues
        # has_any_TE_in_insertion = (ins_match_arr.count("T") > 0)
        # has_any_TE_in_insertion = (ins_match_arr.count("T") > 1) # require at least 2 TE matches to avoid misclassifying WT-like reads as PE-like.
        has_any_TE_in_insertion = (ins_match_arr.count("T") >= num_changes_to_check)
        # has_any_WT_in_insertion = (ins_match_arr.count("W") > 0)

        total_WT_count = match_arr.count("W")
        has_all_WT = (total_WT_count == len(match_arr))
        has_any_WT = (total_WT_count > 0)

        # ins_TE_count = ins_match_arr.count("T")
        # is_imperfect_PE = (ins_TE_count >= num_changes_to_check and ins_TE_count < len(ins_match_arr))

        # define left/right flaps
        is_left_flap = False
        is_right_flap = False
        if len(ins_match_arr) >= num_changes_to_check:
            is_left_flap = all(m == "T" for m in ins_match_arr[:num_changes_to_check])
            is_right_flap = all(m == "T" for m in ins_match_arr[-num_changes_to_check:])
            # is_left_flap = all(m == "T" for m in match_arr[del_end-del_start+1:del_end-del_start+1+num_changes_to_check])
            # is_right_flap = all(m == "T" for m in match_arr[-num_changes_to_check:])

        if has_all_TE and has_indel:
            df_alleles.at[idx,'Category'] = "PE_Indel"
        elif has_all_TE:
            df_alleles.at[idx,'Category'] = "Perfect_PE"
        elif has_all_WT and has_indel:
            df_alleles.at[idx,'Category'] = "WT_Indel"
        elif has_all_WT:
            df_alleles.at[idx,'Category'] = "WT"
        elif is_left_flap and not is_right_flap:
            df_alleles.at[idx,'Category'] = "Left_Flap"
        elif is_right_flap and not is_left_flap:
            df_alleles.at[idx,'Category'] = "Right_Flap"
        elif has_any_WT and not has_any_TE_in_insertion:
            df_alleles.at[idx,'Category'] = "Imperfect_WT"
        elif has_any_TE:
            df_alleles.at[idx,'Category'] = "Imperfect_PE"
        else:
            df_alleles.at[idx,'Category'] = "Uncategorized"

    return {
        "df_alleles": df_alleles, 
        "pegRNA_cut_points": ref["sgRNA_cut_points"],
        "pegRNA_plot_cut_points": ref["sgRNA_plot_cut_points"],
        "pegRNA_intervals": ref["sgRNA_intervals"],
        "pegRNA_mismatches": ref["sgRNA_mismatches"],
        "pegRNA_names": ref["sgRNA_names"], 
        "bp_changes_arr": bp_changes_arr, 
        "del_start": del_start, 
        "del_end": del_end, 
        "ins_start": ins_start, 
        "ins_end": ins_end, 
        "ins_region_len": ins_region_len
    }

def analyze_collapsed_categorized_alleles(
    # crispresso_output_folder, 
    df_alleles_final, 
    twinpe_8cat_results_folder, 
    # wt_seq=None, 
    # twin_seq=None, 
    # wt_aln_seq_a=None, 
    # twin_aln_seq_a=None, 
    # comp_ref_seq_a=None, 
    # wt_aln_seq_b=None, 
    # twin_aln_seq_b=None, 
    # comp_ref_seq_b=None, 
    bp_changes_arr_a=None, 
    del_start_a=None, 
    del_end_a=None, 
    ins_start_a=None, 
    ins_end_a=None, 
    ins_region_len_a=None, 
    pegRNA_intervals_a=None, 
    bp_changes_arr_b=None, 
    del_start_b=None, 
    del_end_b=None, 
    ins_start_b=None, 
    ins_end_b=None, 
    ins_region_len_b=None, 
    pegRNA_intervals_b=None, 
    # num_changes_to_check=2, 
    ignore_extraspacer_deletions=False, 
    # produce_png=False
):
    if ins_region_len_a != ins_region_len_b:
        raise ValueError("Insertion region lengths for A and B do not match.")
    
    # Only plot ins region
    edit_counts = [0] * ins_region_len_a
    from_left_all_edit_counts = [0] * ins_region_len_a
    from_right_all_edit_counts = [0] * ins_region_len_a
    cat_perfect_pe_count = 0
    cat_perfect_pe_count_arr = [0] * ins_region_len_a
    cat_pe_indels_count = 0
    cat_pe_indels_count_arr = [0] * ins_region_len_a
    cat_left_flap_count = 0
    cat_left_flap_count_arr = [0] * ins_region_len_a
    cat_right_flap_count = 0
    cat_right_flap_count_arr = [0] * ins_region_len_a
    cat_imperfect_pe_count = 0
    cat_imperfect_pe_count_arr = [0] * ins_region_len_a
    cat_imperfect_wt_count = 0
    # cat_imperfect_wt_count_arr = [0] * ins_region_len_a
    cat_wt_indel_count = 0
    # cat_wt_indel_count_arr = [0] * ins_region_len_a
    cat_wt_count = 0
    # cat_wt_count_arr = [0] * ins_region_len_a
    cat_uncategorized_count = 0
    # cat_uncategorized_count_arr = [0] * ins_region_len_a
    from_left_all_edit_counts_with_indels = [0] * ins_region_len_a
    # from_left_all_edit_counts_no_indels = [0] * ins_region_len_a
    from_right_all_edit_counts_with_indels = [0] * ins_region_len_a
    # from_right_all_edit_counts_no_indels = [0] * ins_region_len_a
    edit_counts_with_indels = [0] * ins_region_len_a
    # edit_counts_no_indels = [0] * ins_region_len_a
    # deletion_counts = [0] * ins_region_len # check if fix needed
    # insertion_counts = [0] * ins_region_len_a
    # substitution_counts = [0] * ins_region_len_a
    allele_counts = {} # e.g. TTRTT > 100
    # allele_categories = {} # e.g. TTRTT > Left_flap
    perfect_T_count = 0

    if len(bp_changes_arr_a) != len(bp_changes_arr_b):
        raise ValueError("Base changes array lengths for A and B do not match.")

    # Plot full ins del region
    full_edit_counts = [0] * len(bp_changes_arr_a)
    full_deletion_counts = [0] * len(bp_changes_arr_a)
    full_insertion_counts = [0] * len(bp_changes_arr_a)
    full_substitution_counts = [0] * len(bp_changes_arr_a)
    # full_edit_counts_with_indels = [0] * len(bp_changes_arr_a)
    # full_edit_counts_no_indels = [0] * len(bp_changes_arr_a)
    # full_perfect_T_count = 0

    total_alleles = 0
    total_alleles_reads = 0
    # total_alleles_deletions_reads = 0

    # has_indel_in_match_arr_count = 0
    # has_any_indel_byproduct_count = 0

    # iterate all alleles in input allele table
    for idx, allele in df_alleles_final.iterrows():
        comp_aln_seq_read = allele.Reference_Sequence
        read_aln_seq = allele.Aligned_Sequence

        total_alleles += 1
        total_alleles_reads += allele["#Reads"]

        # Map this reads aligned bases to Compound ref aligned bases 
        read_map = get_refpos_values(comp_aln_seq_read, read_aln_seq)

        if allele['allele_source'] == 'A':
            # Determine reads base matching for insertion/deletion region
            match_arr, ins_match_arr, del_match_arr = get_read_match_array(bp_changes_arr_a, read_map, del_start_a, del_end_a, ins_start_a, ins_end_a)
            has_indel = indel_checking(allele.all_insertion_left_positions, allele.all_deletion_positions, del_start_a, del_end_a, ins_start_a, ins_end_a, ignore_extraspacer_deletions, pegRNA_intervals_a)
        else:
            match_arr, ins_match_arr, del_match_arr = get_read_match_array(bp_changes_arr_b, read_map, del_start_b, del_end_b, ins_start_b, ins_end_b)
            has_indel = indel_checking(allele.all_insertion_left_positions, allele.all_deletion_positions, del_start_b, del_end_b, ins_start_b, ins_end_b, ignore_extraspacer_deletions, pegRNA_intervals_b)

        # TODO: Unused - use for indel tri-barplot?
        # if has_indel:
        #     total_alleles_deletions_reads += allele['#Reads']

        # Use compound variant A for both A and B alleles
        match_arr_str = "\t".join(del_match_arr) + "\t".join(ins_match_arr) + "\t" + str(has_indel)

        if match_arr_str not in allele_counts:
            allele_counts[match_arr_str] = 0
        allele_counts[match_arr_str] += allele['#Reads']

        if allele['Category'] == "PE_Indel":
            cat_pe_indels_count += allele['#Reads']
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                # Redundant check
                # if match == 'T':
                cat_pe_indels_count_arr[pos_idx] += allele['#Reads']
        
        elif allele['Category'] == "Perfect_PE":
            cat_perfect_pe_count += allele['#Reads']
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                # Redundant check
                # if match == 'T':
                cat_perfect_pe_count_arr[pos_idx] += allele['#Reads']

        elif allele['Category'] == "WT_Indel":
            cat_wt_indel_count += allele['#Reads']
            # This can't be true by definition
            # for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
            #     if match == 'T':
            #         cat_wt_indel_count_arr[pos_idx] += allele['#Reads']

        elif allele['Category'] == "WT":
            cat_wt_count += allele['#Reads']
            # this can't be true by definition
            # for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
            #     if match == 'T':
            #         cat_wt_count_arr[pos_idx] += allele['#Reads']

        elif allele['Category'] == "Left_Flap":
            cat_left_flap_count += allele['#Reads']
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match == 'T':
                    cat_left_flap_count_arr[pos_idx] += allele['#Reads']

        elif allele['Category'] == "Right_Flap":
            cat_right_flap_count += allele['#Reads']
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match == 'T':
                    cat_right_flap_count_arr[pos_idx] += allele['#Reads']

        elif allele['Category'] == "Imperfect_WT":
            cat_imperfect_wt_count += allele['#Reads']
            # The cat_imperfect_wt_count_arr is currently not used for plotting
            # for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
            #     if match == 'T':
            #         cat_imperfect_wt_count_arr[pos_idx] += allele['#Reads']

        elif allele['Category'] == "Imperfect_PE":
            cat_imperfect_pe_count += allele['#Reads']
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match == 'T':
                    cat_imperfect_pe_count_arr[pos_idx] += allele['#Reads']

        else:
            cat_uncategorized_count += allele['#Reads']
            # This cat_uncategorized_count_arr is currently not used for plotting
            # for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
            #     if match == 'T':
            #         cat_uncategorized_count_arr[pos_idx] += allele['#Reads']

        # TODO: Removed this due to dual compound conflicts - Is this still necessary?
        # if (
        #     match_arr_str in allele_categories
        #     and allele_categories[match_arr_str] != allele['Category']
        # ):
        #     raise Exception(
        #         f"Conflicting categories for {match_arr_str} {allele['Category']} vs {allele_categories[match_arr_str]}"
        #     )
        # allele_categories[match_arr_str] = allele['Category']

        for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
            if match == 'T':
                from_left_all_edit_counts[pos_idx] += allele['#Reads']

                if has_indel:
                    from_left_all_edit_counts_with_indels[pos_idx] += allele['#Reads']
                # else:
                #     from_left_all_edit_counts_no_indels[pos_idx] += allele['#Reads']

            else:
                break

        for pos_idx, match in zip(reversed(range(len(ins_match_arr))), reversed(ins_match_arr)):
            if match == 'T':
                from_right_all_edit_counts[pos_idx] += allele['#Reads']

                if has_indel:
                    from_right_all_edit_counts_with_indels[pos_idx] += allele['#Reads']
                # else:
                #     from_right_all_edit_counts_no_indels[pos_idx] += allele['#Reads']
            else:
                break

        for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
            if match == 'T':
                edit_counts[pos_idx] += allele['#Reads']
                if has_indel:
                    edit_counts_with_indels[pos_idx] += allele['#Reads']
            #     else:
            #         edit_counts_no_indels[pos_idx] += allele['#Reads']
            # if match == 'D':
            #     deletion_counts[pos_idx] += allele['#Reads']
            # if match == "I":
            #     insertion_counts[pos_idx] += allele["#Reads"]
            # if match == "S":
            #     substitution_counts[pos_idx] += allele["#Reads"]
        if ins_match_arr == ['T']*len(ins_match_arr):
            perfect_T_count += allele['#Reads']

        for pos_idx, match in zip(range(len(bp_changes_arr_a)), del_match_arr+ins_match_arr): # match_arr):
            if match == 'T':
                full_edit_counts[pos_idx] += allele['#Reads']
                # if has_indel:
                #     full_edit_counts_with_indels[pos_idx] += allele['#Reads']
                # else:
                #     full_edit_counts_no_indels[pos_idx] += allele['#Reads']
            if match == 'D':
                full_deletion_counts[pos_idx] += allele['#Reads']
            if match == "I":
                full_insertion_counts[pos_idx] += allele["#Reads"]
            if match == "S":
                full_substitution_counts[pos_idx] += allele["#Reads"]
        # if match_arr == ['T']*len(match_arr):
        #     full_perfect_T_count += allele['#Reads']

    folder_category_counts = {
        "WT": cat_wt_count,
        "WT Indel": cat_wt_indel_count,
        "Imperfect WT": cat_imperfect_wt_count,
        "Left Flap": cat_left_flap_count,
        "Right Flap": cat_right_flap_count,
        "Perfect PE": cat_perfect_pe_count,
        "Imperfect PE": cat_imperfect_pe_count,
        "PE Indel": cat_pe_indels_count,
        "Uncategorized": cat_uncategorized_count,
    }

    with open(
        os.path.join(twinpe_8cat_results_folder, "c4.top_alleles_by_category.txt"), "w"
    ) as fout:
        for c in df_alleles_final['Category'].unique():
            fc = c.replace("_", " ")            
            fout.write(f"Category: {c}, Total reads: {folder_category_counts[fc]}\n")
            for i, row in df_alleles_final[df_alleles_final['Category'] == c].sort_values(by='#Reads', ascending=False).head(50).iterrows():
                fout.write(f"Read {i} count: {row['#Reads']}\n")
                fout.write(f"{row['Aligned_Sequence']}\n")
                fout.write(f"{row['Reference_Sequence']}\n")
            fout.write("\n")

    imperfect_from_left_all_edit_counts = [x-perfect_T_count for x in from_left_all_edit_counts] # these are reads that start with a twin-edited, but are not completely twin-edited. Note that from_left_all_edit_counts includes perfect matches
    imperfect_from_right_all_edit_counts = [x-perfect_T_count for x in from_right_all_edit_counts]
    # number of perfect_PE alleles at each position - now excludes those with indels
    # should this take indels into account? - what exactly are these two plots that use this trying to show?
    perfect_edit_counts = [perfect_T_count] * ins_region_len_a
    # full_perfect_edit_counts = [full_perfect_T_count] * len(bp_changes_arr_a)
    total_counts = [total_alleles_reads] * ins_region_len_a
    full_total_counts = [total_alleles_reads] * len(bp_changes_arr_a)

    with open(twinpe_8cat_results_folder + "/c3.allele_counts.txt", "w") as fout:
        sorted_allele_counts = sorted(
            allele_counts.keys(), key=lambda x: allele_counts[x], reverse=True
        )
        fout.write("\t".join([str(x) for x in bp_changes_arr_a]) + "\thas_indel\tcount\n")
        for allele_str in sorted_allele_counts:
            fout.write(allele_str + "\t" + str(allele_counts[allele_str]) + "\n")

    with open(twinpe_8cat_results_folder + "/c2.arrays.txt", "w") as fout:
        fout.write("Class\t" + "\t".join([str(x) for x in bp_changes_arr_a]) + "\n")
        fout.write("total_counts\t" + "\t".join([str(x) for x in total_counts]) + "\n")
        fout.write("edit_counts\t" + "\t".join([str(x) for x in edit_counts]) + "\n")
        fout.write(
            "perfect_edit_counts\t"
            + "\t".join([str(x) for x in perfect_edit_counts])
            + "\n"
        )
        fout.write(
            "from_left_all_edit_counts\t"
            + "\t".join([str(x) for x in from_left_all_edit_counts])
            + "\n"
        )
        fout.write(
            "from_right_all_edit_counts\t"
            + "\t".join([str(x) for x in from_right_all_edit_counts])
            + "\n"
        )
        fout.write(
            "imperfect_from_left_all_edit_counts\t"
            + "\t".join([str(x) for x in imperfect_from_left_all_edit_counts])
            + "\n"
        )
        fout.write(
            "imperfect_from_right_all_edit_counts\t"
            + "\t".join([str(x) for x in imperfect_from_right_all_edit_counts])
            + "\n"
        )
        # fout.write(
        #     "deletion_counts\t" + "\t".join([str(x) for x in deletion_counts]) + "\n"
        # )
        # fout.write(
        #     "insertion_counts\t" + "\t".join([str(x) for x in insertion_counts]) + "\n"
        # )
        # fout.write(
        #     "substitution_counts\t" + "\t".join([str(x) for x in substitution_counts]) + "\n"
        # )

    with open(twinpe_8cat_results_folder + "/c1.counts.txt", "w") as fout:
        fout.write(
            "\t".join(
                [
                    "Perfect_PE",
                    "PE_indels",
                    "Imperfect_PE",
                    "Left_flap",
                    "Right_flap",
                    "Imperfect_WT",
                    "WT_indels",
                    "WT",
                    "Uncategorized",
                ]
            )
            + "\n"
        )
        fout.write(
            "\t".join(
                [
                    str(x)
                    for x in [
                        cat_perfect_pe_count,
                        cat_pe_indels_count,
                        cat_imperfect_pe_count,
                        cat_left_flap_count,
                        cat_right_flap_count,
                        cat_imperfect_wt_count,
                        cat_wt_indel_count,
                        cat_wt_count,
                        cat_uncategorized_count,
                    ]
                ]
            )
            + "\n"
        )
    # No longer set through this function
    # min_frequency = crispresso2_info["running_info"][
    #     "args"
    # ].min_frequency_alleles_around_cut_to_plot
    # max_n_rows = crispresso2_info["running_info"][
    #     "args"
    # ].max_rows_alleles_around_cut_to_plot

    df_alleles_final.to_csv(
        os.path.join(twinpe_8cat_results_folder, "c8.detailed_allele_table_with_categories.csv"), 
                     index=False
    )

    return {
        # input to multiple functions
        "edit_counts": edit_counts, 
        "from_left_all_edit_counts": from_left_all_edit_counts, 
        "from_right_all_edit_counts": from_right_all_edit_counts, 
        "perfect_edit_counts": perfect_edit_counts,
        # inputs to plot_successful_twin_edit_counts_by_category only
        "folder_category_counts": folder_category_counts, 
        "cat_perfect_pe_count_arr": cat_perfect_pe_count_arr,
        "cat_left_flap_count_arr": cat_left_flap_count_arr,
        "cat_right_flap_count_arr": cat_right_flap_count_arr, 
        "cat_imperfect_pe_count_arr": cat_imperfect_pe_count_arr,
        "cat_pe_indels_count_arr": cat_pe_indels_count_arr,
        # inputs to plot_total_read_counts only 
        "total_counts": total_counts,
        # inputs to plot_edit_read_counts only
        # inputs to plot_edit_read_counts_with_indels only 
        "edit_counts_with_indels":  edit_counts_with_indels,
        "from_right_all_edit_counts_with_indels": from_right_all_edit_counts_with_indels,
        "from_left_all_edit_counts_with_indels": from_left_all_edit_counts_with_indels, 
        # inputs to plot_editing_summary only
        "full_deletion_counts": full_deletion_counts,
        "full_insertion_counts": full_insertion_counts,
        "full_substitution_counts": full_substitution_counts, 
        "full_edit_counts": full_edit_counts, 
        "full_total_counts": full_total_counts

        # Not used
        # "df_alleles": df_alleles, 
        # "bp_changes_arr": bp_changes_arr, 
        # "cat_imperfect_wt_count_arr": cat_imperfect_wt_count_arr, 
        # "cat_wt_indel_count_arr": cat_wt_indel_count_arr,
        # "cat_wt_count_arr": cat_wt_count_arr,
        # "cat_uncategorized_count_arr": cat_uncategorized_count_arr,
        # "deletion_counts": deletion_counts,
        # "insertion_counts": insertion_counts,
        # "substitution_counts": substitution_counts, 
        # "del_start": del_start,
        # "del_end": del_end, 
        # "ins_start": ins_start, 
        # "ins_end": ins_end, 
        # "wt_seq": wt_seq,
        # "wt_aln_seq": wt_aln_seq,
        # "twin_aln_seq": twin_aln_seq, 
        # "min_frequency": min_frequency,
        # "max_n_rows": max_n_rows, 
        # "pegRNA_cut_points": pegRNA_cut_points,
        # "pegRNA_plot_cut_points": pegRNA_plot_cut_points,
        # "pegRNA_intervals": pegRNA_intervals,
        # "pegRNA_mismatches": pegRNA_mismatches,
        # "pegRNA_names": pegRNA_names, 
    }


def get_allele_df_keys(df):
    # Create sequence_key and take the lexicographically smaller of the forward and reverse complement
    df = df.copy()
    df['sequence_key_fw'] = df['Aligned_Sequence'].str.replace('-', '', regex=False)
    df['sequence_key_rc'] = df['sequence_key_fw'].apply(reverse_complement)
    df['sequence_key'] = df[['sequence_key_fw', 'sequence_key_rc']].min(axis=1)

    return df


def collapse_allele_categories(df_alleles_a, df_alleles_b):
    # Prepare sequence_key for both dataframes
    a = get_allele_df_keys(df_alleles_a)
    b = get_allele_df_keys(df_alleles_b)

    # Check for duplicates in sequence_key within each dataframe
    dupes_a = a.loc[a['sequence_key'].duplicated(), 'sequence_key']
    if len(dupes_a):
        raise ValueError(f"Duplicate sequence_key in df_alleles_a: {dupes_a.unique()}")

    dupes_b = b.loc[b['sequence_key'].duplicated(), 'sequence_key']
    if len(dupes_b):
        raise ValueError(f"Duplicate sequence_key in df_alleles_b: {dupes_b.unique()}")

    # Merge only the needed columns from B
    b_subset = b[['sequence_key', 'Category', 'Aligned_Reference_Scores']]
    merged = a.merge(b_subset, on='sequence_key', how='left', suffixes=('', '_B'), validate='one_to_one')

    # Convert reference scores to numeric for comparison
    ref_a = pd.to_numeric(merged['Aligned_Reference_Scores'], errors='raise')
    ref_b = pd.to_numeric(merged['Aligned_Reference_Scores_B'], errors='raise')

    # Decide whether to keep original A category or use B category
    use_b = (merged['Category'] != merged['Category_B']) & (ref_b > ref_a)

    # Update Category if B wins
    merged['allele_source'] = np.where(use_b, 'B', 'A')
    merged['Category'] = np.where(use_b, merged['Category_B'], merged['Category'])

    # Drop helper columns from B
    merged = merged.drop(columns=['sequence_key_fw', 'sequence_key_rc', 'sequence_key', 'Category_B', 'Aligned_Reference_Scores_B'])

    return merged.reset_index(drop=True)


def setBarMatplotlibDefaults():
    matplotlib.use("AGG")
    matplotlib.rcParams["font.sans-serif"] = [
        "Arial",
        "Liberation Sans",
        "Bitstream Vera Sans",
    ]
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["axes.facecolor"] = "white"
    plt.ioff()


def plot_reads_input_summary_barplot(crispresso_output_folder, fig_root=None, produce_png=False):

    crispresso_mapping_statistics_file = os.path.join(crispresso_output_folder, 'CRISPResso_mapping_statistics.txt')
    read_stats = pd.read_csv(crispresso_mapping_statistics_file, sep="\t")
    counts = [read_stats['READS IN INPUTS'][0], read_stats['READS AFTER PREPROCESSING'][0], read_stats['READS ALIGNED'][0]]
    labels = ["Reads Input", "Reads After Preprocessing", "Reads Aligned"]
    total = read_stats['READS IN INPUTS'][0]

    sorted_pairs = sorted(zip(labels, counts), key=lambda x: x[1], reverse=True)
    sorted_labels = [p[0] for p in sorted_pairs]
    sorted_values = [p[1] for p in sorted_pairs]

    percent_labels = [f"{lab}\n({val/total*100:.1f}%)" for lab, val in sorted_pairs]
    # percent_labels = [f"{lab}\n({val:}%)" for lab, val in sorted_pairs]
    # percent_labels = [f"{lab}\n({val:,})" if val == total else f"{lab}\n({val:,} | {val/total*100:.1f}%)" for lab, val in sorted_pairs]

    # use default mpl color cycle
    # colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    # bar_colors = [colors[i % len(colors)] for i in range(len(sorted_values))]

    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)

    # plot bar chart
    bars = ax.bar(sorted_labels, sorted_values, 
                color="lightgrey", edgecolor="white", linewidth=0.6)

    # add percentage labels
    for bar in bars:
        height = bar.get_height()
        # percentage = (height / total) * 100
        # percetange_label = f"{percentage:.1f}%"
            
        ax.text(bar.get_x() + bar.get_width()/2.,
                height,
                f"{height:,}", 
                # percetange_label,
                ha='center',
                va='bottom',
                fontsize=10,
                color='black'
            )
        
    ax.set_ylabel("Counts", fontsize=10)

    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(percent_labels, rotation=0, ha="center")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.minorticks_on()
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    ax.tick_params(axis="x", which="minor", bottom=False)

    plt.tight_layout()

    plt.savefig(fig_root + "/a1.Reads_input_summary.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a1.Reads_input_summary.png", bbox_inches='tight')


def plot_category_stacked_summary_barplot(crispresso_output_folder, counts_dict, fig_root=None, produce_png=False, category_colors=None):

    aligned_labels = list(counts_dict.keys())
    aligned_counts = list(counts_dict.values())
    crispresso_mapping_statistics_file = os.path.join(crispresso_output_folder, 'CRISPResso_mapping_statistics.txt')
    read_stats = pd.read_csv(crispresso_mapping_statistics_file, sep="\t")
    unaligned_count = read_stats['READS AFTER PREPROCESSING'][0] - read_stats['READS ALIGNED'][0]

    sorted_counts_labels = sorted(zip(aligned_labels, aligned_counts), key=lambda x: x[1], reverse=True)
    total = sum(aligned_counts)
    legend_labels = [f"{lab} ({val/total*100:.1f}%)" for lab, val in sorted_counts_labels]

    fig, axes = plt.subplots(1, 2, figsize=(3, 6), dpi=300, sharey=True, gridspec_kw={'wspace': 0}) # dpi=300

    # plot stacked bar
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

    # plot second bar
    bar2 = axes[1].bar([0], unaligned_count, color='lightgrey')
    # axes[1].bar([0], [unaligned_count], color='orange', edgecolor='white', linewidth=0.3)

    axes[1].text([0][0], unaligned_count, f"{unaligned_count:,}", ha='center', va='bottom', fontsize=8)

    handles, labels = axes[0].get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    handles.append(bar2[0])
    # labels.append(f"Unaligned Reads")
    # labels.append(f"Unaligned ({unaligned_count/ (total) * 100:.1f}%)")
    axes[0].legend(
        handles, labels,
        bbox_to_anchor=(2.03, 0.5),
        loc="center left",
        borderaxespad=0.25,
        fontsize=8
    )
    # axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axes[0].set_xticks([])
    axes[0].set_xlabel("Aligned Reads", fontsize=8)
    axes[0].set_ylabel("Counts", fontsize=8)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].minorticks_on()
    axes[0].tick_params(axis="y", which="minor", length=3, width=0.8)
    axes[0].tick_params(axis="y", which="major", length=5, width=1)
    axes[0].tick_params(axis="y", labelsize=8)

    axes[1].set_xticks([])
    axes[1].set_xlabel("Unaligned Reads", fontsize=8)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[1].minorticks_on()
    axes[1].tick_params(axis="y", which="minor", length=3, width=0.8, left=False)
    axes[1].tick_params(axis="y", which="major", length=5, width=1, left=False)

    # axes[0].set_title("Aligned Reads")
    # axes[1].set_title("Unaligned Reads")

    # ax.legend(bbox_to_anchor=(1, .5), loc="center left", borderaxespad=0.25)

    # plt.tight_layout()
    plt.savefig(fig_root + "/a2.Category_summary_stacked.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a2.Category_summary_stacked.png", bbox_inches='tight')


def plot_category_summary_barplot(counts_dict, fig_root=None, produce_png=False, category_colors=None):
    labels = list(counts_dict.keys())
    counts = list(counts_dict.values())

    
    sorted_pairs = sorted(zip(labels, counts), key=lambda x: x[1], reverse=True)
    sorted_labels = [p[0] for p in sorted_pairs]
    sorted_values = [p[1] for p in sorted_pairs]

    total = sum(sorted_values)

    percent_labels = [f"{lab}\n({val/total*100:.1f}%)" for lab, val in sorted_pairs]

    # use default mpl color cycle
    # colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    # bar_colors = [colors[i % len(colors)] for i in range(len(sorted_values))]
    colors = [category_colors.get(lab, "black") for lab in sorted_labels]

    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)

    # plot bar chart
    bars = ax.bar(sorted_labels, sorted_values, 
                color=colors, edgecolor="white", linewidth=0.6)
    
    # add value labels
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.,
                height,
                f"{height:,}",
                ha='center',
                va='bottom',
                fontsize=8,
            )

    ax.set_ylabel("Counts", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)

    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(percent_labels, rotation=0, ha="center", fontsize=7)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.minorticks_on()
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    ax.tick_params(axis="x", which="minor", bottom=False)

    plt.tight_layout()

    plt.savefig(fig_root + "/a3.Category_summary.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a3.Category_summary.png", bbox_inches='tight')


# def plot_stacked_summary_barplot(counts_dict, fig_root=None, produce_png=False):
#     labels = list(counts_dict.keys())
#     counts = list(counts_dict.values())

#     sorted_counts_labels = sorted(zip(labels, counts), key=lambda x: x[1], reverse=True)
#     total = sum(counts)
#     legend_labels = [f"{lab} ({val/total*100:.1f}%)" for lab, val in sorted_counts_labels]
#     fig, ax = plt.subplots(figsize=(3.5, 6), dpi=300) # dpi=300

#     x = [0]
#     bottom = 0
#     # plot stacked bar
#     for (lab, val), legend_label in zip(sorted_counts_labels, legend_labels):
#         ax.bar(
#             x, val, 
#             bottom=bottom, 
#             label=legend_label, 
#             edgecolor='white',
#             linewidth=.3
#         )
#         bottom += val

#     ax.set_xticks([])
#     ax.set_ylabel("Counts", fontsize=12)

#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)

#     ax.minorticks_on()
#     ax.tick_params(axis="y", which="minor", length=3, width=0.8)
#     ax.tick_params(axis="y", which="major", length=5, width=1)


#     handles, labels = ax.get_legend_handles_labels()
#     ax.legend(
#         handles[::-1], labels[::-1],
#         bbox_to_anchor=(1, .5),
#         loc="center left",
#         borderaxespad=0.25
#     )
#     # ax.legend(bbox_to_anchor=(1, .5), loc="center left", borderaxespad=0.25)

#     plt.tight_layout()

#     plt.savefig(fig_root + "/Summary_barplot_stacked.pdf")
#     if produce_png:
#         plt.savefig(fig_root + "/Summary_barplot_stacked.png")


def plot_successful_twin_edit_counts_by_category(
    # bp_changes_arr,
    edit_counts,
    cat_perfect_pe_count_arr,
    cat_left_flap_count_arr,
    cat_right_flap_count_arr,
    cat_imperfect_pe_count_arr,
    # cat_imperfect_wt_count_arr,
    cat_pe_indels_count_arr,
    # cat_wt_indel_count_arr,
    # cat_wt_count_arr,
    # cat_uncategorized_count_arr,
    ins_start, 
    ins_end,
    fig_root=None,
    produce_png=False, 
    category_colors=None,
):

    # Indices for the x-axis
    indices = np.arange(len(edit_counts))

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 6))

    # This is actually the only stacked bar plot - none of the others are stacked
    # Stacked bar plot
    ax.bar(indices, cat_perfect_pe_count_arr, label="Perfect PE", color=category_colors["Perfect PE"])
    bottom_so_far = cat_perfect_pe_count_arr

    ax.bar(indices, cat_left_flap_count_arr, label="Left Flap", color=category_colors["Left Flap"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_left_flap_count_arr)]

    ax.bar(indices, cat_right_flap_count_arr, label="Right Flap", color=category_colors["Right Flap"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_right_flap_count_arr)]

    ax.bar(indices, cat_imperfect_pe_count_arr, label="Imperfect PE", color=category_colors["Imperfect PE"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_imperfect_pe_count_arr)]
    
    ax.bar(indices, cat_pe_indels_count_arr, label="PE Indel", color=category_colors["PE Indel"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_pe_indels_count_arr)]

    # By definition these cannot have any 'successful twin edits' so not including
    # ax.bar(indices, cat_imperfect_wt_count_arr, label="Imperfect WT", bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_imperfect_wt_count_arr)]

    # ax.bar(indices, cat_wt_indel_count_arr, label="WT Indels", bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_wt_indel_count_arr)]

    # ax.bar(indices, cat_wt_count_arr, label="WT", bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_wt_count_arr)]

    # ax.bar(indices, cat_uncategorized_count_arr, label="Uncategorized", bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, cat_uncategorized_count_arr)]

    # Adding labels and title
    ax.set_xlabel("Compound Reference Index")
    ax.set_ylabel("Read Counts")
    # ax.set_title("Successful Twin Edit Counts by Category")
    ax.set_title("Stacked Per-Position Twin Prime Editing Across Inserted Region by Category")

    ax.set_xticks(indices)
    ax.set_xticklabels([str(x) for x in range(ins_start, ins_end+1)], rotation=45, ha="center")
    # Dynamic? >35 x-tick = font 8 & rotate 45
    ax.tick_params(axis="x", labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    plt.tight_layout()
    
    plt.savefig(fig_root + "/a4.TPE_counts_by_category_stacked.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a4.TPE_counts_by_category_stacked.png", bbox_inches='tight')


def plot_total_read_counts(
    # bp_changes_arr,
    total_counts,
    edit_counts,
    from_right_all_edit_counts,
    from_left_all_edit_counts,
    perfect_edit_counts,
    # deletion_counts,
    # insertion_counts,
    # substitution_counts,
    ins_start, 
    ins_end,
    fig_root=None,
    produce_png=False, 
    category_colors=None
):

    # Indices for the x-axis
    indices = np.arange(len(edit_counts))

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 6))

    # total number of alleles at each position
    ax.bar(indices, total_counts, label="Total Reads", color=category_colors["PE Indel"], alpha=1.0)
    # number of alleles with "T"s at each position
    ax.bar(indices, edit_counts, label="Total TPEs", color=category_colors["Imperfect PE"], alpha=1.0)
    # number of alleles with continuous "T"s from right at each position
    ax.bar(indices, from_right_all_edit_counts, color=category_colors["Right Flap"], label="Continuous TPEs From Right", alpha=1.0)
    # number of alleles with continuous "T"s from left at each position
    ax.bar(indices, from_left_all_edit_counts, color=category_colors["Left Flap"], label="Continuous TPEs From Left", alpha=0.75)
    # this perfect_PE seems out of place for this plot?
    # number of perfect_PE alleles at each position - changed above to support indels - is this wanted?
    ax.bar(indices, perfect_edit_counts, label="Perfect TPEs", color=category_colors["Perfect PE"], alpha=1.0)
    # # number of alleles with "D"s at each position
    # ax.bar(indices, deletion_counts, label="Deletions")
    # # number of alleles with "I"s at each position
    # ax.bar(indices, insertion_counts, label="Insertions")
    # # number of alleles with "S"s at each position
    # ax.bar(indices, substitution_counts, label="Substitutions")

    # Adding labels and title
    ax.set_xlabel('Compound Reference Index')
    ax.set_ylabel('Read Counts')
    # ax.set_title("Edit Counts")
    ax.set_title("Per-Position Twin Prime Editing Across Inserted Region - All Reads")

    ax.set_xticks(indices)
    ax.set_xticklabels([str(x) for x in range(ins_start, ins_end+1)], rotation=45, ha="center")
    # Dynamic? >35 x-tick = font 8 & rotate 45
    ax.tick_params(axis="x", labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    plt.tight_layout()

    plt.savefig(fig_root + "/a5.TPE_profile_all_reads.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a5.TPE_profile_all_reads.png", bbox_inches='tight')


def plot_edit_read_counts(
    # bp_changes_arr,
    edit_counts,
    from_right_all_edit_counts,
    from_left_all_edit_counts,
    perfect_edit_counts,
    # deletion_counts,
    # insertion_counts,
    # substitution_counts,
    ins_start, 
    ins_end,
    fig_root=None,
    produce_png=False, 
    category_colors=None
):

    # Indices for the x-axis
    indices = np.arange(len(edit_counts))

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 6))

    # Stacked bar plot
    ax.bar(indices, edit_counts, label="Total TPEs", color=category_colors["Imperfect PE"], alpha=1.0)
    ax.bar(indices, from_right_all_edit_counts, label="Continuous TPEs From Right", color=category_colors["Right Flap"], alpha=1.0)
    ax.bar(indices, from_left_all_edit_counts, label="Continuous TPEs From Left", color=category_colors["Left Flap"], alpha=0.75)
    ax.bar(indices, perfect_edit_counts, label="Perfect TPEs", color=category_colors["Perfect PE"], alpha=1.0)
    # ax.bar(indices, deletion_counts, label="Deletions")
    # ax.bar(indices, insertion_counts, label="Insertions")
    # ax.bar(indices, substitution_counts, label="Substitutions")

    # Adding labels and title
    ax.set_xlabel("Compound Reference Index")
    ax.set_ylabel("Read Counts")
    # ax.set_title("Successful Twin Edit Counts")
    ax.set_title("Per-Position Twin Prime Editing Across Inserted Region - Edited Reads")

    ax.set_xticks(indices)
    ax.set_xticklabels([str(x) for x in range(ins_start, ins_end+1)], rotation=45, ha="center")
    # Dynamic? >35 x-tick = font 8 & rotate 45
    ax.tick_params(axis="x", labelsize=8)
    ax.legend(loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig(fig_root + "/a6.TPE_profile_edited_reads.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a6.TPE_profile_edited_reads.png", bbox_inches='tight')


def plot_edit_read_counts_with_indels(
    edit_counts,
    # bp_changes_arr,
    edit_counts_with_indels,
    from_right_all_edit_counts,
    from_right_all_edit_counts_with_indels,
    from_left_all_edit_counts,
    from_left_all_edit_counts_with_indels,
    ins_start, 
    ins_end,
    fig_root=None,
    produce_png=False, 
    category_colors=None
):

    # Indices for the x-axis
    indices = np.arange(len(edit_counts))

    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 18))

    # Stacked bar plot
    # number of "T"s at each position
    ax1.bar(indices, edit_counts, label="Total TPEs", color=category_colors["Perfect PE"])
    # number of "T"s at each position for match_arr's that have indels
    ax1.bar(
        indices, edit_counts_with_indels, label="TPEs with Indels", color=category_colors["WT"]
    )

    # Adding labels and title
    ax1.set_xlabel('Compound Reference Index')
    ax1.set_ylabel('Read Counts')
    ax1.set_title('Per-Position Twin Prime Editing With and Without Indels Across Inserted Region')

    ax1.set_xticks(indices)
    ax1.set_xticklabels([str(x) for x in range(ins_start, ins_end+1)], rotation=45, ha="center")
    ax1.tick_params(axis="x", labelsize=8)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.legend()

    # Stacked bar plot
    # number of continuous "T"s from left at each position
    ax2.bar(indices, from_left_all_edit_counts, label="Continuous TPEs From Left", color=category_colors["Left Flap"])
    # number of continuous "T"s from left at each position for match_arr's that have indels
    ax2.bar(
        indices,
        from_left_all_edit_counts_with_indels,
        label="Continuous TPEs From Left with Indels", 
        color=category_colors["WT"]
    )

    # Adding labels and title
    ax2.set_xlabel('Compound Reference Index')
    ax2.set_ylabel('Read Counts')
    # ax2.set_title('Twin Edit From Left Counts')

    ax2.set_xticks(indices)
    ax2.set_xticklabels([str(x) for x in range(ins_start, ins_end+1)], rotation=45, ha="center")
    ax2.tick_params(axis="x", labelsize=8)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.legend()

    # Stacked bar plot
    # number of continuous "T"s from right at each position
    ax3.bar(indices, from_right_all_edit_counts, label="Continuous TPEs From Right", color=category_colors["Right Flap"])
    # number of continuous "T"s from right at each position for match_arr's that have indels
    ax3.bar(
        indices,
        from_right_all_edit_counts_with_indels,
        label="Continuous TPEs From Right with Indels", 
        color=category_colors["WT"]
    )

    # Adding labels and title
    ax3.set_xlabel('Compound Reference Index')
    ax3.set_ylabel('Read Counts')
    # ax3.set_title('Twin Edit From Right Counts')

    ax3.set_xticks(indices)
    ax3.set_xticklabels([str(x) for x in range(ins_start, ins_end+1)], rotation=45, ha="center")
    # Dynamic? >35 x-tick = font 8 & rotate 45
    ax3.tick_params(axis="x", labelsize=8)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.legend()
    plt.tight_layout()

    plt.savefig(fig_root + "/a7.TPEs_vs_indels.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a7.TPEs_vs_indels.png", bbox_inches='tight')


def plot_editing_summary(
        full_deletion_counts, 
        full_insertion_counts, 
        full_substitution_counts, 
        full_edit_counts, 
        full_total_counts, 
        # ins_start, 
        ins_end, 
        del_start, 
        # del_end, 
        fig_root=None, 
        produce_png=False, 
        category_colors=None
    ):
    x = np.arange(len(full_deletion_counts))
    w = 0.27
    x_labels = [str(x) for x in range(del_start, ins_end+1)]

    misedit_counts = np.array(full_deletion_counts) + np.array(full_insertion_counts) + np.array(full_substitution_counts)
    unedited_counts = np.array(full_total_counts) - np.array(full_edit_counts) - np.array(misedit_counts)
    misedit_counts = misedit_counts.tolist()
    unedited_counts = unedited_counts.tolist()

    plt.figure(figsize=(16,4), dpi=300)

    plt.bar(x - w, full_edit_counts, width=w, label="Edited", color=category_colors["Perfect PE"])
    plt.bar(x, unedited_counts, width=w, label="Unedited", color=category_colors["Left Flap"])
    plt.bar(x + w, misedit_counts, width=w, label="Misedited", color=category_colors["Right Flap"])

    # Show all x-axis labels
    plt.xticks(x, x_labels, rotation=45, ha="center", fontsize=7)
    plt.minorticks_on()
    plt.tick_params(axis="y", which="minor")
    plt.tick_params(axis="y", which="major")
    plt.tick_params(axis="x", which="minor", bottom=False)

    # plt.xlabel("TwinPE Edit Index")
    plt.xlabel("Compound Reference Index")
    plt.ylabel("Counts")
    plt.title("Per-Position Editing Across Deleted and Inserted Regions")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()

    plt.savefig(fig_root + "/a8.Edit_type_summary.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a8.Edit_type_summary.png", bbox_inches='tight')


def plot_nonprogrammed_edit_counts(full_deletion_counts, full_insertion_counts, full_substitution_counts, ins_start, ins_end, del_start, del_end, fig_root=None, produce_png=False):
    x = np.arange(len(full_deletion_counts))
    w = 0.27
    x_labels = [str(x) for x in range(del_start, ins_end+1)]


    plt.figure(figsize=(16,4), dpi=300)

    plt.bar(x - w, full_deletion_counts, width=w, label="Deletions")
    plt.bar(x, full_insertion_counts, width=w, label="Insertions")
    plt.bar(x + w, full_substitution_counts, width=w, label="Substitutions")

    # Show all x-axis labels
    plt.xticks(x, x_labels, rotation=45, ha="center", fontsize=7)
    plt.minorticks_on()
    plt.tick_params(axis="y", which="minor")
    plt.tick_params(axis="y", which="major")
    plt.tick_params(axis="x", which="minor", bottom=False)

    # plt.xlabel("TwinPE Edit Index")
    plt.xlabel("Compound Reference Index")
    plt.ylabel("Counts")
    plt.title("Non-programmed Edits")
    plt.legend()
    plt.tight_layout()

    plt.savefig(fig_root + "/a9.Nonprogrammed_edits.pdf", bbox_inches='tight')
    if produce_png:
        plt.savefig(fig_root + "/a9.Nonprogrammed_edits.png", bbox_inches='tight')


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
    if ref_aln_seq[0] == "-":
        aln_index = 0
        read_start_bases = ""
        while aln_index < len(ref_aln_seq) and ref_aln_seq[aln_index] == "-":
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
        if ref_base == "-":
            refpos_dict[last_nongap_ref_pos] += read_base
        else:
            refpos_dict[ref_pos] += read_base
            last_nongap_ref_pos = ref_pos
            ref_pos += 1
    return refpos_dict


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
                    extend_left_non_gap=None):         # per-sgRNA left extension in non-gap bases
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
        display_start = max(0, this_sgRNA_start - left_extra_cols)

        # Draw as multiple rectangles over non-gap runs only
        runs = _non_gap_runs(ref_row_seq, display_start, this_sgRNA_end)
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


### Allele plot
# We need to customize the seaborn heatmap class and function
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


def prep_alleles_table(
    df_alleles,
    reference_seq,
    ref_aln_seq_region,
    twin_aln_seq_region,
    MAX_N_ROWS,
    MIN_FREQUENCY,
    pegRNA_intervals
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
        ]
        to_append = np.array([{}] * len(idx), dtype=object)
        to_append[idxs_sub] = {"weight": "bold", "color": "black", "size": 16}
        per_element_annot_kws.append(to_append)

    if twin_aln_seq_region is not None and len(twin_aln_seq_region) > 0:
        # Append two blank rows to the end of outputs
        for _ in range(num_blanks):
            blank_row = " " * len(twin_aln_seq_region)
            X.append(seq_to_numbers(blank_row))
            annot.append(list(blank_row))
            y_labels.append("")
            is_reference.append(False)
            to_append = np.array([{"color": "white"}] * len(blank_row), dtype=object)
            per_element_annot_kws.append(to_append)

        # Append twin_aln_seq_region to the end of outputs
        X.append(seq_to_numbers(twin_aln_seq_region.upper()))
        annot.append(list(twin_aln_seq_region))
        y_labels.append("TwinPE Reference")
        # if len(pegRNA_intervals) > 1:
        #     y_labels.append("TwinPE Reference")
        # else:
        #     y_labels.append("PE Reference")
        is_reference.append(False)
        # # detect substitutions (non-gap mismatches) between twin_aln_seq_region and ref_aln_seq_region to style them in bold/black
        # idxs_sub = [
        #     i_sub
        #     for i_sub in range(len(twin_aln_seq_region))
        #     if (ref_aln_seq_region[i_sub] != twin_aln_seq_region[i_sub])
        #     and (ref_aln_seq_region[i_sub] != "-")
        #     and (twin_aln_seq_region[i_sub] != "-")
        # ]
        to_append = np.array([{}] * len(twin_aln_seq_region), dtype=object)
        # to_append[idxs_sub] = {"weight": "bold", "color": "black", "size": 16}
        per_element_annot_kws.append(to_append)
        # # track insertions in twin_aln_seq_region relative to ref_aln_seq_region
        # for p in re_find_indels.finditer(ref_aln_seq_region):
        #     insertion_dict[idx_row + num_blanks].append((p.start(), p.end()))

    return X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference


def plot_alleles_heatmap(
        reference_aln_seq,
        # twin_aln_seq=twin_aln_seq_region,
        X,
        annot,
        y_labels,
        insertion_dict,
        per_element_annot_kws,
        fig_filename_root=None,
        custom_colors=None,
        SAVE_ALSO_PNG=False,
        plot_cut_point=True,
        cut_point_ind=None,
        sgRNA_intervals=None,
        sgRNA_names=None,
        sgRNA_mismatches=None,
        extend_left_non_gap=None,
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
    -SAVE_ALSO_PNG: whether to write png file as well
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

    if N_ROWS < 1:
        fig, ax = plt.subplots()
        fig.text(0.5, 0.5, 'No Alleles', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
        ax.set_clip_on(False)

        if fig_filename_root is None:
            plt.show()
        else:
            fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight')
            if SAVE_ALSO_PNG:
                fig.savefig(fig_filename_root+'.png', bbox_inches='tight')
        plt.close(fig)
        return

    sgRNA_rows = []
    num_sgRNA_rows = 0

    if sgRNA_intervals and len(sgRNA_intervals) > 0:
        sgRNA_rows = get_rows_for_sgRNA_annotation(sgRNA_intervals, plot_nuc_len)
        num_sgRNA_rows = max(sgRNA_rows) + 1
        fig=plt.figure(figsize=(plot_nuc_len*0.3, (N_ROWS+1 + num_sgRNA_rows)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+2, N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+2, N_COLUMNS)
        #ax_hm_ref heatmap for the reference
        ax_hm_ref=plt.subplot(gs1[0:1,:])
        ax_hm=plt.subplot(gs2[2:,:])
    else:
        fig=plt.figure(figsize=(plot_nuc_len*0.3, (N_ROWS+1)*0.6))
        gs1 = gridspec.GridSpec(N_ROWS+1, N_COLUMNS)
        gs2 = gridspec.GridSpec(N_ROWS+1, N_COLUMNS)
        #ax_hm_ref heatmap for the reference
        ax_hm_ref=plt.subplot(gs1[0,:])
        ax_hm=plt.subplot(gs2[1:,:])


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

        # Normalize extend_left_non_gap to a list the length of sgRNA_intervals
        if extend_left_non_gap is None:
            extend = [0]*len(sgRNA_intervals)
        elif isinstance(extend_left_non_gap, int):
            extend = [extend_left_non_gap]*len(sgRNA_intervals)
        elif isinstance(extend_left_non_gap, dict):
            # allow mapping by index or (if available) name
            extend = [0]*len(sgRNA_intervals)
            for i in range(len(sgRNA_intervals)):
                extend[i] = extend_left_non_gap.get(
                    sgRNA_names[i] if sgRNA_names and i < len(sgRNA_names) else i,
                    extend_left_non_gap.get(i, 0)
                )
        else:
            extend = list(extend_left_non_gap)
            if len(extend) < len(sgRNA_intervals):
                extend += [0]*(len(sgRNA_intervals)-len(extend))
            else:
                extend = extend[:len(sgRNA_intervals)]

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
            extend_left_non_gap=extend
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
    #### End of New method ####
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

    if fig_filename_root is None:
        plt.show()
    else:
        fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
        if SAVE_ALSO_PNG:
            fig.savefig(fig_filename_root+'.png', bbox_inches='tight', bbox_extra_artists=(lgd,))
    plt.close(fig)


# def get_dataframe_allele_region(
#     df_alleles, pegRNA_intervals, ref_seq, ref_aln_seq, twin_aln_seq, cut_points
# ):
#     """
#     Extract region of interest for plotting allele tables based on pegRNA positions on reference sequence.
#     Returns modified reference sequence and cut points.
#     """
#     #### Previous new method ####
#     # support single or dual pegRNAs
#     # if df_alleles.shape[0] == 0:
#     #     return df_alleles
#     # if len(pegRNA_intervals) == 2:
#     #     left_end = max(0, pegRNA_intervals[0][0] - 20)
#     #     right_end = min(len(ref_seq), pegRNA_intervals[1][1] + 20)
#     #     cut_points[0] = cut_points[0] - left_end
#     #     cut_points[1] = cut_points[1] - left_end
#     # elif len(pegRNA_intervals) == 1:
#     #     left_end = max(0, pegRNA_intervals[0][0] - 80)
#     #     right_end = min(len(ref_seq), pegRNA_intervals[0][1] + 80)
#     #     cut_points = pegRNA_intervals[0][1] - 3
#     #### End of Previous new method ####

#     #### New method for single and double guide RNAs - don't need to account for single guide RNAs
#     if df_alleles.shape[0] == 0:
#         # Return empty frame and pass-through sequences, with None for coords
#         return df_alleles, ref_seq, ref_aln_seq, twin_aln_seq, cut_points, pegRNA_intervals

#     # Default padding
#     if len(pegRNA_intervals) == 2:
#         # left_end  = max(0, pegRNA_intervals[0][0] - 20)
#         left_end  = max(0, pegRNA_intervals[0][0] - 50)
#         right_end = min(len(ref_seq), pegRNA_intervals[1][1] + 80)
#     elif len(pegRNA_intervals) == 1:
#         # left_end  = max(0, pegRNA_intervals[0][0] - 80)
#         left_end  = 0
#         # right_end = min(len(ref_seq), pegRNA_intervals[0][1] + 80)
#         right_end = len(ref_seq)-1
#     else:
#         # Fallback: show full sequence
#         left_end, right_end = 0, len(ref_seq)

#     # Shift cut points into window coords (if provided)
#     cut_points_window = None
#     if cut_points is not None:
#         try:
#             cut_points_window = [cp - left_end for cp in cut_points]
#         except TypeError:
#             # Single integer
#             cut_points_window = [cut_points - left_end]

#     # Shift sgRNA intervals into window coords
#     pegRNA_intervals_window = []
#     if pegRNA_intervals is not None:
#         for (s, e) in pegRNA_intervals:
#             pegRNA_intervals_window.append((s - left_end, e - left_end))
#     ##### END OF NEW METHOD ######

#     #### Original method ####
#     # if df_alleles.shape[0] == 0:
#     #     return df_alleles

#     # left_end = max(0, pegRNA_intervals[0][0] - 20)
#     # right_end = min(len(ref_seq), pegRNA_intervals[1][1] + 20)

#     # cut_points[0] = cut_points[0] - left_end
#     # cut_points[1] = cut_points[1] - left_end
#     #### End of Original method ####


#     # don’t overwrite original
#     df = df_alleles.copy()

#     # group first
#     df = (
#         df.groupby(['Aligned_Sequence', 'Reference_Sequence', 'Read_Status',
#                     'n_deleted', 'n_inserted', 'n_mutated'], dropna=False)
#           .sum()
#           .reset_index()
#     )
#     # then slice
#     df['Aligned_Sequence']   = df['Aligned_Sequence'].str.slice(left_end, right_end)
#     df['Reference_Sequence'] = df['Reference_Sequence'].str.slice(left_end, right_end)

#     df = df.set_index('Aligned_Sequence')

#     # Caused bug by slicing before grouping
#     # # Slice the two sequence columns only
#     # df["Aligned_Sequence"] = df["Aligned_Sequence"].str.slice(left_end, right_end)
#     # df["Reference_Sequence"] = df["Reference_Sequence"].str.slice(left_end, right_end)

#     # # Aggregate counts by sequence + metadata
#     # df = (
#     #     df.groupby(
#     #         [
#     #             "Aligned_Sequence",
#     #             "Reference_Sequence",
#     #             "Read_Status",
#     #             "n_deleted",
#     #             "n_inserted",
#     #             "n_mutated",
#     #         ],
#     #         dropna=False,
#     #     )
#     #     .sum()
#     #     .reset_index()
#     #     .set_index("Aligned_Sequence")
#     # )

#     df.sort_values(
#         by=["#Reads", "Aligned_Sequence", "Reference_Sequence"],
#         ascending=[False, True, True],
#         inplace=True,
#     )

#     df["Unedited"] = df["Read_Status"].eq("UNMODIFIED")

#     return (
#         df,
#         ref_seq[left_end:right_end],
#         ref_aln_seq[left_end:right_end],
#         twin_aln_seq[left_end:right_end],
#         cut_points,
#         pegRNA_intervals_window
#     )


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
    df = (
        df.groupby(
            ['Aligned_Sequence', 'Reference_Sequence', 'Read_Status',
             'n_deleted', 'n_inserted', 'n_mutated'],
            dropna=False
        ).sum().reset_index()
    )

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
        by=["#Reads", "Aligned_Sequence", "Reference_Sequence"],
        ascending=[False, True, True],
        inplace=True,
    )
    df["Unedited"] = df["Read_Status"].eq("UNMODIFIED")

    return (
        df,
        ref_seq,                 # keep raw ref_seq unchanged
        ref_aln_seq_trim,        # windowed+trimmed aligned reference
        twin_aln_seq_trim,       # windowed+trimmed aligned twin
        cut_points_window if cut_points_window is not None else cut_points,
        pegRNA_intervals_window
    )


def setAlleleMatplotlibDefaults():
    font = {"size": 22}
    matplotlib.rc("font", **font)
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    sns.set(style="white", font_scale=2.2)


# def plot_categorical_allele_tables(
#     min_frequency,
#     max_n_rows,
#     df_alleles,
#     ref_seq,
#     ref_aln_seq,
#     twin_aln_seq,
#     pegRNA_cut_points,
#     pegRNA_plot_cut_points,
#     pegRNA_intervals,
#     pegRNA_mismatches,
#     pegRNA_names,
#     spacer_info,
#     fig_root=None,
#     produce_png=False,
# ):
#     """ """

#     # need to implement custom names in args
#     if pegRNA_names == ["", ""]:
#         pegRNA_names = ["pegRNA", "pegRNA"]

#     if spacer_info['spacer_a_index_wt'] < spacer_info['spacer_b_index_wt']:
#         b_extra = spacer_info['spacer_b_num_bases_removed']
#     else:
#         b_extra = spacer_info['spacer_a_num_bases_removed']

#     pegRNA_cut_points[1] = pegRNA_cut_points[1] - b_extra

#     # allele_table_plotted_count = 0
#     # for cat in df_alleles["Category"].unique():
#     for cat, name in [("Perfect_PE", "b1.Perfect_PE"), 
#                       ("PE_Indel", "b2.PE_Indel"),
#                       ("Imperfect_PE", "b3.Imperfect_PE"), 
#                       ("Left_Flap", "b4.Left_Flap"), 
#                       ("Right_Flap", "b5.Right_Flap"), 
#                       ("Imperfect_WT", "b6.Imperfect_WT"), 
#                       ("WT_Indel", "b7.WT_Indel"), 
#                       ("WT", "b8.WT"), 
#                       ("Uncategorized", "b9.Uncategorized")]:
#         if len(df_alleles[df_alleles["Category"] == cat]) > 0:
#             df_alleles_cat = df_alleles[df_alleles["Category"] == cat]
#             (
#                 df_alleles_around_region,
#                 ref_seq_region,
#                 ref_aln_seq_region,
#                 twin_aln_seq_region,
#                 pegRNA_cut_points,
#                 pegRNA_intervals_region
#             ) = get_dataframe_allele_region(
#                 df_alleles_cat,
#                 pegRNA_intervals,
#                 ref_seq,
#                 ref_aln_seq,
#                 twin_aln_seq,
#                 pegRNA_cut_points,
#             )

#             X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = (
#                 prep_alleles_table(
#                     df_alleles_around_region,
#                     ref_seq_region,
#                     ref_aln_seq_region,
#                     twin_aln_seq_region,
#                     max_n_rows,
#                     min_frequency,
#                     pegRNA_intervals
#                 )
#             )

#             plot_alleles_heatmap(
#                 # reference_seq=ref_seq_region,
#                 reference_aln_seq=ref_aln_seq_region,
#                 X=X,
#                 annot=annot,
#                 y_labels=y_labels,
#                 insertion_dict=insertion_dict,
#                 per_element_annot_kws=per_element_annot_kws,
#                 fig_filename_root=fig_root + f"/{name}",
#                 SAVE_ALSO_PNG=produce_png,
#                 plot_cut_point=pegRNA_plot_cut_points,
#                 cut_point_ind=pegRNA_cut_points,
#                 sgRNA_intervals=pegRNA_intervals_region,
#                 sgRNA_names=pegRNA_names,
#                 sgRNA_mismatches=pegRNA_mismatches,
#                 category=cat,
#                 extend_left_non_gap=[0, b_extra]   # per-sgRNA non-gap extension
#             )
#             # allele_table_plotted_count += 1


def plot_categorical_allele_tables(
    min_frequency,
    max_n_rows,
    df_alleles,
    ref_seq,
    ref_aln_seq,
    twin_aln_seq,
    pegRNA_cut_points,
    pegRNA_plot_cut_points,
    pegRNA_intervals,
    pegRNA_mismatches,
    pegRNA_names,
    spacer_info,
    fig_root=None,
    produce_png=False, 
    plot_full_reads=False
):
    """ """
    # Ensure names
    if pegRNA_names == ["", ""]:
        pegRNA_names = ["pegRNA", "pegRNA"]

    # Compute b_extra once
    if spacer_info['spacer_a_index_wt'] < spacer_info['spacer_b_index_wt']:
        b_extra = spacer_info['spacer_b_num_bases_removed']
    else:
        b_extra = spacer_info['spacer_a_num_bases_removed']

    # Do NOT mutate pegRNA_cut_points in place; make ordered, adjusted copies per plot
    # Build a single ordering for intervals and keep related arrays consistent
    order = None
    if pegRNA_intervals and len(pegRNA_intervals) == 2:
        order = np.argsort([x[0] for x in pegRNA_intervals])
        intervals_ordered = [pegRNA_intervals[i] for i in order]
        cut_points_ordered = [pegRNA_cut_points[i] for i in order]
        plot_cut_points_ordered = [pegRNA_plot_cut_points[i] for i in order]
        mismatches_ordered = [pegRNA_mismatches[i] for i in order] if pegRNA_mismatches else pegRNA_mismatches
        names_ordered = [pegRNA_names[i] for i in order] if pegRNA_names else pegRNA_names
    else:
        intervals_ordered = pegRNA_intervals
        cut_points_ordered = pegRNA_cut_points if isinstance(pegRNA_cut_points, list) else [pegRNA_cut_points]
        plot_cut_points_ordered = pegRNA_plot_cut_points
        mismatches_ordered = pegRNA_mismatches
        names_ordered = pegRNA_names

    # Adjust second cut point (post-ordering)
    if cut_points_ordered and len(cut_points_ordered) == 2:
        cut_points_adjusted = [cut_points_ordered[0], cut_points_ordered[1] - b_extra]
    else:
        cut_points_adjusted = cut_points_ordered

    for cat, name in [
        ("Perfect_PE", "b1.Perfect_PE"),
        ("PE_Indel", "b2.PE_Indel"),
        ("Left_Flap", "b3.Left_Flap"),
        ("Right_Flap", "b4.Right_Flap"),
        ("Imperfect_PE", "b5.Imperfect_PE"),
        ("Imperfect_WT", "b6.Imperfect_WT"),
        ("WT_Indel", "b7.WT_Indel"),
        ("WT", "b8.WT"),
        ("Uncategorized", "b9.Uncategorized"),
    ]:
        if len(df_alleles[df_alleles["Category"] == cat]) == 0:
            continue

        df_alleles_cat = df_alleles[df_alleles["Category"] == cat]

        # (
        #     df_alleles_around_region,
        #     ref_seq_region,
        #     ref_aln_seq_region,
        #     twin_aln_seq_region,
        #     cut_points_window,
        #     pegRNA_intervals_region
        # ) = get_dataframe_allele_region(
        #     df_alleles_cat,
        #     intervals_ordered,
        #     ref_seq,
        #     ref_aln_seq,
        #     twin_aln_seq,
        #     cut_points_adjusted,
        # )


        (
            df_alleles_around_region,
            ref_seq_region,
            ref_aln_seq_region,
            twin_aln_seq_region,
            cut_points_window,
            pegRNA_intervals_region
        ) = get_dataframe_allele_region(
            df_alleles_cat,
            intervals_ordered,
            ref_seq,
            ref_aln_seq,
            twin_aln_seq,
            cut_points_adjusted,
            window_by_intervals=not plot_full_reads, 
            left_pad=20,
            right_pad=20
        )


        X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference = prep_alleles_table(
            df_alleles_around_region,
            ref_seq_region,
            ref_aln_seq_region,
            twin_aln_seq_region,
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
            SAVE_ALSO_PNG=produce_png,
            plot_cut_point=plot_cut_points_ordered,
            cut_point_ind=cut_points_window,
            sgRNA_intervals=pegRNA_intervals_region,
            sgRNA_names=names_ordered,
            sgRNA_mismatches=mismatches_ordered,
            category=cat,
            extend_left_non_gap=[0, b_extra] if pegRNA_intervals and len(pegRNA_intervals) == 2 else None
        )


if __name__ == "__main__":
    main()
