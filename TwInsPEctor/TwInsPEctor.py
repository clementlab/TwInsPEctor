import os
import re
import sys
import html
import shlex
import shutil
import zipfile
import argparse
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

# Recent CRISPresso2 update causes global style changes that break plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib import colors as colors_mpl
import matplotlib.ticker
import matplotlib.path as mpath
# snapshot current rcParams before importing CRISPResso2
_RC_BEFORE_CRISP = matplotlib.rcParams.copy()
from CRISPResso2 import CRISPRessoShared
# restore rcParams to the snapshot to undo CRISPResso2 global style changes
matplotlib.rcParams.update(_RC_BEFORE_CRISP)

CATEGORY_COLORS = { 
    "Perfect TPE":    "#244d6e",  # Deep ocean navy
    "Dual Flap":      "#6491ac",  # Steel denim blue
    "Flap A":         "#868686",  # Neutral mid-gray
    "Flap B":         "#cccccc",  # Light silver-gray
    "Flap A Hybrid":  "#8e7a68",  # Brown-taupe
    "Flap B Hybrid":  "#dad0bc",  # Almond cream
    "Imperfect TPE":  "#b4c9cb",  # Soft slate blue
    "Null":           "#ab9da5",  # Muted slate mauve
    "Imperfect WT":   "#e09e85",  # Warm terracotta
    "WT":             "#8c1c24",  # Deep brick red
}

CATEGORY_ORDER = [
    "Perfect TPE",
    "Dual Flap",
    "Flap A", 
    "Flap B", 
    "Flap A Hybrid", 
    "Flap B Hybrid", 
    "Imperfect TPE", 
    "Null", 
    "Imperfect WT",
    "WT"
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


def get_rt_templates(rt_templates_arg):
    rt_templates = [s.strip() for s in rt_templates_arg.split(",")]
    if len(rt_templates) != 2:
        raise ValueError("pegRNA reverse transcriptase templates must be entered as two comma-separated sequences")
    rt_template_a, rt_template_b = rt_templates
    rt_template_a = rt_template_a.upper()
    rt_template_b = rt_template_b.upper()

    return rt_template_a, rt_template_b


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    return ''.join(complement.get(base, base) for base in reversed(seq.upper()))


def get_replacement_base_changes(comp_ref_seq=None, wt_aln_seq=None, tpe_aln_seq=None):
    bp_changes_arr = []
    for idx in range(len(comp_ref_seq)):
        wt_base = wt_aln_seq[idx]
        twin_base = tpe_aln_seq[idx]
        if wt_base != '-' and twin_base == '-':
            bp_changes_arr.append((idx, wt_base, twin_base))
        elif wt_base == '-' and twin_base != '-':
            bp_changes_arr.append((idx, wt_base, twin_base))
        elif wt_base != twin_base:
            raise ValueError('Substitution detected in Replacement mode.')

    return bp_changes_arr


def get_recoding_base_changes(wt_seq=None, tpe_seq=None, composite_wt=None, composite_tpe=None):
    std_bp_changes_arr = []
    for idx, (wt_base, twin_base) in enumerate(zip(wt_seq, tpe_seq)):
        if wt_base != twin_base:
            std_bp_changes_arr.append((idx, wt_base, twin_base))

    comp_bp_changes_arr = []
    for idx, (wt_base, twin_base) in enumerate(zip(composite_wt, composite_tpe)):
        if wt_base != twin_base:
            comp_bp_changes_arr.append((idx, wt_base, twin_base))
    
    return std_bp_changes_arr, comp_bp_changes_arr



def get_template_overlap(rt_template_a, rt_template_a_start_tpe, rt_template_b, rt_template_b_start_tpe, tpe_ref, wt_ref=None, composite_tpe=None, composite_wt=None, recoding_mode=False):
    if (
        rt_template_a is None 
        or rt_template_b is None 
        or rt_template_a_start_tpe is None 
        or rt_template_b_start_tpe is None
    ):
        return None

    rt_template_a_end_tpe = rt_template_a_start_tpe + len(rt_template_a)
    rt_template_b_end_tpe = rt_template_b_start_tpe + len(rt_template_b)

    overlap_start = max(rt_template_a_start_tpe, rt_template_b_start_tpe)
    overlap_end = min(rt_template_a_end_tpe, rt_template_b_end_tpe)

    if overlap_start < overlap_end:
        overlap_seq = tpe_ref[overlap_start:overlap_end]
        overlap_length = len(overlap_seq)
    else:
        overlap_seq = None
        overlap_length = 0

    bp_changes_a = []
    bp_changes_b = []
    bp_changes_overlap = []

    if recoding_mode:
        std_bp_changes_arr, _ = get_recoding_base_changes(wt_ref, tpe_ref, composite_wt, composite_tpe)
        for change in std_bp_changes_arr:
            pos, wt_base, edited_base = change
            
            in_a = rt_template_a_start_tpe <= pos < rt_template_a_end_tpe
            
            in_b = rt_template_b_start_tpe <= pos < rt_template_b_end_tpe

            if in_a:
                bp_changes_a.append(change)
            if in_b:
                bp_changes_b.append(change)
            if in_a and in_b:
                bp_changes_overlap.append(change)

    result = {
        "overlap_start_tpe": overlap_start if overlap_length > 0 else None,
        "overlap_end_tpe": overlap_end if overlap_length > 0 else None,
        "overlap_length": overlap_length,
        "overlap_sequence": overlap_seq, 
        "inserted_sequence": tpe_ref[rt_template_a_start_tpe:rt_template_b_end_tpe],
    }

    if recoding_mode:
        result.update({
            "rt_template_a_base_change_coverage": bp_changes_a,
            "rt_template_b_base_change_coverage": bp_changes_b,
            "rt_overlap_base_change_coverage": bp_changes_overlap, 
            "std_bp_changes_arr": std_bp_changes_arr
        })

    return result


def analyze_references(wt_seq, tpe_seq, spacer_a, spacer_b, rt_template_a, rt_template_b, cleavage_offset_a, cleavage_offset_b, output_root, recoding_mode=False):
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

    # Find shared starting and ending bases in inserted_seq and deleted_seq for replacement mode, also needed for homology adjustment in both modes
    num_bases_shared_start = 0
    num_bases_shared_end = 0
    num_bases_shared_start_for_homology_adj = 0
    num_bases_shared_end_for_homology_adj = 0
    for i in range(min(len(tpe_inserted_seq), len(wt_deleted_seq))):
        if tpe_inserted_seq[i] == wt_deleted_seq[i]:
            num_bases_shared_start_for_homology_adj += 1
            if not recoding_mode:
                num_bases_shared_start += 1
        else:
            break

    for i in range(1, min(len(tpe_inserted_seq), len(wt_deleted_seq)) + 1):
        if tpe_inserted_seq[-i] == wt_deleted_seq[-i]:
            num_bases_shared_end_for_homology_adj += 1
            if not recoding_mode:
                num_bases_shared_end += 1
        else:
            break

    # Optional template overlap analysis
    # Find rt_template A in TPE sequence
    if rt_template_a:
        if tpe_seq.find(rt_template_a) != -1:
            rt_template_a_start_tpe = tpe_seq.find(rt_template_a)
            is_rt_template_a_rc = False
        elif tpe_seq.find(reverse_complement(rt_template_a)) != -1:
            rt_template_a_start_tpe = tpe_seq.find(reverse_complement(rt_template_a))
            is_rt_template_a_rc = True
        else:
            print("Warning: Could not find pegRNA reverse transcriptase template A in TPE sequence")

    # Find rt_template B in TPE sequence
    if rt_template_b:
        if tpe_seq.find(rt_template_b) != -1:
            rt_template_b_start_tpe = tpe_seq.find(rt_template_b)
            is_rt_template_b_rc = False
        elif tpe_seq.find(reverse_complement(rt_template_b)) != -1:
            rt_template_b_start_tpe = tpe_seq.find(reverse_complement(rt_template_b))
            is_rt_template_b_rc = True
        else:
            print("Warning: Could not find pegRNA reverse transcriptase template B in TPE sequence")
  
    # Get rt template A and B overlap info
    overlap_info = None
    rt_template_b_start_inserted_seq = None
    rt_template_b_start_bp_changes_arr = None
    std_bp_changes_arr_len = None
    rt_template_a_base_change_coverage_len = None
    if rt_template_a and rt_template_b and recoding_mode:
        overlap_info = get_template_overlap(rt_template_a, rt_template_a_start_tpe, rt_template_b, rt_template_b_start_tpe, tpe_seq, wt_seq, composite_tpe, composite_wt, recoding_mode)
        rt_template_b_start_bp_changes_arr = len(overlap_info["std_bp_changes_arr"]) - len(overlap_info["rt_template_b_base_change_coverage"])
        std_bp_changes_arr_len = len(overlap_info["std_bp_changes_arr"])
        rt_template_a_base_change_coverage_len = len(overlap_info["rt_template_a_base_change_coverage"])
    elif rt_template_a and rt_template_b and not recoding_mode:
        overlap_info = get_template_overlap(rt_template_a, rt_template_a_start_tpe, rt_template_b, rt_template_b_start_tpe, tpe_seq)
        rt_template_b_start_inserted_seq = len(overlap_info["inserted_sequence"]) - len(rt_template_b)

    # Check that inserted sequence based on cleavage offsets matches that based on rt templates
    if rt_template_a and rt_template_b and overlap_info["inserted_sequence"] != tpe_inserted_seq:
        print("Warning: Inserted sequence based on cleavage offsets (default: -3) does not match inserted sequence based on provided rt templates")

    with open(
        os.path.join(output_root, "d8.reference_sequences.txt"), "w"
    ) as fout:
        fout.write(f"@Sequence inputs\n\n")
        fout.write(f">Wildtype reference sequence\n{wt_seq}\n\n")
        fout.write(f">TwinPE reference sequence\n{tpe_seq}\n\n")
        fout.write(f">pegRNA spacer a sequence\n{spacer_a}\n\n")
        fout.write(f">pegRNA spacer b sequence\n{spacer_b}\n\n\n")
        fout.write(f"@Composite A alignments\n\n")
        fout.write(f">Wildtype reference sequence alignment\n{wt_aln_seq_a}\n\n")
        fout.write(f">Composite A reference sequence\n{composite_a_ref_seq}\n\n")
        fout.write(f">TwinPE reference sequence alignment\n{tpe_aln_seq_a}\n\n\n")
        fout.write(f"@Composite B alignments\n\n")
        fout.write(f">Wildtype reference sequence alignment\n{wt_aln_seq_b}\n\n")
        fout.write(f">Composite B reference sequence\n{composite_b_ref_seq}\n\n")
        fout.write(f">TwinPE reference sequence alignment\n{tpe_aln_seq_b}\n\n\n")
        if recoding_mode:
            fout.write(f"@Recoding mode base change reference sequences\n\n")
            fout.write(f">Composite WT reference sequence\n{composite_wt}\n\n")
            fout.write(f">Composite TPE reference sequence\n{composite_tpe}\n\n\n")
        fout.write(f"@Modified pegRNA spacer sequences\n\n")
        fout.write(f">pegRNA spacer a for twinPE and composite b reference sequences\n{spacer_a_tpe}\n\n")
        fout.write(f">pegRNA spacer b for twinPE and composite a reference sequences\n{spacer_b_tpe}\n\n")

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
        "wt_seq_replacement_bp_changes": wt_seq_replacement_bp_changes, 
        "num_bases_shared_start": num_bases_shared_start,
        "num_bases_shared_end": num_bases_shared_end, 
        "num_bases_shared_start_for_homology_adj": num_bases_shared_start_for_homology_adj, 
        "num_bases_shared_end_for_homology_adj": num_bases_shared_end_for_homology_adj, 
        "rt_template_a": rt_template_a, 
        "rt_template_b": rt_template_b, 
        "rt_template_a_length": len(rt_template_a) if rt_template_a else None, 
        "rt_template_b_length": len(rt_template_b) if rt_template_b else None, 
        "rt_overlap_length": overlap_info["overlap_length"] if overlap_info else None, 
        "rt_overlap_sequence": overlap_info["overlap_sequence"] if overlap_info else None, 
        "rt_overlap_start_tpe": overlap_info["overlap_start_tpe"] if overlap_info else None,
        "rt_overlap_end_tpe": overlap_info["overlap_end_tpe"] if overlap_info else None,
        "rt_template_a_start_tpe": rt_template_a_start_tpe if rt_template_a else None, 
        "rt_template_b_start_tpe": rt_template_b_start_tpe if rt_template_b else None, 
        "is_rt_template_a_rc": is_rt_template_a_rc if rt_template_a else None, 
        "is_rt_template_b_rc": is_rt_template_b_rc if rt_template_b else None, 
        "rt_template_b_start_inserted_seq": rt_template_b_start_inserted_seq, 
        "rt_template_a_base_change_coverage_len": rt_template_a_base_change_coverage_len, 
        "rt_template_b_start_bp_changes_arr": rt_template_b_start_bp_changes_arr, 
        "std_bp_changes_arr_len": std_bp_changes_arr_len
    }

    return reference_info


def get_crispresso_command(args, extra_crispresso_args, ref_seq, ref_name, spacer_a, spacer_b, 
                           crispresso_output_folder, twinspector_results_folder,
                           n_processes=1, append=True):
    
    cmd = [
        "CRISPResso",
        "--fastq_r1", args.fastq_r1,
        "--amplicon_seq", str(ref_seq),
        "--amplicon_name", str(ref_name),
        "--guide_seq", f"{spacer_a},{spacer_b}", 
        "--default_min_aln_score", "0", 
        "--output_folder", os.path.dirname(crispresso_output_folder), 
        "--write_detailed_allele_table", 
    ]

    if not any(
        argument == "--n_processes" or argument.startswith("--n_processes=")
        for argument in extra_crispresso_args
    ):
        cmd.extend(["--n_processes", str(n_processes)])

    if args.fastq_r2:
        cmd.extend(["--fastq_r2", args.fastq_r2])
    if args.no_rerun:
        cmd.append("--no_rerun")

    # Append the raw list of extra CRISPResso2 arguments parsed from the CLI string
    if extra_crispresso_args:
        cmd.extend(extra_crispresso_args)

    with open(os.path.join(twinspector_results_folder, "d10.crispresso2_commands.txt"), "a" if append else "w") as fout:
        fout.write(" ".join(cmd) + "\n\n")

    return cmd


def run_crispresso_command(cmd, verbose=False):
    if verbose:
        print("Running CRISPResso2 with command:\n", " ".join(cmd), "\n")
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_crispresso_commands_parallel(crispresso_tasks, n_processes, verbose=False):
    with ProcessPoolExecutor(max_workers=min(n_processes, len(crispresso_tasks))) as executor:
        futures = [executor.submit(run_crispresso_command, cmd, verbose) for cmd in crispresso_tasks]
        for future in as_completed(futures):
            future.result()


def get_allele_df_keys(df):
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


def get_allele_match_array(bp_changes_arr, allele_map, del_start, del_end, ins_start, ins_end, num_bases_shared_start, num_bases_shared_end):
    match_arr = ["0"] * len(bp_changes_arr)

    for ind, (comp_ind, wt_base, tpe_base) in enumerate(bp_changes_arr):
        allele_base = allele_map.get(comp_ind, "")

        if allele_base == wt_base:
            match_arr[ind] = "W"  
        elif allele_base == tpe_base:
            match_arr[ind] = "T"  
        elif len(allele_base) > 1:
            if allele_base[0] == wt_base:
                match_arr[ind] = "WI" 
            elif allele_base[0] == tpe_base:
                match_arr[ind] = "TI"  
            elif allele_base[0] in {"A", "C", "G", "T"}:
                match_arr[ind] = "SI"  
            else:
                match_arr[ind] = "NI"
        elif allele_base in {"A", "C", "G", "T"}:
            match_arr[ind] = "S" 
        elif allele_base == "-":
            match_arr[ind] = "D"  
        else:
            match_arr[ind] = "N"  

    # Split match_arr into insertion and deletion regions while recording the corresponding indices in match_arr.
    del_match_arr = []
    ins_match_arr = []
    del_indices = []
    ins_indices = []

    for ind, (comp_ind, _, _) in enumerate(bp_changes_arr):
        if del_start <= comp_ind <= del_end:
            del_match_arr.append(match_arr[ind])
            del_indices.append(ind)
        else:
            ins_match_arr.append(match_arr[ind])
            ins_indices.append(ind)

    # Mark shared bases at the start
    if num_bases_shared_start > 0:
        # Prepend "=" to match_arr
        for idx in ins_indices[:num_bases_shared_start] + del_indices[:num_bases_shared_start]:
            match_arr[idx] = "=" + match_arr[idx]   
        # Prepend "=" to the split arrays
        ins_match_arr[:num_bases_shared_start] = ["=" + val for val in ins_match_arr[:num_bases_shared_start]]
        del_match_arr[:num_bases_shared_start] = ["=" + val for val in del_match_arr[:num_bases_shared_start]]

    # Mark shared bases at the end
    if num_bases_shared_end > 0:
        # Prepend "=" to match_arr
        for idx in ins_indices[-num_bases_shared_end:] + del_indices[-num_bases_shared_end:]:
            match_arr[idx] = "=" + match_arr[idx]
        # Prepend "=" to the split arrays
        ins_match_arr[-num_bases_shared_end:] = ["=" + val for val in ins_match_arr[-num_bases_shared_end:]]
        del_match_arr[-num_bases_shared_end:] = ["=" + val for val in del_match_arr[-num_bases_shared_end:]]

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

    CATEGORY_PRIORITY = {
        "Perfect_TPE": 1,
        "TPE_Indel": 2,
        "WT": 3,
        "WT_Indel": 4,
        "Imperfect_TPE": 5,
        "Flap_A": 6,
        "Flap_B": 6,
        "Imperfect_WT": 7,
    }

    DEFAULT_PRIORITY = 99

    flap_disagreement = (
        (category_a == "Flap_A" and category_b == "Flap_B")
        or
        (category_a == "Flap_B" and category_b == "Flap_A")
    )

    if flap_disagreement and abs(score_a - score_b) <= flap_score_delta:
        return "Imperfect_TPE", "Composite_A&B"

    if not category_a:
        return category_b, "Composite_B"

    if not category_b:
        return category_a, "Composite_A"

    rank_a = CATEGORY_PRIORITY.get(category_a, DEFAULT_PRIORITY)
    rank_b = CATEGORY_PRIORITY.get(category_b, DEFAULT_PRIORITY)

    if rank_b < rank_a:
        return category_b, "Composite_B"

    if rank_a < rank_b:
        return category_a, "Composite_A"

    if score_b > score_a:
        return category_b, "Composite_B"

    return category_a, "Composite_A"


def categorize_alleles(
        df_merged=None, 
        wt_seq=None, 
        tpe_seq=None, 
        reference_info=None, 
        # pegrna_info=None, 
        min_num_base_edits=2, 
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
    df_merged['TPE_Indel'] = False
    df_merged['Category_comp_b'] = ""
    df_merged['tpe_ins_match_arr'] = ""
    df_merged['comp_a_ins_match_arr'] = ""
    df_merged['comp_b_ins_match_arr'] = ""
    df_merged['wt_del_match_arr'] = ""
    df_merged['comp_a_del_match_arr'] = ""
    df_merged['comp_b_del_match_arr'] = ""
    df_merged['wt_all_ins_pos'] = "[]"
    df_merged['wt_all_del_pos'] = "[]"
    df_merged['wt_all_sub_pos'] = "[]"
    df_merged['tpe_all_ins_pos'] = "[]"
    df_merged['tpe_all_del_pos'] = "[]"
    df_merged['tpe_all_sub_pos'] = "[]"

    # Base changes
    if recoding_mode:
        # Uses composite wt/tpe references (and standard wt/tpe references for plotting only)
        std_bp_changes_arr, comp_bp_changes_arr = get_recoding_base_changes(wt_seq, tpe_seq, composite_wt, composite_tpe)
        bp_changes_arrs = {"std_bp_changes_arr": std_bp_changes_arr, "comp_bp_changes_arr": comp_bp_changes_arr}
    else:
        # Uses composite a/b references (and altered standard wt/tpe references that do not directly compare bases for plotting only)
        comp_a_bp_changes_arr = get_replacement_base_changes(comp_a_ref_seq, wt_aln_seq_comp_a, tpe_aln_seq_comp_a)
        comp_b_bp_changes_arr = get_replacement_base_changes(comp_b_ref_seq, wt_aln_seq_comp_b, tpe_aln_seq_comp_b)
        wt_bp_changes_arr = get_replacement_base_changes(comp_ref_seq=wt_seq, wt_aln_seq=wt_seq, tpe_aln_seq=reference_info["tpe_seq_replacement_bp_changes"])
        tpe_bp_changes_arr = get_replacement_base_changes(comp_ref_seq=tpe_seq, wt_aln_seq=reference_info["wt_seq_replacement_bp_changes"], tpe_aln_seq=tpe_seq)
        bp_changes_arrs = {"comp_a_bp_changes_arr": comp_a_bp_changes_arr, "comp_b_bp_changes_arr": comp_b_bp_changes_arr, "wt_bp_changes_arr": wt_bp_changes_arr, "tpe_bp_changes_arr": tpe_bp_changes_arr}

    for idx, allele in df_merged.iterrows():
        # Classify perfect WT alleles and some Imperfect WT alleles using WT reference alignment
        wt_seq_aln_allele = allele.Reference_Sequence_wt
        allele_seq_aln_wt = allele.Aligned_Sequence_wt

        wt_map = get_refpos_values(wt_seq_aln_allele, allele_seq_aln_wt)

        if recoding_mode:
            wt_del_match_arr, _, _ = get_allele_match_array(std_bp_changes_arr, wt_map, del_start=0, del_end=0, ins_start=0, ins_end=0, num_bases_shared_start=reference_info["num_bases_shared_start"], num_bases_shared_end=reference_info["num_bases_shared_end"])           
        else:
            wt_del_match_arr, _, _ = get_allele_match_array(wt_bp_changes_arr, wt_map, del_start=0, del_end=0, ins_start=0, ins_end=0, num_bases_shared_start=reference_info["num_bases_shared_start"], num_bases_shared_end=reference_info["num_bases_shared_end"])           

        wt_all_sub_pos, wt_sub_between_cuts, wt_all_del_pos, wt_del_between_cuts, wt_all_ins_pos, wt_has_substitutions, wt_has_deletions, wt_has_insertions = get_mutations(wt_map, wt_seq, reference_info["cut_points_wt"])

        df_merged.at[idx, 'wt_all_ins_pos'] = str(wt_all_ins_pos)
        df_merged.at[idx, 'wt_all_del_pos'] = str(wt_all_del_pos)
        df_merged.at[idx, 'wt_all_sub_pos'] = str(wt_all_sub_pos)

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
            tpe_ins_match_arr, _, _ = get_allele_match_array(std_bp_changes_arr, tpe_map, del_start=0, del_end=0, ins_start=0, ins_end=0, num_bases_shared_start=reference_info["num_bases_shared_start"], num_bases_shared_end=reference_info["num_bases_shared_end"])           
        else:
            tpe_ins_match_arr, _, _ = get_allele_match_array(tpe_bp_changes_arr, tpe_map, del_start=0, del_end=0, ins_start=0, ins_end=0, num_bases_shared_start=reference_info["num_bases_shared_start"], num_bases_shared_end=reference_info["num_bases_shared_end"])

        tpe_all_sub_pos, tpe_sub_between_cuts, tpe_all_del_pos, tpe_del_between_cuts, tpe_all_ins_pos, tpe_has_substitutions, tpe_has_deletions, tpe_has_insertions = get_mutations(tpe_map, tpe_seq, reference_info["cut_points_tpe"])

        df_merged.at[idx, 'tpe_all_ins_pos'] = str(tpe_all_ins_pos)
        df_merged.at[idx, 'tpe_all_del_pos'] = str(tpe_all_del_pos)
        df_merged.at[idx, 'tpe_all_sub_pos'] = str(tpe_all_sub_pos)

        if ignore_extraspacer_deletions:
            tpe_has_all_extraspacer_deletions = all(x < reference_info["spacer_a_start_tpe"] or x > (reference_info["spacer_b_start_tpe"]+len(reference_info["spacer_b_tpe"])) for x in tpe_all_del_pos)

        if allele['Aligned_Reference_Score_tpe'] == 100.0:
            df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'
        elif not tpe_has_insertions and not tpe_has_deletions and not tpe_sub_between_cuts:
            df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'  # Perfect TPEs with substitutions outside of cut sites
        elif ignore_extraspacer_deletions and not tpe_has_insertions and len(tpe_all_sub_pos) <= 1 and tpe_has_deletions and tpe_has_all_extraspacer_deletions:
            if not tpe_sub_between_cuts:
                df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'  # Perfect TPEs with deletions outside of the spacer regions and 0 or 1 substitutions outside of the cut sites
            else:
                df_merged.at[idx, 'Category_tpe'] = 'Imperfect_TPE'  # Imperfect TPEs with deletions outside of the spacer regions but subtitutions within the cut sites
        elif not tpe_has_insertions and not tpe_has_deletions and len(tpe_all_sub_pos) <= 1:
            if not tpe_sub_between_cuts:
                df_merged.at[idx, 'Category_tpe'] = 'Perfect_TPE'  # Perfect TPEs with 0 or 1 substitutions outside of the cut sites
            # Permits a Flap allele with a single substitution withint min_num_base_edits adjacent to the opposing cut site to be classified as Imperfect TPE
            # else:
            #     df_merged.at[idx, 'Category_tpe'] = 'Imperfect_TPE'  # Imperfect TPEs with 0 or 1 substitutions between cut sites

        # Classify all alleles using Composite A reference alignment
        comp_a_seq_aln_allele = allele.Reference_Sequence_comp_a
        allele_seq_aln_comp_a = allele.Aligned_Sequence_comp_a

        comp_a_map = get_refpos_values(comp_a_seq_aln_allele, allele_seq_aln_comp_a)

        if recoding_mode:
            comp_a_match_arr, comp_a_ins_match_arr, comp_a_del_match_arr = get_allele_match_array(comp_bp_changes_arr, comp_a_map, reference_info["composite_a_del_start"], reference_info["composite_a_del_end"], reference_info["composite_a_ins_start"], reference_info["composite_a_ins_end"], reference_info["num_bases_shared_start"], reference_info["num_bases_shared_end"])
        else:
            comp_a_match_arr, comp_a_ins_match_arr, comp_a_del_match_arr = get_allele_match_array(comp_a_bp_changes_arr, comp_a_map, reference_info["composite_a_del_start"], reference_info["composite_a_del_end"], reference_info["composite_a_ins_start"], reference_info["composite_a_ins_end"], reference_info["num_bases_shared_start"], reference_info["num_bases_shared_end"])

        comp_a_all_sub_pos, comp_a_sub_between_cuts, comp_a_all_del_pos, comp_a_del_between_cuts, comp_a_all_ins_pos, comp_a_has_substitutions, comp_a_has_deletions, comp_a_has_insertions = get_mutations(comp_a_map, comp_a_ref_seq, reference_info["cut_points_composite_a"])

        comp_a_has_indel = check_indel_positions(comp_a_all_ins_pos, comp_a_all_del_pos, reference_info["composite_a_del_start"], reference_info["composite_a_del_end"], reference_info["composite_a_ins_start"], reference_info["composite_a_ins_end"], ignore_extraspacer_deletions, reference_info["pegRNA_intervals_composite_a"])

        comp_a_total_TPE_count = comp_a_match_arr.count("T") + comp_a_match_arr.count("TI")
        comp_a_total_shared_base_TPE_count = comp_a_match_arr.count("=T") + comp_a_match_arr.count("=TI")
        comp_a_has_all_TPE = (comp_a_total_TPE_count + comp_a_total_shared_base_TPE_count == len(comp_a_match_arr))
        comp_a_has_any_TPE = (comp_a_total_TPE_count >= min_num_base_edits)
        comp_a_has_any_TPE_in_insertion = (comp_a_ins_match_arr.count("T") + comp_a_ins_match_arr.count("TI") >= min_num_base_edits)

        comp_a_total_WT_count = comp_a_match_arr.count("W") + comp_a_match_arr.count("WI")
        comp_a_total_shared_base_WT_count = comp_a_match_arr.count("=W") + comp_a_match_arr.count("=WI")
        comp_a_has_all_WT = (comp_a_total_WT_count + comp_a_total_shared_base_WT_count == len(comp_a_match_arr))
        comp_a_has_any_WT = (comp_a_total_WT_count > 0)

        comp_a_has_flap_a = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])
        comp_a_has_flap_b = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

        if comp_a_has_all_TPE and comp_a_has_indel:
            df_merged.at[idx,'Category_comp_a'] = "TPE_Indel"
            df_merged.at[idx,'TPE_Indel'] = True
        elif comp_a_has_all_TPE:
            df_merged.at[idx,'Category_comp_a'] = "Perfect_TPE"
        elif comp_a_has_all_WT and comp_a_has_indel:
            df_merged.at[idx,'Category_comp_a'] = "WT_Indel"
        elif comp_a_has_all_WT:
            df_merged.at[idx,'Category_comp_a'] = "WT"
        elif comp_a_has_flap_a and not comp_a_has_flap_b:
            df_merged.at[idx,'Category_comp_a'] = "Flap_A"
        elif comp_a_has_flap_b and not comp_a_has_flap_a:
            df_merged.at[idx,'Category_comp_a'] = "Flap_B"
        elif comp_a_has_any_WT and not comp_a_has_any_TPE_in_insertion:
            df_merged.at[idx,'Category_comp_a'] = "Imperfect_WT"
        elif comp_a_has_any_TPE:
            df_merged.at[idx,'Category_comp_a'] = "Imperfect_TPE"
        elif not comp_a_has_flap_a and not comp_a_has_flap_b and not comp_a_has_any_TPE:
            df_merged.at[idx,'Category_comp_a'] = "Imperfect_WT"
        else:
            df_merged.at[idx,'Category_comp_a'] = "Uncategorized"

        # Classify all alleles using Composite B reference alignment
        comp_b_seq_aln_allele = allele.Reference_Sequence_comp_b
        allele_seq_aln_comp_b = allele.Aligned_Sequence_comp_b

        comp_b_map = get_refpos_values(comp_b_seq_aln_allele, allele_seq_aln_comp_b)

        if recoding_mode:
            comp_b_match_arr, comp_b_ins_match_arr, comp_b_del_match_arr = get_allele_match_array(comp_bp_changes_arr, comp_b_map, reference_info["composite_b_del_start"], reference_info["composite_b_del_end"], reference_info["composite_b_ins_start"], reference_info["composite_b_ins_end"], reference_info["num_bases_shared_start"], reference_info["num_bases_shared_end"])
        else:
            comp_b_match_arr, comp_b_ins_match_arr, comp_b_del_match_arr = get_allele_match_array(comp_b_bp_changes_arr, comp_b_map, reference_info["composite_b_del_start"], reference_info["composite_b_del_end"], reference_info["composite_b_ins_start"], reference_info["composite_b_ins_end"], reference_info["num_bases_shared_start"], reference_info["num_bases_shared_end"])

        comp_b_all_sub_pos, comp_b_sub_between_cuts, comp_b_all_del_pos, comp_b_del_between_cuts, comp_b_all_ins_pos, comp_b_has_substitutions, comp_b_has_deletions, comp_b_has_insertions = get_mutations(comp_b_map, comp_b_ref_seq, reference_info["cut_points_composite_b"])

        comp_b_has_indel = check_indel_positions(comp_b_all_ins_pos, comp_b_all_del_pos, reference_info["composite_b_del_start"], reference_info["composite_b_del_end"], reference_info["composite_b_ins_start"], reference_info["composite_b_ins_end"], ignore_extraspacer_deletions, reference_info["pegRNA_intervals_composite_b"])

        comp_b_total_TPE_count = comp_b_match_arr.count("T") + comp_b_match_arr.count("TI")
        comp_b_total_shared_base_TPE_count = comp_b_match_arr.count("=T") + comp_b_match_arr.count("=TI")
        comp_b_has_all_TPE = (comp_b_total_TPE_count + comp_b_total_shared_base_TPE_count == len(comp_b_match_arr))
        comp_b_has_any_TPE = (comp_b_total_TPE_count >= min_num_base_edits)
        comp_b_has_any_TPE_in_insertion = (comp_b_ins_match_arr.count("T") + comp_b_ins_match_arr.count("TI") >= min_num_base_edits)

        comp_b_total_WT_count = comp_b_match_arr.count("W") + comp_b_match_arr.count("WI")
        comp_b_total_shared_base_WT_count = comp_b_match_arr.count("=W") + comp_b_match_arr.count("=WI")
        comp_b_has_all_WT = (comp_b_total_WT_count + comp_b_total_shared_base_WT_count == len(comp_b_match_arr))
        comp_b_has_any_WT = (comp_b_total_WT_count > 0)

        comp_b_has_flap_a = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])
        comp_b_has_flap_b = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

        if comp_b_has_all_TPE and comp_b_has_indel:
            df_merged.at[idx,'Category_comp_b'] = "TPE_Indel"
            df_merged.at[idx,'TPE_Indel'] = True
        elif comp_b_has_all_TPE:
            df_merged.at[idx,'Category_comp_b'] = "Perfect_TPE"
        elif comp_b_has_all_WT and comp_b_has_indel:
            df_merged.at[idx,'Category_comp_b'] = "WT_Indel"
        elif comp_b_has_all_WT:
            df_merged.at[idx,'Category_comp_b'] = "WT"
        elif comp_b_has_flap_a and not comp_b_has_flap_b:
            df_merged.at[idx,'Category_comp_b'] = "Flap_A"
        elif comp_b_has_flap_b and not comp_b_has_flap_a:
            df_merged.at[idx,'Category_comp_b'] = "Flap_B"
        elif comp_b_has_any_WT and not comp_b_has_any_TPE_in_insertion:
            df_merged.at[idx,'Category_comp_b'] = "Imperfect_WT"
        elif comp_b_has_any_TPE:
            df_merged.at[idx,'Category_comp_b'] = "Imperfect_TPE"
        elif not comp_b_has_flap_a and not comp_b_has_flap_b and not comp_b_has_any_TPE:
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

    # Extend categories
    # Split out dual flap outcomes from Imperfect TPEs
    df_merged["dual_flap"] = None
    for idx, allele in df_merged.iterrows():

        if allele.Category_final == "Imperfect_TPE":
            has_tpe_flap_a = None
            has_tpe_flap_b = None

            if allele.Classified_by == "TPE":
                tpe_ins_match_arr = allele.tpe_ins_match_arr

                has_tpe_flap_a = all(base in {"T", "TI"} for base in tpe_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])
                has_tpe_flap_b = all(base in {"T", "TI"} for base in tpe_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            elif allele.Classified_by == "Composite_A":
                comp_a_ins_match_arr = allele.comp_a_ins_match_arr

                has_tpe_flap_a = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])
                has_tpe_flap_b = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            elif allele.Classified_by == "Composite_B":
                comp_b_ins_match_arr = allele.comp_b_ins_match_arr

                has_tpe_flap_a = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])
                has_tpe_flap_b = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            elif allele.Classified_by == "Composite_A&B":
                comp_a_ins_match_arr = allele.comp_a_ins_match_arr
                comp_b_ins_match_arr = allele.comp_b_ins_match_arr

                has_tpe_flap_a = all(base in {"T", "TI"} for base in comp_b_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])
                has_tpe_flap_b = all(base in {"T", "TI"} for base in comp_a_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            if has_tpe_flap_a is True and has_tpe_flap_b is True:
                df_merged.at[idx, "dual_flap"] = True
            elif has_tpe_flap_a is False and has_tpe_flap_b is False:
                df_merged.at[idx, "dual_flap"] = False

    # Split out hybrid flap outcomes from flap A category
    df_merged["flap_a_hybrid"] = None
    for idx, allele in df_merged.iterrows():

        if allele.Category_final == "Flap_A":
            has_wt_flap_b = None

            if allele.Classified_by == "TPE":
                tpe_ins_match_arr = allele.tpe_ins_match_arr
                has_wt_flap_b = all(base in {"W", "WI"} for base in tpe_ins_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            elif allele.Classified_by == "Composite_A":
                comp_a_del_match_arr = allele.comp_a_del_match_arr
                has_wt_flap_b = all(base in {"W", "WI"} for base in comp_a_del_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            elif allele.Classified_by == "Composite_B":
                comp_b_del_match_arr = allele.comp_b_del_match_arr
                has_wt_flap_b = all(base in {"W", "WI"} for base in comp_b_del_match_arr[-(reference_info["num_bases_shared_end"] + min_num_base_edits):None if reference_info["num_bases_shared_end"] == 0 else -reference_info["num_bases_shared_end"]])

            if has_wt_flap_b is True:
                df_merged.at[idx, "flap_a_hybrid"] = True
            elif has_wt_flap_b is False:
                df_merged.at[idx, "flap_a_hybrid"] = False
    
    # Split out hybrid flap outcomes from flap B category
    df_merged["flap_b_hybrid"] = None
    for idx, allele in df_merged.iterrows():

        if allele.Category_final == "Flap_B":
            has_wt_flap_a = None

            if allele.Classified_by == "TPE":
                tpe_ins_match_arr = allele.tpe_ins_match_arr
                has_wt_flap_a = all(base in {"W", "WI"} for base in tpe_ins_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])

            elif allele.Classified_by == "Composite_A":
                comp_a_del_match_arr = allele.comp_a_del_match_arr
                has_wt_flap_a = all(base in {"W", "WI"} for base in comp_a_del_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])

            elif allele.Classified_by == "Composite_B":
                comp_b_del_match_arr = allele.comp_b_del_match_arr
                has_wt_flap_a = all(base in {"W", "WI"} for base in comp_b_del_match_arr[reference_info["num_bases_shared_start"]:reference_info["num_bases_shared_start"] + min_num_base_edits])

            if has_wt_flap_a is True:
                df_merged.at[idx, "flap_b_hybrid"] = True
            elif has_wt_flap_a is False:
                df_merged.at[idx, "flap_b_hybrid"] = False

    # Split out alleles lacking any sequence (WT or TPE) between nick sites from Imperfect WTs
    df_merged["null"] = None
    for idx, allele in df_merged.iterrows():

        has_null = None
        if allele.Classified_by == "WT":
            wt_del_match_arr = allele.wt_del_match_arr
            has_null = all(base not in {"W", "WI", "T", "TI"} for base in wt_del_match_arr)

        if allele.Classified_by == "TPE":
            wt_del_match_arr = allele.wt_del_match_arr
            has_null = all(base not in {"W", "WI", "T", "TI"} for base in tpe_ins_match_arr)

        elif allele.Classified_by == "Composite_A":
            comp_a_del_match_arr = allele.comp_a_del_match_arr
            comp_a_ins_match_arr = allele.comp_a_ins_match_arr
            has_null = all(base in {"T", "TI", "=T", "=TI"} for base in comp_a_del_match_arr) and all(base in {"W", "WI", "=W", "=WI"} for base in comp_a_ins_match_arr)

        elif allele.Classified_by == "Composite_B":
            comp_b_del_match_arr = allele.comp_b_del_match_arr
            comp_b_ins_match_arr = allele.comp_b_ins_match_arr
            has_null = all(base in {"T", "TI", "=T", "=TI"} for base in comp_b_del_match_arr) and all(base in {"W", "WI", "=W", "=WI"} for base in comp_b_ins_match_arr)

        elif allele.Classified_by == "Composite_A&B":
            comp_a_del_match_arr = allele.comp_a_del_match_arr
            comp_a_ins_match_arr = allele.comp_a_ins_match_arr
            comp_b_del_match_arr = allele.comp_b_del_match_arr
            comp_b_ins_match_arr = allele.comp_b_ins_match_arr
            has_null = all(base in {"T", "TI", "=T", "=TI"} for base in comp_a_del_match_arr) and all(base in {"W", "WI", "=W", "=WI"} for base in comp_a_ins_match_arr) and all(base in {"T", "TI", "=T", "=TI"} for base in comp_b_del_match_arr) and all(base in {"W", "WI", "=W", "=WI"} for base in comp_b_ins_match_arr)
        
        if has_null is True:
            df_merged.at[idx, "null"] = True
        elif has_null is False:
            df_merged.at[idx, "null"] = False

    # Incorporate extended categories into summary counts
    for idx, allele in df_merged.iterrows():
        # Keep cat name as imperfect tpe
        # if allele["Category_final"] == "Imperfect_TPE" and allele["dual_flap"] == False and allele["null"] == False:
        #     df_merged.at[idx, "Category_final"] = "Imperfect_TPE"
        if allele["null"] == True:
            df_merged.at[idx, "Category_final"] = "Null"
        elif allele["dual_flap"] == True:
            df_merged.at[idx, "Category_final"] = "Dual_Flap"
        elif allele["flap_a_hybrid"] == True:
            df_merged.at[idx, "Category_final"] = "Flap_A_Hybrid"
        elif allele["flap_b_hybrid"] == True:
            df_merged.at[idx, "Category_final"] = "Flap_B_Hybrid"
        # Indel categories leading to confusion - best to look at mutation plots for more useful breakdown
        if allele["Category_final"] == "TPE_Indel":
            df_merged.at[idx, "Category_final"] = "Dual_Flap"
        if allele["Category_final"] == "WT_Indel":
            df_merged.at[idx, "Category_final"] = "Imperfect_WT"


    return df_merged, bp_changes_arrs


def filter_out_of_bounds_positions(all_ins_pos, all_del_pos, all_sub_pos, ref_seq_len, pegRNA_intervals=None, ignore_extraspacer_deletions=False):
    
    def process_positions(pos_string, ins=False):
        stripped_string = pos_string.strip('[]').strip()
        if not stripped_string:
            positions = []
        else:
            positions = [int(val) for val in stripped_string.split(',')]
            
        # Filter out positions beyond the reference length (specific to sample issues)
        if ins:
            return [pos for pos in positions if pos < ref_seq_len-1]
        else:
            return [pos for pos in positions if pos < ref_seq_len]
        
    filtered_ins = process_positions(all_ins_pos, ins=True)
    filtered_del = process_positions(all_del_pos)
    filtered_sub = process_positions(all_sub_pos)
    
    # Remove deletions that fall outside the pegRNA window (specific to sample issues)
    if ignore_extraspacer_deletions and pegRNA_intervals is not None:
        min_bound = pegRNA_intervals[0][0]
        max_bound = pegRNA_intervals[1][1]
        
        filtered_del = [pos for pos in filtered_del if min_bound <= pos <= max_bound]
        
    return filtered_ins, filtered_del, filtered_sub


def tally_mutation_positions(mutation_positions_list):
    """
    Tally the occurrences of each mutation position in a list of mutation position lists.
    mutation_positions_list: List of tuples containing actual lists instead of strings.
    Example: [ ([1, 3], [], [7, 7]), ([], [4, 5], [8]) ]
    """
    counts = {"ins": {}, "del": {}, "sub": {}}

    for ins_list, del_list, sub_list in mutation_positions_list:
        for pos in ins_list:
            counts["ins"][pos] = counts["ins"].get(pos, 0) + 1
        for pos in del_list:
            counts["del"][pos] = counts["del"].get(pos, 0) + 1
        for pos in sub_list:
            counts["sub"][pos] = counts["sub"].get(pos, 0) + 1

    return counts


def combine_mutation_position_dicts(dicts_to_combine):
    master_ins = Counter()
    master_del = Counter()
    master_sub = Counter()

    for d in dicts_to_combine:
        master_ins.update(d.get("ins", {}))
        master_del.update(d.get("del", {}))
        master_sub.update(d.get("sub", {}))

    return {
        "ins": dict(master_ins), 
        "del": dict(master_del), 
        "sub": dict(master_sub)
    }


def process_mutation_category(df_categorized, cat_name, prefix, is_wt_mode, wt_len, tpe_len, peg_intervals, ignore_del):
    """Filters data for a specific category, maps positional mutations, and tallies counts."""
    df_cat = df_categorized[df_categorized["Category_final"] == cat_name].copy()
    
    for col in ["ins", "del", "sub", "mutation_subtype"]:
        df_cat[col] = None
        
    positions = []
    seq_len = wt_len if is_wt_mode else tpe_len
    peg_int = peg_intervals["wt"] if is_wt_mode else peg_intervals["tpe"]
    
    for idx, allele in df_cat.iterrows():
        ins_pos = allele.wt_all_ins_pos if is_wt_mode else allele.tpe_all_ins_pos
        del_pos = allele.wt_all_del_pos if is_wt_mode else allele.tpe_all_del_pos
        sub_pos = allele.wt_all_sub_pos if is_wt_mode else allele.tpe_all_sub_pos

        f_ins, f_del, f_sub = filter_out_of_bounds_positions(
            ins_pos, del_pos, sub_pos, seq_len, peg_int, ignore_del
        )

        has_ins, has_del, has_sub = bool(f_ins), bool(f_del), bool(f_sub)
        
        if has_ins: df_cat.at[idx, "ins"] = True
        if has_del: df_cat.at[idx, "del"] = True
        if has_sub: df_cat.at[idx, "sub"] = True

        subs = []
        if has_ins: subs.append("Ins")
        if has_del: subs.append("Del")
        if has_sub: subs.append("Sub")

        if subs:
            df_cat.at[idx, "mutation_subtype"] = f"{prefix} + " + " + ".join(subs)

        positions.append((f_ins, f_del, f_sub))

    return df_cat, tally_mutation_positions(positions)


def classify_flap(df, rt_template, prefix, is_flap_b=False):
    """Determines if a flap is 'Full' or 'Partial' based on the template match array."""
    df["flap_subtype"] = None
    valid_bases = {"T", "TI", "=T", "=TI"}

    if isinstance(rt_template, tuple):
        dual_flap = True
        rt_len_a = len(rt_template[0])
        rt_len_b = len(rt_template[1])
    else:
        dual_flap = False
        rt_len = len(rt_template)
    
    for idx, allele in df.iterrows():
        if dual_flap:
            if allele.Classified_by == "Composite_A&B":
                arr_a = allele.comp_a_ins_match_arr
                arr_b = allele.comp_b_ins_match_arr

                if not arr_a or not arr_b:
                    continue

                seg_a = arr_b[:rt_len_a]
                seg_b = arr_a[-rt_len_b:]
            else:
                if allele.Classified_by == "TPE":
                    arr = allele.tpe_ins_match_arr
                elif allele.Classified_by == "Composite_A":
                    arr = allele.comp_a_ins_match_arr
                elif allele.Classified_by == "Composite_B":
                    arr = allele.comp_b_ins_match_arr
                else:
                    arr = []

                if not arr:
                    continue

                seg_a = arr[:rt_len_a]
                seg_b = arr[-rt_len_b:]

            is_full = (
                len(seg_a) == rt_len_a
                and len(seg_b) == rt_len_b
                and all(base in valid_bases for base in seg_a)
                and all(base in valid_bases for base in seg_b)
            )

        else:
            if allele.Classified_by == "TPE":
                arr = allele.tpe_ins_match_arr
            elif allele.Classified_by == "Composite_A":
                arr = allele.comp_a_ins_match_arr
            elif allele.Classified_by == "Composite_B":
                arr = allele.comp_b_ins_match_arr
            else:
                arr = []

            if not arr:
                continue

            segment = arr[-rt_len:] if is_flap_b else arr[:rt_len]
            is_full = (
                len(segment) == rt_len
                and all(base in valid_bases for base in segment)
            )

        df.at[idx, "flap_subtype"] = (
            f"Full {prefix}" if is_full else f"Partial {prefix}"
        )
        
    return df


def mutation_analysis(wt_seq_len, tpe_seq_len, df_categorized=None, twinspector_results_folder=None, rt_template_a=None, rt_template_b=None, pegRNA_intervals=None, ignore_extraspacer_deletions=False, recoding_mode=False):
    
    # Configuration for all mutation categories (Category_final string, Prefix for dicts, is_wt_mode boolean)
    category_configs = [
        ("Perfect_TPE", "Perfect TPE", False),
        # ("TPE_Indel", "TPE Indel", False),
        ("Dual_Flap", "Dual Flap", False),
        ("Flap_A", "Flap A", False),
        ("Flap_A_Hybrid", "Flap A Hybrid", False),
        ("Flap_B", "Flap B", False),
        ("Flap_B_Hybrid", "Flap B Hybrid", False),
        ("Imperfect_TPE", "Imperfect TPE", False),
        ("Null", "Null", False),
        ("Imperfect_WT", "Imperfect WT", True),
        # ("WT_Indel", "WT Indel", True),
        ("WT", "WT", True),
    ]

    dfs_processed = {}
    tally_dicts = {}
    ins_del_sub_counts_dict = {}
    mutations_subtype_counts_dict = {}

    # Process all categories through a single loop
    for cat_name, prefix, is_wt in category_configs:
        df, tally = process_mutation_category(
            df_categorized, cat_name, prefix, is_wt, 
            wt_seq_len, tpe_seq_len, pegRNA_intervals, ignore_extraspacer_deletions
        )
        dfs_processed[cat_name] = df
        tally_dicts[cat_name] = tally

        # Aggregate counts iteratively
        total_reads = df["#Reads"].sum()
        ins_reads = df[df["ins"] == True]["#Reads"].sum()
        del_reads = df[df["del"] == True]["#Reads"].sum()
        sub_reads = df[df["sub"] == True]["#Reads"].sum()

        ins_del_sub_counts_dict.update({
            f"{prefix} Total": total_reads,
            f"{prefix} + Ins": ins_reads,
            f"{prefix} + Del": del_reads,
            f"{prefix} + Sub": sub_reads
        })
        
        mutations_subtype_counts_dict.update({
            f"{prefix} Total": total_reads,
            f"{prefix} + Ins": df[df["mutation_subtype"] == f"{prefix} + Ins"]["#Reads"].sum(),
            f"{prefix} + Del": df[df["mutation_subtype"] == f"{prefix} + Del"]["#Reads"].sum(),
            f"{prefix} + Sub": df[df["mutation_subtype"] == f"{prefix} + Sub"]["#Reads"].sum(),
            f"{prefix} + Ins + Del": df[df["mutation_subtype"] == f"{prefix} + Ins + Del"]["#Reads"].sum(),
            f"{prefix} + Ins + Sub": df[df["mutation_subtype"] == f"{prefix} + Ins + Sub"]["#Reads"].sum(),
            f"{prefix} + Del + Sub": df[df["mutation_subtype"] == f"{prefix} + Del + Sub"]["#Reads"].sum(),
            f"{prefix} + Ins + Del + Sub": df[df["mutation_subtype"] == f"{prefix} + Ins + Del + Sub"]["#Reads"].sum(),
        })

    # Combine tally dictionaries
    tpe_categories = ["Perfect_TPE", "Dual_Flap", "Flap_A", "Flap_B", "Flap_A_Hybrid", "Flap_B_Hybrid", "Imperfect_TPE", "Null"]  # "TPE_Indel", 
    wt_categories = ["Imperfect_WT", "WT"]  # "WT_Indel", 
    
    tpe_dicts = [tally_dicts[cat] for cat in tpe_categories]
    wt_dicts = [tally_dicts[cat] for cat in wt_categories]

    if recoding_mode:
        all_tpe_aligned = None
        all_wt_aligned = None
        all_combined = combine_mutation_position_dicts(tpe_dicts + wt_dicts)
    else:
        all_tpe_aligned = combine_mutation_position_dicts(tpe_dicts)
        all_wt_aligned = combine_mutation_position_dicts(wt_dicts)
        all_combined = None

    # Flap Analysis
    flap_completion_counts_dict = {}
    if rt_template_a and rt_template_b:
        flap_configs = [
            ("Dual_Flap", (rt_template_a, rt_template_b), False),
            ("Flap_A", rt_template_a, False),
            ("Flap_A_Hybrid", rt_template_a, False),
            ("Flap_B", rt_template_b, True),
            ("Flap_B_Hybrid", rt_template_b, True)
        ]
        
        for flap_cat, rt_temp, is_b in flap_configs:
            df = classify_flap(dfs_processed[flap_cat], rt_temp, flap_cat, is_flap_b=is_b)
            
            total = df["#Reads"].sum()
            full_reads = df[df["flap_subtype"] == f"Full {flap_cat}"]["#Reads"].sum()
            partial_reads = df[df["flap_subtype"] == f"Partial {flap_cat}"]["#Reads"].sum()
            
            flap_completion_counts_dict.update({
                f"{flap_cat} Total": total,
                f"Full {flap_cat}": full_reads,
                f"Partial {flap_cat}": partial_reads
            })

    with open(f"{twinspector_results_folder}/d4.ins_del_sub_counts.txt", "w") as fout:
        fout.write("\t".join(ins_del_sub_counts_dict.keys()) + "\n")
        fout.write("\t".join(map(str, ins_del_sub_counts_dict.values())) + "\n")

    with open(f"{twinspector_results_folder}/d5.mutations_subtype_counts.txt", "w") as fout:
        fout.write("\t".join(mutations_subtype_counts_dict.keys()) + "\n")
        fout.write("\t".join(map(str, mutations_subtype_counts_dict.values())) + "\n")

    if rt_template_a and rt_template_b:
        with open(f"{twinspector_results_folder}/d3.flap_completion_counts.txt", "w") as fout:
            fout.write("\t".join(flap_completion_counts_dict.keys()) + "\n")
            fout.write("\t".join(map(str, flap_completion_counts_dict.values())) + "\n")

    with open(f"{twinspector_results_folder}/d6.ins_del_sub_positions.txt", "w") as fout:
        fout.write("alignment\tcategory\tmut_type\tposition\tcount\n")
        alignment_groups = [
            ("TPE", tpe_categories),
            ("WT", wt_categories)
        ]
        for alignment_label, categories in alignment_groups:
            for cat in categories:
                cat_data = tally_dicts[cat][0] if isinstance(tally_dicts[cat], list) else tally_dicts[cat]
                for mut_type, positions in cat_data.items():
                    for position, count in positions.items():
                        fout.write(f"{alignment_label}\t{cat}\t{mut_type}\t{position}\t{count}\n")

    return {
        "ins_del_sub_counts_dict": ins_del_sub_counts_dict,
        "mutations_subtype_counts_dict": mutations_subtype_counts_dict,
        "perfect_tpe_mutation_position_counts_dict": tally_dicts["Perfect_TPE"], 
        # "tpe_indel_mutation_position_counts_dict": tally_dicts["TPE_Indel"], 
        "dual_flap_mutation_position_counts_dict": tally_dicts["Dual_Flap"], 
        "flap_a_only_mutation_position_counts_dict": tally_dicts["Flap_A"], 
        "flap_b_only_mutation_position_counts_dict": tally_dicts["Flap_B"], 
        "flap_a_hybrid_mutation_position_counts_dict": tally_dicts["Flap_A_Hybrid"], 
        "flap_b_hybrid_mutation_position_counts_dict": tally_dicts["Flap_B_Hybrid"], 
        "null_mutation_position_counts_dict": tally_dicts["Null"],
        "aberrant_tpe_mutation_position_counts_dict": tally_dicts["Imperfect_TPE"], 
        "imperfect_wt_mutation_position_counts_dict": tally_dicts["Imperfect_WT"], 
        # "wt_indel_mutation_position_counts_dict": tally_dicts["WT_Indel"], 
        "wt_mutation_position_counts_dict": tally_dicts["WT"], 
        "all_wt_aligned_mutation_position_counts_dict": all_wt_aligned, 
        "all_tpe_aligned_mutation_position_counts_dict": all_tpe_aligned, 
        "all_mutation_position_counts_dict": all_combined, 
        "flap_completion_counts_dict": flap_completion_counts_dict,
    }


#### Plotting Functions ####
def get_plotting_stats(df=None, reference_info=None, bp_changes_arrs=None,twinspector_results_folder=None, recoding_mode=False, max_n_alleles_to_write=50):
          
    # Stats for summary barplots
    perfect_tpe_count = df[df["Category_final"] == "Perfect_TPE"]["#Reads"].sum()
    # tpe_indel_count = df[df["Category_final"] == "TPE_Indel"]["#Reads"].sum()
    dual_flap_count = df[df["Category_final"] == "Dual_Flap"]["#Reads"].sum()
    flap_a_only_count = df[df["Category_final"] == "Flap_A"]["#Reads"].sum()
    flap_b_only_count = df[df["Category_final"] == "Flap_B"]["#Reads"].sum()
    flap_a_hybrid_count = df[df["Category_final"] == "Flap_A_Hybrid"]["#Reads"].sum()
    flap_b_hybrid_count = df[df["Category_final"] == "Flap_B_Hybrid"]["#Reads"].sum()
    aberrant_tpe_count = df[df["Category_final"] == "Imperfect_TPE"]["#Reads"].sum()
    null_count = df[df["Category_final"] == "Null"]["#Reads"].sum()
    imperfect_wt_count = df[df["Category_final"] == "Imperfect_WT"]["#Reads"].sum()
    # wt_indel_count = df[df["Category_final"] == "WT_Indel"]["#Reads"].sum()
    wt_count = df[df["Category_final"] == "WT"]["#Reads"].sum()
    uncategorized_count = df[df["Category_final"] == "Uncategorized"]["#Reads"].sum()

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

    # Build arrays for contiguous flap base position integration/removal
    total_read_bases = 0
    all_base_integration_counts_arr = [0] * ins_region_len
    from_flap_b_contiguous_base_integration_counts = [0] * ins_region_len
    from_flap_a_contiguous_base_integration_counts = [0] * ins_region_len
    all_base_removal_counts_arr = [0] * del_region_len
    from_flap_b_contiguous_base_removal_counts = [0] * del_region_len
    from_flap_a_contiguous_base_removal_counts = [0] * del_region_len

    # build arrays for base position integration by category
    cat_perfect_tpe_base_integration_counts_arr = [0] * ins_region_len
    cat_tpe_indel_base_integration_counts_arr = [0] * ins_region_len
    cat_dual_flap_base_integration_counts_arr = [0] * ins_region_len
    cat_flap_a_only_base_integration_counts_arr = [0] * ins_region_len
    cat_flap_b_only_base_integration_counts_arr = [0] * ins_region_len
    cat_flap_a_hybrid_base_integration_counts_arr = [0] * ins_region_len
    cat_flap_b_hybrid_base_integration_counts_arr = [0] * ins_region_len
    cat_aberrant_tpe_base_integration_counts_arr = [0] * ins_region_len
    cat_null_base_integration_counts_arr = [0] * ins_region_len
    cat_imperfect_wt_base_integration_counts_arr = [0] * ins_region_len
    # cat_wt_indel_base_integration_counts_arr = [0] * ins_region_len
    cat_wt_base_integration_counts_arr = [0] * ins_region_len
    cat_uncategorized_base_integration_counts_arr = [0] * ins_region_len

    for idx, allele in df.iterrows():

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

        # Build arrays for contiguous flap base position integration
        if ins_match_arr:
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match in {"T", "TI", "=T", "=TI"}:
                    all_base_integration_counts_arr[pos_idx] += allele['#Reads']

            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if match in {"T", "TI", "=T", "=TI"}:
                    from_flap_a_contiguous_base_integration_counts[pos_idx] += allele['#Reads']
                else:
                    break

            for pos_idx, match in zip(reversed(range(len(ins_match_arr))), reversed(ins_match_arr)):
                if match in {"T", "TI", "=T", "=TI"}:
                    from_flap_b_contiguous_base_integration_counts[pos_idx] += allele['#Reads']
                else:
                    break

        # Build arrays for contiguous flap base position removal
        if del_match_arr:    
            for pos_idx, match in zip(range(len(del_match_arr)), del_match_arr):
                if match in {"T", "TI", "=T", "=TI"}:
                    all_base_removal_counts_arr[pos_idx] += allele['#Reads']

            for pos_idx, match in zip(range(len(del_match_arr)), del_match_arr):
                if match in {"T", "TI", "=T", "=TI"}:
                    from_flap_a_contiguous_base_removal_counts[pos_idx] += allele['#Reads']
                else:
                    break

            for pos_idx, match in zip(reversed(range(len(del_match_arr))), reversed(del_match_arr)):
                if match in {"T", "TI", "=T", "=TI"}:
                    from_flap_b_contiguous_base_removal_counts[pos_idx] += allele['#Reads']
                else:
                    break

        # Build arrays for base position edits by category
        if ins_match_arr:
            for pos_idx, match in zip(range(len(ins_match_arr)), ins_match_arr):
                if allele.Category_final == "Perfect_TPE" and match in {"T", "TI", "=T", "=TI"}:
                    cat_perfect_tpe_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.TPE_Indel == True and match in {"T", "TI", "=T", "=TI"}:
                    cat_tpe_indel_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Dual_Flap" and match in {"T", "TI", "=T", "=TI"}:
                    cat_dual_flap_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Flap_A" and match in {"T", "TI", "=T", "=TI"}:
                    cat_flap_a_only_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Flap_B" and match in {"T", "TI", "=T", "=TI"}:
                    cat_flap_b_only_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Flap_A_Hybrid" and match in {"T", "TI", "=T", "=TI"}:
                    cat_flap_a_hybrid_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Flap_B_Hybrid" and match in {"T", "TI", "=T", "=TI"}:
                    cat_flap_b_hybrid_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Imperfect_TPE" and match in {"T", "TI", "=T", "=TI"}:
                    cat_aberrant_tpe_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Null" and match in {"T", "TI", "=T", "=TI"}:
                    cat_null_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Imperfect_WT" and match in {"T", "TI", "=T", "=TI"}:
                    cat_imperfect_wt_base_integration_counts_arr[pos_idx] += allele['#Reads']

                # if allele.Category_final == "WT_Indel" and match in {"T", "TI", "=T", "=TI"}:
                #     cat_wt_indel_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "WT" and match in {"T", "TI", "=T", "=TI"}:
                    cat_wt_base_integration_counts_arr[pos_idx] += allele['#Reads']

                if allele.Category_final == "Uncategorized" and match in {"T", "TI", "=T", "=TI"}:
                    cat_uncategorized_base_integration_counts_arr[pos_idx] += allele['#Reads']

        total_read_bases += allele["#Reads"]

    total_read_bases_ins_region_arr = [total_read_bases] * ins_region_len
    total_read_bases_del_region_arr = [total_read_bases] * del_region_len
    # 
    perfect_base_removal_counts_arr = [perfect_tpe_count] * del_region_len

    # Add reads classified by TPE reference as these have fully deleted the wt sequence but do not have a del_match_arr and thus are not included
    classified_by_tpe_ref_count = df[df['Classified_by'] == 'TPE']['#Reads'].sum()
    all_base_removal_counts_arr = [x + classified_by_tpe_ref_count for x in all_base_removal_counts_arr]
    from_flap_a_contiguous_base_removal_counts = [x + classified_by_tpe_ref_count for x in from_flap_a_contiguous_base_removal_counts]
    from_flap_b_contiguous_base_removal_counts = [x + classified_by_tpe_ref_count for x in from_flap_b_contiguous_base_removal_counts]

    # count arrays for base edits by category plot
    category_count_arrs = {
        "total_read_bases_ins_region_arr": total_read_bases_ins_region_arr, 
        "Perfect TPE": cat_perfect_tpe_base_integration_counts_arr,
        "TPE Indel": cat_tpe_indel_base_integration_counts_arr, 
        "Dual Flap": cat_dual_flap_base_integration_counts_arr, 
        "Flap A": cat_flap_a_only_base_integration_counts_arr,
        "Flap B": cat_flap_b_only_base_integration_counts_arr, 
        "Flap A Hybrid": cat_flap_a_hybrid_base_integration_counts_arr, 
        "Flap B Hybrid": cat_flap_b_hybrid_base_integration_counts_arr, 
        "Imperfect TPE": cat_aberrant_tpe_base_integration_counts_arr, 
        "Null": cat_null_base_integration_counts_arr, 
        "Imperfect WT": cat_imperfect_wt_base_integration_counts_arr,
        # "WT Indel": cat_wt_indel_base_integration_counts_arr,
        "WT": cat_wt_base_integration_counts_arr,
        "Uncategorized": cat_uncategorized_base_integration_counts_arr
    }

    category_counts = {
        "Perfect TPE": perfect_tpe_count,
        # "TPE Indel": tpe_indel_count,
        "Dual Flap": dual_flap_count,
        "Flap A": flap_a_only_count,
        "Flap B": flap_b_only_count,
        "Flap A Hybrid": flap_a_hybrid_count,
        "Flap B Hybrid": flap_b_hybrid_count,
        "Imperfect TPE": aberrant_tpe_count,
        "Null": null_count,
        "Imperfect WT": imperfect_wt_count,
        # "WT Indel": wt_indel_count,
        "WT": wt_count,
        "Uncategorized": uncategorized_count,
    }

    stats = {
        "category_counts": category_counts,
        "total_read_bases_ins_region_arr": total_read_bases_ins_region_arr, 
        "all_base_integration_counts_arr": all_base_integration_counts_arr, 
        "from_flap_a_contiguous_base_integration_counts": from_flap_a_contiguous_base_integration_counts, 
        "from_flap_b_contiguous_base_integration_counts": from_flap_b_contiguous_base_integration_counts,
        "perfect_tpe_counts_arr": cat_perfect_tpe_base_integration_counts_arr, 
        "tpe_indel_counts_arr": cat_aberrant_tpe_base_integration_counts_arr, 
        # "tpe_indel_counts_arr": cat_tpe_indel_base_integration_counts_arr, 
        "dual_flap_counts_arr": cat_dual_flap_base_integration_counts_arr, 
        "inserted_seq": inserted_seq, 
        "total_read_bases_del_region_arr": total_read_bases_del_region_arr,
        "perfect_base_removal_counts_arr": perfect_base_removal_counts_arr,
        "all_base_removal_counts_arr": all_base_removal_counts_arr,
        "from_flap_a_contiguous_base_removal_counts": from_flap_a_contiguous_base_removal_counts, 
        "from_flap_b_contiguous_base_removal_counts": from_flap_b_contiguous_base_removal_counts,
        "deleted_seq": deleted_seq, 
        "category_count_arrs": category_count_arrs
    }
    
    # Write to files
    with open(twinspector_results_folder + "/d1.category_counts.txt", "w") as fout:
        fout.write("\t".join(CATEGORY_ORDER + ["Uncategorized"]) + "\n")
        fout.write("\t".join([str(category_counts[cat]) for cat in CATEGORY_ORDER + ["Uncategorized"]]) + "\n")

    with open(twinspector_results_folder + "/d7.top_alleles_by_category.txt", "w") as fout:
        for c in CATEGORY_ORDER[::-1] + ["Uncategorized"]:
            fc = c.replace(" ", "_")
            if fc not in df['Category_final'].values:
                continue
            fout.write(f"Category: {c}, Total Reads: {category_counts[c]}\n\n")
            for idx, row in df[df['Category_final'] == fc].sort_values(by='#Reads', ascending=False).head(max_n_alleles_to_write).iterrows():
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
            fout.write("\n\n")

    with open(twinspector_results_folder + "/d2.base_counts.txt", "w") as fout:
        fout.write("inserted_tpe_bases\t" + "\t".join([str(x) for x in inserted_seq]) + "\n")
        # fout.write("all_base_integration_counts\t" + '\t'.join([str(x) for x in all_base_integration_counts_arr]) + "\n")
        fout.write("from_flap_b_contiguous_base_integration_counts\t" + '\t'.join([str(x) for x in from_flap_b_contiguous_base_integration_counts]) + "\n")
        fout.write("from_flap_a_contiguous_base_integration_counts\t" + '\t'.join([str(x) for x in from_flap_a_contiguous_base_integration_counts]) + "\n")
        fout.write("full_tpe_+_mutations\t" + '\t'.join([str(x) for x in cat_tpe_indel_base_integration_counts_arr]) + "\n")
        fout.write("perfectly_edited_base_counts\t" + '\t'.join([str(x) for x in cat_perfect_tpe_base_integration_counts_arr]) + "\n")
        fout.write("total_read_counts\t" + '\t'.join([str(x) for x in total_read_bases_ins_region_arr]) + "\n\n")
        fout.write("removed_wt_bases\t" + '\t'.join([str(x) for x in deleted_seq]) + "\n")
        # fout.write("all_base_removal_counts_arr\t" + '\t'.join([str(x) for x in all_base_removal_counts_arr]) + "\n")
        fout.write("from_flap_b_contiguous_base_removal_counts\t" + '\t'.join([str(x) for x in from_flap_b_contiguous_base_removal_counts]) + "\n")
        fout.write("from_flap_a_contiguous_base_removal_counts\t" + '\t'.join([str(x) for x in from_flap_a_contiguous_base_removal_counts]) + "\n")
        fout.write("perfectly_edited_base_counts\t" + '\t'.join([str(x) for x in perfect_base_removal_counts_arr]) + "\n")
        fout.write("total_read_counts\t" + '\t'.join([str(x) for x in total_read_bases_del_region_arr]) + "\n")

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


def save_plot(filename, plot_formats, fig_root=None, fig=None, **kwargs):
    for ext in plot_formats:
        save_kwargs = kwargs.copy()
        if ext == "png":
            save_kwargs.setdefault("dpi", 300)
        if fig:
            fig.savefig(f"{fig_root}/{filename}.{ext}", **save_kwargs)
        else:
            plt.savefig(f"{fig_root}/{filename}.{ext}", **save_kwargs)
    if fig:
        plt.close(fig)
    else:
        plt.close()


#### Summary barplots ####
def plot_reads_input_summary_barplot(crispresso_output_folder, counts_dict, fig_root=None, plot_formats=False):

    crispresso_mapping_statistics_file = os.path.join(crispresso_output_folder, 'CRISPResso_mapping_statistics.txt')
    read_stats = pd.read_csv(crispresso_mapping_statistics_file, sep="\t")
    # Update Reads Aligned count for post-CRISPResso homology filtering
    num_input = read_stats['READS IN INPUTS'][0]
    num_after_preprocessing = read_stats['READS AFTER PREPROCESSING'][0]
    num_analyzed = sum(counts_dict.values())
    # num_discarded = read_stats['READS IN INPUTS'][0] - sum(counts_dict.values())
    counts = [num_input, num_after_preprocessing, num_analyzed] #, num_discarded]  # read_stats['READS ALIGNED'][0]]
    labels = ["Input", "After Preprocessing", "Analyzed"]  #, "Discarded"]
    total = read_stats['READS IN INPUTS'][0]

    sorted_pairs = sorted(zip(labels, counts), key=lambda x: x[1], reverse=True)
    sorted_labels = [p[0] for p in sorted_pairs]
    sorted_values = [p[1] for p in sorted_pairs]

    percent_labels = [f"{lab}\n({val/total*100:.1f}%)" for lab, val in sorted_pairs]

    width = 0.6

    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

    bars = ax.bar(sorted_labels, sorted_values, width,  
                color="lightgrey", edgecolor="white", linewidth=0.4)

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
        
    ax.set_ylabel("Reads", fontsize=10)

    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(percent_labels, rotation=0, ha="center")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ax.minorticks_on()
    # ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    # ax.tick_params(axis="x", which="minor", bottom=False)

    plt.tight_layout()

    save_plot("a1.Reads_input", plot_formats, fig_root, fig, bbox_inches="tight")


def plot_categorized_stacked_barplot(counts_dict, fig_root=None, plot_formats=False, category_colors=None, categorized=False):

    total = sum(counts_dict.values())
    if categorized:
        aligned_labels = list(counts_dict.keys())
        aligned_counts = list(counts_dict.values())
        plot_order = CATEGORY_ORDER
    else:
        perfect_tpe = counts_dict["Perfect TPE"]
        wt = counts_dict["WT"]
        aligned_labels = ["Perfect TPE", "Other", "WT"]
        aligned_counts = [perfect_tpe, total - perfect_tpe - wt, wt]
        plot_order = ["Perfect TPE", "Other", "WT"]

    count_dict = dict(zip(aligned_labels, aligned_counts))
    sorted_counts_labels = [
        (label, count_dict.get(label, 0))
        for label in plot_order
    ]

    legend_labels = [f"{lab} ({val/total*100:.1f}%)" for lab, val in sorted_counts_labels]
    
    fig, ax = plt.subplots(figsize=(1, 6), dpi=300)

    x = [0]
    bottom = 0
    for (lab, val), legend_label in zip(sorted_counts_labels, legend_labels):
        color = category_colors.get(lab, "#5f3353")
        ax.bar(
            x, val,   
            bottom=bottom, 
            label=legend_label, 
            color=color, 
            edgecolor='white',
            linewidth=.2
        )
        bottom += val

    ax.text(x[0], total, f"{total:,}", ha='center', va='bottom', fontsize=8)

    handles, labels = ax.get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    ax.legend(
        handles, labels,
        bbox_to_anchor=(1.05, 0.5),
        loc="center left",
        borderaxespad=0.25,
        fontsize=8
    )

    def double_label_formatter(y, pos):
        pct = (y / total * 100) if total > 0 else 0
        return f"{int(y):,} ({pct:.1f}%)"

    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(double_label_formatter))

    ax.set_xticks([0], labels=["Analyzed"], fontsize=8)
    ax.set_ylabel("Reads", fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.minorticks_on()
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    ax.tick_params(axis="y", labelsize=8)

    file_name = "a3.Outcomes_categorized_stacked" if categorized else "a2.Outcomes_stacked"

    save_plot(file_name, plot_formats, fig_root, fig, bbox_inches="tight")


def plot_categorized_barplot(counts_dict, fig_root=None, plot_formats=False, category_colors=None, categorized=False):

    total = sum(counts_dict.values())

    if categorized:
        plot_order = CATEGORY_ORDER
        sorted_pairs = [
            (label, counts_dict.get(label, 0))
            for label in reversed(plot_order)
            if label in counts_dict
        ]
    else:
        perfect_tpe = counts_dict.get("Perfect TPE", 0)
        wt = counts_dict.get("WT", 0)
        other_count = total - perfect_tpe - wt
        
        plot_order = ["WT", "Other", "Perfect TPE"]  # Kept in reversed order to match your original logic
        sorted_pairs = [
            ("WT", wt),
            ("Other", other_count),
            ("Perfect TPE", perfect_tpe)
        ]

    sorted_labels = [p[0] for p in sorted_pairs]
    sorted_values = [p[1] for p in sorted_pairs]

    percent_labels = [
        f"{lab}\n({val/total*100:.1f}%)" if total > 0 else f"{lab}\n(0.0%)" 
        for lab, val in sorted_pairs
    ]

    colors = [category_colors.get(lab, "#5f3353") for lab in sorted_labels]

    fig, ax = plt.subplots(figsize=(8 if categorized else 3, 4), dpi=300)

    bars = ax.bar(sorted_labels, sorted_values, 
                color=colors, edgecolor="white", linewidth=0.6)
    
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.,
                height,
                f"{int(height):,}",
                ha='center',
                va='bottom',
                fontsize=8,
            )

    ax.set_ylabel("Reads", fontsize=8)
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

    file_name = "a5.Outcomes_categorized" if categorized else "a4.Outcomes"

    save_plot(file_name, plot_formats, fig_root, fig, bbox_inches="tight")


def plot_summary_barplots(category_counts, crispresso_output_folder, twinspector_results_folder, crispresso_wt, plot_formats):

    setBarMatplotlibDefaults()
  
    plot_reads_input_summary_barplot(
        crispresso_output_folder,
        category_counts,
        fig_root=twinspector_results_folder,
        plot_formats=plot_formats
    )

    plot_categorized_stacked_barplot(
        category_counts, 
        fig_root=twinspector_results_folder, 
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS, 
    )

    plot_categorized_stacked_barplot(
        category_counts, 
        fig_root=twinspector_results_folder, 
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS, 
        categorized=True
    )

    plot_categorized_barplot(
        category_counts,
        fig_root=twinspector_results_folder, 
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS
    )

    plot_categorized_barplot(
        category_counts,
        fig_root=twinspector_results_folder, 
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS, 
        categorized=True
    )


#### Mutation barplots ####
def plot_ins_del_sub_by_category_barplot(counts_dict, category_counts, twinspector_results_folder=None, plot_formats=False, plot_percentages=False):

    labels = list(reversed(CATEGORY_ORDER))
    # labels = list(["WT", "Imperfect\nWT", "Null", "Imperfect\nTPE", "Flap B\nHybrid", "Flap A\nHybrid", "Flap B", "Flap A", "Dual\nFlap", "Perfect\nTPE"])
    
    # Sum counts for each category in category_counts dict
    totals = [category_counts.get(lab, 0) for lab in labels]
    ins_counts = [counts_dict[f"{lab} + Ins"] for lab in labels]
    del_counts = [counts_dict[f"{lab} + Del"] for lab in labels]
    sub_counts = [counts_dict[f"{lab} + Sub"] for lab in labels]

    grand_total = sum(totals)

    labels.append("Total")
    totals.append(grand_total)

    ins_counts.append(sum(ins_counts))
    del_counts.append(sum(del_counts))
    sub_counts.append(sum(sub_counts))

    ins_pcts = [(i / grand_total * 100) if grand_total > 0 else 0 for i in ins_counts]
    del_pcts = [(d / grand_total * 100) if grand_total > 0 else 0 for d in del_counts]
    sub_pcts = [(s / grand_total * 100) if grand_total > 0 else 0 for s in sub_counts]

    # ins_pcts = [(i / t * 100) if t > 0 else 0 for i, t in zip(ins_counts, totals)]
    # del_pcts = [(d / t * 100) if t > 0 else 0 for d, t in zip(del_counts, totals)]
    # sub_pcts = [(s / t * 100) if t > 0 else 0 for s, t in zip(sub_counts, totals)]

    # Conditionally select the data and axis label based on the flag
    if plot_percentages:
        data_ins, data_del, data_sub = ins_pcts, del_pcts, sub_pcts
    else:
        data_ins, data_del, data_sub = ins_counts, del_counts, sub_counts
    y_label = "Reads"

    # count_labels = []
    # for lab, val in zip(labels, totals):
    #     if lab == "All":
    #         cat_share = (val / grand_total * 100) if grand_total > 0 else 0
    #         # count_labels.append(f"{lab}\n{val:,} reads\n({cat_share:.1f}%)")
    #         count_labels.append(f"{lab}\n({cat_share:.1f}%)")
    #     else:
    #         cat_share = (val / grand_total * 100) if grand_total > 0 else 0
    #         cat_share = (val / grand_total * 100) if grand_total > 0 else 0
    #         # count_labels.append(f"{lab}\n{val:,}\n({cat_share:.1f})")
    #         count_labels.append(f"{lab}\n({cat_share:.1f}%)")

    fig, ax = plt.subplots(figsize=(5, 4))

    x = np.arange(len(labels))
    width = 0.29 

    colors = {
        "Ins": "#C29B38",
        "Del": "#C67D3A",
        "Sub": "#6B8E63"
    }

    ax.axvline(x=len(labels) - 1.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    bars_ins = ax.bar(x - width, data_ins, width, color=colors["Ins"], edgecolor="white", linewidth=0.6, label="Ins", alpha=1.0)
    bars_del = ax.bar(x,         data_del, width, color=colors["Del"], edgecolor="white", linewidth=0.6, label="Del", alpha=1.0)
    bars_sub = ax.bar(x + width, data_sub, width, color=colors["Sub"], edgecolor="white", linewidth=0.6, label="Sub", alpha=1.0)

    max_height = max(data_ins + data_del + data_sub) if any(data_ins + data_del + data_sub) else 1
    y_offset = max_height * 0.02

    # Custom formatter function to make percentages primary
    def y_formatter(y, pos):
        if plot_percentages:
            # If plotting percentages, 'y' is already a percentage. Just add the % sign.
            count = (y / 100) * grand_total if grand_total > 0 else 0
            return f"{y:.1f}% ({int(count):,})"
        else:
            # If plotting counts, calculate the percentage and put it first.
            pct = (y / grand_total * 100) if grand_total > 0 else 0
            return f"{pct:.1f}% ({int(y):,})"

    # Apply the formatter to the y-axis
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(y_formatter))

    ax.set_ylabel(y_label, fontsize=10)
    # ax.spines['left'].set_visible(False)
    # ax.tick_params(axis="y", left=False, labelleft=False)
    # ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.tick_params(axis="y", which="major", labelsize=8)
    ax.tick_params(axis="y", which="minor", left=True, length=2, color="black")

    ax.set_xticks(x)

    # Map the original dictionary keys to your desired multi-line display labels
    label_mapping = {
        "WT": "WT",
        "Imperfect WT": "Imp WT",
        "Null": "Null",
        "Imperfect TPE": "Imp TPE",
        "Flap B Hybrid": "Flap B Hyb",
        "Flap A Hybrid": "Flap A Hyb",
        "Flap B": "Flap B",
        "Flap A": "Flap A",
        "Dual Flap": "Dual Flap",
        "Perfect TPE": "Perfect TPE",
        "Total": "Total"
    }
    # Generate the new list of labels solely for the x-axis
    display_labels = [label_mapping.get(lab, lab) for lab in labels]

    # Apply the formatted display labels
    ax.set_xticklabels(display_labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=8)

    # ax.set_xticklabels(count_labels, rotation=70, ha="center", fontsize=8)
    # ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=8)
    ax.set_xlim(-0.5, len(labels) - 0.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis="x", which="minor", bottom=False)
    
    # ax.set_ylim(0, max_height * 1.15)
    # Conditionally set the y-limit so it scales correctly for both counts and percentages
    if plot_percentages:
        # ax.set_ylim(0, 30)
        # ax.set_ylim(0, 25)
        # ax.set_ylim(0, max_height * 1.15)
        ax.set_ylim(0, max_height)
    else:
        ax.set_ylim(0, 42500) # Keeps your original manual limit for raw counts

    # fig.suptitle(f"({grand_total:,} reads)", fontsize=10, y=0.87, x=0.55)
    
    handles, labels_legend = ax.get_legend_handles_labels()
    # ax.legend(handles[::-1], labels_legend[::-1], fontsize=10, frameon=False, loc="upper right", bbox_to_anchor=(1.015, 1.0))
    # ax.legend(handles[::-1], labels_legend[::-1], fontsize=10, frameon=False, loc="upper right", bbox_to_anchor=(1.1, 1.0))
    # ax.legend(handles[::-1], labels_legend[::-1], fontsize=10, frameon=False, loc="upper left")
    ax.legend(handles, labels_legend, fontsize=8, frameon=False, loc="upper left")

    plt.tight_layout()

    file_name = "c1.Ins_del_sub_by_category"

    save_plot(file_name, plot_formats, twinspector_results_folder, fig, bbox_inches="tight")


def plot_ins_del_sub_combinations_by_category_stacked_barplot(counts_dict, twinspector_results_folder=None, plot_formats=False):

    expected_labels = list(reversed(CATEGORY_ORDER))
    
    valid_labels = [lab for lab in expected_labels if f"{lab} Total" in counts_dict]
    
    totals = [counts_dict.get(f"{lab} Total", 0) for lab in valid_labels]

    plot_labels = list(valid_labels)
    plot_labels.append("Total")
    totals.append(sum(totals))

    label_mapping = {
        "WT": "WT",
        "Imperfect WT": "Imp WT",
        "Null": "Null",
        "Imperfect TPE": "Imp TPE",
        "Flap B Hybrid": "Flap B Hyb",
        "Flap A Hybrid": "Flap A Hyb",
        "Flap B": "Flap B",
        "Flap A": "Flap A",
        "Dual Flap": "Dual Flap",
        "Perfect TPE": "Perfect TPE",
        "All": "All"
    }
    display_labels = [label_mapping.get(lab, lab) for lab in plot_labels]

    count_labels = [f"{lab}\n({val:,})" for lab, val in zip(display_labels, totals)]
    # count_labels = [f"{lab}\n{val:,} reads\n({val/grand_total*100:.1f}%)" for lab, val in zip(plot_labels, totals)]

    combos = [
        " + Ins", " + Del", " + Sub", 
        " + Ins + Del", " + Ins + Sub", " + Del + Sub", " + Ins + Del + Sub"
    ]
    
    legend_labels = [c.replace(" + ", "", 1) for c in combos]

    colors = [
        "#C29B38",
        "#C67D3A",
        "#6B8E63",
        "#a68c9e",
        "#ccc0b0",
        "#7d7168",
        "#444d4f"
    ]

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    
    x = np.arange(len(plot_labels))
    width = 0.8

    ax.axvline(x=len(plot_labels) - 1.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    bottoms = np.zeros(len(plot_labels))

    for i, combo in enumerate(combos):
        raw_counts = [counts_dict.get(f"{lab}{combo}", 0) for lab in valid_labels]
        
        raw_counts.append(sum(raw_counts))

        pcts = [(count / t * 100) if t > 0 else 0 for count, t in zip(raw_counts, totals)]
        
        bars = ax.bar(x, pcts, width, bottom=bottoms, 
                      color=colors[i], edgecolor="white", linewidth=0.2, 
                      label=legend_labels[i], alpha=1.0)
        
        for j, bar in enumerate(bars):
            height = bar.get_height()
            if height > 4.0:  
                y_pos = bottoms[j] + (height / 2)
                ax.text(bar.get_x() + bar.get_width()/2.,
                        y_pos,
                        f"{height:.1f}%",
                        ha='center',
                        va='center',
                        fontsize=7,
                        color='white' if i in [5, 6] else 'black')
        
        bottoms += pcts

    ax.set_ylabel("Category reads (%)", fontsize=12)
    
    # 1. Increase y-limit to 110 to make room for the totals above the bars
    ax.set_ylim(0, 100) 
    
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.tick_params(axis="y", which="major", length=5, width=1, labelsize=10)
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)

    # 2. Just use the category names for the x-axis
    ax.set_xticks(x)
    ax.set_xticklabels(display_labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
    
    ax.set_xlim(-0.5, len(plot_labels) - 0.5)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis="x", which="minor", bottom=False)
    
    # 3. Add the total read counts just above the 100% mark
    for pos, total in zip(x, totals):
        ax.text(pos, 100.5, f"n={total:,}", 
                ha='center', va='bottom', 
                fontsize=8, color='black', 
                rotation=45) # Rotated to prevent large numbers from overlapping
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], fontsize=10, title_fontsize=12, frameon=False, loc="upper left", bbox_to_anchor=(0.99, 1.0))

    plt.tight_layout()

    file_name = "c2.Ins_del_sub_combinations_by_category"

    save_plot(file_name, plot_formats, twinspector_results_folder, fig, bbox_inches="tight")


def plot_ins_del_sub_positions_barplot(counts_dict, total_reads, ref_seq_len, twinspector_results_folder=None, plot_formats=False, num=None, cat=None, vlines=None, recoding_mode=None):

    display_total_reads = total_reads 
    if total_reads == 0:
        total_reads = 1

    colors = {
        "Ins": "#C29B38",
        "Del": "#C67D3A",
        "Sub": "#6B8E63"
    }

    # hardcoded_ylims = [(0, 0.6), (0, 5.5), (0, 1.3)]
    
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7.5, 5.5), dpi=300, sharex=True)
    
    mutation_types = [("ins", "Ins"), ("del", "Del"), ("sub", "Sub")]
    
    for i, (ax, (mut_key, label)) in enumerate(zip(axes, mutation_types)):
        mut_data = counts_dict.get(mut_key, {})
        
        # Build legend handles list for this subplot
        legend_handles = [patches.Patch(color=colors[label], label=label)]
        
        # Draw vertical grey dashed lines and add 'Nick sites' to legend
        if vlines:
            if ax == axes[0]:
                nick_line_handle = matplotlib.lines.Line2D([], [], color='gray', linestyle='--', linewidth=1, alpha=0.7, label='Nick sites')
                legend_handles.append(nick_line_handle)
            
            for x_pos in vlines:
                ax.axvline(x=x_pos, color='gray', linestyle='--', linewidth=1, alpha=0.7, zorder=6)

        # Draw the updated legend with both handles
        # ax.legend(handles=legend_handles, fontsize=12, frameon=False, loc="upper left", bbox_to_anchor=(.75, .9))
        # ax.legend(handles=legend_handles, fontsize=12, frameon=False, loc="upper left", bbox_to_anchor=(.72, .9))
        # ax.legend(handles=legend_handles, fontsize=12, frameon=False, loc="upper left", bbox_to_anchor=(.82, 1.25))

        if mut_data:
            positions = sorted(mut_data.keys())
            counts = [mut_data[pos] for pos in positions]
            percentages = [(count / total_reads) * 100 for count in counts]
            
            ax.bar(positions, percentages, width=1.0, color=colors[label], 
                   linewidth=0, alpha=1.0)
            
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=3))
            ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        else:
            ax.text(0.5, 0.5, f"No {label} detected", ha='center', va='center', 
                    transform=ax.transAxes, color='gray', fontsize=12, fontstyle='italic')
            
            # # Lock the Y-axis to a default 0-1 scale so the text sits nicely in the middle
            # ax.set_ylim(0, 1)
            # Remove y-ticks for the empty plot to keep it clean
            ax.set_yticks([])

        # ax.set_ylim(hardcoded_ylims[i])

        ax.set_ylabel("Reads (%)", fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ax.tick_params(axis="y", which="major", labelsize=10)
        ax.tick_params(axis="y", which="minor", left=True) 
        ax.tick_params(axis="x", which="minor", bottom=False)

    # Bottom plot X-axis styling
    axes[-1].set_xlabel(
        "WT Reference Position" if "WT" in cat else 
        "Reference Position" if recoding_mode else
        "TPE Reference Position", fontsize=12, labelpad=10
    )
    axes[-1].tick_params(axis="x", labelsize=10)

    # Force x-axis to span reference length
    axes[-1].set_xlim(-0.5, ref_seq_len + 0.5)

    cat_title = cat.replace("_", " ")
    fig.suptitle(f"{cat_title} ({display_total_reads:,} reads)", fontsize=12, y=0.94)
    # fig.suptitle("24bp", fontsize=16, y=0.94)

    plt.tight_layout(rect=[0, 0, 1, 0.98])  # Leaves top 4% of canvas for the title
    # plt.tight_layout()

    file_name = f"c{num}.{cat}.ins_del_sub_positions"

    save_plot(file_name, plot_formats, twinspector_results_folder, fig, bbox_inches="tight")

    plt.close(fig)



def plot_mutation_barplots(mutation_dicts, category_counts, cut_points, wt_seq_len, tpe_seq_len, twinspector_results_folder=None, plot_formats=False, recoding_mode=False):

    plot_ins_del_sub_by_category_barplot(
        mutation_dicts["ins_del_sub_counts_dict"], 
        category_counts, 
        twinspector_results_folder=twinspector_results_folder, 
        plot_formats=plot_formats, 
        plot_percentages=True
    )

    plot_ins_del_sub_combinations_by_category_stacked_barplot(
        mutation_dicts["mutations_subtype_counts_dict"], 
        twinspector_results_folder=twinspector_results_folder, 
        plot_formats=plot_formats
    )

    all_tpe_keys = ["Perfect TPE", "Dual Flap", "Flap A", "Flap B", "Flap A Hybrid", "Flap B Hybrid", "Imperfect TPE", "Null"]  # "TPE Indel", 
    all_wt_keys = ["Imperfect WT", "WT"]  # "WT Indel", 

    if recoding_mode:
        all_cats = {
            "all_mutation_position_counts_dict": ("All", 1, sum(category_counts.values()))
        }
    else:
        all_cats = {
            "all_tpe_aligned_mutation_position_counts_dict": ("All_TPE", 1, sum(category_counts[k] for k in all_tpe_keys)), "all_wt_aligned_mutation_position_counts_dict": ("All_WT", 0, sum(category_counts[k] for k in all_wt_keys))
        }

    other_cats = {
        "perfect_tpe_mutation_position_counts_dict": ("Perfect_TPE", 1, category_counts["Perfect TPE"]),  # "tpe_indel_mutation_position_counts_dict": ("TPE_Indel", 1, category_counts["TPE Indel"]), 
        "dual_flap_mutation_position_counts_dict": ("Dual_Flap", 1, category_counts["Dual Flap"]), "flap_a_only_mutation_position_counts_dict": ("Flap_A", 1, category_counts["Flap A"]), 
        "flap_b_only_mutation_position_counts_dict": ("Flap_B", 1, category_counts["Flap B"]), "flap_a_hybrid_mutation_position_counts_dict": ("Flap_A_Hybrid", 1, category_counts["Flap A Hybrid"]), 
        "flap_b_hybrid_mutation_position_counts_dict": ("Flap_B_Hybrid", 1, category_counts["Flap B Hybrid"]), "aberrant_tpe_mutation_position_counts_dict": ("Imperfect_TPE", 1, category_counts["Imperfect TPE"]), 
        "null_mutation_position_counts_dict": ("Null", 1, category_counts["Null"]), "imperfect_wt_mutation_position_counts_dict": ("Imperfect_WT", 0, category_counts["Imperfect WT"]), 
        "wt_mutation_position_counts_dict": ("WT", 0, category_counts["WT"]),  # "wt_indel_mutation_position_counts_dict": ("WT_Indel", 0, category_counts["WT Indel"]), 
    }

    mutation_position_dicts = {**all_cats, **other_cats}

    for idx, (k, (cat, cps, tr)) in enumerate(mutation_position_dicts.items(), start=3):
        plot_ins_del_sub_positions_barplot(
            mutation_dicts[k], 
            total_reads=tr, 
            ref_seq_len=wt_seq_len if cps == 0 else tpe_seq_len,
            twinspector_results_folder=twinspector_results_folder,
            plot_formats=plot_formats,
            num=str(idx), 
            cat=cat, 
            vlines=cut_points[cps], 
            recoding_mode=recoding_mode
        )


#### Per-base barplots ####
def plot_base_integration_by_category(
    total_counts, 
    category_count_arrs, 
    insert_sequence, 
    rt_template_a_span_end_idx=None, 
    rt_template_b_span_start_idx=None, 
    recoding_mode=False, 
    title=None, 
    fig_root=None, 
    plot_formats=False, 
    category_colors=None, 
):
    has_spans = (rt_template_a_span_end_idx is not None) or (rt_template_b_span_start_idx is not None)

    n = len(insert_sequence)
    indices = np.arange(n)

    physical_bar_width = 0.25 
    physical_gap_width = 0.02 

    width_per_base = physical_bar_width + physical_gap_width
    bar_width = physical_bar_width / width_per_base

    spine_gap = bar_width / 2 + 0.15
    data_range = (n - 1 + spine_gap) - (-spine_gap)
    axes_width = data_range * width_per_base

    min_fig_width = 13
    fig_height = 7.5 if has_spans else 6.25
    y_label_space = 1.0  
    right_margin = 0.5   

    natural_fig_width = axes_width + y_label_space + right_margin
    fig_width = max(min_fig_width, natural_fig_width)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

    extra_padding = fig_width - natural_fig_width
    axes_left = y_label_space + (extra_padding / 2)

    axes_bottom = 0.30 if has_spans else 0.24
    axes_height = 0.66 if has_spans else 0.69

    ax = fig.add_axes([axes_left / fig_width, axes_bottom, axes_width / fig_width, axes_height])
    if title:
        fig.suptitle(title,fontsize=24,y=0.97,ha="center")

    total_reads = max(total_counts)

    perfect_tpe_pct = np.array(category_count_arrs["Perfect TPE"]) / total_reads * 100
    # tpe_indel_pct = np.array(category_count_arrs["TPE Indel"]) / total_reads * 100
    dual_flap_pct = np.array(category_count_arrs["Dual Flap"]) / total_reads * 100
    flap_a_only_pct = np.array(category_count_arrs["Flap A"]) / total_reads * 100
    flap_b_only_pct = np.array(category_count_arrs["Flap B"]) / total_reads * 100
    flap_a_hybrid_pct = np.array(category_count_arrs["Flap A Hybrid"]) / total_reads * 100
    flap_b_hybrid_pct = np.array(category_count_arrs["Flap B Hybrid"]) / total_reads * 100
    aberrant_tpe_pct = np.array(category_count_arrs["Imperfect TPE"]) / total_reads * 100
    null_pct = np.array(category_count_arrs["Null"]) / total_reads * 100
    imperfect_wt_pct = np.array(category_count_arrs["Imperfect WT"]) / total_reads * 100
    # wt_indel_pct = np.array(category_count_arrs["WT Indel"]) / total_reads * 100
    wt_pct = np.array(category_count_arrs["WT"]) / total_reads * 100

    # Stacked bar plot
    ax.bar(indices, perfect_tpe_pct, width=bar_width, label="Perfect TPE", color=category_colors["Perfect TPE"])
    bottom_so_far = perfect_tpe_pct

    # ax.bar(indices, tpe_indel_pct, width=bar_width, label="TPE Indel", color=category_colors["TPE Indel"], bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, tpe_indel_pct)]

    ax.bar(indices, dual_flap_pct, width=bar_width, label="Dual Flap", color=category_colors["Dual Flap"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, dual_flap_pct)]

    ax.bar(indices, flap_a_only_pct, width=bar_width, label="Flap A", color=category_colors["Flap A"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, flap_a_only_pct)]

    ax.bar(indices, flap_b_only_pct, width=bar_width, label="Flap B", color=category_colors["Flap B"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, flap_b_only_pct)]

    ax.bar(indices, flap_a_hybrid_pct, width=bar_width, label="Flap A Hybrid", color=category_colors["Flap A Hybrid"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, flap_a_hybrid_pct)]

    ax.bar(indices, flap_b_hybrid_pct, width=bar_width, label="Flap B Hybrid", color=category_colors["Flap B Hybrid"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, flap_b_hybrid_pct)]

    ax.bar(indices, aberrant_tpe_pct, width=bar_width, label="Imperfect TPE", color=category_colors["Imperfect TPE"], bottom=bottom_so_far)
    bottom_so_far = [x + y for x, y in zip(bottom_so_far, aberrant_tpe_pct)]

    # ax.bar(indices, category_count_arrs["Imperfect WT"], label="Imperfect WT", color=category_colors["Imperfect WT"], bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, category_count_arrs["Imperfect WT"])]

    # ax.bar(indices, category_count_arrs["WT Indel"], label="WT Indel", color=category_colors["WT Indel"], bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, category_count_arrs["WT Indel"])]

    # ax.bar(indices, category_count_arrs["WT"], label="WT", color=category_colors["WT"], bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, category_count_arrs["WT"])]

    # ax.bar(indices, category_count_arrs["Uncategorized"], label="Uncategorized", bottom=bottom_so_far)
    # bottom_so_far = [x + y for x, y in zip(bottom_so_far, category_count_arrs["Uncategorized"])]

    plot_max = max(bottom_so_far) * 1.01

    ax.set_ylabel("Reads (%)", fontsize=20)
    ax.set_xlim(-spine_gap, n - 1 + spine_gap)
    ax.set_ylim(0, plot_max)

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.tick_params(axis='x', length=0)

    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', length=3, width=0.8)
    ax.tick_params(axis='y', which='major', labelsize=20, length=6, width=1)

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
            fontsize=16,
            clip_on=False
        )

    # Span calculations 
    region_text_offset = 0.01 * max_height
    label_text = "Programmed Base Changes" if recoding_mode else "Programmed Sequence"
    
    label_y = y_base - region_text_offset
    ax.text(
        (n - 1) / 2,
        label_y,
        label_text,
        ha="center",
        va="top",
        fontsize=20,
        color="black",
        clip_on=False
    )

    if has_spans:
        region_gap_inches = 0.05
        region_height_inches = 0.24
        arrow_tip_inches = 0.10

        label_height_inches = 0.28 
        
        label_clearance_data = label_height_inches / ax_height_in * max_height
        region_gap_data = region_gap_inches / ax_height_in * max_height
        region_height = region_height_inches / ax_height_in * max_height
        
        ax_pos = ax.get_position()
        ax_width_in = fig_width * ax_pos.width
        arrow_tip_data = (arrow_tip_inches / ax_width_in) * data_range

        y_region_1 = label_y - label_clearance_data - region_height
        y_region_2 = y_region_1 - region_gap_data - region_height

        Path = mpath.Path

        # Flap A span
        if rt_template_a_span_end_idx is not None:
            x1_start = -bar_width / 2
            x1_end = rt_template_a_span_end_idx + bar_width / 2

            span1_verts = [
                (x1_start, y_region_1 + region_height),
                (x1_end - arrow_tip_data, y_region_1 + region_height),
                (x1_end, y_region_1 + region_height / 2),
                (x1_end - arrow_tip_data, y_region_1),
                (x1_start, y_region_1),
                (x1_start, y_region_1 + region_height)
            ]
            span1_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
            path1 = Path(span1_verts, span1_codes)
            ax.add_patch(patches.PathPatch(path1, facecolor="#9e9e9e", edgecolor="none", alpha=0.6, clip_on=False))

            ax.text(
                (x1_start + bar_width * 0.5),
                y_region_1 + region_height / 2,
                # f"Flap A ({rt_template_a_span_end_idx + 1} bp)", 
                "Flap A RTT", 
                ha="left",
                va="center",
                fontsize=17,
                color="black",
                clip_on=False
            )

        # Flap B span
        if rt_template_b_span_start_idx is not None:
            x2_start = rt_template_b_span_start_idx - bar_width / 2
            x2_end = (n - 1) + bar_width / 2

            span2_verts = [
                (x2_end, y_region_2 + region_height),
                (x2_start + arrow_tip_data, y_region_2 + region_height),
                (x2_start, y_region_2 + region_height / 2),
                (x2_start + arrow_tip_data, y_region_2),
                (x2_end, y_region_2),
                (x2_end, y_region_2 + region_height)
            ]
            span2_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
            path2 = Path(span2_verts, span2_codes)
            ax.add_patch(patches.PathPatch(path2, facecolor="#d9d9d9", edgecolor="none", alpha=0.6, clip_on=False))

            ax.text(
                (x2_end - bar_width * 0.5),
                y_region_2 + region_height / 2,
                # f"Flap B ({n - rt_template_b_span_start_idx} bp)", 
                "Flap B RTT", 
                ha="right",
                va="center",
                fontsize=17,
                color="black",
                clip_on=False
            )

    # Legend
    legend_fontsize = 20
    points_per_inch = 72
    perfect_handle_length = (physical_bar_width * points_per_inch) / legend_fontsize

    fig.legend(
        loc="upper center",
        bbox_to_anchor=(0.50, 0.14 if has_spans else 0.16),
        ncol=4 if n > 33 else 3,
        frameon=False, 
        fontsize=legend_fontsize,
        handlelength=perfect_handle_length,
    )

    save_plot("b2.3'_base_integration_by_category", plot_formats, fig_root, fig)


def plot_flap_integration(
    total_counts,
    edit_counts,
    from_right_all_edit_counts,
    from_left_all_edit_counts, 
    tpe_indel_counts, 
    perfect_edit_counts, 
    insert_sequence,
    rt_template_a_span_end_idx=None, 
    rt_template_b_span_start_idx=None, 
    show_total_reads=False, 
    recoding_mode=False, 
    title=None, 
    fig_root=None,
    plot_formats=False, 
    category_colors=None
):
    has_spans = (rt_template_a_span_end_idx is not None) or (rt_template_b_span_start_idx is not None)

    n = len(insert_sequence)
    indices = np.arange(n)

    physical_bar_width = 0.25 
    physical_gap_width = 0.02 

    width_per_base = physical_bar_width + physical_gap_width
    bar_width = physical_bar_width / width_per_base

    spine_gap = bar_width / 2 + 0.15
    data_range = (n - 1 + spine_gap) - (-spine_gap)
    axes_width = data_range * width_per_base

    min_fig_width = 17.5 if show_total_reads else 15
    fig_height = 7 if has_spans else 6.25
    y_label_space = 1.0  
    right_margin = 0.5   

    natural_fig_width = axes_width + y_label_space + right_margin
    fig_width = max(min_fig_width, natural_fig_width)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

    extra_padding = fig_width - natural_fig_width
    axes_left = y_label_space + (extra_padding / 2)

    # Change spacing
    axes_bottom = 0.28 if has_spans else 0.22
    axes_height = 0.64 if has_spans else 0.67

    ax = fig.add_axes([axes_left / fig_width, axes_bottom, axes_width / fig_width, axes_height])

    if title:
        fig.suptitle(title,fontsize=24,y=0.97,ha="center")

    total_reads = max(total_counts)
    if from_left_all_edit_counts[-1] == from_right_all_edit_counts[0]:
        tpe_indel_counts = np.array(perfect_edit_counts) + from_right_all_edit_counts[0]
    # total_counts_pct = np.array(total_counts) / total_reads * 100
    edit_counts_pct = np.array(edit_counts) / total_reads * 100
    from_right_all_edit_counts_pct = np.array(from_right_all_edit_counts) / total_reads * 100
    from_left_all_edit_counts_pct = np.array(from_left_all_edit_counts) / total_reads * 100
    tpe_indel_perfect_tpe_counts = np.array(tpe_indel_counts) - np.array(perfect_edit_counts)
    tpe_indel_counts_pct = np.array(tpe_indel_perfect_tpe_counts) / total_reads * 100
    perfect_edit_counts_pct = np.array(perfect_edit_counts) / total_reads * 100

    # Mask Flap A to only plot up to rt_template_a_span_end_idx
    from_left_masked = np.array(from_left_all_edit_counts_pct, dtype=float)
    # if rt_template_a_span_end_idx is not None:
    #     from_left_masked[rt_template_a_span_end_idx + 1:] = np.nan

    # Mask Flap B to only plot from rt_template_b_span_start_idx onward
    from_right_masked = np.array(from_right_all_edit_counts_pct, dtype=float)
    # if rt_template_b_span_start_idx is not None:
    #     from_right_masked[:rt_template_b_span_start_idx] = np.nan

    if show_total_reads:
        # ax.bar(indices, total_counts_pct, width=bar_width, label="Base Absent", color=category_colors["WT"], alpha=1.0)
        # plot_max = 100
        ax.bar(indices, edit_counts_pct, width=bar_width, label="Base Present", color=category_colors["Dual Flap"], alpha=1.0)
        plot_max = min(100, max(edit_counts_pct) * 1.01)
    else:
        # plot_max = min(100, max(edit_counts_pct) * 1.05)
        plot_max = min(100, max(max(from_left_all_edit_counts_pct) * 1.01, max(from_right_all_edit_counts_pct) * 1.01))

    ax.bar(indices, from_left_masked, width=bar_width, color="#9e9e9e", label=f"Contiguous from Flap A", alpha=1.0)
    ax.bar(indices, from_right_masked, width=bar_width, color="#d9d9d9", label=f"Contiguous from Flap B", alpha=0.75)
    ax.bar(indices, tpe_indel_counts_pct, width=bar_width, label=f"{'All TPE Bases' if recoding_mode else 'Full TPE Sequence'} + Mutations", color=category_colors["Dual Flap"], alpha=1.0)
    ax.bar(indices, perfect_edit_counts_pct, width=bar_width, label="Perfect TPE", color=category_colors["Perfect TPE"], alpha=1.0)

    ax.set_ylabel("Reads (%)", fontsize=20)
    ax.set_xlim(-spine_gap, n - 1 + spine_gap)
    ax.set_ylim(0, plot_max)

    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([])
    ax.tick_params(axis='x', length=0)

    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', length=3, width=0.8)
    ax.tick_params(axis='y', which='major', labelsize=20, length=6, width=1)

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
            fontsize=16,
            clip_on=False
        )

    # Span calculations 
    region_text_offset = 0.01 * max_height
    label_text = "Programmed Base Changes" if recoding_mode else "Programmed Sequence"
    
    label_y = y_base - region_text_offset
    ax.text(
        (n - 1) / 2,
        label_y,
        label_text,
        ha="center",
        va="top",
        fontsize=20,
        color="black",
        clip_on=False
    )

    if has_spans:
        region_gap_inches = 0.05
        region_height_inches = 0.24
        arrow_tip_inches = 0.10

        label_height_inches = 0.28 
        
        label_clearance_data = label_height_inches / ax_height_in * max_height
        region_gap_data = region_gap_inches / ax_height_in * max_height
        region_height = region_height_inches / ax_height_in * max_height
        
        ax_pos = ax.get_position()
        ax_width_in = fig_width * ax_pos.width
        arrow_tip_data = (arrow_tip_inches / ax_width_in) * data_range

        y_region_1 = label_y - label_clearance_data - region_height
        y_region_2 = y_region_1 - region_gap_data - region_height

        Path = mpath.Path

        # Flap A span
        if rt_template_a_span_end_idx is not None:
            x1_start = -bar_width / 2
            x1_end = rt_template_a_span_end_idx + bar_width / 2

            span1_verts = [
                (x1_start, y_region_1 + region_height),
                (x1_end - arrow_tip_data, y_region_1 + region_height),
                (x1_end, y_region_1 + region_height / 2),
                (x1_end - arrow_tip_data, y_region_1),
                (x1_start, y_region_1),
                (x1_start, y_region_1 + region_height)
            ]
            span1_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
            path1 = Path(span1_verts, span1_codes)
            ax.add_patch(patches.PathPatch(path1, facecolor="#9e9e9e", edgecolor="none", alpha=0.6, clip_on=False))

            ax.text(
                (x1_start + bar_width * 0.5),
                y_region_1 + region_height / 2,
                # f"Flap A ({rt_template_a_span_end_idx + 1} bp)", 
                "Flap A RTT", 
                ha="left",
                va="center",
                fontsize=17,
                color="black",
                clip_on=False
            )

        # Flap B span
        if rt_template_b_span_start_idx is not None:
            x2_start = rt_template_b_span_start_idx - bar_width / 2
            x2_end = (n - 1) + bar_width / 2

            span2_verts = [
                (x2_end, y_region_2 + region_height),
                (x2_start + arrow_tip_data, y_region_2 + region_height),
                (x2_start, y_region_2 + region_height / 2),
                (x2_start + arrow_tip_data, y_region_2),
                (x2_end, y_region_2),
                (x2_end, y_region_2 + region_height)
            ]
            span2_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
            path2 = Path(span2_verts, span2_codes)
            ax.add_patch(patches.PathPatch(path2, facecolor="#d9d9d9", edgecolor="none", alpha=0.6, clip_on=False))

            ax.text(
                (x2_end - bar_width * 0.5),
                y_region_2 + region_height / 2,
                # f"Flap B ({n - rt_template_b_span_start_idx} bp)", 
                "Flap B RTT", 
                ha="right",
                va="center",
                fontsize=17,
                color="black",
                clip_on=False
            )

    # Legend
    legend_fontsize = 20
    points_per_inch = 72
    perfect_handle_length = (physical_bar_width * points_per_inch) / legend_fontsize

    fig.legend(
            loc="upper center",
            bbox_to_anchor=(0.50, 0.09 if has_spans else 0.11),
            ncol=5 if n > 33 else 3,
            frameon=False,
            fontsize=legend_fontsize,
            handlelength=perfect_handle_length,
        )

    # file_name = "b2.3'_flap_integration_all_edits" if show_total_reads else "b1.3'_flap_integration"
    file_name = "b1.3'_flap_integration"

    save_plot(file_name, plot_formats, fig_root, fig)


def plot_flap_removal(
    total_counts,
    edit_counts,
    from_right_all_edit_counts_del_region,
    from_left_all_edit_counts_del_region,
    perfect_edit_counts_del_region,
    deleted_sequence,
    rt_template_a_span_end_idx=None,
    rt_template_b_span_start_idx=None,
    show_total_reads=True,
    recoding_mode=False,
    title=None,
    fig_root=None,
    plot_formats=False,
    category_colors=None,
):
    has_spans = (rt_template_a_span_end_idx is not None) or (rt_template_b_span_start_idx is not None)

    n = len(deleted_sequence)
    indices = np.arange(n)

    physical_bar_width = 0.25 
    physical_gap_width = 0.02 

    width_per_base = physical_bar_width + physical_gap_width
    bar_width = physical_bar_width / width_per_base

    spine_gap = bar_width / 2 + 0.15
    data_range = (n - 1 + spine_gap) - (-spine_gap)
    axes_width = data_range * width_per_base

    min_fig_width = 17 if show_total_reads else 13
    fig_height = 7 if has_spans else 6.25
    y_label_space = 1.0  
    right_margin = 0.5   

    natural_fig_width = axes_width + y_label_space + right_margin
    fig_width = max(min_fig_width, natural_fig_width)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)

    extra_padding = fig_width - natural_fig_width
    axes_left = y_label_space + (extra_padding / 2)

    # Change spacing
    axes_bottom = 0.28 if has_spans else 0.22
    axes_height = 0.64 if has_spans else 0.67

    ax = fig.add_axes([axes_left / fig_width, axes_bottom, axes_width / fig_width, axes_height])

    if title:
        fig.suptitle(title, fontsize=20, y=0.97, ha="center")

    total_reads = max(total_counts)

    total_counts_pct = np.array(total_counts) / total_reads * 100
    edit_counts_pct = np.array(edit_counts) / total_reads * 100
    from_right_all_edit_counts_pct = np.array(from_right_all_edit_counts_del_region) / total_reads * 100
    from_left_all_edit_counts_pct = np.array(from_left_all_edit_counts_del_region) / total_reads * 100
    perfect_edit_counts_pct = np.array(perfect_edit_counts_del_region) / total_reads * 100

    # Mask Flap A to only plot up to rt_template_a_span_end_idx
    from_left_masked = np.array(from_left_all_edit_counts_pct, dtype=float)
    if rt_template_a_span_end_idx is not None:
        from_left_masked[rt_template_a_span_end_idx + 1:] = np.nan

    # Mask Flap B to only plot from rt_template_b_span_start_idx onward
    from_right_masked = np.array(from_right_all_edit_counts_pct, dtype=float)
    if rt_template_b_span_start_idx is not None:
        from_right_masked[:rt_template_b_span_start_idx] = np.nan

    # NaN-safe plot_max computation matching plot_flap_integration math
    if show_total_reads:
        ax.bar(indices, edit_counts_pct, width=bar_width, label="Base Removed", color=category_colors["Dual Flap"], alpha=1.0)
        peak_val = np.nanmax(np.nan_to_num(edit_counts_pct, nan=0.0))
        plot_max = min(100.0, max(1.0, peak_val * 1.03))
    else:
        valid_left = np.nan_to_num(from_left_masked, nan=0.0)
        valid_right = np.nan_to_num(from_right_masked, nan=0.0)
        valid_perfect = np.nan_to_num(perfect_edit_counts_pct, nan=0.0)

        peak_val = max(np.max(valid_left), np.max(valid_right), np.max(valid_perfect))
        plot_max = min(100.0, max(1.0, peak_val * 1.03))
    
    ax.bar(indices, from_left_masked, width=bar_width, color="#9e9e9e", label="Contiguous 5' Flap A Removal", alpha=1.0)
    ax.bar(indices, from_right_masked, width=bar_width, color="#d9d9d9", label="Contiguous 5' Flap B Removal", alpha=0.75)
    ax.bar(indices, perfect_edit_counts_pct, width=bar_width, label="Perfect TPE", color=category_colors["Perfect TPE"], alpha=1.0)

    ax.set_ylabel("Reads (%)", fontsize=16)

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

    # Span calculations
    region_text_offset = 0.01 * max_height
    label_text = "Wild-type Bases" if recoding_mode else "Wild-type Sequence"

    label_y = y_base - region_text_offset
    ax.text(
        (n - 1) / 2,
        label_y,
        label_text,
        ha="center",
        va="top",
        fontsize=16,
        color="black",
        clip_on=False,
    )

    if has_spans:
        region_gap_inches = 0.05
        region_height_inches = 0.24
        arrow_tip_inches = 0.10
        label_height_inches = 0.28

        label_clearance_data = label_height_inches / ax_height_in * max_height
        region_gap_data = region_gap_inches / ax_height_in * max_height
        region_height = region_height_inches / ax_height_in * max_height

        ax_width_in = fig_width * ax_pos.width
        arrow_tip_data = (arrow_tip_inches / ax_width_in) * data_range

        y_region_1 = label_y - label_clearance_data - region_height
        y_region_2 = y_region_1 - region_gap_data - region_height

        Path = mpath.Path

        # Flap A span
        if rt_template_a_span_end_idx is not None:
            x1_start = -bar_width / 2
            x1_end = rt_template_a_span_end_idx + bar_width / 2

            span1_verts = [
                (x1_start, y_region_1 + region_height),
                (x1_end - arrow_tip_data, y_region_1 + region_height),
                (x1_end, y_region_1 + region_height / 2),
                (x1_end - arrow_tip_data, y_region_1),
                (x1_start, y_region_1),
                (x1_start, y_region_1 + region_height),
            ]
            span1_codes = [
                Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
            ]
            path1 = Path(span1_verts, span1_codes)
            ax.add_patch(
                patches.PathPatch(
                    path1,
                    facecolor="#9e9e9e",
                    edgecolor="none",
                    alpha=0.6,
                    clip_on=False,
                )
            )

            padding_x = bar_width * 0.5
            ax.text(
                x1_start + padding_x,
                y_region_1 + region_height / 2,
                "WT Flap A",
                ha="left",
                va="center",
                fontsize=16,
                color="black",
                clip_on=False,
            )

        # Flap B span
        if rt_template_b_span_start_idx is not None:
            x2_start = rt_template_b_span_start_idx - bar_width / 2
            x2_end = (n - 1) + bar_width / 2

            span2_verts = [
                (x2_end, y_region_2 + region_height),
                (x2_start + arrow_tip_data, y_region_2 + region_height),
                (x2_start, y_region_2 + region_height / 2),
                (x2_start + arrow_tip_data, y_region_2),
                (x2_end, y_region_2),
                (x2_end, y_region_2 + region_height),
            ]
            span2_codes = [
                Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
            ]
            path2 = Path(span2_verts, span2_codes)
            ax.add_patch(
                patches.PathPatch(
                    path2,
                    facecolor="#d9d9d9",
                    edgecolor="none",
                    alpha=0.6,
                    clip_on=False,
                )
            )

            padding_x = bar_width * 0.5
            ax.text(
                x2_end - padding_x,
                y_region_2 + region_height / 2,
                "WT Flap B",
                ha="right",
                va="center",
                fontsize=16,
                color="black",
                clip_on=False,
            )

    # Legend
    legend_fontsize = 16
    points_per_inch = 72
    perfect_handle_length = (physical_bar_width * points_per_inch) / legend_fontsize

    fig.legend(
        loc="upper center",
        bbox_to_anchor=(0.50, 0.09),
        ncol=5 if n > 33 else 3,
        frameon=False,
        fontsize=legend_fontsize,
        handlelength=perfect_handle_length,
    )

    # file_name = "b5.5'_flap_removal_all_edits" if show_total_reads else "b4.5'_flap_removal"
    file_name = "b3.5'_flap_removal"

    save_plot(file_name, plot_formats, fig_root, fig)


def plot_flap_distribution(counts_dict, fig_root=None, plot_formats=False):
    flap_types = [
        ("Dual Flap", "Dual_Flap"), 
        ("Flap A\nOnly",   "Flap_A"),
        ("Flap A\nHybrid", "Flap_A_Hybrid"),
        ("Flap B\nOnly",   "Flap_B"),
        ("Flap B\nHybrid", "Flap_B_Hybrid"),
    ]

    labels = [ft[0] for ft in flap_types]
    full_vals = [counts_dict.get(f"Full {ft[1]}", 0) for ft in flap_types]
    partial_vals = [counts_dict.get(f"Partial {ft[1]}", 0) for ft in flap_types]
    
    totals = [f + p for f, p in zip(full_vals, partial_vals)]
    max_total = max(totals) if totals else 1
    
    visibility_threshold = max_total * 0.02
    
    fig, ax = plt.subplots(figsize=(4, 5), dpi=300)
    
    # Plot stacked bars
    bars_full = ax.bar(range(len(labels)), full_vals, color='#2d5a44', edgecolor="white", linewidth=0.2, label='Full\nFlap')
    bars_partial = ax.bar(range(len(labels)), partial_vals, bottom=full_vals, color='#739e82', edgecolor="white", linewidth=0.2, label='Partial\nFlap')
    
    # Add total count on top and segment percentages inside
    for i, total in enumerate(totals):
        if total == 0:
            continue

        # Total count above bar
        # ax.text(i, total + (max_total * 0.015), f"{int(total):,}", ha='center', va='bottom', fontsize=8)
        ax.text(i, total, f"n={int(total):,}", ha='center', va='bottom', fontsize=7, rotation=45)

        # Percentage within Full segment (only if physically tall enough)
        if full_vals[i] > visibility_threshold:
            full_pct = (full_vals[i] / total) * 100
            ax.text(i, full_vals[i] / 2, f"{full_pct:.1f}%", ha='center', va='center', fontsize=6.5, color='white')

        # Percentage within Partial segment (only if physically tall enough)
        if partial_vals[i] > visibility_threshold:
            partial_pct = (partial_vals[i] / total) * 100
            ax.text(i, full_vals[i] + (partial_vals[i] / 2), f"{partial_pct:.1f}%", ha='center', va='center', fontsize=6.5, color='white')

    ax.set_ylabel("Reads", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=7)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.minorticks_on()
    ax.tick_params(axis="y", which="minor", length=3, width=0.8)
    ax.tick_params(axis="y", which="major", length=5, width=1)
    ax.tick_params(axis="x", which="minor", bottom=False)
    
    ax.legend(fontsize=8, frameon=False, bbox_to_anchor=(.95, 1), loc='upper left')

    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], legend_labels[::-1], fontsize=8, frameon=False, bbox_to_anchor=(.95, 1), loc='upper left')

    plt.tight_layout()
    file_name = "b4.3'_flap_completion"
    save_plot(file_name, plot_formats, fig_root, fig, bbox_inches="tight")


def plot_per_base_pos_barplots(plotting_info, reference_info=None, twinspector_results_folder=None, plot_formats=False, recoding_mode=False, rt_templates=False, flap_data=None):

    setBarMatplotlibDefaults()

    if recoding_mode and rt_templates:
        rt_template_a_span_end_idx_integration = reference_info["rt_template_a_base_change_coverage_len"]-1
        rt_template_b_span_start_idx_integration = reference_info["rt_template_b_start_bp_changes_arr"]
        rt_template_a_span_end_idx_removal = reference_info["std_bp_changes_arr_len"]-1
        rt_template_b_span_start_idx_removal = 0
    elif not recoding_mode and rt_templates:
        rt_template_a_span_end_idx_integration = reference_info["rt_template_a_length"]-1
        rt_template_b_span_start_idx_integration = reference_info["rt_template_b_start_inserted_seq"]
        rt_template_a_span_end_idx_removal = reference_info["del_region_len"]-1
        rt_template_b_span_start_idx_removal = 0
    else:
        rt_template_a_span_end_idx_integration = None
        rt_template_b_span_start_idx_integration = None
        rt_template_a_span_end_idx_removal = None
        rt_template_b_span_start_idx_removal = None  

    plot_flap_integration(
        total_counts=plotting_info["total_read_bases_ins_region_arr"], 
        edit_counts=plotting_info["all_base_integration_counts_arr"],
        from_right_all_edit_counts=plotting_info["from_flap_b_contiguous_base_integration_counts"],
        from_left_all_edit_counts=plotting_info["from_flap_a_contiguous_base_integration_counts"],
        tpe_indel_counts=plotting_info["tpe_indel_counts_arr"],
        # tpe_indel_counts=plotting_info["dual_flap_counts_arr"],
        perfect_edit_counts=plotting_info["perfect_tpe_counts_arr"],
        insert_sequence=plotting_info["inserted_seq"], 
        rt_template_a_span_end_idx=rt_template_a_span_end_idx_integration, 
        rt_template_b_span_start_idx=rt_template_b_span_start_idx_integration, 
        show_total_reads=False, 
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder,
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS 
    )

    # plot_flap_integration(
    #     total_counts=plotting_info["total_read_bases_ins_region_arr"],  
    #     edit_counts=plotting_info["all_base_integration_counts_arr"], 
    #     from_right_all_edit_counts=plotting_info["from_flap_b_contiguous_base_integration_counts"], 
    #     from_left_all_edit_counts=plotting_info["from_flap_a_contiguous_base_integration_counts"], 
    #     tpe_indel_counts=plotting_info["tpe_indel_counts_arr"],
    #     perfect_edit_counts=plotting_info["perfect_tpe_counts_arr"], 
    #     insert_sequence=plotting_info["inserted_seq"], 
    #     rt_template_a_span_end_idx=rt_template_a_span_end_idx_integration, 
    #     rt_template_b_span_start_idx=rt_template_b_span_start_idx_integration, 
    #     show_total_reads=True,
    #     recoding_mode=recoding_mode, 
    #     fig_root=twinspector_results_folder,
    #     plot_formats=plot_formats, 
    #     category_colors=CATEGORY_COLORS 
    # )

    plot_base_integration_by_category(
        total_counts=plotting_info["total_read_bases_ins_region_arr"],  
        category_count_arrs=plotting_info["category_count_arrs"], 
        insert_sequence=plotting_info["inserted_seq"], 
        rt_template_a_span_end_idx=rt_template_a_span_end_idx_integration, 
        rt_template_b_span_start_idx=rt_template_b_span_start_idx_integration, 
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder, 
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS
    )  
   
    plot_flap_removal(
        total_counts=plotting_info["total_read_bases_del_region_arr"], 
        edit_counts=plotting_info["all_base_removal_counts_arr"],
        from_right_all_edit_counts_del_region=plotting_info["from_flap_b_contiguous_base_removal_counts"],
        from_left_all_edit_counts_del_region=plotting_info["from_flap_a_contiguous_base_removal_counts"],
        perfect_edit_counts_del_region=plotting_info["perfect_base_removal_counts_arr"],
        deleted_sequence=plotting_info["deleted_seq"], 
        rt_template_a_span_end_idx=rt_template_a_span_end_idx_removal,
        rt_template_b_span_start_idx=rt_template_b_span_start_idx_removal, 
        show_total_reads=False, 
        recoding_mode=recoding_mode, 
        fig_root=twinspector_results_folder,
        plot_formats=plot_formats, 
        category_colors=CATEGORY_COLORS
    )

    # plot_flap_removal(
    #     total_counts=plotting_info["total_read_bases_del_region_arr"], 
    #     edit_counts=plotting_info["all_base_removal_counts_arr"],
    #     from_right_all_edit_counts_del_region=plotting_info["from_flap_b_contiguous_base_removal_counts"],
    #     from_left_all_edit_counts_del_region=plotting_info["from_flap_a_contiguous_base_removal_counts"],
    #     perfect_edit_counts_del_region=plotting_info["perfect_base_removal_counts_arr"],
    #     deleted_sequence=plotting_info["deleted_seq"], 
    #     rt_template_a_span_end_idx=rt_template_a_span_end_idx_removal,
    #     rt_template_b_span_start_idx=rt_template_b_span_start_idx_removal, 
    #     show_total_reads=True, 
    #     recoding_mode=recoding_mode, 
    #     fig_root=twinspector_results_folder,
    #     plot_formats=plot_formats, 
    #     category_colors=CATEGORY_COLORS
    # )

    if flap_data:
        plot_flap_distribution(
            counts_dict=flap_data, 
            fig_root=twinspector_results_folder, 
            plot_formats=plot_formats
        )


def adjust_microhomology_alignment(
    df,
    motifs,
    comp,
    cut_points=None,
    junction_index=None,
):
    """
    Shift only the requested shared edge motifs across the composite WT/TPE
    junction, while restricting the affected alignment span to cut_points.
    """
    df_out = df.copy()

    if cut_points is None or junction_index is None:
        return df_out

    edit_start = int(cut_points[0]) + 1
    edit_end = int(cut_points[1])

    col_aligned = next(
        (c for c in df.columns if "Aligned_Sequence" in c and comp in c),
        None,
    )
    col_ref = next(
        (c for c in df.columns if "Reference_Sequence" in c and comp in c),
        None,
    )

    if col_aligned is None or col_ref is None:
        print(f"Warning: Could not find alignment columns for '{comp}'")
        return df_out

    def get_target_direction(category, reference_name):
        category = str(category).replace("_", " ")

        if "comp_a" in str(reference_name) or reference_name == "a":
            if category == "Perfect TPE":
                return "push_right"
            if category == "WT":
                return "push_left"

        if "comp_b" in str(reference_name) or reference_name == "b":
            if category == "Perfect TPE":
                return "push_left"
            if category == "WT":
                return "push_right"

        return None

    for idx, row in df_out.iterrows():
        aligned_seq = row.get(col_aligned)
        ref_seq = row.get(col_ref)
        category = row.get("Category_final")

        if pd.isna(aligned_seq) or pd.isna(ref_seq) or pd.isna(category):
            continue

        direction = get_target_direction(category, comp)
        if direction is None:
            continue

        new_aligned = str(aligned_seq)

        for motif_str in sorted(motifs, key=len, reverse=True):
            motif = re.escape(motif_str)

            if direction == "push_right":
                pattern = f"({motif})(-+)"

                def replacer(match):
                    match_start, match_end = match.span()
                    motif_text = match.group(1)
                    gaps = match.group(2)

                    # push_right
                    is_at_junction = (
                        match_end - len(motif_text) == junction_index
                    )
                    is_in_edit_window = (
                        edit_start <= match_start
                        and match_end - 1 <= edit_end
                    )
                    matches_reference = (
                        ref_seq[
                            match_end - len(motif_text):match_end
                        ].upper()
                        == motif_text.upper()
                    )

                    if (
                        is_at_junction
                        and is_in_edit_window
                        and matches_reference
                    ):
                        return gaps + motif_text

                    return match.group(0)

            else:
                pattern = f"(-+)({motif})"

                def replacer(match):
                    match_start, match_end = match.span()
                    gaps = match.group(1)
                    motif_text = match.group(2)

                    # push_left
                    is_at_junction = (
                        match_start + len(motif_text) == junction_index
                    )
                    is_in_edit_window = (
                        edit_start <= match_start
                        and match_end - 1 <= edit_end
                    )
                    matches_reference = (
                        ref_seq[
                            match_start:match_start + len(motif_text)
                        ].upper()
                        == motif_text.upper()
                    )

                    if (
                        is_at_junction
                        and is_in_edit_window
                        and matches_reference
                    ):
                        return motif_text + gaps

                    return match.group(0)

            new_aligned = re.sub(pattern, replacer, new_aligned)

        if new_aligned != aligned_seq:
            df_out.at[idx, col_aligned] = new_aligned

    return df_out


#### Allele tables ####
def get_dataframe_allele_region(
    df_alleles,
    pegRNA_intervals,
    ref_seq,
    ref_aln_seq,
    tpe_aln_seq=None,
    cut_points=None,
    window_by_intervals=False,
    left_pad=20,
    right_pad=20, 
    collapse_displayed_alleles=False,
    tpe_cats=None,
    tpe_freqs=None,
):
    """
    Return aligned sequences trimmed so all arrays match.
    If window_by_intervals is True, first trim to a window spanning pegRNA intervals
    with left_pad bases before the first interval start and right_pad bases after
    the last interval end; then harmonize to the shortest length within that window.
    """
    if df_alleles.shape[0] == 0:
        return df_alleles, ref_seq, ref_aln_seq, tpe_aln_seq, cut_points, pegRNA_intervals

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
    if tpe_aln_seq:
        tpe_aln_seq_win = tpe_aln_seq[start_idx:stop_idx_excl] if isinstance(tpe_aln_seq, str) else tpe_aln_seq

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

    # Slice allele strings to the window (if present)
    if 'Aligned_Sequence' in df.columns:
        df['Aligned_Sequence'] = df['Aligned_Sequence'].astype(str).str.slice(start_idx, stop_idx_excl)
    if 'Reference_Sequence' in df.columns:
        df['Reference_Sequence'] = df['Reference_Sequence'].astype(str).str.slice(start_idx, stop_idx_excl)

    # Determine smallest length across aligned strings in the window
    min_len_candidates = []
    if isinstance(ref_aln_seq_win, str):
        min_len_candidates.append(len(ref_aln_seq_win))
    if tpe_aln_seq and isinstance(tpe_aln_seq_win, str):
        min_len_candidates.append(len(tpe_aln_seq_win))
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
    tpe_aln_seq_trim = None
    if tpe_aln_seq:
        tpe_aln_seq_trim = tpe_aln_seq_win[:min_len] if isinstance(tpe_aln_seq_win, str) else tpe_aln_seq_win
    if 'Aligned_Sequence' in df.columns:
        df['Aligned_Sequence'] = df['Aligned_Sequence'].astype(str).str.slice(0, min_len)
    if 'Reference_Sequence' in df.columns:
        df['Reference_Sequence'] = df['Reference_Sequence'].astype(str).str.slice(0, min_len)

    if collapse_displayed_alleles:
        df["#Reads"] = pd.to_numeric(df["#Reads"], errors="raise")
        df["%Reads"] = pd.to_numeric(df["%Reads"], errors="raise")

        # Base grouping columns
        group_cols = ["Aligned_Sequence"]
        
        # Base aggregation dictionary
        agg_dict = {
            "#Reads": "sum",
            "%Reads": "sum",
            "Reference_Sequence": "first",
            "sequence_key": "first",
        }
        
        # If Category_final exists, add it to the GROUPING keys so it doesn't merge across categories
        if "Category_final" in df.columns:
            group_cols.append("Category_final")

        df = (
            df.groupby(group_cols, as_index=False, sort=False)
            .agg(agg_dict)
        )

    # Adjust cut points from global to window coords, then clamp to min_len
    cut_points_window = None
    if cut_points is not None:
        try:
            cps = list(cut_points)
        except TypeError:
            cps = [cut_points]
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

    if tpe_cats is not None and tpe_freqs is not None and "Category_final" in df.columns:
        # Ensure %Reads and freqs are numeric for accurate >= comparison
        df["%Reads"] = pd.to_numeric(df["%Reads"], errors="coerce")
        tpe_freqs_num = [float(x) for x in tpe_freqs]
        
        # Map thresholds
        threshold_map = dict(zip(tpe_cats, tpe_freqs_num))
        row_thresholds = df["Category_final"].map(threshold_map)
        
        # Filter dataframe
        df = df[df["%Reads"] >= row_thresholds].copy()
        
        # Create Ordered Categorical column
        df["Category_final"] = pd.Categorical(
            df["Category_final"], 
            categories=tpe_cats, 
            ordered=True
        )
        
        # Sort by custom Category order, then %Reads (decreasing)
        df.sort_values(
            by=["Category_final", "%Reads"],
            ascending=[True, False],
            inplace=True
        )
    else:
        # Fallback to original sorting if no filters are provided
        sort_cols = ["#Reads", "sequence_key"] if "sequence_key" in df.columns else ["#Reads"]
        sort_asc = [False, True] if "sequence_key" in df.columns else [False]
        
        df.sort_values(
            by=sort_cols,
            ascending=sort_asc,
            inplace=True,
        )

    return (
        df,
        ref_seq,                 
        ref_aln_seq_trim,        
        tpe_aln_seq_trim,       
        cut_points_window if cut_points_window is not None else cut_points,
        pegRNA_intervals_window
    )


def prep_alleles_table(
    df_alleles,
    reference_seq,
    ref_aln_seq_region,
    alternate_aln_seq_region,
    MAX_N_ROWS,
    MIN_FREQUENCY,
    pegRNA_intervals, 
    alternate_label="TwinPE Reference"
):
    """
    Prepare matrices and metadata required to render an allele heatmap.
    """
    dna_to_numbers = {"-": 0, "A": 1, "T": 2, "C": 3, "G": 4, "N": 5, " ": 6}
    seq_to_numbers = lambda seq: [dna_to_numbers[x] for x in seq]
    X = []
    annot = []
    y_labels = []
    insertion_dict = defaultdict(list)
    per_element_annot_kws = []
    is_reference = []
    category_headers = []  # Tracks (row_index, category_name)
    num_blanks = 2

    # Regex to find contiguous gap runs ('-') in the reference alignment
    re_find_indels = re.compile("(-*-)")
    idx_row = 0
    prev_cat = None

    for idx, row in df_alleles[df_alleles["%Reads"] >= MIN_FREQUENCY][
        :MAX_N_ROWS
    ].iterrows():
        
        # Category spacer and header logic
        if "Category_final" in row:
            curr_cat = row["Category_final"]
            if pd.notna(curr_cat):
                if prev_cat is None:
                    # Record the start of the first category
                    category_headers.append((idx_row, curr_cat))
                    prev_cat = curr_cat
                elif curr_cat != prev_cat:
                    # Inject blank rows for new categories
                    for _ in range(num_blanks):
                        blank_row = " " * len(idx)
                        X.append(seq_to_numbers(blank_row))
                        annot.append(list(blank_row))
                        y_labels.append("")
                        is_reference.append(False)
                        
                        to_append = np.array([{"color": "white"}] * len(blank_row), dtype=object)
                        per_element_annot_kws.append(to_append)
                        
                        idx_row += 1
                    
                    # Record the start of the new category after blank rows
                    category_headers.append((idx_row, curr_cat))
                    prev_cat = curr_cat

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

    if alternate_aln_seq_region is not None and len(alternate_aln_seq_region) > 0:
        for _ in range(num_blanks):
            blank_row = " " * len(alternate_aln_seq_region)
            X.append(seq_to_numbers(blank_row))
            annot.append(list(blank_row))
            y_labels.append("")
            is_reference.append(False)

            to_append = np.array(
                [{"color": "white"}] * len(blank_row),
                dtype=object,
            )
            per_element_annot_kws.append(to_append)

        X.append(seq_to_numbers(alternate_aln_seq_region.upper()))
        annot.append(list(alternate_aln_seq_region))
        y_labels.append(alternate_label)
        is_reference.append(False)

        to_append = np.array(
            [{}] * len(alternate_aln_seq_region),
            dtype=object,
        )
        per_element_annot_kws.append(to_append)

    return X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference, category_headers


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
        reference_label="WT Reference", 
        has_alternate_reference=False,
        fig_filename=None, 
        fig_root=None,
        custom_colors=None,
        plot_formats=False,
        plot_cut_point=True,
        cut_point_ind=None,
        sgRNA_intervals=None,
        sgRNA_names=None,
        sgRNA_mismatches=None,
        extend_left_non_gap=None, 
        extend_right_non_gap=None, 
        category=None, 
        category_headers=None, 
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
    # -fig_filename_root: figure filename to plot (not including '.pdf' or '.png'). If None, plots are shown interactively.
    -fig_filename: figure filename to plot (including '.pdf' or '.png'). If None, plots are shown interactively.
    -fig_root: figure filename root to plot (not including '.pdf' or '.png'). If None, plots are shown interactively.
    -custom_colors: dict of colors to plot (e.g. colors['A'] = (1,0,0,0.4) # red,blue,green,alpha )
    -plot_formats: 
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

    # Draw category section titles
    if category_headers:
        N_ROWS = len(annot)
        # Find the middle of the X-axis 
        N_COLS = len(annot[0]) if N_ROWS > 0 else 0
        x_center = (N_COLS / 2) - 4
        
        for r, cat_name in category_headers:
            display_name = str(cat_name).replace("_", " ")
            y_top_of_block = N_ROWS - r
            
            if r == 0:
                # Top category
                y_text = y_top_of_block + 0.15
                va_align = "bottom"
            else:
                # Subsequent categories
                y_text = y_top_of_block + 0.15
                va_align = "bottom"
                
            ax_hm.text(
                x_center,
                y_text, 
                display_name, 
                ha="center",
                va=va_align, 
                fontsize=24, 
                # fontweight="bold", 
                color="black",
                clip_on=False
            )

    # place ticks at row centers
    ax_hm.yaxis.tick_right()
    ax_hm.set_yticks(np.arange(N_ROWS) + 0.5)
    
    # apply labels in reverse order
    reversed_labels = y_labels[::-1]
    ax_hm.set_yticklabels(reversed_labels, rotation=True, va='center')
    
    # Dynamically turn off tick marks for ALL blank spacer rows
    for i, tick in enumerate(ax_hm.yaxis.get_major_ticks()):
        if i < len(reversed_labels) and reversed_labels[i] == "":
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)
            
    ax_hm.xaxis.set_ticks([])

    # Right-side labels and extended grey pegRNA spacerB span by caller
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

    #create boxes for ins
    for idx, lss in insertion_dict.items():
        for ls in lss:
            ax_hm.add_patch(patches.Rectangle((ls[0], N_ROWS-idx-1), ls[1]-ls[0], 1, linewidth=3, edgecolor='r', fill=False))

    # cut point vertical line at correct position and for WT and TwinPE refs rows
    if plot_cut_point:
        if cut_point_ind is None:
            cut_point_ind = [plot_nuc_len / 2]

        def _is_blank_row(chars):
            return all(c == " " for c in chars)

        # Convert base indices to boundary positions (draw line AFTER the specified base).
        raw_points = cut_point_ind if isinstance(cut_point_ind, (list, tuple, np.ndarray)) else [cut_point_ind]
        xs = [cp + 1 for cp in raw_points if cp is not None]

        # Map out contiguous Y-bound segments for blocks of data rows
        ymins = []
        ymaxs = []
        block_start = None
        
        for i in range(len(annot)):
            if not _is_blank_row(annot[i]):
                if block_start is None:
                    block_start = i  # Start of a new contiguous block
            else:
                if block_start is not None:
                    ymins.append(block_start)
                    ymaxs.append(i)  # End of the contiguous block
                    block_start = None
                    
        # Close the final block if the heatmap doesn't end on a blank row
        if block_start is not None:
            ymins.append(block_start)
            ymaxs.append(len(annot))

        # Draw segmented lines that skip the whitespaces
        if ymins:
            for x in xs:
                # Duplicate the x-coordinate so it matches the length of ymins/ymaxs
                x_coords = [x] * len(ymins)
                ax_hm.vlines(x_coords, ymins, ymaxs, linestyles="dashed", colors="black")

        # on WT Reference (ax_hm_ref): single row
        ax_hm_ref.vlines(xs, 0, 1, linestyles="dashed", colors="black")

    ax_hm_ref.yaxis.tick_right()
    ax_hm_ref.xaxis.set_ticks([])
    ax_hm_ref.yaxis.set_ticklabels(
        [reference_label],
        rotation=True,
        va="center",
    )
    gs2.update(left=0, right=1, hspace=0.05, wspace=0, top=1*(((N_ROWS)*1.13))/(N_ROWS))
    gs1.update(left=0, right=1, hspace=0.05, wspace=0,)

    sns.set_context(rc={'axes.facecolor':'white','lines.markeredgewidth': 1,'mathtext.fontset' : 'stix','text.usetex':True,'text.latex.unicode':True} )

    proxies = [matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='black',
                    mec='none', marker=r'$\mathbf{{{}}}$'.format('bold'), ms=16),
               matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='r', marker='s', ms=6, markeredgewidth=2),
              matplotlib.lines.Line2D([0], [0], linestyle='none', mfc='none',
                    mec='black', marker='_', ms=2,)]
    descriptions=['Substitutions', 'Insertions', 'Deletions']

    if plot_cut_point:
        proxies.append(
              matplotlib.lines.Line2D([0], [1], linestyle='--', c='black', ms=6))
        descriptions.append('Nick site')

    if category:
        category = category.replace("_", r"\ ")
        proxies.append(patches.Patch(color='none'))
        descriptions.append(rf"$\bf{{{category}}}$")

    #ax_hm_ref.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='center', bbox_to_anchor=(0.5, 4),ncol=1)
    lgd = ax_hm.legend(proxies, descriptions, numpoints=1, markerscale=2, loc='upper center', bbox_to_anchor=(0.5, 0), ncol=1, fancybox=True, shadow=False)

    save_plot(fig_filename, plot_formats, fig=fig, fig_root=fig_root, bbox_inches='tight', bbox_extra_artists=(lgd,))

    # if fig_filename_root is None:
    #     plt.show()
    # else:
    # fig.savefig(fig_filename_root+'.png', bbox_inches='tight', bbox_extra_artists=(lgd,), dpi=300)
    # if plot_formats:
    #     fig.savefig(fig_filename_root+'.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
    # plt.close(fig)


def plot_categorical_ref_allele_tables(
    df_alleles, 
    reference_info, 
    ref_seq, 
    tpe_seq, 
    ref_aln_seq,
    tpe_aln_seq, 
    no_alignment_adjustments, 
    min_frequency,
    max_n_rows, 
    plot_full_reads=False, 
    fig_root=None,
    plot_formats=False, 
    ref_type=None, 
    collapse_displayed_alleles=False
):

    same_length = len(ref_seq) == len(tpe_seq)
    if ref_type == "wt":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_wt", None)
        cut_points = reference_info.get("cut_points_wt", None)
        extend_left_non_gap = [0, 0]
        extend_right_non_gap = [0, 0]
        primary_seq = ref_seq
        primary_aln_seq = ref_aln_seq
        alternate_aln_seq = tpe_aln_seq if same_length else None
        primary_label = "WT Reference"
        alternate_label = "TwinPE Reference"
    elif ref_type == "tpe":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_tpe", None)
        cut_points = reference_info.get("cut_points_tpe", None)
        extend_left_non_gap = [0, 0]
        extend_right_non_gap = [0, 0]
        if same_length:
            # Two-reference tables always show WT first, TPE second.
            primary_seq = ref_seq
            primary_aln_seq = ref_aln_seq
            alternate_aln_seq = tpe_aln_seq
            primary_label = "WT Reference"
            alternate_label = "TwinPE Reference"
        else:
            # TPE-only table.
            primary_seq = tpe_seq
            primary_aln_seq = tpe_aln_seq
            alternate_aln_seq = None
            primary_label = "TwinPE Reference"
            alternate_label = "WT Reference"
    elif ref_type == "comp_a":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_composite_a", None)
        cut_points = reference_info.get("cut_points_composite_a", None)
        extend_left_non_gap = [0, abs(reference_info['cleavage_offset_b'])]
        extend_right_non_gap = [0, 0]
        primary_seq = ref_seq
        primary_aln_seq = ref_aln_seq
        alternate_aln_seq = tpe_aln_seq
        primary_label = "WT Reference"
        alternate_label = "TwinPE Reference"
    elif ref_type == "comp_b":
        pegRNA_intervals = reference_info.get("pegRNA_intervals_composite_b", None)
        cut_points = reference_info.get("cut_points_composite_b", None)
        extend_left_non_gap = [0, 0]
        extend_right_non_gap = [abs(reference_info['cleavage_offset_a']), 0]
        primary_seq = ref_seq
        primary_aln_seq = ref_aln_seq
        alternate_aln_seq = tpe_aln_seq
        primary_label = "WT Reference"
        alternate_label = "TwinPE Reference"

    for cat, name_comp_a, name_comp_b, name_tpe, name_wt in [
        ("Perfect_TPE", "e1.Perfect_TPE.aligned_to_composite_a", "e1.Perfect_TPE.aligned_to_composite_b", "e1.Perfect_TPE.aligned_to_tpe", "e1.Perfect_TPE.aligned_to_wt"), 
        ("Dual_Flap", "e2.Dual_Flap.aligned_to_composite_a", "e2.Dual_Flap.aligned_to_composite_b", "e2.Dual_Flap.aligned_to_tpe", "e2.Dual_Flap.aligned_to_wt"),
        ("Flap_A", "e3.Flap_A.aligned_to_composite_a", "e3.Flap_A.aligned_to_composite_b", "e3.Flap_A.aligned_to_tpe", "e3.Flap_A.aligned_to_wt"),
        ("Flap_B", "e4.Flap_B.aligned_to_composite_a", "e4.Flap_B.aligned_to_composite_b", "e4.Flap_B.aligned_to_tpe", "e4.Flap_B.aligned_to_wt"),
        ("Flap_A_Hybrid", "e5.Flap_A_Hybrid.aligned_to_composite_a", "e5.Flap_A_Hybrid.aligned_to_composite_b", "e5.Flap_A_Hybrid.aligned_to_tpe", "e5.Flap_A_Hybrid.aligned_to_wt"),
        ("Flap_B_Hybrid", "e6.Flap_B_Hybrid.aligned_to_composite_a", "e6.Flap_B_Hybrid.aligned_to_composite_b", "e6.Flap_B_Hybrid.aligned_to_tpe", "e6.Flap_B_Hybrid.aligned_to_wt"),
        ("Imperfect_TPE", "e7.Imperfect_TPE.aligned_to_composite_a", "e7.Imperfect_TPE.aligned_to_composite_b", "e7.Imperfect_TPE.aligned_to_tpe", "e7.Imperfect_TPE.aligned_to_wt"), 
        ("Null", "e8.Null.aligned_to_composite_a", "e8.Null.aligned_to_composite_b", "e8.Null.aligned_to_tpe", "e8.Null.aligned_to_wt"),
        ("Imperfect_WT", "e9.Imperfect_WT.aligned_to_composite_a", "e9.Imperfect_WT.aligned_to_composite_b", "e9.Imperfect_WT.aligned_to_tpe", "e9.Imperfect_WT.aligned_to_wt"),
        ("WT", "e10.WT.aligned_to_composite_a", "e10.WT.aligned_to_composite_b", "e10.WT.aligned_to_tpe", "e10.WT.aligned_to_wt"),
        ("Uncategorized", "e11.Uncategorized.aligned_to_composite_a", "e11.Uncategorized.aligned_to_composite_b", "e11.Uncategorized.aligned_to_tpe", "e11.Uncategorized.aligned_to_wt"),
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

        # Adjust homologies
        motifs_to_fix = []
        if not no_alignment_adjustments:
            cut_points_comp_a = reference_info["cut_points_composite_a"]
            cut_points_comp_b = reference_info["cut_points_composite_b"]
            if cat in ["Perfect_TPE", "WT"] and ref_type in ["comp_a", "comp_b"]:
                comp_a_junction = reference_info["composite_a_ins_start"]
                comp_b_junction = reference_info["composite_b_del_start"]
                if reference_info["num_bases_shared_start_for_homology_adj"]:
                    motifs_to_fix.append(reference_info["inserted_seq"][:reference_info["num_bases_shared_start_for_homology_adj"]])
                if reference_info["num_bases_shared_end_for_homology_adj"]:
                    motifs_to_fix.append(reference_info["inserted_seq"][-reference_info["num_bases_shared_end_for_homology_adj"]:])
                if motifs_to_fix:
                    junction_index = comp_a_junction if ref_type == "comp_a" else comp_b_junction
                    df_adjusted = adjust_microhomology_alignment(
                        df_alleles_cat, 
                        motifs_to_fix,
                        comp=ref_type, 
                        cut_points=cut_points_comp_a if ref_type == "comp_a" else cut_points_comp_b,
                        junction_index=junction_index
                    )

        df_alleles_around_region, ref_seq_region, ref_aln_seq_region, alternate_aln_seq_region, cut_points_window, pegRNA_intervals_region = get_dataframe_allele_region(
            df_adjusted if motifs_to_fix else df_alleles_cat,
            pegRNA_intervals,
            primary_seq,
            primary_aln_seq,
            alternate_aln_seq,
            cut_points,
            window_by_intervals=not plot_full_reads,
            left_pad=6,
            right_pad=6, 
            collapse_displayed_alleles=collapse_displayed_alleles
        )

        X, annot, y_labels, insertion_dict, per_element_annot_kws, is_reference, _ = prep_alleles_table(
            df_alleles_around_region,
            ref_seq_region,
            ref_aln_seq_region,
            alternate_aln_seq_region,
            max_n_rows,
            min_frequency,
            pegRNA_intervals_region,
            alternate_label=alternate_label,
        )

        has_alternate_reference = (
            alternate_aln_seq_region is not None
            and len(alternate_aln_seq_region) > 0
        )

        plot_alleles_heatmap(
            reference_aln_seq=ref_aln_seq_region,
            X=X,
            annot=annot,
            y_labels=y_labels,
            insertion_dict=insertion_dict,
            per_element_annot_kws=per_element_annot_kws, 
            reference_label=primary_label, 
            has_alternate_reference=has_alternate_reference,
            fig_filename=name, 
            fig_root=fig_root,
            plot_formats=plot_formats,
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
    df_alleles, reference_info, args, twinspector_results_folder, ref_type, ref_aln_seq, tpe_aln_seq, plot_formats = task

    setAlleleMatplotlibDefaults()

    plot_categorical_ref_allele_tables(
        df_alleles=df_alleles,
        reference_info=reference_info,
        ref_seq=args.wt_seq,
        tpe_seq=args.tpe_seq,
        ref_aln_seq=ref_aln_seq,
        tpe_aln_seq=tpe_aln_seq, 
        no_alignment_adjustments=args.no_alignment_adjustments,
        min_frequency=args.min_frequency_alleles,
        max_n_rows=args.max_n_rows,
        plot_full_reads=args.plot_full_reads,
        fig_root=twinspector_results_folder,
        ref_type=ref_type, 
        plot_formats=plot_formats, 
        collapse_displayed_alleles=args.collapse_displayed_alleles
    )
    return ref_type


# Dynamic allele-table plotting
def plot_ref_allele_tables(args, df_categorized, reference_info, twinspector_results_folder, plot_formats, n_processes=1):
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
        (df_categorized[cols], reference_info, args, twinspector_results_folder, ref_type, ref_aln_seq, tpe_aln_seq, plot_formats)
        for ref_type, cols, ref_aln_seq, tpe_aln_seq in ref_specs
    ]

    if n_processes == 1:
        setAlleleMatplotlibDefaults()
        for task in tasks:
            plot_one_allele_table(task)
        return

    with ProcessPoolExecutor(max_workers=min(n_processes, len(tasks))) as executor:
        futures = [executor.submit(plot_one_allele_table, task) for task in tasks]
        for future in as_completed(futures):
            future.result()


def sort_alleles_by_category(df, category_order):
    result = df.copy()
    result["category_order"] = pd.Categorical(
        result["Category_final"],
        categories=category_order,
        ordered=True,
    )
    return (
        result
        .sort_values(
            ["category_order", "%Reads", "sequence_key"],
            ascending=[True, False, True],
        )
        .drop(columns="category_order")
    )


def allele_table_summary_figures(
    df_categorized,
    reference_info,
    args,
    ref_type,
    ref_aln_seq,
    tpe_aln_seq,
    output_folder,
    plot_formats,
):

    all_summary_categories = list(reversed([
        cat.replace(" ", "_")
        for cat in CATEGORY_ORDER
    ]))

    categories_by_reference = {
        "wt": [
            "WT",
            "Imperfect_WT",
            "Null",
        ],
        "tpe": [
            "Perfect_TPE",
            "Dual_Flap",
            "Flap_A",
            "Flap_B",
            "Imperfect_TPE",
            "Null",
        ],
        "comp_a": all_summary_categories,
        "comp_b": all_summary_categories,
    }

    category_order = categories_by_reference[ref_type]

    columns = {
        "wt": ["sequence_key", "#Reads", "%Reads",
               "Aligned_Sequence_wt", "Reference_Sequence_wt",
               "Category_final", "Classified_by"],
        "tpe": ["sequence_key", "#Reads", "%Reads",
                "Aligned_Sequence_tpe", "Reference_Sequence_tpe",
                "Category_final", "Classified_by"],
        "comp_a": ["sequence_key", "#Reads", "%Reads",
                   "Aligned_Sequence_comp_a", "Reference_Sequence_comp_a",
                   "Category_final", "Classified_by"],
        "comp_b": ["sequence_key", "#Reads", "%Reads",
                   "Aligned_Sequence_comp_b", "Reference_Sequence_comp_b",
                   "Category_final", "Classified_by"],
    }

    reference_config = {
        "wt": {
            "ref_seq": args.wt_seq,
            "ref_aln_seq": args.wt_seq,
            "alternate_aln_seq": None,
            "label": "WT Reference",
            "intervals_key": "pegRNA_intervals_wt",
            "cut_points_key": "cut_points_wt",
            "extend_left": [0, 0],
            "extend_right": [0, 0],
        },
        "tpe": {
            "ref_seq": args.tpe_seq,
            "ref_aln_seq": args.tpe_seq,
            "alternate_aln_seq": None,
            "label": "TwinPE Reference",
            "intervals_key": "pegRNA_intervals_tpe",
            "cut_points_key": "cut_points_tpe",
            "extend_left": [0, 0],
            "extend_right": [0, 0],
        },
        "comp_a": {
            "ref_seq": args.wt_seq,
            "ref_aln_seq": reference_info["wt_aln_seq_a"],
            "alternate_aln_seq": reference_info["tpe_aln_seq_a"],
            "label": "WT Reference",
            "intervals_key": "pegRNA_intervals_composite_a",
            "cut_points_key": "cut_points_composite_a",
            "extend_left": [0, abs(reference_info["cleavage_offset_b"])],
            "extend_right": [0, 0],
        },
        "comp_b": {
            "ref_seq": args.wt_seq,
            "ref_aln_seq": reference_info["wt_aln_seq_b"],
            "alternate_aln_seq": reference_info["tpe_aln_seq_b"],
            "label": "WT Reference",
            "intervals_key": "pegRNA_intervals_composite_b",
            "cut_points_key": "cut_points_composite_b",
            "extend_left": [0, 0],
            "extend_right": [abs(reference_info["cleavage_offset_a"]), 0],
        },
    }
    
    config = reference_config[ref_type]
    pegRNA_intervals = reference_info[config["intervals_key"]]
    cut_points = reference_info[config["cut_points_key"]]

    df_selected = df_categorized[columns[ref_type]].copy()

    (
        df_region,
        ref_seq_region,
        ref_aln_seq_region,
        alternate_aln_seq_region,
        cut_points_window,
        pegRNA_intervals_region,
    ) = get_dataframe_allele_region(
        df_selected,
        pegRNA_intervals,
        config["ref_seq"],
        config["ref_aln_seq"],
        config["alternate_aln_seq"],
        cut_points,
        window_by_intervals=not args.plot_full_reads,
        left_pad=6,
        right_pad=6,
        collapse_displayed_alleles=args.collapse_displayed_alleles,
    )
  
    df_region = df_region[
        df_region["Category_final"].isin(category_order)
        & (df_region["%Reads"] >= args.min_frequency_alleles)
    ] 

    df_region = sort_alleles_by_category(df_region, category_order)

    df_region = (
        df_region
        .groupby("Category_final", sort=False, observed=True)
        .head(10)
        # .head(args.max_n_rows)
    )

    (
        X,
        annot,
        y_labels,
        insertion_dict,
        per_element_annot_kws,
        is_reference,
        category_headers,
    ) = prep_alleles_table(
        df_region,
        ref_seq_region,
        ref_aln_seq_region,
        alternate_aln_seq_region,
        MAX_N_ROWS=len(df_region),
        MIN_FREQUENCY=0,
        pegRNA_intervals=pegRNA_intervals_region,
        alternate_label="TwinPE Reference",
    )

    filename_by_reference = {
        "comp_a": "a6.Top_alleles.aligned_to_comp_a",
        "comp_b": "a7.Top_alleles.aligned_to_comp_b",
        "tpe": "a8.Top_alleles.aligned_to_tpe",
        "wt": "a9.Top_alleles.aligned_to_wt",
    }

    plot_alleles_heatmap(
        reference_aln_seq=ref_aln_seq_region,
        X=X,
        annot=annot,
        y_labels=y_labels,
        insertion_dict=insertion_dict,
        per_element_annot_kws=per_element_annot_kws,
        reference_label=config["label"],
        fig_filename=filename_by_reference[ref_type],
        fig_root=output_folder,
        plot_formats=plot_formats,
        plot_cut_point=True,
        cut_point_ind=cut_points_window,
        sgRNA_intervals=pegRNA_intervals_region,
        sgRNA_names=["pegRNA a", "pegRNA b"],
        sgRNA_mismatches=[[], []],
        category_headers=category_headers, 
        extend_left_non_gap=config["extend_left"],
        extend_right_non_gap=config["extend_right"],      
    )


def plot_summarized_allele_tables(
    args,
    df_categorized,
    reference_info,
    output_folder,
    plot_formats,
    n_processes=1,
):

    setAlleleMatplotlibDefaults()

    ref_specs = [
        ("comp_a", reference_info["wt_aln_seq_a"], reference_info["tpe_aln_seq_a"]),
        ("comp_b", reference_info["wt_aln_seq_b"], reference_info["tpe_aln_seq_b"]),
        ("tpe", reference_info["tpe_aln_seq_a"], None),
        ("wt", reference_info["wt_aln_seq_a"], None),
    ]

    tasks = [
        (df_categorized, reference_info, args, ref_type,
         ref_aln_seq, tpe_aln_seq, output_folder, plot_formats)
        for ref_type, ref_aln_seq, tpe_aln_seq in ref_specs
    ]

    if n_processes == 1:
        setAlleleMatplotlibDefaults()
        for task in tasks:
            allele_table_summary_figures(*task)
        return

    with ProcessPoolExecutor(max_workers=min(n_processes, len(tasks))) as executor:
        futures = [
            executor.submit(allele_table_summary_figures, *task)
            for task in tasks
        ]
        for future in as_completed(futures):
            future.result()


def generate_report(twinspector_results_folder, parent_folder, html_filename="TwInsPEctor_report.html"):
    """
    Bundle every plot and table inside the structured `twinspector_results_folder` 
    into a single HTML report. The HTML report is saved directly inside `twinspector_results_folder`.
    """
    base_dir = Path(twinspector_results_folder)
    if not base_dir.is_dir():
        raise ValueError(f"{twinspector_results_folder} is not a valid directory")

    captions = {
        "Reads input": "Figure 1: Description for reads input goes ...",
        "Outcomes": "Figure 2: Description for uncategorized outcomes goes ...",
        "Outcomes stacked": "Figure 3: Description for stacked uncategorized outcomes...",
        "Categorized": "Figure 4: Description for categorized outcomes...",
        "Categorized stacked": "Figure 5: Description for stacked categorized outcomes...",
        "Aligned to Composite A": "Description for alignments to Composite A...",
        "3' Flap integration": "Description detailing the 3' flap integration...",
        "By category": "Description of the data broken down by category...",
    }

    # Helper for natural sorting
    def natural_sort_key(s):
        return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

    # Glob for files with the updated directory structure
    a_files = sorted(base_dir.glob("a*.png"), key=lambda p: natural_sort_key(p.name))
    b_files = sorted((base_dir / "b.base_plots").glob("b*.png"), key=lambda p: natural_sort_key(p.name)) if (base_dir / "b.base_plots").is_dir() else []
    c_files = sorted((base_dir / "c.mutation_plots").glob("c*.png"), key=lambda p: natural_sort_key(p.name)) if (base_dir / "c.mutation_plots").is_dir() else []
    d_files = sorted((base_dir / "d.text_files").glob("d*.txt"), key=lambda p: natural_sort_key(p.name)) if (base_dir / "d.text_files").is_dir() else []
    e_files = sorted((base_dir / "e.allele_tables").glob("e*.png"), key=lambda p: natural_sort_key(p.name)) if (base_dir / "e.allele_tables").is_dir() else []

    cards_html = []
    uid = 0

    # Helper to generate standard Tabbed Cards
    def build_tabbed_card(title, files, get_label_fn, uid_prefix, img_max_width="85%", img_max_height="none"):
        nonlocal uid
        if not files: return ""
        
        buttons = []
        panels = []
        
        for idx, filepath in enumerate(files):
            is_first = (idx == 0)
            active = "active" if is_first else ""
            show_active = "show active" if is_first else ""
            aria_selected = "true" if is_first else "false"
            
            panel_id = f"panel-{uid_prefix}-{idx}"
            
            # Use the provided labeling function to get the tab name
            label_text = html.escape(get_label_fn(filepath.name))
            
            # Lookup the caption, or use a default placeholder
            caption_text = captions.get(get_label_fn(filepath.name), "[ Add figure description for this plot here ]")
            
            # Since the HTML is inside base_dir, paths are simply relative to base_dir.
            rel_path = filepath.relative_to(base_dir)
            img_path = html.escape(str(rel_path.as_posix()))
            
            # Href now points directly to the PNG for zooming
            href_path = img_path 
            
            buttons.append(
                f'''<li class="nav-item">
                    <button class="nav-link {active}" data-bs-toggle="tab" data-bs-target="#{panel_id}" role="tab" aria-controls="{panel_id}" aria-selected="{aria_selected}">{label_text}</button>
                </li>'''
            )
            
            panels.append(f"""
                <div class="tab-pane fade {show_active}" id="{panel_id}" role="tabpanel">
                    <div class="d-flex flex-column align-items-center">
                        <a href="{href_path}" data-fancybox="gallery" style="width: 100%; text-align: center;">
                            <img src="{img_path}" class="report-img" style="max-width: {img_max_width}; max-height: {img_max_height}; width: auto; height: auto; margin: auto;">
                        </a>
                        <div class="caption-text mt-3">
                            {caption_text}
                        </div>
                    </div>
                </div>
            """)
            
        uid += 1
        return f"""
        <div class='card mb-3 breakinpage'>
            <div class='card-header text-center'>
                <h5>{title}</h5>
                <ul class="nav nav-tabs justify-content-center" id="tab-{uid_prefix}" role="tablist">
                    {''.join(buttons)}
                </ul>
            </div>
            <div class='card-body text-center'>
                <div class="tab-content" id="tabContent-{uid_prefix}">
                    {''.join(panels)}
                </div>
            </div>
        </div>
        """

    # Build Summary Cards (A files)
    a_files_1_5 = [f for f in a_files if re.match(r'^a[1-5]\.', f.name)]
    a_files_6_plus = [f for f in a_files if re.match(r'^a[6-9]\.', f.name)]

    def get_a1_label(name):
        norm = name.split('.')[1].replace('_', ' ').lower()
        mapping = {
            "reads input": "Reads input",
            "outcomes": "Outcomes",
            "outcomes stacked": "Outcomes stacked",
            "outcomes categorized": "Categorized",
            "outcome categorized stacked": "Categorized stacked",
            "outcomes categorized stacked": "Categorized stacked" # Catch alternative spelling
        }
        return mapping.get(norm, norm.capitalize())
        
    def get_align_label(name):
        match = re.search(r'\.aligned_to_(.+?)\.png', name)
        if match:
            norm = match.group(1).lower()
            mapping = {
                "comp_a": "Aligned to Composite A",
                "composite_a": "Aligned to Composite A",
                "comp_b": "Aligned to Composite B",
                "composite_b": "Aligned to Composite B",
                "tpe": "Aligned to TwinPE",
                "wt": "Aligned to WT"
            }
            return mapping.get(norm, f"Aligned to {norm.upper()}")
        return name

    if a_files_1_5:
        cards_html.append(build_tabbed_card("Results Summary", a_files_1_5, get_a1_label, "A1", img_max_height="60vh"))
        
    if a_files_6_plus:
        cards_html.append(build_tabbed_card("Most Frequent Alleles", a_files_6_plus, get_align_label, "A2", img_max_width="100%"))

    # Build Base Plots Card (B files)
    def get_b_label(name):
        norm = name.split('.')[1].replace('_', ' ').lower()
        mapping = {
            "3' flap integration": "3' Flap integration",
            "3 flap integration": "3' Flap integration",
            "3' base integration by category": "By category",
            "3 base integration by category": "By category",
            "3' flap completion": "3' Flap completion",
            "3 flap completion": "3' Flap completion",
            "5' flap removal": "5' Flap removal",
        }
        return mapping.get(norm, norm.capitalize())
    
    if b_files:
        cards_html.append(build_tabbed_card("Contiguous Editing", b_files, get_b_label, "B", img_max_height="60vh"))

    # Build Mutation Plots Card (C files)
    def get_c_label(name):
            match = re.search(r'c\d+\.(.+?)(?:\.ins_del|\.png)', name)
            return match.group(1).replace('_', ' ') if match else name

    if c_files:
        cards_html.append(build_tabbed_card("Insertions, Deletions, Substitutions", c_files, get_c_label, "C", img_max_height="60vh"))

    # Build Allele Tables Cards (E files)
    if e_files:
        allele_groups = {}
        for f in e_files:
            m = re.match(r"^(?P<stem>e\d+\..+?)\.aligned_to_(?P<target>wt|tpe|composite_a|composite_b)\.png$", f.name)
            if m:
                stem = m.group("stem")
                allele_groups.setdefault(stem, []).append(f)
                
        for stem in sorted(allele_groups.keys(), key=natural_sort_key):
            stem_files = allele_groups[stem]
            stem_title = stem.split('.')[1].replace('_', ' ') + " Alleles"
            
            # We reuse the get_align_label function to prefix tabs with "Aligned to "
            cards_html.append(build_tabbed_card(stem_title, stem_files, get_align_label, f"E_{stem.replace('.', '_')}"))

    # Build Text Data Files List (D files)
    if d_files:
        file_links = []
        for f in d_files:
            clean_name = f.name.split('.')[1].replace('_', ' ').title()
            rel_path = f.relative_to(base_dir).as_posix()
            href_path = html.escape(rel_path)
            file_links.append(f"<div class='py-1'><i class='far fa-file-alt text-muted me-2'></i><a href='{href_path}' target='_blank' class='text-decoration-none'>{clean_name}</a></div>")
            
        cards_html.append(f"""
        <div class='card mb-3 breakinpage'>
            <div class='card-header text-center'>
                <h5>Text Reports & Data Files</h5>
            </div>
            <div class='card-body d-flex flex-column align-items-center'>
                <div class="text-start" style="min-width: 250px;">
                    {''.join(file_links)}
                </div>
            </div>
        </div>
        """)

    # HTML Template
    html_out = f"""<!DOCTYPE html>
<html lang="en">
  <head>
    <title>TwInsPEctor Report</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link href='https://fonts.googleapis.com/css?family=Montserrat:300,400,500,600' rel='stylesheet' type='text/css'>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.1.3/css/bootstrap.min.css" integrity="sha512-GQGU0fMMi238uA+a/bdWJfpUGKUkBdgfFdgBm72SUQ6BeyWjoY/ton0tEjH+OSH9iP4Dfh+7HM0I9f5eR0L/4w==" crossorigin="anonymous" referrerpolicy="no-referrer" />    
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/fancybox/3.3.5/jquery.fancybox.min.css" />
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.0.13/css/all.css" integrity="sha384-DNOHZ68U8hZfKXOrtjWvjxusGo9WQnrNx2sqG0tfsghAvtVlRW3tvkXWZh58N9jp" crossorigin="anonymous">
    
    <style>
      body {{ font-family: 'Montserrat', "Segoe UI", Helvetica, Arial, sans-serif !important; background-color: #f4f6f9; color: #333; }}
      
      /* Header & Logo Styles */
      .wrapper {{ text-align: center; margin: 1rem 0; }}
      .logo {{ font-size: 70px; font-weight: 350; letter-spacing: 1px; display: inline-flex; align-items: baseline; }}
      .twinpe {{ color: black; }}
      .spector {{ color: gray; font-weight: 300; }}
      .magnifier {{ width: 80px; height: 95px; margin: 0 -8px; overflow: visible; }}
      
      /* Subtitle & Path Sizing */
      .subtitle {{ font-size: 18px; letter-spacing: 2px; color: gray; margin-top: -20px; }}
      .report-path {{ font-size: 14px; margin-top: 20px; word-break: break-all; color: #888; }}

      /* Card Customization */
      .card {{ border: none; border-radius: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.04); overflow: hidden; }}
      .card-header {{ background-color: #ffffff; border-bottom: none; padding: 1.5rem 1.5rem 0 1.5rem; }}
      
      /* UPDATED: Smaller Panel Titles */
      .card-header h5 {{ font-weight: 600; color: #333; text-transform: uppercase; letter-spacing: 1px; font-size: 0.95rem; margin-bottom: 1rem; }}
      .card-body {{ background-color: #ffffff; padding: 2rem; }}

      /* Clean Tabs */
      .nav-tabs {{ border-bottom: 2px solid #edf1f5; gap: 8px; }}
      
      /* UPDATED: Smaller Tab Buttons (Figure Names) */
      .nav-tabs .nav-link {{ border: none; border-bottom: 3px solid transparent; color: #777; font-weight: 500; padding: 0.5rem 1rem; font-size: 0.85rem; transition: all 0.2s ease; background: transparent; }}
      .nav-tabs .nav-link:hover {{ color: #2166ac; border-bottom-color: #c2d5e6; }}
      .nav-tabs .nav-link.active {{ border: none; border-bottom: 3px solid #2166ac; color: #2166ac; background: transparent; font-weight: 600; }}

      /* Image & Text Formatting */
      .report-img {{ border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.03); border: 1px solid #f0f0f0; transition: transform 0.2s; }}
      .report-img:hover {{ transform: scale(1.005); }}
      
      /* NEW: Left-aligned, readable caption text */
      .caption-text {{ font-weight: 400; color: #555; font-size: 0.9rem; text-align: left; max-width: 85%; width: 100%; line-height: 1.5; }}
      
      a {{ color: #2166ac; }}
      a:hover {{ color: #16477a; }}

      @media print {{
        .tab-content > .tab-pane {{ display: block !important; opacity: 1 !important; visibility: visible !important; margin-bottom: 2rem !important; }}
        .nav-tabs {{ display:none !important; visibility:hidden !important; }}
        .breakinpage {{ clear: both; page-break-before: always !important; display: block; }}
        .card {{ box-shadow: none !important; border: 1px solid #eee !important; }}
      }}
    </style>

    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/fancybox/3.3.5/jquery.fancybox.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
  </head>

  <body>
    <div class="container pb-5">
      
      <!-- Centered Header Block -->
      <div class="wrapper">
        <div class="logo">
          <span class="twinpe">TwIn</span>
          <span class="spector">s</span>
          <svg class="magnifier" viewBox="0 0 94 120">
            <line x1="19" y1="59" x2="19" y2="125" stroke="black" stroke-width="5" stroke-linecap="round"/>
            <circle cx="50" cy="60" r="27" stroke="black" stroke-width="5" fill="none"/>
            <circle cx="55" cy="55" r="16" fill="gray" opacity="0.72"/>
            <circle cx="55" cy="56" r="6" fill="black"/>
            <path d="M24 27 C37 16, 47 13, 57 19 S75 30, 80 28 C83 28, 85 32, 81 34 C71 41, 62 31, 53 25 S40 19, 26 26 Z" fill="black"/>
          </svg>
          <span class="twinpe">E</span>
          <span class="spector">ctor</span>
        </div>
        <div class="subtitle">TWIN PRIME EDITING ANALYSIS</div>
        <div class="report-path">{html.escape(parent_folder)}</div>
      </div>

      <!-- Report Cards Container -->
      <div class="row justify-content-center">
        <div class="col-lg-10 col-md-12">
            <div id='jumbotron_content'>
                {''.join(cards_html)}
                <div align="center" class='p-4'>
                  <button class='btn btn-outline-secondary hidden-print' onclick='window.print();' style="border-radius: 20px; padding: 8px 24px;">
                    <i class="fas fa-print me-2"></i> Print Report
                  </button>
                </div>
            </div>
        </div>
      </div>

    </div>
  </body>
</html>
"""

    output_filepath = base_dir / html_filename
    with open(output_filepath, 'w') as f:
        f.write(html_out)


def parse_args():
    prog = 'TwInsPEctor'
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Analyze amplicon sequencing data from twin prime editing experiments. TwInsPEctor aligns reads to standard and composite references using CRISPResso2, classifies allelic outcomes into 10 categories, and provides detailed visualizations of the editing results.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            f'Example: {prog}'
            "--fastq_r1 <fastq R1 file> --fastq_r2 <fastq R2 file> "
            "--wt_seq <full wild-type amplicon sequence> --tpe_seq <full twin prime edited amplicon sequence> "
            "--peg_spacers <pegRNA-a spacer sequence>,<pegRNA-b spacer sequence>\n"
            "--rt_templates: <RT template A>,<RT template B>\n\n"
            "----Category Definitions----\n"
            "Perfect TPE: Complete incorporation of the programmed sequence without unintended substitutions between the nick sites and without unintended insertions and deletions anywhere. Single base substitutions outside of the nick sites are permitted."
            "Dual Flap: Meets both Flap A and Flap B criteria.\n"
            "Flap A: At least N consecutive bases originating from the start of the pegRNA-a encoded sequence. Unintended substitutions, insertions, and deletions are permitted. Default: N=3 (sequence replacement), N=2 (sequence recoding)\n"
            "Flap B: At least N consecutive bases originating from the start of the pegRNA-b encoded sequence. Unintended substitutions, insertions, and deletions are permitted.\n"
            "Flap A Hybrid: Meets Flap A criteria and retains wild-type (5`) Flap B.\n"
            "Flap B Hybrid: Meets Flap B criteria and retains wild-type (5`) Flap A.\n"
            "Imperfect TPE: Contains partials of the programmed sequence that meet neither Flap A nor Flap B criteria. Residual wild-type sequence and unintended substitutions, insertions, and deletions are permitted.\n"
            "Null: Lacks any detectable wild-type and programmed sequence between the two pegRNA nick sites.\n"
            "Imperfect WT: Partial wild-type sequence with none of the programmed sequence. Unintended substitutions, insertions, and deletions are permitted.\n"
            "Imperfect WT: Partial wild-type sequence with none of the programmed sequence. Unintended substitutions, insertions, and deletions are permitted.\n"
            "WT: Retention of the complete wild-type sequence without unintended substitutions between the nick sites and without unintended insertions and deletions anywhere. Single base substitutions outside of the nick sites are permitted.\n"
        )
    )

    parser.add_argument("-r1", "--fastq_r1", type=str, required=True, help="Path to FASTQ R1 file.")
    parser.add_argument("-r2", "--fastq_r2", type=str, required=False, help="Path to FASTQ R2 file for paired-end data.")
    parser.add_argument("-w", "--wt_seq", type=str, required=True, help="Full wild-type reference amplicon sequence including spacers.")
    parser.add_argument("-t", "--tpe_seq", type=str, required=True, help="Full Twin prime edited reference amplicon sequence with 5' & 3' ends identical to wildtype reference amplicon.")
    parser.add_argument("-g", "--peg_spacers", type=str, required=True, help="Comma-separated pegRNA spacer sequences: <spacer A>,<spacer B>. Should include bases immediately adjacent to but not including the PAM sequence (usually 20nt 5' of NGG).")
    parser.add_argument("-rt", "--rt_templates", type=str, required=False, default=None, help="Comma-separated pegRNA reverse transcriptase templates: <RT template A>,<RT template B>. Informs flap analysis and plotting.")
    parser.add_argument("-o", "--output_root", type=str, default=None, help="Root output folder for CRISPResso2 and TwInsPEctor results. If not provided, a folder will be created in the current working directory based on the input fastq file names.")
    parser.add_argument("-rcm", "--recoding_mode", action="store_true", help="Run in recoding mode if the wild-type and twin prime edited sequences are the same length and should be evaluated as having only base substitutions.")
    parser.add_argument("-ne", "--min_num_base_edits", type=int, default=None, help="Minimum number of base changes required for a read to be considered edited. Default is 3 for replacment mode and 2 for recoding mode.")
    parser.add_argument("-dmas", "--default_min_aln_score", type=int, default=30, help="Default minimum homology score for a read to align to the compound reference amplicon")
    parser.add_argument("-pfr", "--plot_full_reads", action="store_true", help="Allele tables will display full read sequences.")
    parser.add_argument("-ncda", "--no_collapse_displayed_alleles", action='store_false', dest="collapse_displayed_alleles", default=True, help="Do not combine alleles that become identical in the displayed allele table window.")
    parser.add_argument("-naa", "--no_alignment_adjustments", action="store_true", default=False, help="Do not visually adjust homologies in the allele tables. This does not affect allele classification.")
    parser.add_argument("-ied", "--ignore_extraspacer_deletions", action="store_true", help="Classification ignores deletions occurring beyond the spacers (outside edit window).")
    parser.add_argument("-nf", "--no_figures", action="store_true", help="Skip all figures if only text outputs are desired.")
    parser.add_argument("-nsf", "--no_summary_figures", action="store_true", help="Skip summary barplots if they are not desired.")
    parser.add_argument("-nbf", "--no_per_base_figures", action="store_true", help="Skip per-base barplots if they are not desired.")
    parser.add_argument("-nmf", "--no_mutation_figures", action="store_true", help="Skip mutation barplots if they are not desired.")
    parser.add_argument("-pet", "--plot_extended_tables", action="store_true", help="Generates separate allele tables for each category. Use --max_n_rows and --min_frequency_alleles to control how many alleles are displayed in each table.")
    parser.add_argument("-pdf", "--save_pdf", action="store_true", help="Only save PDF versions of all plots.")
    parser.add_argument("-npng", "--no_save_png", action="store_true", help="Do not save PNG versions of all plots.")
    parser.add_argument("-mfa", "--min_frequency_alleles", type=float, default=0.1, help="Minimum percent read frequency required to report an allele in the alleles tables.")
    parser.add_argument("-mnr", "--max_n_rows", type=int, default=25, help="Maximum number of allele rows to display in the allele tables by category.")
    parser.add_argument("-mna", "--max_n_alleles_to_write", type=int, default=50, help="Maximum number of alleles per category to write to the f7 text file.")
    parser.add_argument("-nrr", "--no_rerun", action="store_true", help="Don't rerun CRISPResso2 if a run using the same parameters has already been finished.")
    parser.add_argument("-kco", "--keep_crispresso_outputs", action="store_true", help="Don't delete CRISPResso2 output folders after analysis.")
    parser.add_argument("--crispresso_args", type=str, default="", help='Additional arguments to pass to CRISPResso2 (wrapped in quotes, use equal sign); Example: --crispresso_args="--trim_sequences"; do not use --n_processes.')
    parser.add_argument("-coa", "--cleavage_offset_a", type=int, default=-3, help="Cleavage offset for pegRNA spacer A (default: -3).")
    parser.add_argument("-cob", "--cleavage_offset_b", type=int, default=-3, help="Cleavage offset for pegRNA spacer B (default: -3).")
    parser.add_argument("-p", "--n_processes", "--n_threads", dest="n_processes", type=int, default=1, metavar="N", help="Total process budget. CRISPResso2 divides it across up to four reference runs; allele-table plotting uses up to N processes (default: 1).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose CRISPResso2 output.")
    parser.add_argument("-V", "--version", action="version", version="%(prog)s 1.0.0")

    args = parser.parse_args()

    if args.n_processes < 1:
        parser.error("--n_processes/--n_threads must be at least 1")

    command_str = shlex.join(sys.argv)

    # By default, set to 3 for replacement mode, 2 for recoding mode
    if args.min_num_base_edits is None:
        args.min_num_base_edits = 2 if args.recoding_mode else 3

    args.wt_seq = args.wt_seq.upper()
    args.tpe_seq = args.tpe_seq.upper()

    extra_crispresso_args = shlex.split(args.crispresso_args)

    if any(
        argument in ("--n_processes", "-p") or 
        argument.startswith(("--n_processes=", "-p=")) or 
        (argument.startswith("-p") and argument[2:].isdigit())
        for argument in extra_crispresso_args
    ):
        parser.error("Do not include --n_processes in --crispresso_args; use TwInsPEctor's --n_processes/--n_threads total process budget instead.")

    return args, extra_crispresso_args, command_str


def main():
    print("\nStarting TwInsPEctor...")
    args, extra_crispresso_args, command_str = parse_args()

    if args.save_pdf:
        plot_formats = ["png", "pdf"]
    elif args.no_save_png:
        plot_formats = ["pdf"]
    else:
        plot_formats = ["png"]

    parent_folder, crispresso_wt, crispresso_tpe, crispresso_composite_a, crispresso_composite_b, twinspector_results_folder = get_folder_names(args)
    text_output_dir = os.path.join(twinspector_results_folder, "d.text_files")
    os.makedirs(twinspector_results_folder, exist_ok=True)
    os.makedirs(crispresso_wt, exist_ok=True)
    os.makedirs(crispresso_tpe, exist_ok=True)
    os.makedirs(crispresso_composite_a, exist_ok=True)
    os.makedirs(crispresso_composite_b, exist_ok=True)
    os.makedirs(text_output_dir, exist_ok=True)

    with open(os.path.join(text_output_dir, "d9.twinspector_command.txt"), "w") as fout:
        fout.write(f"{command_str}\n")

    print("Analyzing reference inputs...")
    spacer_a, spacer_b = get_spacer_seqs(args.peg_spacers)
    if args.rt_templates:
        rt_template_a, rt_template_b = get_rt_templates(args.rt_templates)
    else:
        rt_template_a, rt_template_b = None, None
    reference_info = analyze_references(args.wt_seq, args.tpe_seq, spacer_a, spacer_b, rt_template_a, rt_template_b, cleavage_offset_a=args.cleavage_offset_a, cleavage_offset_b=args.cleavage_offset_b, output_root=text_output_dir, recoding_mode=args.recoding_mode)

    crispresso_jobs = min(4, args.n_processes)
    crispresso_processes = max(1, args.n_processes // crispresso_jobs)

    crispresso_cmd_wt = get_crispresso_command(args=args, extra_crispresso_args=extra_crispresso_args, ref_seq=args.wt_seq, ref_name="WT", spacer_a=reference_info["spacer_a_wt"], spacer_b=reference_info["spacer_b_wt"], crispresso_output_folder=crispresso_wt, twinspector_results_folder=text_output_dir, n_processes=crispresso_processes, append=False)
    crispresso_cmd_tpe = get_crispresso_command(args=args,extra_crispresso_args=extra_crispresso_args, ref_seq=args.tpe_seq, ref_name="TPE", spacer_a=reference_info["spacer_a_tpe"], spacer_b=reference_info["spacer_b_tpe"], crispresso_output_folder=crispresso_tpe, twinspector_results_folder=text_output_dir, n_processes=crispresso_processes)
    crispresso_cmd_composite_a = get_crispresso_command(args=args, extra_crispresso_args=extra_crispresso_args, ref_seq=reference_info["composite_a_ref_seq"], ref_name="Composite_A", spacer_a=reference_info["spacer_a_composite_a"], spacer_b=reference_info["spacer_b_composite_a"], crispresso_output_folder=crispresso_composite_a, twinspector_results_folder=text_output_dir, n_processes=crispresso_processes)
    crispresso_cmd_composite_b = get_crispresso_command(args=args, extra_crispresso_args=extra_crispresso_args, ref_seq=reference_info["composite_b_ref_seq"], ref_name="Composite_B", spacer_a=reference_info["spacer_a_composite_b"], spacer_b=reference_info["spacer_b_composite_b"], crispresso_output_folder=crispresso_composite_b, twinspector_results_folder=text_output_dir, n_processes=crispresso_processes)

    if args.n_processes == 1:
        print("Running CRISPResso2...")
        run_crispresso_command(crispresso_cmd_wt, verbose=args.verbose)
        run_crispresso_command(crispresso_cmd_tpe, verbose=args.verbose)
        run_crispresso_command(crispresso_cmd_composite_a, verbose=args.verbose)
        run_crispresso_command(crispresso_cmd_composite_b, verbose=args.verbose)
    else:
        print("Running CRISPResso2 in parallel...")
        crispresso_tasks = [crispresso_cmd_wt, crispresso_cmd_tpe, crispresso_cmd_composite_a, crispresso_cmd_composite_b]
        run_crispresso_commands_parallel(crispresso_tasks, n_processes=crispresso_jobs, verbose=args.verbose)

    print("Analyzing CRISPResso2 outputs...")
    df_merged = merge_crispresso_allele_tables(crispresso_wt=crispresso_wt, crispresso_tpe=crispresso_tpe, crispresso_composite_a=crispresso_composite_a, crispresso_composite_b=crispresso_composite_b)
    df_categorized, bp_changes_arrs = categorize_alleles(df_merged=df_merged, wt_seq=args.wt_seq, tpe_seq=args.tpe_seq, reference_info=reference_info, min_num_base_edits=args.min_num_base_edits, ignore_extraspacer_deletions=args.ignore_extraspacer_deletions, default_min_aln_score=args.default_min_aln_score, recoding_mode=args.recoding_mode)
    plotting_info = get_plotting_stats(df=df_categorized, reference_info=reference_info, bp_changes_arrs=bp_changes_arrs, twinspector_results_folder=text_output_dir, recoding_mode=args.recoding_mode, max_n_alleles_to_write=args.max_n_alleles_to_write)

    mutation_dicts = mutation_analysis(wt_seq_len=len(args.wt_seq), tpe_seq_len=len(args.tpe_seq), df_categorized=df_categorized, twinspector_results_folder=text_output_dir, rt_template_a=rt_template_a, rt_template_b=rt_template_b, pegRNA_intervals={"wt": reference_info["pegRNA_intervals_wt"], "tpe": reference_info["pegRNA_intervals_tpe"]}, ignore_extraspacer_deletions=args.ignore_extraspacer_deletions, recoding_mode=args.recoding_mode)

    if not args.no_figures:
        if not args.no_summary_figures:
            print("Generating summary bar plots...")
            plot_summary_barplots(category_counts=plotting_info["category_counts"], crispresso_output_folder=crispresso_wt, twinspector_results_folder=twinspector_results_folder, crispresso_wt=crispresso_wt, plot_formats=plot_formats)

        if not args.no_per_base_figures:
            print("Generating base bar plots...")
            base_plot_dir = os.path.join(twinspector_results_folder, "b.base_plots")
            os.makedirs(base_plot_dir, exist_ok=True)
            plot_per_base_pos_barplots(plotting_info=plotting_info, reference_info=reference_info, twinspector_results_folder=base_plot_dir, plot_formats=plot_formats, recoding_mode=args.recoding_mode, rt_templates=args.rt_templates, flap_data=mutation_dicts["flap_completion_counts_dict"])

        if not args.no_mutation_figures:
            print("Generating mutation bar plots...")
            mutation_plot_dir = os.path.join(twinspector_results_folder, "c.mutation_plots")
            os.makedirs(mutation_plot_dir, exist_ok=True)
            plot_mutation_barplots(mutation_dicts=mutation_dicts, category_counts=plotting_info["category_counts"], cut_points=[reference_info["cut_points_wt"], reference_info["cut_points_tpe"]], wt_seq_len=len(args.wt_seq), tpe_seq_len=len(args.tpe_seq), twinspector_results_folder=mutation_plot_dir, plot_formats=plot_formats, recoding_mode=args.recoding_mode)

        if not args.no_summary_figures:
            print("Generating summarized allele tables...")
            plot_summarized_allele_tables(
                args=args,
                df_categorized=df_categorized,
                reference_info=reference_info,
                output_folder=twinspector_results_folder,
                plot_formats=plot_formats,
                n_processes=args.n_processes,
            )

        if args.plot_extended_tables:
            allele_table_dir = os.path.join(twinspector_results_folder, "e.allele_tables")
            os.makedirs(allele_table_dir, exist_ok=True)
            print("Generating extended allele tables by category...")
            plot_ref_allele_tables(
                args=args,
                df_categorized=df_categorized,
                reference_info=reference_info,
                twinspector_results_folder=allele_table_dir,
                plot_formats=plot_formats,
                n_processes=args.n_processes,
            )

        if "png" in plot_formats:
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
