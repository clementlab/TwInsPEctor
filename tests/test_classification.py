"""
Synthetic twin prime editing (twinPE) classification tests for TwInsPEctor.

Covers all 10 TwInsPEctor allele categories across three twinPE architectures:
  1. Pure insertion at nick A (05_abdalla locus: all 29 bp of WT between nick
     sites are shared with the 3' end of the TPE edit window).
  2. Standard sequence replacement (demo/demo.md locus: 116 bp WT replaced by
     50 bp programmed sequence).
  3. Pure insertion at nick B (mirror case: entire WT edit window shared with
     the 5' start of the TPE edit window).

Usage:
  # Run unit tests directly (uses built-in affine-gap aligner + categorize_alleles)
  python tests/test_classification.py

  # Or with pytest
  pytest tests/test_classification.py

  # Optional: write synthetic FASTQ + run script for end-to-end CRISPResso2 testing
  python tests/test_classification.py --write-fastq /tmp/twinspector_test_fastq
"""

import argparse
import gzip
import os
import sys
import tempfile
import unittest

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from TwInsPEctor.TwInsPEctor import analyze_references, categorize_alleles, reverse_complement


# ---------------------------------------------------------------------------
# Pure-Python affine-gap global aligner (mirrors CRISPResso2 Needleman-Wunsch)
# ---------------------------------------------------------------------------
def global_align(read_seq, ref_seq, match=2, mismatch=-2, gap_open=-5, gap_extend=-1):
    """Affine-gap global alignment returning (aligned_read, aligned_ref, pct_identity)."""
    m, n = len(read_seq), len(ref_seq)
    neg_inf = -10**9

    # M[i][j]: aligned read[i-1] to ref[j-1]
    # X[i][j]: read[i-1] aligned to gap in ref (insertion in read)
    # Y[i][j]: gap in read aligned to ref[j-1] (deletion in read)
    M = [[neg_inf] * (n + 1) for _ in range(m + 1)]
    X = [[neg_inf] * (n + 1) for _ in range(m + 1)]
    Y = [[neg_inf] * (n + 1) for _ in range(m + 1)]
    trM = [[0] * (n + 1) for _ in range(m + 1)]
    trX = [[0] * (n + 1) for _ in range(m + 1)]
    trY = [[0] * (n + 1) for _ in range(m + 1)]

    M[0][0] = 0
    for i in range(1, m + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
        trX[i][0] = 1 if i > 1 else 0
    for j in range(1, n + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend
        trY[0][j] = 2 if j > 1 else 0

    for i in range(1, m + 1):
        ri = read_seq[i - 1]
        Mi, Mi1 = M[i], M[i - 1]
        Xi, Xi1 = X[i], X[i - 1]
        Yi, Yi1 = Y[i], Y[i - 1]
        tMi, tXi, tYi = trM[i], trX[i], trY[i]
        for j in range(1, n + 1):
            # X: gap in ref
            o_x = Mi1[j] + gap_open
            e_x = Xi1[j] + gap_extend
            y_x = Yi1[j] + gap_open
            if e_x >= o_x and e_x >= y_x:
                Xi[j], tXi[j] = e_x, 1
            elif o_x >= y_x:
                Xi[j], tXi[j] = o_x, 0
            else:
                Xi[j], tXi[j] = y_x, 2

            # Y: gap in read
            o_y = Mi[j - 1] + gap_open
            e_y = Yi[j - 1] + gap_extend
            x_y = Xi[j - 1] + gap_open
            if e_y >= o_y and e_y >= x_y:
                Yi[j], tYi[j] = e_y, 2
            elif o_y >= x_y:
                Yi[j], tYi[j] = o_y, 0
            else:
                Yi[j], tYi[j] = x_y, 1

            # M: match / mismatch
            sc = match if ri == ref_seq[j - 1] else mismatch
            best_prev = Mi1[j - 1]
            src = 0
            if Xi1[j - 1] > best_prev:
                best_prev = Xi1[j - 1]
                src = 1
            if Yi1[j - 1] > best_prev:
                best_prev = Yi1[j - 1]
                src = 2
            Mi[j] = best_prev + sc
            tMi[j] = src

    # Traceback
    state = 0
    best = M[m][n]
    if X[m][n] > best:
        best, state = X[m][n], 1
    if Y[m][n] > best:
        best, state = Y[m][n], 2

    aln_read, aln_ref = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if state == 0:
            prev = trM[i][j]
            aln_read.append(read_seq[i - 1])
            aln_ref.append(ref_seq[j - 1])
            i -= 1
            j -= 1
            state = prev
        elif state == 1:
            prev = trX[i][j]
            aln_read.append(read_seq[i - 1])
            aln_ref.append("-")
            i -= 1
            state = prev
        else:
            prev = trY[i][j]
            aln_read.append("-")
            aln_ref.append(ref_seq[j - 1])
            j -= 1
            state = prev

    aln_read = "".join(reversed(aln_read))
    aln_ref = "".join(reversed(aln_ref))
    matches = sum(1 for a, b in zip(aln_read, aln_ref) if a == b and a != "-")
    pct_identity = 100.0 * matches / len(aln_ref) if aln_ref else 0.0
    return aln_read, aln_ref, pct_identity


def classify_synthetic_reads(wt_seq, tpe_seq, spacer_a, spacer_b, sample_reads, min_num_base_edits=3):
    """
    Align synthetic sample reads against WT, TPE, Composite A, and Composite B
    and run TwInsPEctor.categorize_alleles().

    sample_reads: list of (label, read_seq, expected_category)
    Returns: dict mapping label -> (observed_category, classified_by)
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_info = analyze_references(
            wt_seq=wt_seq,
            tpe_seq=tpe_seq,
            spacer_a=spacer_a,
            spacer_b=spacer_b,
            rt_template_a=None,
            rt_template_b=None,
            cleavage_offset_a=-3,
            cleavage_offset_b=-3,
            output_root=tmpdir,
            recoding_mode=False,
        )

    comp_a_ref = ref_info["composite_a_ref_seq"]
    comp_b_ref = ref_info["composite_b_ref_seq"]

    rows = []
    label_by_key = {}
    for label, read_seq, _ in sample_reads:
        key = min(read_seq, reverse_complement(read_seq))
        label_by_key[key] = label

        ar_wt, ref_wt, sc_wt = global_align(read_seq, wt_seq)
        ar_tpe, ref_tpe, sc_tpe = global_align(read_seq, tpe_seq)
        ar_ca, ref_ca, sc_ca = global_align(read_seq, comp_a_ref)
        ar_cb, ref_cb, sc_cb = global_align(read_seq, comp_b_ref)

        rows.append({
            "sequence_key": key,
            "#Reads": 10,
            "%Reads": 100.0 / len(sample_reads),
            "Aligned_Sequence_wt": ar_wt,
            "Reference_Sequence_wt": ref_wt,
            "Aligned_Reference_Score_wt": sc_wt,
            "Aligned_Sequence_tpe": ar_tpe,
            "Reference_Sequence_tpe": ref_tpe,
            "Aligned_Reference_Score_tpe": sc_tpe,
            "Aligned_Sequence_comp_a": ar_ca,
            "Reference_Sequence_comp_a": ref_ca,
            "Aligned_Reference_Score_comp_a": sc_ca,
            "Aligned_Sequence_comp_b": ar_cb,
            "Reference_Sequence_comp_b": ref_cb,
            "Aligned_Reference_Score_comp_b": sc_cb,
        })

    df_merged = pd.DataFrame(rows)
    df_cat, _ = categorize_alleles(
        df_merged=df_merged,
        wt_seq=wt_seq,
        tpe_seq=tpe_seq,
        reference_info=ref_info,
        min_num_base_edits=min_num_base_edits,
        ignore_extraspacer_deletions=False,
        default_min_aln_score=30,
        recoding_mode=False,
    )

    results = {}
    for _, row in df_cat.iterrows():
        label = label_by_key[row["sequence_key"]]
        results[label] = (row["Category_final"], row["Classified_by"])
    return results


# ---------------------------------------------------------------------------
# Locus 1: 05_abdalla (pure 42-bp insertion at nick A; all 29 WT bp shared at 3')
# ---------------------------------------------------------------------------
ABDALLA_WT = (
    "AGAGCCAGATGCTGAGCCCGCAGAGGGAGAGGAAGAGGAGGCCCAGCCAGCTGCCACAGATGCCCCCAACGGAGACCTCACGCCCGAGATGATCCTCAGCATGATGGACCGG"
    "TGATGGCGGAGCCTTGTGCGTCGCCAGGA"
    "CTTCTCTGGGCTGTGTTTAAACGGCCCGCATCTTAATTTTTCTCCCTTCTTTCTTTTTTTGGCTTTGGGAAAAGCATCATTTTACCAAACATACCGA"
)
ABDALLA_TPE = (
    "AGAGCCAGATGCTGAGCCCGCAGAGGGAGAGGAAGAGGAGGCCCAGCCAGCTGCCACAGATGCCCCCAACGGAGACCTCACGCCCGAGATGATCCTCAGCATGATGGACCGG"
    "GGAGCAACAACGCTGCTTCGGCCTGGCAGTGGTTCGGGAGGT"
    "TGATGGCGGAGCCTTGTGCGTCGCCAGGA"
    "CTTCTCTGGGCTGTGTTTAAACGGCCCGCATCTTAATTTTTCTCCCTTCTTTCTTTTTTTGGCTTTGGGAAAAGCATCATTTTACCAAACATACCGA"
)
ABDALLA_SPACER_A = "TCAGCATGATGGACCGGTGA"
ABDALLA_SPACER_B = "AACACAGCCCAGAGAAGTCC"

_ABD_PREFIX = ABDALLA_WT[:112]
_ABD_WT_MID = ABDALLA_WT[112:141]      # TGATGGCGGAGCCTTGTGCGTCGCCAGGA (29 bp, shared at 3' end of TPE)
_ABD_INS = ABDALLA_TPE[112:154]        # GGAGCAACAACGCTGCTTCGGCCTGGCAGTGGTTCGGGAGGT (42 bp novel insert)
_ABD_SUFFIX = ABDALLA_WT[141:]

ABDALLA_SAMPLE_READS = [
    ("WT_exact", _ABD_PREFIX + _ABD_WT_MID + _ABD_SUFFIX, "WT"),
    ("WT_sub_outside_nicks", "T" + _ABD_PREFIX[1:] + _ABD_WT_MID + _ABD_SUFFIX, "WT"),
    ("Imperfect_WT_partial_del", _ABD_PREFIX + _ABD_WT_MID[:14] + _ABD_SUFFIX, "Imperfect_WT"),
    ("Imperfect_WT_sub_in_window", _ABD_PREFIX + _ABD_WT_MID[:10] + "A" + _ABD_WT_MID[11:] + _ABD_SUFFIX, "Imperfect_WT"),
    ("Null_complete_deletion", _ABD_PREFIX + _ABD_SUFFIX, "Null"),
    ("Perfect_TPE_exact", _ABD_PREFIX + _ABD_INS + _ABD_WT_MID + _ABD_SUFFIX, "Perfect_TPE"),
    ("Perfect_TPE_sub_outside", "T" + _ABD_PREFIX[1:] + _ABD_INS + _ABD_WT_MID + _ABD_SUFFIX, "Perfect_TPE"),
    ("Dual_Flap_internal_del", _ABD_PREFIX + _ABD_INS[:12] + _ABD_INS[-12:] + _ABD_WT_MID + _ABD_SUFFIX, "Dual_Flap"),
    ("Flap_A_only", _ABD_PREFIX + _ABD_INS[:16] + _ABD_SUFFIX, "Flap_A"),
    ("Flap_A_Hybrid", _ABD_PREFIX + _ABD_INS[:16] + _ABD_WT_MID + _ABD_SUFFIX, "Flap_A_Hybrid"),
    ("Flap_B_only", _ABD_PREFIX + _ABD_INS[-16:] + _ABD_WT_MID + _ABD_SUFFIX, "Flap_B"),
    ("Flap_B_Hybrid", _ABD_PREFIX + _ABD_WT_MID + _ABD_INS[-16:] + _ABD_WT_MID + _ABD_SUFFIX, "Flap_B_Hybrid"),
    ("Imperfect_TPE_middle_only", _ABD_PREFIX + _ABD_INS[10:28] + _ABD_SUFFIX, "Imperfect_TPE"),
]


# ---------------------------------------------------------------------------
# Locus 2: demo/demo.md (standard replacement: 116 bp WT -> 50 bp TPE)
# ---------------------------------------------------------------------------
DEMO_WT = (
    "CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGGAGTCCATTAGCAGAAGTGGCCCTTGGCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC"
)
DEMO_TPE = (
    "CGCCGGGAACTGCCGCTGGCCCCCCACCGCCAGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC"
)
DEMO_SPACER_A = "GCTGGCCCCCCACCGCCCCA"
DEMO_SPACER_B = "GCGACTCCTGGAAGTGGCCA"


def _build_replacement_reads(wt_seq, tpe_seq, spacer_a, spacer_b):
    with tempfile.TemporaryDirectory() as tmpdir:
        info = analyze_references(wt_seq, tpe_seq, spacer_a, spacer_b, None, None, -3, -3, tmpdir)
    nick_a_wt, nick_b_wt = info["cut_points_wt"]
    prefix = wt_seq[: nick_a_wt + 1]
    wt_mid = info["deleted_seq"]
    tpe_mid = info["inserted_seq"]
    suffix = wt_seq[nick_b_wt + 1 :]
    s_start = info["num_bases_shared_start"]
    s_end = info["num_bases_shared_end"]

    return [
        ("WT_exact", prefix + wt_mid + suffix, "WT"),
        ("Imperfect_WT_partial_del", prefix + wt_mid[:25] + suffix, "Imperfect_WT"),
        ("Null_complete_deletion", prefix + suffix, "Null"),
        ("Perfect_TPE_exact", prefix + tpe_mid + suffix, "Perfect_TPE"),
        ("Dual_Flap_internal_del", prefix + tpe_mid[:14] + tpe_mid[-14:] + suffix, "Dual_Flap"),
        ("Flap_A_only", prefix + tpe_mid[: s_start + 14] + suffix, "Flap_A"),
        ("Flap_A_Hybrid", prefix + tpe_mid[: s_start + 14] + wt_mid[-25:] + suffix, "Flap_A_Hybrid"),
        ("Flap_B_only", prefix + tpe_mid[-(s_end + 14) :] + suffix, "Flap_B"),
        ("Flap_B_Hybrid", prefix + wt_mid[:25] + tpe_mid[-(s_end + 14) :] + suffix, "Flap_B_Hybrid"),
        ("Imperfect_TPE_middle_only", prefix + tpe_mid[s_start + 8 : len(tpe_mid) - s_end - 8] + suffix, "Imperfect_TPE"),
    ]


# ---------------------------------------------------------------------------
# Locus 3: Mirror pure insertion at nick B (all 29 WT bp shared at 5' of TPE)
# ---------------------------------------------------------------------------
MIRROR_TPE = _ABD_PREFIX + _ABD_WT_MID + _ABD_INS + _ABD_SUFFIX
MIRROR_SAMPLE_READS = [
    ("WT_exact", _ABD_PREFIX + _ABD_WT_MID + _ABD_SUFFIX, "WT"),
    ("Imperfect_WT_partial_del", _ABD_PREFIX + _ABD_WT_MID[-14:] + _ABD_SUFFIX, "Imperfect_WT"),
    ("Null_complete_deletion", _ABD_PREFIX + _ABD_SUFFIX, "Null"),
    ("Perfect_TPE_exact", _ABD_PREFIX + _ABD_WT_MID + _ABD_INS + _ABD_SUFFIX, "Perfect_TPE"),
    ("Dual_Flap_internal_del", _ABD_PREFIX + _ABD_WT_MID + _ABD_INS[:12] + _ABD_INS[-12:] + _ABD_SUFFIX, "Dual_Flap"),
    ("Flap_A_only", _ABD_PREFIX + _ABD_WT_MID + _ABD_INS[:16] + _ABD_SUFFIX, "Flap_A"),
    ("Flap_A_Hybrid", _ABD_PREFIX + _ABD_WT_MID + _ABD_INS[:16] + _ABD_WT_MID + _ABD_SUFFIX, "Flap_A_Hybrid"),
    ("Flap_B_only", _ABD_PREFIX + _ABD_INS[-16:] + _ABD_SUFFIX, "Flap_B"),
    ("Flap_B_Hybrid", _ABD_PREFIX + _ABD_WT_MID + _ABD_INS[-16:] + _ABD_SUFFIX, "Flap_B_Hybrid"),
    ("Imperfect_TPE_middle_only", _ABD_PREFIX + _ABD_INS[10:28] + _ABD_SUFFIX, "Imperfect_TPE"),
]


class TestTwinPEClassification(unittest.TestCase):
    def _assert_all_classified(self, wt_seq, tpe_seq, spacer_a, spacer_b, sample_reads):
        observed = classify_synthetic_reads(wt_seq, tpe_seq, spacer_a, spacer_b, sample_reads)
        for label, _, expected in sample_reads:
            obs_cat, obs_by = observed[label]
            self.assertEqual(
                obs_cat,
                expected,
                msg=f"{label}: expected '{expected}', got '{obs_cat}' (classified by {obs_by})",
            )

    def test_pure_insertion_at_nick_a_abdalla(self):
        self._assert_all_classified(
            ABDALLA_WT, ABDALLA_TPE, ABDALLA_SPACER_A, ABDALLA_SPACER_B, ABDALLA_SAMPLE_READS
        )

    def test_standard_replacement_demo(self):
        demo_reads = _build_replacement_reads(DEMO_WT, DEMO_TPE, DEMO_SPACER_A, DEMO_SPACER_B)
        self._assert_all_classified(
            DEMO_WT, DEMO_TPE, DEMO_SPACER_A, DEMO_SPACER_B, demo_reads
        )

    def test_pure_insertion_at_nick_b_mirror(self):
        self._assert_all_classified(
            ABDALLA_WT, MIRROR_TPE, ABDALLA_SPACER_A, ABDALLA_SPACER_B, MIRROR_SAMPLE_READS
        )


def write_synthetic_fastq(out_dir, reads_per_allele=10):
    """Write a gzipped FASTQ of the 05_abdalla synthetic alleles and a helper run script."""
    os.makedirs(out_dir, exist_ok=True)
    fq_path = os.path.join(out_dir, "abdalla_synthetic_10cat.fastq.gz")
    with gzip.open(fq_path, "wt") as fout:
        for label, seq, expected in ABDALLA_SAMPLE_READS:
            qual = "I" * len(seq)
            for rep in range(reads_per_allele):
                fout.write(f"@{label}_{expected}_rep{rep}\n{seq}\n+\n{qual}\n")

    sh_path = os.path.join(out_dir, "run_synthetic_test.sh")
    with open(sh_path, "w") as fout:
        fout.write(
            "#!/bin/bash\n"
            f"TwInsPEctor -r1 {fq_path} \\\n"
            f"  -w {ABDALLA_WT} \\\n"
            f"  -t {ABDALLA_TPE} \\\n"
            f"  -g {ABDALLA_SPACER_A},{ABDALLA_SPACER_B} \\\n"
            f"  -o {os.path.join(out_dir, 'TwInsPEctor_synthetic_out')}\n"
        )
    print(f"Wrote synthetic FASTQ ({len(ABDALLA_SAMPLE_READS) * reads_per_allele} reads) to: {fq_path}")
    print(f"Wrote run script to: {sh_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--write-fastq", type=str, default=None)
    known, remaining = parser.parse_known_args()
    if known.write_fastq:
        write_synthetic_fastq(known.write_fastq)
    else:
        unittest.main(argv=[sys.argv[0]] + remaining, verbosity=2)
