# REQUIRED ARGUMENTS:
# -r1 <input fastq file>
# -w  <wildtype sequence>
# -t  <twinPE edited sequence>
# -g  <comma-separated pegRNA spacer sequences>

# OPTIONAL ARGUMENTS:
# -r2   <second input fastq for paired-end data>
# -o    <output directory>
# -ne   <minimum number of programmed bases that must be edited for a read to be classified> [default: 2]
# -dmas <default minimum homology score for CRISPResso2 to align read> [default: 50]
# -mfa  <minimum percent read frequency required to report an allele in the alleles tables>
# -mnr  <maximum number of allele rows to display in the allele tables>

# FLAGS:
# -ied [ignore deletions beyond the spacers - use if reads were truncated and are being misclassified with deletions]
# -nat [used to decrease run time if allele tables are not needed]
# -rcm [run in recoding mode if programmed edit should be treated as base substitutions only]
# -pfr [plot entire read sequences in allele tables]
# -nrr [skip rerunning CRISPResso2 if a run using the same parameters has already been completed]

# EXAMPLE USAGE: EXPT5-PE6c-AAVS1-PBS11-NM-REP1_S73_L001_R1_001.fastq (trimmed front 4 bases)
# Script:
python twin_prime_editing/TwInsPEctor/TwInsPEctor.py -r1 03_nate_practicum/data/scripps_collab/trimmed/EXPT5-PE6c-AAVS1-PBS11-NM-REP1_S73_L001_R1_001.fastq -w CCAGGATCAGTGAAACGCACCAGACAGCCGCGTCAGAGCAGCTCAGGTTCTGGGAGAGGGTAGCGCAGGGTGGCCACTGAGAACCGGGCAGGTCACGCATCCCCCCCTTCCCTCCCACCCCCTGCCAAGCTCTCCCTCCCAGGATCCTCTCTGGCTCCATCGTAAGCAAACCTTAGAGGTTCTGGCAAGGAGAGAGATGGCTCCAGGAAAT -t CCAGGATCAGTGAAACGCACCAGACAGCCGCGTCAGAGCAGCTCAGGTTCTGGGAGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCTTCCTCTCTGGCTCCATCGTAAGCAAACCTTAGAGGTTCTGGCAAGGAGAGAGATGGCTCCAGGAAAT -g AGCTCAGGTTCTGGGAGAGG,GATGGAGCCAGAGAGGATCC -o /uufs/chpc.utah.edu/common/home/u0493285/twinpe/demo_no_tables -dmas 30 -ied -nat
# Conda package:
TwInsPEctor -r1 03_nate_practicum/data/scripps_collab/trimmed/EXPT5-PE6c-AAVS1-PBS11-NM-REP1_S73_L001_R1_001.fastq -w CCAGGATCAGTGAAACGCACCAGACAGCCGCGTCAGAGCAGCTCAGGTTCTGGGAGAGGGTAGCGCAGGGTGGCCACTGAGAACCGGGCAGGTCACGCATCCCCCCCTTCCCTCCCACCCCCTGCCAAGCTCTCCCTCCCAGGATCCTCTCTGGCTCCATCGTAAGCAAACCTTAGAGGTTCTGGCAAGGAGAGAGATGGCTCCAGGAAAT -t CCAGGATCAGTGAAACGCACCAGACAGCCGCGTCAGAGCAGCTCAGGTTCTGGGAGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCTTCCTCTCTGGCTCCATCGTAAGCAAACCTTAGAGGTTCTGGCAAGGAGAGAGATGGCTCCAGGAAAT -g AGCTCAGGTTCTGGGAGAGG,GATGGAGCCAGAGAGGATCC -o /uufs/chpc.utah.edu/common/home/u0493285/twinpe/NM_test -dmas 30 -ied -nat





















# Test fastq:
# python twin_prime_editing/TwInsPEctor/TwInsPEctor.py -r1 twin_prime_editing/demo.fastq.gz -w CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGGAGTCCATTAGCAGAAGTGGCCCTTGGCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC -t CGCCGGGAACTGCCGCTGGCCCCCCACCGCCAGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC -g GCTGGCCCCCCACCGCCCCA GCGACTCCTGGAAGTGGCCA -o /uufs/chpc.utah.edu/common/home/u0493285/twinpe/demo/ -dmas 10 --no_allele_tables --no_rerun
