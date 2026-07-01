import gzip
import argparse


def revcomp(seq):
    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(complement)[::-1]


def filter_paired_fastq(
    r1_in,
    r2_in,
    r1_out,
    r2_out,
    seq_start,
    seq_end,
):
    rc_end = revcomp(seq_end)

    total_pairs = 0
    kept_pairs = 0

    with gzip.open(r1_in, "rt") as f1, \
         gzip.open(r2_in, "rt") as f2, \
         gzip.open(r1_out, "wt") as o1, \
         gzip.open(r2_out, "wt") as o2:

        while True:
            r1 = [f1.readline() for _ in range(4)]
            r2 = [f2.readline() for _ in range(4)]

            if not r1[0]:
                break

            total_pairs += 1

            r1_seq = r1[1].rstrip()
            r2_seq = r2[1].rstrip()

            if (
                r1_seq.startswith(seq_start)
                and
                r2_seq.startswith(rc_end)
            ):
                o1.writelines(r1)
                o2.writelines(r2)
                kept_pairs += 1

    return {
        "total_pairs": total_pairs,
        "kept_pairs": kept_pairs,
        "filtered_pairs": total_pairs - kept_pairs,
    }
def main():
    parser = argparse.ArgumentParser(
        description="Filter paired-end FASTQ.gz files."
    )

    parser.add_argument("--r1", required=True,
                        help="Input R1 FASTQ.gz")
    parser.add_argument("--r2", required=True,
                        help="Input R2 FASTQ.gz")

    parser.add_argument("--r1-out", required=True,
                        help="Output filtered R1 FASTQ.gz")
    parser.add_argument("--r2-out", required=True,
                        help="Output filtered R2 FASTQ.gz")

    parser.add_argument("--seq-start", required=True,
                        help="Required R1 prefix")
    parser.add_argument("--seq-end", required=True,
                        help="Required R1 suffix")

    args = parser.parse_args()

    stats = filter_paired_fastq(
        r1_in=args.r1,
        r2_in=args.r2,
        r1_out=args.r1_out,
        r2_out=args.r2_out,
        seq_start=args.seq_start,
        seq_end=args.seq_end,
    )

    print(
        f"Processed {stats['total_pairs']:,} read pairs\n"
        f"Kept {stats['kept_pairs']:,} pairs "
        f"({100 * stats['kept_pairs'] / stats['total_pairs']:.2f}%)"
    )


if __name__ == "__main__":
    main()
