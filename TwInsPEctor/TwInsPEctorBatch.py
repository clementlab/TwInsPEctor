import argparse
import os
import pandas as pd
import multiprocessing as mp

def main(args):
    input_batch_file = args.batch_file
    if not os.path.isfile(input_batch_file):
        print(f"Error: Batch file '{input_batch_file}' does not exist.")
        return
    
    batch_df = pd.read_csv(input_batch_file, sep='\t')

    for col in batch_df.columns:
        if col not in args.__dict__.keys():
            print(f"Warning: Column '{col}' in batch file is not a recognized argument and will be ignored.")

    required_args = ['fastq_r1', 'wt_seq', 'twin_seq', 'peg_spacers']
    for arg in required_args:
        if arg not in batch_df.columns and getattr(args, arg) is None:
            print(f"Error: Required argument '{arg}' is missing from both command line and batch file.")
            return

    commands = []
    for index, row in batch_df.iterrows():
        command = 'TwInsPEctor'
        for arg_name in args.__dict__.keys():
            arg_value = args.__dict__[arg_name]
            #if arg is provided in file, overwrite the command line argument
            if arg_name in batch_df.columns:
                    arg_value = row[arg_name]
                    if isinstance(arg_value, str) and arg_value.lower() == "true":
                        arg_value = True

            if isinstance(arg_value, bool) and arg_value: # if arg is a boolean
                command += f' --{arg_name.replace("_", "-")}'
            else:
                command += f' --{arg_name.replace("_", "-")} {arg_value}'

            command = command.replace("  ", " ") # remove double spaces if any
    
    #run commands in parallel using multiprocessing
    pool = mp.Pool(mp.cpu_count())
    pool.map(os.system, commands)
    pool.close()
    pool.join() 

    print(f'Finished running TwInsPEctor on {len(commands)} samples.')
    






def parse_args():
    prog = 'TwInsPEctorBatch'
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Analyzes Twin Prime Editing outcomes in batches by running TwInsPEctor on each sample and summarizing the results.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument( '-f', '--batch_file', help='Path to the batch file containing sample information (tab-delimited with columns: sample_name, fastq_r1, fastq_r2, wt_seq, twin_seq, peg_spacers).', required=True)
    parser.add_argument("-w", "--wt_seq", type=str, help="Full wildtype reference amplicon sequence including spacers.", default=None)
    parser.add_argument("-t", "--twin_seq", type=str, help="Full TwinPE reference amplicon sequence with 5' & 3' ends identical to wildtype reference amplicon.", default=None)
    parser.add_argument("-g", "--peg_spacers", type=str, help="Comma-separated pegRNA spacer sequences: <spacer A>,<spacer B>", default=None)
    parser.add_argument("-o", "--output_root", type=str, default=None, help="Root output folder for CRISPResso2 and TwinPE 8cat results. If not provided, a folder will be created in the current working directory based on the input fastq file names.")
    parser.add_argument("-ne", "--num_changes_to_check", type=int, default=2, help="Minimum number of programmed bases that must be edited for read to be classified.")
    parser.add_argument("-rcm", "--recoding_mode", action="store_true", help="Run in recoding mode when there are only base substitutions.")
    parser.add_argument("-dmas", "--default_min_aln_score", type=int, default=30, help="Default minimum homology score for a read to align to the compound reference amplicon")
    parser.add_argument("-pfr", "--plot_full_reads", action="store_true", help="Allele tables will display full read sequences.")
    parser.add_argument("-ied", "--ignore_extraspacer_deletions", action="store_true", help="Classification ignores deletions occurring beyond the spacers (outside edit window).")
    parser.add_argument("-nat", "--no_allele_tables", action="store_true", help="Decrease run time when allele table are not wanted.")
    parser.add_argument("-pp", "--produce_png", action="store_true", help="Produce PNG versions of all plots in addition to PDF versions.")
    parser.add_argument("-mfa", "--min_frequency_alleles", type=float, default=0.0, help="Minimum percent read frequency required to report an allele in the alleles tables.")
    parser.add_argument("-mnr", "--max_n_rows", type=int, default=50, help="Maximum number of allele rows to display in the allele tables.")
    parser.add_argument("-nrr", "--no_rerun", action="store_true", help="Don't rerun CRISPResso2 if a run using the same parameters has already been finished.")
    parser.add_argument("-kco", "--keep_crispresso_outputs", action="store_true", help="Don't delete CRISPResso2 output folders after analysis.")
    parser.add_argument("-ts", "--trim_string", type=str, default=None, help="String to trim reads using fastp with override options within CRISPResso2 before analysis.")
    parser.add_argument("-fp", "--fastp_command", type=str, default=None, help="Command to run fastp for read trimming within CRISPResso2 before analysis.")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-V", "--version", action="version", version="%(prog)s 1.0")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args)
