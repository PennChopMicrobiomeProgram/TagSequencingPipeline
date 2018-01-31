import sys
import os.path
import argparse
import itertools

class FastqRead(object):
    def __init__(self, read):
        self.desc, self.seq, self.qual = read
    
    def __repr__(self):
        return self.desc + "\n" + self.seq + "\n+\n" + self.qual + "\n"

def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return zip(*args)

def parse_fastq(f):
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield desc, seq, qual

def combine_barcodes_main(argv=None):
    p = argparse.ArgumentParser()
    
    # Input
    p.add_argument(
        "--work-dir", required=True,
        help="Directory where the files are present")
    p.add_argument(
        "--I1-fp", required=False,
        default="Undetermined_S0_L001_I1_001.fastq",
        help="I1 FASTQ files")
    p.add_argument(
        "--I2-fp", required=False,
        default="Undetermined_S0_L001_I2_001.fastq",
        help="I2 FASTQ files")
    # Output
    p.add_argument(
        "--out-fp", required=False,
        default="barcodes.fastq",
        help="")
    args = p.parse_args(argv)
    
    I1 = os.path.join(args.work_dir, args.I1_fp)
    I2 = os.path.join(args.work_dir, args.I2_fp)
    I = os.path.join(args.work_dir, args.out_fp)

    with open(I, "w") as f_out, open(I1,'r') as I1_handle, open(I2, 'r') as I2_handle:
        fwds = (FastqRead(x) for x in parse_fastq(I1_handle))
        revs = (FastqRead(x) for x in parse_fastq(I2_handle))
        for fwd, rev in zip(fwds, revs):
            f_out.write("@%s\n%s\n+\n%s\n" % (fwd.desc, fwd.seq+rev.seq, fwd.qual+rev.qual))

if __name__ == "__main__":
    combine_barcodes_main()
