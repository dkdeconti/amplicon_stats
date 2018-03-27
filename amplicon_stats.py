#!/usr/bin/env python

import argparse
from collections import defaultdict
import numpy
import matplotlib.pyplot as plt
import pysam
import matplotlib
import re


def plot_depth(pileups, locus, bams, buf):
    """Plots depth across given loci."""
    chrom, begin, end = locus
    begin = int(begin)
    filenames = [re.sub('\.bam$', '', bam) for bam in bams]
    plot_file = '-'.join(locus) + '.png'
    fig, axarray = plt.subplots(nrows=len(pileups),
                            ncols=1,
                            figsize=(16, 4*len(pileups)),
                            sharex=True)
    if len(pileups) == 1:
        pileup = pileups[0]
        axarray.plot(range(begin - buf, len(pileup) + begin - buf), pileup)
        axarray.axvline(begin, color='r')
        axarray.axvline(end, color='r')
        axarray.set_title(filenames[0])
    else:
        for i, pileup in enumerate(pileups):
            # TODO add a filepath
            axarray[i].plot(range(begin - buf, len(pileup) + begin - buf),
                            pileup)
            axarray[i].axvline(begin, color='r')
            axarray[i].axvline(end, color='r')
            axarray[i].set_title(filenames[i])
    fig.savefig(plot_file)
    plt.close()
    return plot_file

        
def parse_bam_for_pileup(bams, loci, buf):
    """Parses BAM for pileup data into map by locus key."""
    pileups = {}
    for locus in loci:
        chrom = locus[0]
        begin = int(locus[1]) -1
        end = int(locus[2])
        key = tuple(locus)
        pileups[key] = [get_pileup_vector(bam, chrom, begin, end, buf)
                        for bam in bams]
    return pileups


def get_pileup_vector(bam, chrom, begin, end, buf):
    """Get depth vector from bam given pos."""
    depth = [0]*((end + buf) - (begin - buf))
    for p in pysam.AlignmentFile(bam, 'rb').pileup(chrom, begin-buf, end+buf):
        depth[p.reference_pos - begin] = p.nsegments
    return depth
    #return [p.nsegments
    #        for p in pysam.AlignmentFile(bam, 'rb').pileup(chrom,
    #                                                       begin - buf,
    #                                                       end + buf)]


def main():
    """"""
    # arg parsing
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--loci", required=True,
                        nargs=3, metavar=("chr", "begin", "end"),
                        action="append", help="comma-delimited loci [required]")
    parser.add_argument("--buf", metavar="BUFFER", type=int, default=1000,
                        help="buffer about the loci to include in plots")
    parser.add_argument("-b", "--bam", metavar="BAM", nargs='+',
                        required=True, help="BAM file(s) [required]")
    args = parser.parse_args()
    # central dispatch
    print args.bam
    print args.loci
    pileups = parse_bam_for_pileup(args.bam, args.loci, args.buf)
    for locus in args.loci:
        plot_depth(pileups[tuple(locus)], locus, args.bam, args.buf)


if __name__ == "__main__":
    main()
