#!/usr/bin/env python

import argparse
import jinja2
import numpy
import matplotlib.pyplot as plt
import os
import pysam
import shutil
import matplotlib
import re


def create_report(plots_map, stats_map):
    """Inserts files into template HTML."""
    this_dir = os.path.dirname(os.path.realpath(__file__))
    lib_dir = os.path.join(this_dir, "lib")
    template_filepath = os.path.join(this_dir, "template.html")
    # TODO add dynamic report_dir addition
    report_dir = os.getcwd()
    lib_destination = os.path.join(report_dir, "lib")
    report = os.path.join(report_dir, "html_report.html")
    #report = '/'.join([report_dir, "html_report.html"])
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(this_dir))
    #template = env.get_template(template_filepath)
    template = env.get_template("template.html")
    context = {"loci_depths": plots_map,
               "bam_stats": stats_map}
    with open(report, 'w') as outfile:
        outfile.write(template.render(context))
    shutil.copytree(lib_dir, lib_destination)


def plot_depth_across_locus(pileups, locus, bams, buf):
    """Plots depth across given loci."""
    chrom, begin, end = locus
    begin = int(begin)
    end = int(end)
    filenames = [re.sub('\.bam$', '', bam) for bam in bams]
    plot_file = '-'.join(locus) + '.across_whole_locus_depth.png'
    # TODO add a file path
    fig, axarray = plt.subplots(nrows=len(pileups),
                                ncols=1,
                                figsize=(16, 4*len(pileups)),
                                sharex=True,
                                sharey=True)
    if len(pileups) == 1:
        pileup = pileups[0]
        axarray.plot(range(begin - buf, len(pileup) + begin - buf), pileup)
        axarray.axvline(begin, color='r')
        axarray.axvline(end, color='r')
        axarray.set_title(filenames[0])
        axarray.set_ylabel("depth")
    else:
        for i, pileup in enumerate(pileups):
            axarray[i].plot(range(begin - buf, len(pileup) + begin - buf),
                            pileup)
            axarray[i].axvline(begin, color='r')
            axarray[i].axvline(end, color='r')
            axarray[i].set_title(filenames[i])
            axarray[i].set_ylabel("depth")
    plt.xlabel("genomic position")
    fig.savefig(plot_file)
    plt.close()
    return plot_file


def plot_depth_about_locus_ends(pileups, locus, bams, buf):
    """Plots depth about locus endpoints +/- buffer."""
    chrom, begin, end = locus
    begin = int(begin)
    end = int(end)
    filenames = [re.sub('\.bam$', '', bam) for bam in bams]
    plot_file = '-'.join(locus) + '.cut_site_boundary_depth.png'
    fig, axarray = plt.subplots(nrows=len(pileups),
                                ncols=2,
                                figsize=(16, 4*len(pileups)),
                                sharey=True,
                                sharex='col')
    if len(pileups) == 1:
        pileup = pileups[0]
        axarray[0].plot(range(begin - buf, begin + buf), pileup[:2*buf])
        axarray[0].axvline(begin, color='r')
        axarray[0].set_title(filenames[0] + " - 5\' of locus")
        axarray[0].set_xlabel("genomic position about 5' amplicon site")
        axarray[0].set_ylabel("depth")
        axarray[1].plot(range(end - buf, end + buf), pileup[-(2*buf):])
        axarray[1].axvline(end, color='r')
        axarray[1].set_title(filenames[0] + " - 3\' of locus")
        axarray[1].set_xlabel("genomic position about 3' amplicon site")
    else:
        for i, pileup in enumerate(pileups):
            axarray[i, 0].plot(range(begin - buf, begin + buf), pileup[:2*buf])
            axarray[i, 0].axvline(begin, color='r')
            axarray[i, 0].axvspan(begin - buf, begin, color = 'r', alpha = 0.25)
            axarray[i, 0].set_title(filenames[0] + " - 5\' of locus")
            axarray[i, 0].set_ylabel('depth')
            axarray[i, 1].plot(range(end - buf, end + buf), pileup[-(2*buf):])
            axarray[i, 1].axvline(end, color='r')
            axarray[i, 1].axvspan(end, end + buf, color = 'r', alpha = 0.25)
            axarray[i, 1].set_title(filenames[0] + " - 3\' of locus")
        axarray[i, 0].set_xlabel("genomic position about 5' amplicon site")
        axarray[i, 1].set_xlabel("genomic position about 3' amplicon site")
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


def get_stats(bams, loci):
    """Calculates stats for alignments."""
    stats = {}
    for bamfilename in bams:
        samplename = re.sub("\.bam", "", bamfilename)
        bam = pysam.AlignmentFile(bamfilename, 'rb')
        aln_rate = round(int(bam.mapped)/float(bam.mapped+bam.unmapped)*100, 2)
        loci_counts = [bam.count(chrom, int(begin), int(end))
                       for chrom, begin, end in loci]
        loci_rate = sum(loci_counts) / float(bam.mapped) * 100
        loci_covered = 0
        loci_depth = 0
        for chrom, begin, end in loci:
            begin = int(begin)
            end = int(end)
            for base in bam.count_coverage(chrom, begin, end):
                loci_covered += len([b for b in base if b != 0])
                loci_depth += sum(base)
        loci_area = sum([int(end) - int(begin) for chrom, begin, end in loci])
        loci_mean_depth = loci_depth / float(loci_area)
        loci_coverage = round(loci_covered / float(loci_area) * 100, 2)
        stats[bamfilename] = {"samplename": samplename,
                              "aln_rate": str(aln_rate) + "%",
                              "loci_rate": str(loci_rate) + "%",
                              "loci_coverage": str(loci_coverage) + "%",
                              "loci_mean_depth": loci_mean_depth}
    return stats


def main():
    """"""
    # arg parsing
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--loci", required=True,
                        nargs=3, metavar=("chr", "begin", "end"),
                        action="append", 
                        help="comma-delimited loci [required]")
    parser.add_argument("--buf", metavar="BUFFER", type=int, default=1000,
                        help="buffer about the loci to include in plots")
    parser.add_argument("-b", "--bam", metavar="BAM", nargs='+',
                        required=True, help="BAM file(s) [required]")
    args = parser.parse_args()
    # central dispatch
    pileups = parse_bam_for_pileup(args.bam, args.loci, args.buf)
    plots_map = {}
    for locus in args.loci:
        whole = plot_depth_across_locus(pileups[tuple(locus)],
                                        locus,
                                        args.bam,
                                        args.buf)
        sites = plot_depth_about_locus_ends(pileups[tuple(locus)],
                                            locus,
                                            args.bam,
                                            args.buf)
        plots_map['-'.join(locus)] = {"locus": '-'.join(locus),
                                      "whole": whole,
                                      "sites": sites}
    stats_map = get_stats(args.bam, args.loci)
    create_report(plots_map, stats_map)


if __name__ == "__main__":
    main()
