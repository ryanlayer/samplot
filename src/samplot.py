#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
import matplotlib.gridspec as gridspec
matplotlib.use('Agg')
import pylab
import random
import pysam
import os
import argparse

parser=argparse.ArgumentParser()

parser.add_argument("-b",
                  dest="bams",
                  help="Bam file names",
                  nargs="+",
                  required=True)

parser.add_argument("-o",
                  dest="output_file",
                  help="Output file name",
                  #type=lambda e:file_choices(("csv","tab"),e)
                  required=True)

parser.add_argument("-s",
                  dest="start",
                  help="Start range")

parser.add_argument("-e",
                  dest="end",
                  help="End range")

parser.add_argument("-c",
                  dest="chrom",
                  help="Chromosome range")

parser.add_argument("-w",
                  dest="window",
                  type=int,
                  default=1000,
                  help="Window, default(1000)")

parser.add_argument("--embed",
                  dest="embedded_path",
                  help="Embedded path",
                  required=False)
                  
args = parser.parse_args()

vals = []
all_plot_pairs = []
all_plot_reads = []

for bam_file_name in args.bams:
    bam_file = pysam.AlignmentFile(bam_file_name, "rb")
    pairs = {}

    plot_reads = []

    for read in bam_file.fetch(args.chrom,
                               int(args.start) - args.window,
                               int(args.end) + args.window):
        if (read.is_paired):
            if (read.reference_end) :
                if read.query_name not in pairs:
                    pairs[read.query_name] = []
                pairs[read.query_name].append(read)
        else:
            if (read.reference_end) :
                plot_reads.append([read.reference_start, read.reference_end])


    plot_pairs = []
    
    for pair in pairs:
        if len(pairs[pair]) == 2:

            plot_pairs.append([[pairs[pair][0].reference_start,
                                pairs[pair][0].reference_end,
                                not(pairs[pair][0].is_reverse)],
                               [pairs[pair][1].reference_start,
                                pairs[pair][1].reference_end,
                                not(pairs[pair][1].is_reverse)]])

            vals.append(pairs[pair][0].reference_start)
            if pairs[pair][0].reference_end:
                vals.append(pairs[pair][0].reference_end)

            vals.append(pairs[pair][1].reference_start)
            if pairs[pair][1].reference_end:
                vals.append(pairs[pair][1].reference_end)


    plot_pairs.sort(key=lambda x: x[1][1] - x[0][0])
    #plot_pairs.sort(key=lambda x: x[0][0], reverse=True)
    all_plot_pairs.append(plot_pairs)
    all_plot_reads.append(plot_reads)

all_plot_depths = []

for bam_file_name in args.bams:
    plot_depths = []
    bam_file = pysam.AlignmentFile(bam_file_name, "rb")
    for depth in bam_file.pileup(args.chrom,
                               int(args.start) - 1000,
                               int(args.end) + 1000):
        plot_depths.append([depth.reference_pos,depth.nsegments])

    all_plot_depths.append(plot_depths)

range_min = min(vals)
range_max = max(vals)

matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(8,10),dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)

num_ax = len(args.bams)+1
ax_i = 1

main_ax = fig.add_subplot(111)    # The big subplot

gs = gridspec.GridSpec(num_ax, 1, height_ratios = [1] + [10] * (num_ax-1))

for plot_pairs in all_plot_pairs:
    #ax = fig.add_subplot(num_ax,1,ax_i+1)
    ax =  matplotlib.pyplot.subplot(gs[ax_i])
    bam_file.close()

    c = 0
    min_insert_size=3000000
    max_insert_size=0
    for plot_pair in plot_pairs:
        p = [float(plot_pair[0][0] - range_min)/float(range_max - range_min), \
            float(plot_pair[1][1] - range_min)/float(range_max - range_min)]

        r1=[float(plot_pair[0][0] - range_min)/float(range_max - range_min), \
            float(plot_pair[0][1] - range_min)/float(range_max - range_min)]

        r2=[float(plot_pair[1][0] - range_min)/float(range_max - range_min), \
            float(plot_pair[1][1] - range_min)/float(range_max - range_min)]

        insert_size = plot_pair[1][1] - plot_pair[0][0]
        min_insert_size = min(min_insert_size, insert_size)
        max_insert_size = max(max_insert_size, insert_size)

        colors = {
                (True, False): 'black',
                (False, True): 'red',
                (False, False): 'blue',
                (True, True): 'green',
                }
        color = colors[(plot_pair[0][2], plot_pair[1][2])]

        ax.plot(p,[insert_size,insert_size],'-',color=color, alpha=0.25)
        ax.plot(r1, \
                [insert_size,insert_size], \
                '-', \
                color=color, \
                lw=3, 
                alpha=0.25)

        ax.plot(r2,[insert_size,insert_size],'-',color=color,lw=3, alpha=0.25)

        c+=1
    ax.yaxis.set_ticks(np.arange(min_insert_size, max_insert_size, (max_insert_size-min_insert_size)/4))

    ax.set_title(os.path.basename(args.bams[ax_i-1]))
    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    #ax.set_yticklabels([])
    if ax_i+1 < num_ax:
        ax.set_xticklabels([])
    else:
        labels = [int(range_min + l*(range_max-range_min)) \
                for l in ax.xaxis.get_majorticklocs()]
        ax.set_xticklabels(labels)
        ax.set_xlabel('Position on Chromosome ' + args.chrom)
    ax.set_ylabel('Insert size')
    ax_i += 1

#ax = fig.add_subplot(num_ax,1,1)
ax =  matplotlib.pyplot.subplot(gs[0])
r=[float(int(args.start) - range_min)/float(range_max - range_min), \
    float(int(args.end) - range_min)/float(range_max - range_min)]
ax.plot(r,[insert_size,insert_size],'-',color='black',lw=3)
ax.set_xlim([0,1])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
matplotlib.pyplot.tick_params(axis='x',length=0)
matplotlib.pyplot.tick_params(axis='y',length=0)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title(args.chrom + ':' + \
             str(args.start) + '-' + \
             str(args.end))

if args.embedded_path:
    if not os.path.isdir(args.embedded_path):
        os.mkdir(args.embedded_path)
    filename, ext = os.path.splitext(os.path.basename(args.output_file))
    img_filename = args.embedded_path + "/" + filename + ext
    arg_filename = args.embedded_path + "/" + filename + ".args"
    matplotlib.pyplot.savefig(img_filename,bbox_inches='tight')
    
    with open(arg_filename, 'w') as arg_file:
            keys = [str(x) for x in args.__dict__.keys()]
            keys.append("script")
            arg_file.write("#" + "\t".join(keys) + "\n")
            values = args.__dict__.values()
            for i in range(len(values)):
                if type(values[i]) == list:
                    values[i] = ','.join([os.path.basename(v) for v in values[i]])
                else:
                    values[i] = str(values[i])
            values.append("-")
            arg_file.write("\t".join(values) + "\n")
else:
    matplotlib.pyplot.savefig(args.output_file,bbox_inches='tight')
