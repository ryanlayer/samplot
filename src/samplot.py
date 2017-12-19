#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab
import random
import pysam
import os
import re
import statistics
import random
from optparse import OptionParser

marker_size = 4

def calc_query_pos_from_cigar(cigar, strand):
    cigar_ops = [[int(op[0]),op[1]] for op in re.findall('(\d+)([A-Za-z])', cigar)]

    order_ops = cigar_ops
    if not strand: # - strand
        order_ops = order_ops[::-1]

    qs_pos = 0
    qe_pos = 0
    q_len = 0

    for op_position  in range(len(cigar_ops)):
        op_len = cigar_ops[op_position][0]
        op_type = cigar_ops[op_position][1]

        if op_position == 0 and ( op_type == 'H' or op_type == 'S' ):
            qs_pos += op_len
            qe_pos += op_len
            q_len += op_len
        elif op_type == 'H' or op_type == 'S':
            q_len += op_len
        elif op_type == 'M' or op_type == 'I':
            qe_pos += op_len
            q_len += op_len

    return qs_pos, qe_pos

parser = OptionParser()

parser.add_option("-n",
                  dest="titles",
                  help="Plot title (CSV) ");

parser.add_option("-r",
                  dest="reference",
                  help="Reference file for CRAM");

parser.add_option("-b",
                  dest="bams",
                  help="Bam file names (CSV)")

parser.add_option("-o",
                  dest="output_file",
                  help="Output file name")

parser.add_option("-s",
                  dest="start",
                  help="Start range")

parser.add_option("-e",
                  dest="end",
                  help="End range")

parser.add_option("-c",
                  dest="chrom",
                  help="Chromosome range")

parser.add_option("-w",
                  dest="window",
                  type=int,
                  help="Window size (count of bases to include), default(0.5 * len)")

parser.add_option("-d",
                  dest="max_depth",
                  type=int,
                  help="Sampling depth(100 per 1kb)")

parser.add_option("-t",
                  dest="sv_type",
                  help="SV type")

parser.add_option("-T",
                  dest="transcript_file",
                  help="GFF of transcripts")

parser.add_option("-a",
                  dest="print_args",
                  action="store_true",
                  default=False,
                  help="Print commandline arguments")


(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

window = int((int(options.end) - int(options.start))/2)
if options.window:
    window = options.window

mapping_positions = []
all_pairs = []
all_plot_splits = []

# Get pairs and split
for bam_file_name in options.bams.split(','):
    bam_file = None
    if not options.reference:
        bam_file = pysam.AlignmentFile(bam_file_name, "rb")
    else:
        bam_file = pysam.AlignmentFile(bam_file_name, \
                                       "rc", \
                                       reference_filename=options.reference)
    pairs = {}

    plot_reads = []
    plot_splits = []

    for read in bam_file.fetch(options.chrom,
                               int(options.start) - window,
                               int(options.end) + window):

        if not read.is_secondary and read.has_tag('SA'):
            qs_pos, qe_pos = calc_query_pos_from_cigar(read.cigarstring, \
                                                       (not read.is_reverse))

            split=[[read.reference_end, not read.is_reverse, qs_pos]]

            for sa in read.get_tag('SA').split(';'):
                if len(sa) == 0:
                    continue


                pos = int(sa.split(',')[1])
                strand = sa.split(',')[2] == '+'
                qs_pos, qe_pos = \
                        calc_query_pos_from_cigar(read.cigarstring, strand)
                split.append([pos,strand,qs_pos])

            plot_splits.append(split)

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

            mapping_positions.append(pairs[pair][0].reference_start)
            if pairs[pair][0].reference_end:
                mapping_positions.append(pairs[pair][0].reference_end)

            mapping_positions.append(pairs[pair][1].reference_start)
            if pairs[pair][1].reference_end:
                mapping_positions.append(pairs[pair][1].reference_end)

    plot_pairs.sort(key=lambda x: x[1][1] - x[0][0])
    all_pairs.append(plot_pairs)
    all_plot_splits.append(plot_splits)

all_plot_depths = []
max_plot_depth = 0 

#Get coverage
for bam_file_name in options.bams.split(','):
    plot_depths = []
    bam_file = None

    if not options.reference:
        bam_file = pysam.AlignmentFile(bam_file_name, "rb")
    else:
        bam_file = pysam.AlignmentFile(bam_file_name, \
                                       "rc", \
                                       reference_filename=options.reference)

    for depth in bam_file.pileup(options.chrom,
                               int(options.start) - window - 500,
                               int(options.end) + window + 500):

        #if depth.reference_pos > int(options.start) - window and \
                #depth.reference_pos <  int(options.end) + window :

            plot_depths.append([depth.reference_pos,depth.nsegments])


    a = [x[0] for x in plot_depths]

    alld = set(range( int(options.start) - window - 500,int(options.end) + window + 500))
    rem = alld - set(a)
    plot_depths.extend([r, 0] for r in rem)

    max_plot_depth = max(max_plot_depth,max([x[1] for x in plot_depths]))
    all_plot_depths.append(plot_depths)


range_min = min(mapping_positions) if len(mapping_positions)>0 else int(options.start) - window - 500
range_max = max(mapping_positions) if len(mapping_positions)>0 else int(options.end) + window + 500

# Sample +/- pairs in the normal insert size range
if options.max_depth:
    max_depth=options.max_depth
    sampled_plot_pairs = [] 
    
    for plot_pairs in all_pairs:
        sampled_plot_pair = []
        plus_minus_pairs = []

        for pair in plot_pairs:
            if pair[0][2] == True and pair[1][2] == False:
                plus_minus_pairs.append(pair)
            else:
                sampled_plot_pair.append(pair)

        if len(plus_minus_pairs) > max_depth:
            lens = [pair[1][1] - pair[0][0] for pair in plus_minus_pairs]
            mean = statistics.mean(lens)
            stdev = statistics.stdev(lens)

            outside_norm = [pair for pair in plus_minus_pairs \
                        if pair[1][1] - pair[0][0] >= (mean + 4*stdev)]
            inside_norm = [pair for pair in plus_minus_pairs \
                        if pair[1][1] - pair[0][0] < (mean + 2*stdev)]

            sampled_plot_pair += outside_norm
            if len(inside_norm) > max_depth:
                sampled_plot_pair += random.sample(inside_norm, max_depth)
            else:
                sampled_plot_pair += inside_norm
        else:
            sampled_plot_pair+=plus_minus_pairs
        sampled_plot_pairs.append(sampled_plot_pair)
    all_pairs = sampled_plot_pairs

matplotlib.rcParams.update({'font.size': 12})

height = 1.1 + (1.33 * len(options.bams.split(',')))
fig = matplotlib.pyplot.figure(figsize=(8,height),dpi=300)
fig.subplots_adjust( left=0.075, right=.85, bottom=.15, top=0.9)

num_ax = len(options.bams.split(','))+1

if options.transcript_file:
    num_ax+=1

ax_i = 1

ratios = [1] + [10] * len(options.bams.split(','))

if options.transcript_file:
    ratios += [2]

gs = gridspec.GridSpec(num_ax, 1, height_ratios = ratios)

for plot_pairs in all_pairs:
    ax =  matplotlib.pyplot.subplot(gs[ax_i])
    bam_file.close()

    c = 0
    min_insert_size=3000000
    max_insert_size=0

    # Plot pairs
    for plot_pair in plot_pairs:
        p = [float(plot_pair[0][1] - range_min)/float(range_max - range_min), \
            float(plot_pair[1][0] - range_min)/float(range_max - range_min)]

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

        ax.plot(p,\
                [insert_size,insert_size],\
                '-',color=color, \
                alpha=0.25, \
                lw=0.5, \
                marker='s', \
                markersize=marker_size)

        c+=1

    # Plot splits
    for split in all_plot_splits[ax_i - 1]:
        start = split[0]
        end = split[1]
        if start[0] > end[0]:
            end = split[0]
            start = split[1]

        # Do not plot pairs that extend beyond the current range
        if range_min > start[0] or range_max < end[0]:
            continue
            
        p = [float(start[0] - range_min)/float(range_max - range_min), \
             float(end[0] - range_min)/float(range_max - range_min)]

        color = 'black'

        if start[1] != end[1]: #INV
            if start[1] == '+': 
                color = 'green'
            else:
                color = 'blue'
        elif start[1] == '+' and end[1] == '+' and start[2] < end[2]: #DEL
            color = 'black'
        elif start[1] == '-' and end[1] == '-' and start[2] > end[2]: #DEL
            color = 'black'
        elif start[1] == '+' and end[1] == '+' and start[2] > start[2]: #DUP
            color = 'red'
        elif start[1] == '-' and end[1] == '-' and start[2] < end[2]: #DUP
            color = 'red'
        
        insert_size = end[0] - start[0]

        ax.plot(p,\
                [insert_size,insert_size],\
                ':',\
                color=color, \
                alpha=0.25, \
                lw=1, \
                marker='o',
                markersize=marker_size)


    if max_insert_size-min_insert_size != 0: 
        ax.yaxis.set_ticks(np.arange(min_insert_size, \
                     max_insert_size, \
                     max(2,(max_insert_size-min_insert_size)/4)))

    if options.titles and \
            len(options.titles.split(',')) == len(options.bams.split(',')):
        ax.set_title(options.titles.split(',')[ax_i-1], \
                     fontsize=8, loc='left')
    else:
        ax.set_title(os.path.basename(options.bams.split(',')[ax_i-1]), \
                     fontsize=8, loc='left')


    # Plot depth
    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', labelsize=6)
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    if ax_i+1 < num_ax:
        ax.set_xticklabels([])
    else:
        labels = [int(range_min + l*(range_max-range_min)) \
                for l in ax.xaxis.get_majorticklocs()]
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_xlabel('Chromosomal position on ' + options.chrom, fontsize=8)
    ax.set_ylabel('Insert size', fontsize=8)

    cover_x = []
    cover_y = []

    if max([d[1] for d in all_plot_depths[ax_i-1]]) == 0: continue

    for d in all_plot_depths[ax_i - 1]:
        cover_x.append(float(d[0] - range_min)/float(range_max - range_min))
        cover_y.append(d[1])

    ax2 = ax.twinx()
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,max_plot_depth])
    ax2.fill_between(cover_x, \
                     cover_y, \
                     [0] * len(cover_y),
                     color='grey',
                     alpha=0.25)
 
    ax2.set_ylabel('Coverage', fontsize=8)
    ax2.tick_params(axis='y', colors='grey', labelsize=6)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    ax_i += 1

# Plot the variant
ax =  matplotlib.pyplot.subplot(gs[0])
r=[float(int(options.start) - range_min)/float(range_max - range_min), \
    float(int(options.end) - range_min)/float(range_max - range_min)]
ax.plot(r,[0,0],'-',color='black',lw=3)
ax.set_xlim([0,1])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
matplotlib.pyplot.tick_params(axis='x',length=0)
matplotlib.pyplot.tick_params(axis='y',length=0)
ax.set_xticklabels([])
ax.set_yticklabels([])


sv_size = float(options.end) - float(options.start)
sv_size_unit = 'bp'

if sv_size > 1000000:
    sv_size = "{0:0.2f}".format(sv_size/1000000.0)
    sv_size_unit = 'mb'
if sv_size > 1000:
    sv_size = "{0:0.2f}".format(sv_size/1000.0)
    sv_size_unit = 'kb'

sv_title = str(sv_size) + ' ' + sv_size_unit

if options.sv_type:
    sv_title += ' ' + options.sv_type

ax.set_title(sv_title, fontsize=8)


marker_colors = ['black', 'red', 'blue', 'green']
legend_elements = []

for color in marker_colors:
    legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1], \
            color=color,
            linestyle='-',
            lw=1)]

legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1], \
            markerfacecolor="None",
            markeredgecolor='grey',
            color='grey',
            marker='o', \
            markersize=marker_size,
            linestyle=':',
            lw=1)]

legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1], \
            markerfacecolor="None",
            markeredgecolor='grey',
            color='grey',
            marker='s', \
            markersize=marker_size,
            linestyle='-',
            lw=1)]

fig.legend( legend_elements ,
            ["Deletion/Normal",\
             "Duplication", \
             "Inversion", \
             "Inversion", \
             "Split-read", \
             "Pair-end read"], \
            fontsize = 6,
            bbox_to_anchor=(1, 1),
            bbox_transform=matplotlib.pyplot.gcf().transFigure,
            frameon=False)

# Plot sorted/bgziped/tabixed transcript file
if options.transcript_file:
    tbx = pysam.TabixFile(options.transcript_file)

    genes = {}
    transcripts = {}
    cdss = {}

    for row in tbx.fetch(options.chrom, \
                         int(options.start), \
                         int(options.end)):
        A = row.split()
        if A[1] != 'ensembl_havana': continue

        if A[2] == 'gene':
            info =  dict([list(val.split('=')) for val in A[8].split(';')])

            genes[info['Name']] = [A[0],int(A[3]), int(A[4]),info]

        elif A[2] == 'transcript':
            info =  dict([list(val.split('=')) for val in A[8].split(';')])
            
            if info['Parent'] not in transcripts:
                transcripts[info['Parent']] = {}

            transcripts[info['Parent']][info['ID']] = \
                    [A[0],int(A[3]), int(A[4]),info]

        elif A[2] == 'CDS':

            info =  dict([list(val.split('=')) for val in A[8].split(';')])
            
            if info['Parent'] not in cdss:
                cdss[info['Parent']] = {}

            if info['ID'] not in cdss[info['Parent']]:
                cdss[info['Parent']][info['ID']] = []

            cdss[info['Parent']][info['ID']].append([A[0], \
                                                     int(A[3]), \
                                                     int(A[4]), \
                                                     info])

    ax =  matplotlib.pyplot.subplot(gs[num_ax-1])

    t_i = 0
    for gene in genes:
        gene_id = genes[gene][3]['ID']
        if gene_id not in transcripts: continue
        for transcript in transcripts[gene_id]:
            t_start = max(range_min, transcripts[gene_id][transcript][1])
            t_end = min(range_max, transcripts[gene_id][transcript][2])
            r=[float(t_start - range_min)/float(range_max - range_min), \
               float(t_end - range_min)/float(range_max - range_min)]
            ax.plot(r,[t_i,t_i],'-',color='cornflowerblue',lw=1)

            ax.text(r[0],t_i + 0.02,gene,fontsize=6,color='cornflowerblue')

        
            if transcript in cdss:
                for cds in cdss[transcript]:
                    for exon in cdss[transcript][cds]:
                        e_start = max(range_min,exon[1])
                        e_end = min(range_max,exon[2])
                        r=[float(e_start - range_min)/float(range_max - range_min), \
                            float(e_end - range_min)/float(range_max - range_min)]
                        ax.plot(r,[t_i,t_i],'-',color='cornflowerblue',lw=4)

                t_i += 1

    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    labels = [int(range_min + l*(range_max-range_min)) \
            for l in ax.xaxis.get_majorticklocs()]
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_xlabel('Chromosomal position on ' + options.chrom, fontsize=8)

if len(options.bams.split(',')) > 3:
    matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(options.output_file)


if options.print_args:
    print ('#' + '\t'.join([ 'titles',
                            'reference',
                            'bams',
                            'output_file',
                            'start',
                            'end',
                            'chrom',
                            'window',
                            'max_depth',
                            'sv_type',
                            'transcript_file']))
    print ('\t'.join([ options.titles if options.titles else 'None',
                      options.reference if options.reference else 'None',
                      options.bams,
                      options.output_file,
                      options.start,
                      options.end,
                      options.chrom,
                      str(options.window),
                      str(options.max_depth) if options.max_depth else 'None',
                      options.sv_type,
                      options.transcript_file if options.transcript_file else 'None']))
