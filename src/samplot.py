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
    cigar_ops = [[int(op[0]), op[1]] for op in re.findall('(\d+)([A-Za-z])', \
                  cigar)]

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

def get_depth(chrom, start, end, window, bam_files, reference):
    all_plot_depths = []
    max_plot_depth = 0 

    #Get coverage
    for bam_file_name in bam_files:
        plot_depths = []
        bam_file = None

        if not reference:
            bam_file = pysam.AlignmentFile(bam_file_name, "rb")
        else:
            bam_file = pysam.AlignmentFile(bam_file_name, \
                                           "rc", \
                                           reference_filename=reference)

        last = None
        for depth in bam_file.pileup(chrom,
                                     int(start) - window - 500,
                                     int(end) + window + 500):
            if last != None:
                if depth.reference_pos > last +1 :
                    for i in range(last + 1, depth.reference_pos - 1):
                        plot_depths.append([i,0])

            plot_depths.append([depth.reference_pos,depth.nsegments])
            last = depth.reference_pos

        a = [x[0] for x in plot_depths]

        alld = set(range( int(start) - window - 500,
                          int(end) + window + 500))
        rem = alld - set(a)
        plot_depths.extend([r, 0] for r in rem)

        max_plot_depth = max(max_plot_depth,max([x[1] for x in plot_depths]))
        all_plot_depths.append(plot_depths)

    return max_plot_depth, all_plot_depths

def get_pairs_and_splits(chrom, start, end, window, bam_files, reference):
    mapping_positions = []
    all_pairs = []
    all_plot_splits = []

    # Get pairs and split
    for bam_file_name in bam_files:
        bam_file = None
        if not reference:
            bam_file = pysam.AlignmentFile(bam_file_name, "rb")
        else:
            bam_file = pysam.AlignmentFile(bam_file_name, \
                                           "rc", \
                                           reference_filename=reference)
        pairs = {}

        plot_reads = []
        plot_splits = []

        for read in bam_file.fetch(chrom,
                                   max(0,int(start) - window),
                                   int(end) + window):

            if not read.is_secondary and read.has_tag('SA'):
                qs_pos, qe_pos = \
                    calc_query_pos_from_cigar(read.cigarstring, \
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
                    plot_reads.append([read.reference_start,
                                       read.reference_end])

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

    return [mapping_positions, all_pairs, all_plot_splits]

def sample_normal(max_depth, all_pairs):
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
                        if pair[1][1] - pair[0][0] >= (mean + options.z*stdev)]
            inside_norm = [pair for pair in plus_minus_pairs \
                        if pair[1][1] - pair[0][0] < (mean + options.z*stdev)]

            sampled_plot_pair += outside_norm
            if len(inside_norm) > max_depth:
                sampled_plot_pair += random.sample(inside_norm, max_depth)
            else:
                sampled_plot_pair += inside_norm
        else:
            sampled_plot_pair+=plus_minus_pairs

        sampled_plot_pairs.append(sampled_plot_pair)

    return sampled_plot_pairs


parser = OptionParser()

parser.add_option("-n",
                  dest="titles",
                  help="Plot title (CSV) ");

parser.add_option("-r",
                  dest="reference",
                  help="Reference file for CRAM");

parser.add_option("-z",
                  dest="z",
                  type=int,
                  default=4,
                  help="Number of stdevs from the mean (default 4)");

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

parser.add_option("-A",
                  dest="annotation_file",
                  help="bed.gz tabixed file of transcripts")

parser.add_option("-a",
                  dest="print_args",
                  action="store_true",
                  default=False,
                  help="Print commandline arguments")

parser.add_option("-H",
                  dest="plot_height",
                  type=int,
                  help="Plot height")

parser.add_option("-W",
                  dest="plot_width",
                  type=int,
                  help="Plot width")


(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

if not options.bams:
    parser.error('BAMSs not given')
    
if not options.start:
    parser.error('SV start not given')
    
if not options.end:
    parser.error('SV end not given')

if not options.chrom:
    parser.error('SV chrom not given')

if not options.sv_type:
    parser.error('SV sv_type not given')

plot_height = 5
plot_width = 8

if options.plot_height:
    plot_height = options.plot_height
else:
    plot_height = 2 + len(options.bams.split(','))

if options.plot_width:
    plot_width = options.plot_width

window = int((int(options.end) - int(options.start))/2)
if options.window:
    window = options.window

#get depths from each bam file
max_plot_depth, all_plot_depths =  get_depth(options.chrom,
                                             options.start,
                                             options.end,
                                             window,
                                             options.bams.split(','),
                                             options.reference)

#get pairs and split from each bam file
mapping_positions, all_pairs, all_plot_splits = \
        get_pairs_and_splits(options.chrom,
                             options.start,
                             options.end,
                             window,
                             options.bams.split(','),
                             options.reference)

#get the x-axis range
range_min = min(mapping_positions) \
        if len(mapping_positions) > 0 \
        else int(options.start) - window - 500
range_max = max(mapping_positions) \
        if len(mapping_positions) > 0 \
        else int(options.end) + window + 500

# Sample +/- pairs in the normal insert size range
if options.max_depth:
    all_pairs = sample_normal(options.max_depth, all_pairs)

matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(plot_width, plot_height), dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)

# give one axis to display the SV then one for each bam
num_ax = len(options.bams.split(','))+1

# add another if a annotation file is given
if options.transcript_file:
    num_ax+=1

if options.annotation_file:
    num_ax+=len(options.annotation_file.split(','))

# set the relative sizes for each
ratios = [1] + [10] * len(options.bams.split(','))
if options.annotation_file:
    ratios += [1]*len(options.annotation_file.split(','))
if options.transcript_file:
    ratios += [2]


# set the color of the pair based on the two strands
colors = { (True, False): 'black', # DEL
           (False, True): 'red',   # DUP
           (False, False): 'blue', # INV
           (True, True): 'green' } # INV

gs = gridspec.GridSpec(num_ax, 1, height_ratios = ratios)

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

# make SV title 
sv_size = float(options.end) - float(options.start)
sv_size_unit = 'bp'

if sv_size > 1000000:
    sv_size = "{0:0.2f}".format(sv_size/1000000.0)
    sv_size_unit = 'mb'
elif sv_size > 1000:
    sv_size = "{0:0.2f}".format(sv_size/1000.0)
    sv_size_unit = 'kb'

sv_title = str(sv_size) + ' ' + sv_size_unit + ' ' + options.sv_type
ax.set_title(sv_title, fontsize=8)

sample_axs = []
ax_i = 1
for plot_pairs in all_pairs:
    ax =  matplotlib.pyplot.subplot(gs[ax_i])
    sample_axs.append(ax)

    min_insert_size=3000000
    max_insert_size=0

    # Plot pairs
    for plot_pair in plot_pairs:
        p = [float(plot_pair[0][0] - range_min)/float(range_max - range_min), \
            float(plot_pair[1][1] - range_min)/float(range_max - range_min)]

        # y value is the insert size
        insert_size = plot_pair[1][1] - plot_pair[0][0]

        # use this to scale the y-axis
        min_insert_size = min(min_insert_size, insert_size)
        max_insert_size = max(max_insert_size, insert_size)

        color = colors[(plot_pair[0][2], plot_pair[1][2])]

        # plot the individual pair
        ax.plot(p,\
                [insert_size,insert_size],\
                '-',color=color, \
                alpha=0.25, \
                lw=0.5, \
                marker='s', \
                markersize=marker_size)

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

        # For a given SV, the orientation of the pairs and split do not match
        # so we cannot use the colors dict here
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
        
        # y value is the insert size
        insert_size = end[0] - start[0]

        # plot the individual split
        ax.plot(p,\
                [insert_size,insert_size],\
                ':',\
                color=color, \
                alpha=0.25, \
                lw=1, \
                marker='o',
                markersize=marker_size)


    # try to only have y-lables, which does work if there are too few 
    # observed insert sizes
    if max_insert_size-min_insert_size != 0: 
        ax.yaxis.set_ticks(np.arange(min_insert_size, \
                                     max_insert_size, \
                                     max(2,
                                         (max_insert_size-min_insert_size)/4)))

    # set the axis title to be either one passed in or filename
    if options.titles and \
            len(options.titles.split(',')) == len(options.bams.split(',')):
        ax.set_title(options.titles.split(',')[ax_i-1], \
                     fontsize=8, loc='left')
    else:
        ax.set_title(os.path.basename(options.bams.split(',')[ax_i-1]), \
                     fontsize=8, loc='left')
    
    # set axis parameters
    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', labelsize=6)
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    # only set x-axis labels on the last sample
    if ax_i+1 < num_ax:
        ax.set_xticklabels([])
    else:
        labels = [int(range_min + l*(range_max-range_min)) \
                for l in ax.xaxis.get_majorticklocs()]
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_xlabel('Chromosomal position on ' + options.chrom, fontsize=8)
    ax.set_ylabel('Insert size', fontsize=8)

    # Plot depth
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
 
    # set axis parameters
    ax2.set_ylabel('Coverage', fontsize=8)
    ax2.tick_params(axis='y', colors='grey', labelsize=6)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    ax_i += 1


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

ax.legend( legend_elements ,
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

if options.annotation_file:
    a_i = 0
    for annotation_file in options.annotation_file.split(','):
        tbx = pysam.TabixFile(annotation_file)

        # fetch and parse data from the tabixed gff file
        itr = None
        try:
            itr = tbx.fetch(options.chrom, \
                            int(options.start), \
                            int(options.end))
        except ValueError:
            # try and account for chr/no chr prefix
            chrom = options.chrom
            if chrom[:3] == 'chr':
                chrom = chrom[3:]
            else:
                chrom = 'chr' + chrom

            try:
                itr = tbx.fetch(chrom,
                                int(options.start), \
                                int(options.end))
            except ValueError:
                sys.exit('Error: Could not fetch ' + \
                        options.chrom + ':' + options.start + '-' + \
                        options.end + \
                        'from ' + options.annotation_file)
         
        ax =  matplotlib.pyplot.subplot(gs[
                1 + len(options.bams.split(',')) + a_i])

        for row in itr:
            A = row.rstrip().split()
            t_start = max(range_min, int(A[1]))
            t_end = min(range_max, int(A[2]))

            r=[float(t_start - range_min)/float(range_max - range_min), \
               float(t_end - range_min)/float(range_max - range_min)]

            if len(A) > 3 :
                try:
                    v = float(A[3])
                    ax.plot(r,[0,0],'-',color=str(1-v),lw=5)
                except ValueError:
                    ax.plot(r,[0,0],'-',color='black',lw=5)
                    ax.text(r[0],0 + 0.1,A[3],fontsize=6,color='black')
            else:
                ax.plot(r,[0,0],'-',color='black',lw=5)

        # set axis parameters
        ax.set_xlim([0,1])
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title(os.path.basename(annotation_file), \
                         fontsize=8, loc='left')

        matplotlib.pyplot.tick_params(axis='x',length=0)
        matplotlib.pyplot.tick_params(axis='y',length=0)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        a_i+=1


# Plot sorted/bgziped/tabixed transcript file
if options.transcript_file:
    tbx = pysam.TabixFile(options.transcript_file)

    genes = {}
    transcripts = {}
    cdss = {}

    # fetch and parse data from the tabixed gff file
    itr = None
    try:
        itr = tbx.fetch(options.chrom, \
                        int(options.start), \
                        int(options.end))
    except ValueError:
        # try and account for chr/no chr prefix
        chrom = options.chrom
        if chrom[:3] == 'chr':
            chrom = chrom[3:]
        else:
            chrom = 'chr' + chrom

        try:
            itr = tbx.fetch(chrom,
                            int(options.start), \
                            int(options.end))
        except ValueError:
            sys.exit('Error: Could not fetch ' + \
                    options.chrom + ':' + options.start + '-' + options.end + \
                    'from ' + options.transcript_file)
     
    for row in itr:
        A = row.split()

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
    ax =  matplotlib.pyplot.subplot(gs[-1])

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

    # set axis parameters
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
    ax.set_title(os.path.basename(options.transcript_file), \
                         fontsize=8, loc='left')

# save
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(options.output_file)

# give sv-plaudit output
if options.print_args:
    import json
    args_filename = os.path.splitext(options.output_file)[0] + ".json"
    args_info = {
        'titles': options.titles if options.titles else 'None',
        'reference': options.reference if options.reference else 'None',
        'bams': options.bams.split(','),
        'output_file': options.output_file,
        'start': options.start, 
        'end': options.end,
        'chrom': options.chrom,
        'window': options.window,
        'max_depth': options.max_depth if options.max_depth else 'None',
        'sv_type': options.sv_type,
        'transcript_file': options.transcript_file if options.transcript_file else 'None'
    }
    with open(args_filename, 'w') as outfile:
        json.dump(args_info, outfile)
