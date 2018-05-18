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
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from optparse import OptionParser

marker_size = 4

#{{{def calc_query_pos_from_cigar(cigar, strand):
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
#}}}

#{{{ def sample_normal(max_depth, all_pairs):
def sample_normal(max_depth, all_pairs, z):
    all_sampled_pairs = [] 

    for pairs in all_pairs:
        sampled_pairs = {}
        plus_minus_pairs = {}

        for read_name in pairs:
            pair = pairs[read_name]
            if len(pair) != 2: 
                continue
            if pair[0][2] == True and pair[1][2] == False:
                plus_minus_pairs[read_name] = pair
            else:
                sampled_pairs[read_name] = pair

        if len(plus_minus_pairs) > max_depth:
            lens = [pair[1][1] - pair[0][0] \
                    for pair in plus_minus_pairs.values()]
            mean = statistics.mean(lens)
            stdev = statistics.stdev(lens)

            inside_norm = {}

            for read_name in pairs:
                pair = pairs[read_name]
                if len(pair) != 2: 
                    continue
                if pair[1][1] - pair[0][0] >= mean + options.z*stdev:
                    sampled_pairs[read_name] = pair 
                else:
                    inside_norm[read_name] = pair 

            if len(inside_norm) > max_depth:
                for read_name in random.sample(inside_norm.keys(), max_depth):
                    sampled_pairs[read_name] = inside_norm[read_name]
            else:
                for read_name in inside_norm:
                    sampled_pairs[read_name] = inside_norm[read_name]
        else:
            for read_name in plus_minus_pairs:
                sampled_pairs[read_name] = plus_minus_pairs[read_name]

        all_sampled_pairs.append(sampled_pairs)

    return all_sampled_pairs
#}}}

#{{{ def add_coverage(read, coverage):
def add_coverage(read, coverage):
    curr_pos = read.reference_start
    if not read.cigartuples: return

    for op,length in read.cigartuples:
        if op == 0:
            for pos in range(curr_pos, curr_pos+length+1):
                if pos not in coverage:
                    coverage[pos] = 0
                coverage[pos] += 1
            curr_pos += length
        elif op == 1:
            curr_pos = curr_pos
        elif op == 2:
            curr_pos += length
        elif op == 3:
            curr_pos = length
        elif op == 4:
            curr_pos = curr_pos
        elif op == 5:
            curr_pos = curr_pos
        else:
            curr_pos += length
#}}}

#{{{def add_pair(read, pairs):
def add_pair(read, pairs):
    if read.is_unmapped: return
    if not (read.is_paired): return
    if read.is_secondary: return
    if read.is_supplementary: return

    if read.query_name not in pairs:
        pairs[read.query_name] = []

    pairs[read.query_name].append( [read.reference_start, \
                                    read.reference_end, \
                                    not(read.is_reverse)] )
#}}}

#{{{def add_split(read, splits):
def add_split(read, splits, bam_file):
    if not read.is_secondary and not read.is_supplementary and read.has_tag('SA'):
        qs_pos, qe_pos = \
            calc_query_pos_from_cigar(read.cigarstring, \
                                      (not read.is_reverse))

        splits[read.query_name]=[[read.reference_start,
                                  read.reference_end, \
                                  not read.is_reverse, \
                                  qs_pos]]

        for sa in read.get_tag('SA').split(';'):
            if len(sa) == 0:
                continue
            chrm = sa.split(',')[0]
            if chrm != bam_file.get_reference_name(read.reference_id):
                continue

            pos = int(sa.split(',')[1])
            strand = sa.split(',')[2] == '+'
            cigar = sa.split(',')[3]
            qs_pos, qe_pos = \
                    calc_query_pos_from_cigar(cigar, strand)
            splits[read.query_name].append([pos,pos+qe_pos,strand,qs_pos])

        if len(splits[read.query_name]) == 1:
            del splits[read.query_name]
        else:
            splits[read.query_name].sort(key=lambda x:x[2])
#}}}

#{{{class Alignment:
class Alignment:
    def __init__(self,start,end,strand,query_position):
        self.start = start
        self.end = end
        self.strand = strand
        self.query_position = query_position
    def __str__(self):
        return ','.join([str(x) for x in [self.start, 
                                          self.end, 
                                          self.strand, 
                                          self.query_position]])
#}}}

#{{{class LongRead:
class LongRead:
    def __init__(self,start,end,alignments):
        self.start = start
        self.end = end
        self.alignments = alignments
    def __str__(self):
        return ','.join([str(x) for x in [self.start, 
                                          self.end, 
                                          len(self.alignments)]])
#}}}
    
#{{{def get_alignments_from_cigar(curr_pos, cigartuples, reverse=False):
def get_alignments_from_cigar(curr_pos, strand, cigartuples, reverse=False):
    alignments = []
    q_pos = 0
    if reverse:
        cigartuples = cigartuples[::-1]
    for op,length in cigartuples:
        if op == 0: #M
            alignments.append(Alignment(curr_pos,
                                        curr_pos+length,
                                        strand,
                                        q_pos))
            curr_pos += length
            q_pos += length
        elif op == 1: #I
            q_pos += length
        elif op == 2: #D
            curr_pos += length
        elif op == 3: #N
            curr_pos += length
        elif op == 4: #S
            q_pos += length
    return alignments
#}}}

#{{{def get_cigartuples_from_string(cigarstring):
def get_cigartuples_from_string(cigarstring):
    cigartuples = []

    #pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    #M   BAM_CMATCH  0
    #I   BAM_CINS    1
    #D   BAM_CDEL    2
    #N   BAM_CREF_SKIP   3
    #S   BAM_CSOFT_CLIP  4
    #H   BAM_CHARD_CLIP  5
    #P   BAM_CPAD    6
    #=   BAM_CEQUAL  7
    #X   BAM_CDIFF   8
    #B   BAM_CBACK   9

    cigar_map = { 'M' : 0,
                  'I' : 1,
                  'D' : 2,
                  'N' : 3,
                  'S' : 4,
                  'H' : 5,
                  'P' : 6,
                  '=' : 7,
                  'X' : 8,
                  'B' : 9}

    for match in re.findall(r'(\d+)([A-Z]{1})', cigarstring):
        l = int(match[0])
        o = match[1]
        cigartuples.append((cigar_map[o], l))
    
    return cigartuples
#}}}

#{{{def merge_alignments(min_gap, alignments):
def merge_alignments(min_gap, alignments):

    merged_alignments = []

    for alignment in alignments:
        if len(merged_alignments) == 0:
            merged_alignments.append(alignment)
        else:
            if alignment.start < merged_alignments[-1].end + min_gap:
                merged_alignments[-1].end = alignment.start
            else:
                merged_alignments.append(alignment)
    return merged_alignments
#}}}

#{{{def add_long_reads(read, long_reads, range_min, range_max):
def add_long_reads(read, long_reads, range_min, range_max):

    if read.is_supplementary or read.is_secondary: return

    alignments = get_alignments_from_cigar(read.pos,
                                           not read.is_reverse,
                                           read.cigartuples)
    min_gap = 100
    merged_alignments = merge_alignments(min_gap, alignments)

    read_strand = not read.is_reverse
    read_chrom = read.reference_name

    if read.has_tag('SA'):
        for sa in read.get_tag('SA').split(';'):
            if len(sa) == 0: continue
            rname,pos,strand,cigar,mapq,nm = sa.split(',')
            if rname != read_chrom: continue
            sa_pos = int(pos)
            sa_strand = strand == '+'
            strand_match = read_strand != sa_strand
            sa_cigartuples = get_cigartuples_from_string(cigar)
            sa_alignments = get_alignments_from_cigar(sa_pos, \
                                                      sa_strand, \
                                                      sa_cigartuples, \
                                                      reverse=strand_match)
            sa_merged_alignments =  merge_alignments(min_gap, sa_alignments)
            if (len(sa_merged_alignments) > 0):
                merged_alignments += sa_merged_alignments


    if read.query_name not in long_reads:
        long_reads[read.query_name] = []

    long_reads[read.query_name].append(LongRead(read.reference_start,
                                                read.reference_end,
                                                merged_alignments))
#}}}

#{{{def get_long_read_plan(read_name, long_reads, splits):
def get_long_read_plan(read_name, long_reads, range_min, range_max):
    alignments = []

    for long_read in long_reads[read_name]:
        for alignment in long_read.alignments:
            alignments.append(alignment)

    alignments.sort(key=lambda x: x.query_position)

    gaps = []

    curr = alignments[0]

    primary_strand = alignments[0].strand

    steps = []
    steps.append( [ [curr.start,curr.end], 'ALIGN' ] )

    for i in range(1,len(alignments)):
        last = alignments[i-1]
        curr = alignments[i]

        # figure out what the event is

        # INV
        if (curr.strand != last.strand):
            gap = abs(curr.end - last.end)

            if (curr.strand != primary_strand):
                steps.append( [ [curr.end,last.end], 'INVIN' ] )
            else:
                steps.append( [ [curr.start,last.start], 'INVOUT' ] )

            steps.append( [ [curr.start,curr.end], 'ALIGN' ] )
        # DUP
        elif (curr.start  < last.end):
            gap = abs(last.end - curr.start)

            steps.append( [ [last.end,curr.start], 'DUP' ] )
            steps.append( [ [curr.start,curr.end], 'ALIGN' ] )

        # DEL
        else:
            gap = abs(curr.start - last.end)

            steps.append( [ [last.end,curr.start], 'DEL' ] )
            steps.append( [ [curr.start,curr.end], 'ALIGN' ] )

        if not (min(curr.start,last.end) < range_min or \
                max(curr.start,last.end) > range_max ):
            gaps.append(gap)


    max_gap = 0
    if len(gaps) > 0:
        max_gap = max(gaps)

    plan = [max_gap, steps]

    return plan 
#}}}

#{{{def get_long_read_max_gap(read_name, long_reads):
def get_long_read_max_gap(read_name, long_reads):
    alignments = []
    for long_read in long_reads[read_name]:
        for alignment in long_read.alignments:
            alignments.append(alignment)
    alignments.sort(key=lambda x: x.query_position)

    #find biggest gap
    max_gap = 0
    for i in range(1,len(alignments)):
        curr_gap = abs(alignments[i].start - alignments[i-1].end)
        max_gap = max(max_gap,curr_gap)
    return max_gap
#}}}

#{{{def plot_variant(start, end, sv_type, ax, range_min, range_max):
def plot_variant(start, end, sv_type, ax, range_min, range_max):
    r=[float(int(start) - range_min)/float(range_max - range_min), \
        float(int(end) - range_min)/float(range_max - range_min)]
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

    ## make SV title 
    sv_size = float(end) - float(start)
    sv_size_unit = 'bp'

    if sv_size > 1000000:
        sv_size = "{0:0.2f}".format(sv_size/1000000.0)
        sv_size_unit = 'mb'
    elif sv_size > 1000:
        sv_size = "{0:0.2f}".format(sv_size/1000.0)
        sv_size_unit = 'kb'

    sv_title = str(sv_size) + ' ' + sv_size_unit + ' ' + sv_type
    ax.set_title(sv_title, fontsize=8)
#}}}

#{{{ def plot_pairs(all_pairs
def plot_pairs(pairs,
               ax,
               range_min,
               range_max,
               curr_min_insert_size,
               curr_max_insert_size):

    ## set the color of the pair based on the two strands
    colors = { (True, False): 'black', # DEL
               (False, True): 'red',   # DUP
               (False, False): 'blue', # INV
               (True, True): 'green' } # INV

    for read_name in pairs:
        pair = pairs[read_name]

        if len(pair) != 2: continue
        p = [float(pair[0][0] - range_min)/float(range_max - range_min), \
             float(pair[1][1] - range_min)/float(range_max - range_min)]

        # y value is the insert size
        insert_size = pair[1][1] - pair[0][0]
        # use this to scale the y-axis
        if not curr_min_insert_size or curr_min_insert_size < insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size > insert_size:
            curr_max_insert_size = insert_size

        color = colors[(pair[0][2], pair[1][2])]

        # plot the individual pair
        ax.plot(p,\
                [insert_size,insert_size],\
                '-',color=color, \
                alpha=0.25, \
                lw=0.5, \
                marker='s', \
                markersize=marker_size)
    return [curr_min_insert_size, curr_max_insert_size]
#}}}

#{{{ def plot_splits(all_splits,
def plot_splits(splits,
                ax,
                range_min,
                range_max,
                curr_min_insert_size,
                curr_max_insert_size):

    for read_name in splits:
        split = splits[read_name]
        start = split[0]
        end = split[1]
        if start[0] > end[0]:
            end = split[0]
            start = split[1]

        # Do not plot pairs that extend beyond the current range
        if range_min > start[1] or range_max < end[0]:
            continue
            
        p = [float(start[1] - range_min)/float(range_max - range_min), \
             float(end[0] - range_min)/float(range_max - range_min)]

        # For a given SV, the orientation of the pairs and split do not match
        # so we cannot use the colors dict here
        color = 'black'
        if start[2] != end[2]: #INV
            if start[2] == '+': 
                color = 'green'
            else:
                color = 'blue'
        elif start[2] == True and end[2] == True and start[3] < end[3]: #DEL
            color = 'black'
        elif start[2] == False and end[2] == False and start[3] > end[3]: #DEL
            color = 'black'
        elif start[2] == True and end[2] == True and start[3] > end[3]: #DUP
            color = 'red'
        elif start[2] == False and end[2] == False and start[3] < end[3]: #DUP
            color = 'red'

        # y value is the insert size
        insert_size = abs(end[0] - start[1] - 1)
        if read_name in gap_registry:
            insert_size = gap_registry[read_name]

        # use this to scale the y-axis
        if not curr_min_insert_size or curr_min_insert_size < insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size > insert_size:
            curr_max_insert_size = insert_size

        ax.plot(p,\
                [insert_size,insert_size],\
                ':',\
                color=color, \
                alpha=0.25, \
                lw=1, \
                marker='o',
                markersize=marker_size)
    return [curr_min_insert_size, curr_max_insert_size]
#}}}

#{{{ def plot_long_reads(all_long_reads,
def plot_long_reads(long_reads,
                    ax,
                    range_min,
                    range_max,
                    curr_min_insert_size,
                    curr_max_insert_size):

    Path = mpath.Path

    colors = { 'ALIGN' : 'orange',
               'DEL' : 'black',
               'INVIN' : 'green',
               'INVOUT' : 'blue',
               'DUP' : 'red' }

    for read_name in long_reads:
        long_read_plan = get_long_read_plan(read_name,
                                            long_reads,
                                            range_min,
                                            range_max)
        max_gap = long_read_plan[0]
        steps = long_read_plan[1]
        for step in steps:
            step_cords = step[0]
            step_type = step[1]

            p = [float(step_cords[0]-range_min)/float(range_max - range_min),
                 float(step_cords[1]-range_min)/float(range_max - range_min)]


            if step_type == 'ALIGN':
                ax.plot(p,
                        [max_gap,max_gap],
                        '-', 
                        color=colors[step_type],
                        alpha=0.25,
                        lw=1)
            else:
                x1 = float(step_cords[0]-range_min)/float(range_max-range_min)
                x2 = float(step_cords[1]-range_min)/float(range_max-range_min)

                pp = mpatches.PathPatch(
                        Path([ (x1, max_gap),
                               (x1, max_gap*1.1),
                               (x2, max_gap*1.1),
                               (x2, max_gap)],
                        [ Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                        fc="none",
                        color=colors[step_type],
                        alpha=0.25,
                        lw=1,
                        ls=':')
                ax.add_patch(pp)

    return [curr_min_insert_size, curr_max_insert_size]
#}}}

#{{{parser = OptionParser()
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
                  help="Max number of normal pairs to plot")

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

parser.add_option("--long_read",
                  dest="long_read",
                  type=int,
                  default=1000,
                  help="Min length of a read to be a long-read (default 1000)")

parser.add_option("--common_insert_size",
                  dest="common_insert_size",
                  action="store_true",
                  default=False,
                  help="Set common insert size for all plots")

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

plot_height = 5
plot_width = 8

if options.plot_height:
    plot_height = options.plot_height
else:
    plot_height = 2 + len(options.bams.split(','))

if options.plot_width:
    plot_width = options.plot_width

# if an SV type is given, thn exapand the window around its bounds
if options.sv_type:
    window = int((int(options.end) - int(options.start))/2)
else:
    window = 0
if options.window:
    window = options.window
#}}}

all_pairs = []
all_splits = []
all_coverages = []
all_long_reads = []

max_plot_depth = 0

range_min = max(0,int(options.start) - window)
range_max = int(options.end) + window
min_insert_size = None
max_insert_size = None

bam_files = options.bams.split(',')

#{{{ read data from bams/crams
for bam_file_name in bam_files:
    bam_file = None
    if not options.reference:
        bam_file = pysam.AlignmentFile(bam_file_name, "rb")
    else:
        bam_file = pysam.AlignmentFile(bam_file_name, \
                                       "rc", \
                                       reference_filename=options.reference)

    pairs = {}
    splits = {}
    long_reads = {}
    coverage = {}

    for read in bam_file.fetch(options.chrom,
                               range_min, 
                               range_max):
        if read.query_length >= options.long_read:
            add_long_reads(read, long_reads, range_min, range_max)
        else:
            add_pair(read, pairs)
            add_split(read, splits, bam_file)
        add_coverage(read, coverage)

    pair_insert_sizes = [x[1][1]-x[0][0] for x in pairs.values() if len(x)==2]

    split_insert_sizes = [] 
    for read_name in splits:
        last = splits[read_name][0][0]
        for i in range(1, len(splits[read_name])):
            curr = splits[read_name][i][0]
            if curr >= range_min and curr <= range_max and \
                last >= range_min and last <= range_max:
                split_insert_sizes.append(abs(curr - last))
            last = curr

    long_read_gap_sizes = [] 
    for read_name in long_reads:
        long_read_gap_sizes.append(get_long_read_max_gap(read_name, long_reads))

    insert_sizes = pair_insert_sizes + split_insert_sizes + long_read_gap_sizes

    if not min_insert_size:
        min_insert_size = min(insert_sizes)
    else: 
        min_insert_size = min(min(insert_sizes), min_insert_size)

    if not max_insert_size:
        max_insert_size = max(insert_sizes)
    else: 
        max_insert_size = max(max(insert_sizes), min_insert_size)

    max_plot_depth = max(max_plot_depth, max(coverage.values()))
    all_coverages.append(coverage)
    all_pairs.append(pairs)
    all_splits.append(splits)
    all_long_reads.append(long_reads)
#}}}


# Sample +/- pairs in the normal insert size range
if options.max_depth:
    all_pairs = sample_normal(options.max_depth, all_pairs, options.z)

#{{{ set up sub plots
matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(plot_width, plot_height), dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)

# give one axis to display each sample
num_ax = len(options.bams.split(','))

# add another if we are displaying the SV
if options.sv_type:
    num_ax+=1

# add another if a annotation file is given
if options.transcript_file:
    num_ax+=1

if options.annotation_file:
    num_ax+=len(options.annotation_file.split(','))

# set the relative sizes for each
if options.sv_type:
    ratios = [4] + [10] * len(options.bams.split(','))
else:
    ratios = [10] * len(options.bams.split(','))

if options.annotation_file:
    ratios += [1]*len(options.annotation_file.split(','))
if options.transcript_file:
    ratios += [2]

gs = gridspec.GridSpec(num_ax, 1, height_ratios = ratios)
#}}}

ax_i = 0

if options.sv_type:
    ax = matplotlib.pyplot.subplot(gs[ax_i])
    plot_variant(options.start, options.end, options.sv_type, ax, range_min, range_max)
    ax_i += 1

#{{{ legend
marker_colors = ['black', 'red', 'blue', 'green', 'orange']
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
             "Aligned long read", \
             "Split-read", \
             "Paired-end read"], \
            loc=1,
            fontsize = 6,
            frameon=False)
#}}}

gap_registry = {}

# Plot each sample
for i in range(len(bam_files)):
    ax =  matplotlib.pyplot.subplot(gs[ax_i])

    curr_min_insert_size = None
    curr_max_insert_size = None

    curr_min_insert_size, 
    curr_max_insert_size = plot_long_reads(all_long_reads[i],
                                           ax,
                                           range_min,
                                           range_max,
                                           curr_min_insert_size,
                                           curr_max_insert_size)
    curr_min_insert_size, 
    curr_max_insert_size = plot_pairs(all_pairs[i],
                                      ax,
                                      range_min,
                                      range_max,
                                      curr_min_insert_size,
                                      curr_max_insert_size)

    curr_min_insert_size, 
    curr_max_insert_size = plot_splits(all_splits[i],
                                       ax,
                                       range_min,
                                       range_max,
                                       curr_min_insert_size,
                                       curr_max_insert_size)

    #{{{ try to only have y-lables, which does work if there are too few 
    # observed insert sizes
#    if options.common_insert_size:
#        ax.yaxis.set_ticks(np.arange(min_insert_size, \
#                                     max_insert_size, \
#                                     max(2,
#                                         (max_insert_size-min_insert_size)/4)))
#    else:
#        if max_insert_size-min_insert_size != 0: 
#            ax.yaxis.set_ticks(np.arange(curr_min_insert_size, \
#                                         curr_max_insert_size, \
#                                         max(2,
#                                             (curr_max_insert_size-curr_min_insert_size)/4)))
#    #}}}

    #{{{ set the axis title to be either one passed in or filename
    if options.titles and \
            len(options.titles.split(',')) == len(options.bams.split(',')):
        ax.set_title(options.titles.split(',')[ax_i-1], \
                     fontsize=8, loc='left')
    else:
        ax.set_title(os.path.basename(options.bams.split(',')[ax_i-1]), \
                     fontsize=8, loc='left')
    #}}}
    
    #{{{ set axis parameters
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
    #}}}

    #{{{ Plot coverage
    cover_x = []
    cover_y = []

    #if max([d[1] for d in all_depths[i]]) == 0: continue

    coverage = all_coverages[i]

    for pos in range(range_min,range_max+1):
        if pos in coverage:
            cover_x.append(float(pos-range_min)/float(range_max - range_min))
            cover_y.append(coverage[pos])
        else:
            cover_x.append(float(pos-range_min)/float(range_max - range_min))
            cover_y.append(0)

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
    #}}}

    ax_i += 1

#{{{ Plot annoation files
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
         
        #ax =  matplotlib.pyplot.subplot(gs[
                #1 + len(options.bams.split(',')) + a_i])
        ax =  matplotlib.pyplot.subplot(gs[a_i])
        a_i += 1

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
#}}}

#{{{ Plot sorted/bgziped/tabixed transcript file
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
#}}}

# save
matplotlib.rcParams['agg.path.chunksize'] = 100000
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(options.output_file)

#{{{ give sv-plaudit output
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
#}}}
