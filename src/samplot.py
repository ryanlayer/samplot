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
from matplotlib.offsetbox import AnchoredText

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
def sample_normal(max_depth, pairs, z):
    sampled_pairs = {}
    plus_minus_pairs = {}

    for read_name in pairs:
        pair = pairs[read_name]
        if len(pair) != 2: 
            continue
        if pair[0].strand == True and pair[1].strand == False:
            plus_minus_pairs[read_name] = pair
        else:
            sampled_pairs[read_name] = pair

    if len(plus_minus_pairs) > max_depth:
        lens = [pair[1].end - pair[0].start \
                for pair in plus_minus_pairs.values()]
        mean = statistics.mean(lens)
        stdev = statistics.stdev(lens)

        inside_norm = {}

        for read_name in pairs:
            pair = pairs[read_name]
            if len(pair) != 2: 
                continue
            if pair[1].end - pair[0].start >= mean + options.z*stdev:
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

    return sampled_pairs
#}}}

#{{{ def add_coverage(read, coverage):
def add_coverage(read, coverage):
    hp = 0

    if read.has_tag('HP'):
        hp = int(read.get_tag('HP'))

    if hp not in coverage:
        coverage[hp] = {}

    curr_pos = read.reference_start
    if not read.cigartuples: return

    for op,length in read.cigartuples:
        if op == 0:
            for pos in range(curr_pos, curr_pos+length+1):
                if pos not in coverage[hp]:
                    coverage[hp][pos] = 0
                coverage[hp][pos] += 1
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

#{{{ class PairEnd:
class PairEnd:
    def __init__(self, read):
        self.start = read.reference_start
        self.end = read.reference_end
        self.strand = not(read.is_reverse)
        self.MI = None

        if read.has_tag('MI'):
            self.MI = int(read.get_tag('MI'))

        self.HP = 0

        if read.has_tag('HP'):
            self.HP =  int(read.get_tag('HP'))
#}}}

#{{{def add_pair_end(read, pairs):
def add_pair_end(read, pairs, linked_reads):
    if read.is_unmapped: return
    if not (read.is_paired): return
    if read.is_secondary: return
    if read.is_supplementary: return

    pe = PairEnd(read) 

    if pe.HP not in pairs:
        pairs[pe.HP] = {}

    if read.query_name not in pairs[pe.HP]:
        pairs[pe.HP][read.query_name] = []

    if pe.MI:
        if pe.HP not in linked_reads:
            linked_reads[pe.HP] = {}

        if pe.MI not in linked_reads[pe.HP]:
            linked_reads[pe.HP][pe.MI] = [[],[]]
        linked_reads[pe.HP][pe.MI][0].append(read.query_name)

    pairs[pe.HP][read.query_name].append( pe )
    pairs[pe.HP][read.query_name].sort(key=lambda x:x.start)
#}}}

#{{{ class SplitRead:
class SplitRead:
    def __init__(self, start,end,strand,query_pos):
        self.start = start
        self.end = end
        self.strand = strand
        self.query_pos = query_pos
        self.MI = None

        if read.has_tag('MI'):
            self.MI = int(read.get_tag('MI'))

        self.HP = 0

        if read.has_tag('HP'):
            self.HP =  int(read.get_tag('HP'))
#}}}

#{{{def add_split(read, splits):
def add_split(read, splits, bam_file, linked_reads):
    if not read.is_secondary and \
       not read.is_supplementary and \
       read.has_tag('SA'):
        qs_pos, qe_pos = \
            calc_query_pos_from_cigar(read.cigarstring, \
                                      (not read.is_reverse))

        sr = SplitRead(read.reference_start,
                       read.reference_end,
                       not(read.is_reverse),
                       qs_pos)

        if sr.MI:
            if sr.HP not in linked_reads:
                linked_read[sr.HP] = {}
            if sr.MI not in linked_reads[sr.HP]:
                linked_reads[sr.HP][sr.MI] = [[],[]]
            linked_reads[sr.HP][sr.MI][1].append(read.query_name)

        if sr.HP not in splits:
            splits[sr.HP] = {}

        splits[sr.HP][read.query_name]=[sr]

        
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
            splits[sr.HP][read.query_name].append(SplitRead(pos,
                                                  pos+qe_pos,
                                                  strand,
                                                  qs_pos))

        if len(splits[sr.HP][read.query_name]) == 1:
            del splits[sr.HP][read.query_name]
        else:
            splits[sr.HP][read.query_name].sort(key=lambda x:x.start)
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

    hp = 0 

    if read.has_tag('HP'):
        hp = int(read.get_tag('HP'))

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

    if hp not in long_reads:
        long_reads[hp] = {}

    if read.query_name not in long_reads[hp]:
        long_reads[hp][read.query_name] = []

    long_reads[hp][read.query_name].append(LongRead(read.reference_start,
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
    #matplotlib.pyplot.tick_params(axis='x',length=0)
    #matplotlib.pyplot.tick_params(axis='y',length=0)
    ax.tick_params(axis='x',length=0)
    ax.tick_params(axis='y',length=0)
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

#{{{ def plot_pair(pair, y, ax, range_min, range_max):
def plot_pair(pair, y, ax, range_min, range_max):
    colors = { (True, False): 'black', # DEL
               (False, True): 'red',   # DUP
               (False, False): 'blue', # INV
               (True, True): 'blue' } # INV
               #(True, True): 'green' } # INV

    if pair[0].end < range_min  or pair[1].start > range_max:
        return

    p = [float(pair[0].start - range_min)/float(range_max - range_min), \
         float(pair[1].end - range_min)/float(range_max - range_min)]


    color = colors[(pair[0].strand, pair[1].strand)]

    # plot the individual pair
    ax.plot(p,\
            [y,y],\
            '-',color=color, \
            alpha=0.25, \
            lw=0.5, \
            marker='s', \
            markersize=marker_size)
#}}}

#{{{ def plot_pairs(all_pairs
def plot_pairs(pairs,
               ax,
               range_min,
               range_max,
               curr_min_insert_size,
               curr_max_insert_size):

    ## set the color of the pair based on the two strands
    for read_name in pairs:
        pair = pairs[read_name]

        if len(pair) != 2: continue

        if pair[0].MI or pair[1].MI: continue

        # y value is the insert size
        insert_size = pair[1].end - pair[0].start
        # use this to scale the y-axis
        if not curr_min_insert_size or curr_min_insert_size < insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size > insert_size:
            curr_max_insert_size = insert_size

        plot_pair(pair, insert_size, ax, range_min, range_max)

    return [curr_min_insert_size, curr_max_insert_size]
#}}}

#{{{ def plot_linked_reads(pairs,
def plot_linked_reads(pairs,
                      splits,
                      linked_reads,
                      ax,
                      range_min,
                      range_max,
                      curr_min_insert_size,
                      curr_max_insert_size):

    for linked_read in linked_reads:
        linked_pairs = []
        for name in linked_reads[linked_read][0]:
            if name in pairs and len(pairs[name]) == 2:
                linked_pairs.append(pairs[name])

        linked_splits = []
        for name in linked_reads[linked_read][1]:
            if name in splits:
                linked_splits.append(splits[name])

        if len(linked_pairs) + len(linked_splits) == 0: continue

        ends = []
        for linked_pair in linked_pairs:
            for pair in linked_pair:
                ends.append(pair.start)
                ends.append(pair.end)
        for linked_split in linked_splits:
            for split in linked_split:
                ends.append(split.start)
                ends.append(split.end)

        min_end = min(ends)
        max_end = max(ends)

        gap_sizes = []
        alignments = []
        for linked_pair in linked_pairs:

            alignments.append([linked_pair[0].start,linked_pair[1].end])

            if linked_pair[1].end > range_max or \
               linked_pair[0].start < range_min:
                continue
            gap_sizes.append(abs(linked_pair[1].end - linked_pair[0].start))

        for linked_split in linked_splits:
            alignments.append([linked_split[0].end,linked_split[1].start])
            if linked_split[1].start > range_max or \
               linked_split[0].end < range_min:
                continue
            gap_sizes.append(abs(linked_split[1].start - linked_split[0].end))

        if len(gap_sizes) == 0 : continue

        insert_size = max(gap_sizes)

        if not curr_min_insert_size or curr_min_insert_size < insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size > insert_size:
            curr_max_insert_size = insert_size

        alignments.sort(key=lambda x: x[0])
        compressed_alignments = []
        compressed_alignments.append(alignments[0])

        for alignment in alignments[1:]:
            if alignment[0] < compressed_alignments[-1][1]:
                compressed_alignments[-1][1] = \
                        max(compressed_alignments[-1][1],
                            alignment[1])
            else:
                compressed_alignments.append(alignment)

        alignments=compressed_alignments

        start = alignments[0][0]
        end = alignments[-1][1]
        p = [float(start - range_min)/float(range_max - range_min), \
             float(end - range_min)/float(range_max - range_min)]

        ax.plot(p,
                [insert_size,insert_size],\
                '-',color='green', \
                alpha=0.75, \
                lw=0.25)

#        last = alignments[0]
#        for curr in alignments[1:]:
#            start = last[1]
#            end = curr[0]
#            p = [float(start - range_min)/float(range_max - range_min), \
#                 float(end - range_min)/float(range_max - range_min)]
#
#            ax.plot(p,
#                    [insert_size,insert_size],\
#                    '-',color='purple', \
#                    alpha=0.05, \
#                    lw=4)
#            last = curr
#
        for name in linked_reads[linked_read][0]:
            if name in pairs and len(pairs[name]) == 2:
                pair = pairs[name]
                plot_pair(pair, insert_size, ax, range_min, range_max)

        for name in linked_reads[linked_read][1]:
            if name in splits:
                split = splits[name]
                plot_split(split, insert_size, ax, range_min, range_max)

    return [curr_min_insert_size, curr_max_insert_size]
#}}}

#{{{def plot_split(split, y, ax, range_min, range_max):
def plot_split(split, y, ax, range_min, range_max):
    start = split[0]
    end = split[1]
    if start.start > end.end:
        end = split[0]
        start = split[1]

    # Do not plot pairs that extend beyond the current range
    if range_min > start.end or range_max < end.start:
        return
         
    p = [float(start.end - range_min)/float(range_max - range_min), \
         float(end.start - range_min)/float(range_max - range_min)]
       
    # For a given SV, the orientation of the pairs and split do not match
    # so we cannot use the colors dict here
    color = 'black'
    if start.strand != end.strand: #INV
        if start.strand == True: 
            #color = 'green'
            color = 'blue'
        else:
            color = 'blue'
    elif start.strand == True and \
         end.strand == True and \
         start.query_pos < end.query_pos: #DEL
        color = 'black'
    elif start.strand == False and \
         end.strand == False and \
         start.query_pos > end.query_pos: #DEL
        color = 'black'
    elif start.strand == True and \
         end.strand == True and \
         start.query_pos > end.query_pos: #DUP
        color = 'red'
    elif start.strand == False and \
         end.strand == False and \
         start.query_pos < end.query_pos: #DUP
        color = 'red'

    ax.plot(p,\
            [y,y],\
            ':',\
            color=color, \
            alpha=0.25, \
            lw=1, \
            marker='o',
            markersize=marker_size)
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

        # y value is the insert size
        start = split[0]
        end = split[1]
        if start.start > end.end:
            end = split[0]
            start = split[1]
        insert_size = abs(end.start - start.end - 1)

        # use this to scale the y-axis
        if not curr_min_insert_size or curr_min_insert_size < insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size > insert_size:
            curr_max_insert_size = insert_size

        plot_split(split, insert_size, ax, range_min, range_max)

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
               #'INVIN' : 'green',
               'INVIN' : 'blue',
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

#{{{def plot_coverage(coverage,
def plot_coverage(coverage,
                  ax,
                  range_min,
                  range_max):
    cover_x = []
    cover_y = []

    for pos in range(range_min,range_max+1):
        if pos in coverage:
            cover_x.append(\
                    float(pos-range_min)/float(range_max - range_min))
            cover_y.append(coverage[pos])
        else:
            cover_x.append(\
                    float(pos-range_min)/float(range_max - range_min))
            cover_y.append(0)

    max_plot_depth = max(cover_y)
    ax2 = ax.twinx()
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,max_plot_depth])
    ax2.fill_between(cover_x, \
                     cover_y, \
                     [0] * len(cover_y),
                     color='grey',
                     alpha=0.25)
 
    # set axis parameters
    #ax2.set_ylabel('Coverage', fontsize=8)
    ax2.tick_params(axis='y', colors='grey', labelsize=6)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    #matplotlib.pyplot.tick_params(axis='x',length=0)
    #matplotlib.pyplot.tick_params(axis='y',length=0)
    ax2.tick_params(axis='x',length=0)
    ax2.tick_params(axis='y',length=0)

    return ax2
#}}}

#{{{def get_pair_insert_sizes(pairs):
def get_pair_insert_sizes(pairs):
    pair_insert_sizes = []

    for hp in pairs:
        for read_name in pairs[hp]:
            if len(pairs[hp][read_name]) == 2:
                first = pairs[hp][read_name][0]
                second = pairs[hp][read_name][1]
                pair_insert_sizes.append(second.end - first.start)
    return pair_insert_sizes
#}}}

#{{{def get_split_insert_size(splits):
def get_split_insert_size(splits):
    split_insert_sizes = []

    for hp in splits:
        for read_name in splits[hp]:
            last = splits[hp][read_name][0].start

            for i in range(1, len(splits[hp][read_name])):
                curr = splits[hp][read_name][i].start
                if curr >= range_min and curr <= range_max and \
                    last >= range_min and last <= range_max:
                    split_insert_sizes.append(abs(curr - last))
                last = curr
    return split_insert_sizes
#}}}

#{{{def get_long_read_gap_sizes(long_reads):
def get_long_read_gap_sizes(long_reads):
    long_read_gap_sizes = [] 

    for hp in long_reads:
        for read_name in long_reads[hp]:
            long_read_gap_sizes.append(\
                    get_long_read_max_gap(read_name, long_reads[hp]))
    return long_read_gap_sizes
#}}}

#{{{parser = OptionParser()
parser = OptionParser()

parser.add_option("--marker_size",
                  dest="marker_size",
                  type=int,
                  default=3,
                  help="Size of marks on pairs and splits (default 3) ");


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

parser.add_option("-j",
                  dest="json_only",
                  action="store_true",
                  default=False,
                  help="Create only the json file, not the image plot")

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
    parser.error('BAMs not given')
    
if not options.start:
    parser.error('SV start not given')
    
if not options.end:
    parser.error('SV end not given')

if not options.chrom:
    parser.error('SV chrom not given')

if not options.json_only:
    plot_height = 5
    plot_width = 8
    marker_size = options.marker_size

    if options.plot_height:
        plot_height = options.plot_height
    else:
        plot_height = 2 + len(options.bams.split(','))

    if options.plot_width:
        plot_width = options.plot_width

    # if an SV type is given, then expand the window around its bounds
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
    all_linked_reads = []

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
        linked_reads = {}

        for read in bam_file.fetch(options.chrom,
                                   max(0,range_min-1000), 
                                   range_max+1000):
            if read.query_length >= options.long_read:
                add_long_reads(read, long_reads, range_min, range_max)
            else:
                add_pair_end(read, pairs, linked_reads)
                add_split(read, splits, bam_file, linked_reads)
            add_coverage(read, coverage)


        pair_insert_sizes = get_pair_insert_sizes(pairs)
        split_insert_sizes = get_split_insert_size(splits)
        long_read_gap_sizes = get_long_read_gap_sizes(long_reads)

        insert_sizes = pair_insert_sizes + \
                       split_insert_sizes + \
                       long_read_gap_sizes

        if not min_insert_size:
            min_insert_size = min(insert_sizes)
        else: 
            min_insert_size = min(min(insert_sizes), min_insert_size)

        if not max_insert_size:
            max_insert_size = max(insert_sizes)
        else: 
            max_insert_size = max(max(insert_sizes), min_insert_size)

        all_coverages.append(coverage)
        all_pairs.append(pairs)
        all_splits.append(splits)
        all_long_reads.append(long_reads)
        all_linked_reads.append(linked_reads)
    #}}}

    #{{{ Sample +/- pairs in the normal insert size range
    if options.max_depth:
        for i in range(len(all_pairs)):
            for hp in all_pairs[i]:
                all_pairs[i][hp] = sample_normal(options.max_depth,
                                                 all_pairs[i][hp],
                                                 options.z)
        #all_pairs = sample_normal(options.max_depth, all_pairs, options.z)
    #}}}

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
    ratios = []
    if options.sv_type:
        ratios = [1] 

    for i in range(len(options.bams.split(','))):
        ratios.append( len(all_coverages[i]) * 3 )

    if options.annotation_file:
        ratios += [1]*len(options.annotation_file.split(','))
    if options.transcript_file:
        ratios += [2]

    gs = gridspec.GridSpec(num_ax, 1, height_ratios = ratios)
    #}}}

    ax_i = 0

    #{{{ plot variant on top
    if options.sv_type:
        ax = matplotlib.pyplot.subplot(gs[ax_i])
        plot_variant(options.start,
                     options.end,
                     options.sv_type,
                     ax,
                     range_min,
                     range_max)
        ax_i += 1
    #}}}

    #{{{ plot legend
    #marker_colors = ['black', 'red', 'blue', 'green', 'orange', 'purple']
    marker_colors = ['black', 'red', 'blue', 'orange', 'green']
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
                 "Aligned long read", \
                 "Linked read", \
                 "Split-read", \
                 "Paired-end read"], \
                loc=1,
                fontsize = 6,
                frameon=False)
    #}}}

    # Plot each sample
    for i in range(len(bam_files)):
        ax =  matplotlib.pyplot.subplot(gs[ax_i])
        inner_axs = gridspec.GridSpecFromSubplotSpec(len(all_coverages[i]), 
                                                     1,
                                                     subplot_spec=gs[ax_i],
                                                     wspace=0.0,
                                                     hspace=0.2)

        hps = sorted(all_coverages[i].keys(), reverse=True)
        axs = {}
        for j in range(len(hps)):
            axs[j] = matplotlib.pyplot.subplot(inner_axs[hps[j]])

        curr_min_insert_size = None
        curr_max_insert_size = None

        cover_axs = {}
        for hp in hps:
            curr_ax = axs[hp]

            curr_pairs = []
            if hp in all_pairs[i]:
                curr_pairs = all_pairs[i][hp]

            curr_splits = []
            if hp in all_splits[i]:
                curr_splits = all_splits[i][hp]

            curr_linked_reads = []
            if hp in all_linked_reads[i]:
                curr_linked_reads = all_linked_reads[i][hp]

            curr_long_reads = []
            if hp in all_long_reads[i]:
                curr_long_reads = all_long_reads[i][hp]

            curr_coverage = []
            if hp in all_coverages[i]:
                curr_coverage = all_coverages[i][hp]

            curr_min_insert_size,curr_max_insert_size = plot_linked_reads(curr_pairs,
                                                     curr_splits,
                                                     curr_linked_reads,
                                                     curr_ax,
                                                     range_min,
                                                     range_max,
                                                     curr_min_insert_size,
                                                     curr_max_insert_size)

            curr_min_insert_size,curr_max_insert_size = plot_long_reads(curr_long_reads,
                                                   curr_ax,
                                                   range_min,
                                                   range_max,
                                                   curr_min_insert_size,
                                                   curr_max_insert_size)
            curr_min_insert_size,curr_max_insert_size = plot_pairs(curr_pairs,
                                              curr_ax,
                                              range_min,
                                              range_max,
                                              curr_min_insert_size,
                                              curr_max_insert_size)

            curr_min_insert_size,curr_max_insert_size = plot_splits(curr_splits,
                                               curr_ax,
                                               range_min,
                                               range_max,
                                               curr_min_insert_size,
                                               curr_max_insert_size)

            cover_ax = plot_coverage(curr_coverage,
                                     curr_ax,
                                     range_min,
                                     range_max)
            cover_axs[hp] = cover_ax
     
        #{{{ set axis parameters
        #set the axis title to be either one passed in or filename
        curr_ax = axs[hps[0]]
        if options.titles and \
                len(options.titles.split(',')) == len(options.bams.split(',')):
            curr_ax.set_title(options.titles.split(',')[ax_i-1], \
                         fontsize=8, loc='left')
        else:
            curr_ax.set_title(os.path.basename(options.bams.split(',')[ax_i-1]), \
                         fontsize=8, loc='left')


        if len(axs) > 1:
            for j in axs:
                curr_ax = axs[j]
                fp = dict(size=8, backgroundcolor='white')
                text = 'HP: '
                if j == 0:
                    text += 'Undef'
                else:
                    text += str(j)
                at = AnchoredText(text,
                                  loc=2,
                                  prop=fp,
                                  borderpad=0,
                                  pad=0,
                                  frameon=False)
                curr_ax.add_artist(at)
        
        for j in hps:
            curr_ax = axs[j]
            curr_ax.set_xlim([0,1])
            curr_ax.spines['top'].set_visible(False)
            curr_ax.spines['bottom'].set_visible(False)
            curr_ax.spines['left'].set_visible(False)
            curr_ax.spines['right'].set_visible(False)
            curr_ax.tick_params(axis='y', labelsize=6)
            curr_ax.tick_params(axis='both', length=0)
            curr_ax.set_xticklabels([])

        if ax_i == num_ax - 1:
            curr_ax = axs[ hps[-1] ]
            labels = [int(range_min + l*(range_max-range_min)) \
                    for l in curr_ax.xaxis.get_majorticklocs()]
            curr_ax.set_xticklabels(labels, fontsize=6)
            curr_ax.set_xlabel('Chromosomal position on ' + \
                    options.chrom, fontsize=8)

        curr_ax = axs[hps[ int(len(hps)/2)    ]]
        curr_ax.set_ylabel('Insert size', fontsize=8)
        cover_ax = cover_axs[hps[ int(len(hps)/2)    ]]
        cover_ax.set_ylabel('Coverage', fontsize=8)
        #}}}

        ax_i += 1

    #{{{ Plot annoation files
    if options.annotation_file:
        #a_i = 0
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
             
            ax =  matplotlib.pyplot.subplot(gs[ax_i])
            ax_i += 1

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

            #matplotlib.pyplot.tick_params(axis='x',length=0)
            #matplotlib.pyplot.tick_params(axis='y',length=0)
            ax.tick_params(axis='x',length=0)
            ax.tick_params(axis='y',length=0)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            #a_i+=1
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

        #matplotlib.pyplot.tick_params(axis='x',length=0)
        #matplotlib.pyplot.tick_params(axis='y',length=0)
        ax.tick_params(axis='x',length=0)
        ax.tick_params(axis='y',length=0)
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
if options.print_args or options.json_only:
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
