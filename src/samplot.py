#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import pylab
import pysam
import os
import re
import statistics
import random
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import argparse
from matplotlib.offsetbox import AnchoredText
import matplotlib.ticker as ticker

COLORS = { 
    'Deletion/Normal': 'black',
    'Duplication': 'red',
    'Inversion': 'blue'
}



READ_TYPES_USED = {
    "Deletion/Normal":False,
    "Duplication":False,
    "Inversion":False,
    "Aligned long read": False,
    "Linked read": False,
    "Split-read": False,
    "Paired-end read": False
}

#pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
CIGAR_MAP = { 
    'M' : 0,
    'I' : 1,
    'D' : 2,
    'N' : 3,
    'S' : 4,
    'H' : 5,
    'P' : 6,    
    '=' : 7,
    'X' : 8,
    'B' : 9
}

def calc_query_pos_from_cigar(cigar, strand):
    """Uses the CIGAR string to determine the query position of a read

    The cigar arg is a string like the following: 86M65S
    The strand arg is a boolean, True for forward strand and False for reverse

    Returns pair of ints for query start, end positions
    """

    cigar_ops = [[int(op[0]), op[1]] for op in re.findall('(\d+)([A-Za-z])', cigar)]

    order_ops = cigar_ops
    if not strand: # - strand
        order_ops = order_ops[::-1]

    qs_pos = 0
    qe_pos = 0
    q_len = 0

    for op_position in range(len(cigar_ops)):
        op_len = cigar_ops[op_position][0]
        op_type = cigar_ops[op_position][1]

        if op_position == 0 and ( op_type == 'H' or op_type == 'S' ):
            qs_pos += op_len
            qe_pos += op_len
            q_len += op_len
        elif op_type == 'H' or op_type == 'S':
            q_len += op_len
        elif op_type == 'M' or op_type == 'I' or op_type == 'X':
            qe_pos += op_len
            q_len += op_len

    return qs_pos, qe_pos

def sample_normal(max_depth, pairs, z):
    """Downsamples paired-end reads 
    
    Selects max_depth reads
    Does not remove discordant pairs, those with insert distance greater than z stdevs from mean

    Returns downsampled pairs list
    """

    sampled_pairs = {}
    plus_minus_pairs = {}

    if max_depth == 0:
        return sampled_pairs

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
            if pair[1].end - pair[0].start >= mean + z*stdev:
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

def add_coverage(read, coverage, minq):
    """Adds a read to the known coverage 
    
    Coverage from Pysam read is added to the coverage list
    Coverage list is a pair of high- and low-quality lists
    Quality is determined by minq, which is min quality
    """

    hp = 0

    if read.has_tag('HP'):
        hp = int(read.get_tag('HP'))

    if hp not in coverage:
        coverage[hp] = {}

    curr_pos = read.reference_start
    if not read.cigartuples: return

    for op,length in read.cigartuples:
        if op in [CIGAR_MAP['M'], CIGAR_MAP['='], CIGAR_MAP['X']]:
            for pos in range(curr_pos, curr_pos+length):
                if pos not in coverage[hp]:
                    coverage[hp][pos] = [0,0]

                #the two coverage tracks are [0] high-quality and [1] low-quality
                if read.mapping_quality > minq:
                    coverage[hp][pos][0] += 1
                else:
                    coverage[hp][pos][1] += 1
            curr_pos += length
        elif op == CIGAR_MAP['I']:
            curr_pos = curr_pos
        elif op == CIGAR_MAP['D']:
            curr_pos += length
        elif op == CIGAR_MAP['N']:
            curr_pos = length
        elif op == CIGAR_MAP['S']:
            curr_pos = curr_pos
        elif op == CIGAR_MAP['H']:
            curr_pos = curr_pos
        else:
            curr_pos += length

class PairedEnd:
    """container of paired-end read info

    Contains start(int), end(int), strand(bool True=forward), MI (int molecular identifier), HP (int haplotype)
    """
    
    def __init__(self, start, end, is_reverse, MI_tag, HP_tag):
        """Create PairedEnd instance

        Genomic interval is defined by start and end integers
        Strand is opposite of is_reverse
        Molecular identifier and Haplotype are integers if present, else False
        """
        self.start = start
        self.end = end
        self.strand = not(is_reverse)
        # molecular identifier - linked reads only
        self.MI = None
        #haplotype - phased reads only
        self.HP = 0

        if MI_tag:
            self.MI = MI_tag
        if HP_tag:
            self.HP = HP_tag
    
    def __repr__(self):
        return ('PairedEnd(%s,%s,%s,%s,%s)' % \
                (self.start, self.end, self.strand,self.MI,self.HP))

def add_pair_end(read, pairs, linked_reads):
    """adds a (mapped, primary, non-supplementary, and paired) read to the pairs list

    Pysam read is added as simpified PairedEnd instance to pairs
    Also added to linked_reads list if there is an associated MI tag
    """

    if read.is_unmapped: return
    if not (read.is_paired): return
    if read.is_secondary: return
    if read.is_supplementary: return

    MI_tag = False
    HP_tag = False

    if read.has_tag('MI'):
        MI_tag = int(read.get_tag('MI'))
    if read.has_tag('HP'):
        HP_tag = int(read.get_tag('HP'))
    
    READ_TYPES_USED["Paired-end read"] = True

    pe = PairedEnd(read.reference_start, 
        read.reference_end, 
        read.is_reverse, 
        MI_tag, 
        HP_tag) 

    if pe.HP not in pairs:
        pairs[pe.HP] = {}

    if read.query_name not in pairs[pe.HP]:
        pairs[pe.HP][read.query_name] = []

    if pe.MI:
        READ_TYPES_USED["Linked read"] = True
        if pe.HP not in linked_reads:
            linked_reads[pe.HP] = {}

        if pe.MI not in linked_reads[pe.HP]:
            linked_reads[pe.HP][pe.MI] = [[],[]]
        linked_reads[pe.HP][pe.MI][0].append(read.query_name)

    pairs[pe.HP][read.query_name].append( pe )
    pairs[pe.HP][read.query_name].sort(key=lambda x:x.start)

class SplitRead:
    """container of split read info

    Contains start(int), end(int), strand(bool True=forward), query position (int), MI (int molecular identifier), HP (int haplotype)
    """
    def __init__(self, start,end,strand,query_pos, MI_tag=None, HP_tag=None):
        """Create SplitRead instance

        Genomic interval is defined by start, end, and query_pos integers
        Strand is opposite of is_reverse
        Molecular identifier and Haplotype are integers if present, else False
        """
        self.start = start
        self.end = end
        self.strand = strand
        self.query_pos = query_pos
        # molecular identifier - linked reads only
        self.MI = None
        #haplotype - phased reads only
        self.HP = 0

        if MI_tag:
            self.MI = MI_tag
        if HP_tag:
            self.HP = HP_tag
    
    def __repr__(self):
        return ('SplitRead(%s,%s,%s,%s,%s,%s)' % \
                (self.start, self.end, self.strand, self.query_pos, self.MI, self.HP))


def add_split(read, splits, bam_file, linked_reads):
    """adds a (primary, non-supplementary) read to the splits list

    Pysam read is added as simpified SplitRead instance to splits
    Also added to linked_reads list if there is an associated MI tag
    """
    if read.is_secondary: return
    if read.is_supplementary: return
    if not read.has_tag('SA'): return

    READ_TYPES_USED["Split-read"] = True
    qs_pos, qe_pos = calc_query_pos_from_cigar(read.cigarstring, (not read.is_reverse))
    
    HP_tag = False 
    MI_tag = False
    if read.has_tag('MI'):
        MI_tag = int(read.get_tag('MI'))

    if read.has_tag('HP'):
        HP_tag = int(read.get_tag('HP'))
    sr = SplitRead(read.reference_start,
                   read.reference_end,
                   not(read.is_reverse),
                   qs_pos,
                   MI_tag,
                   HP_tag)

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
        splits[sr.HP][read.query_name].append(SplitRead(pos,pos+qe_pos, strand, qs_pos))

    if len(splits[sr.HP][read.query_name]) == 1:
        del splits[sr.HP][read.query_name]
    else:
        splits[sr.HP][read.query_name].sort(key=lambda x:x.start)

class Alignment:
    """container of alignment info, from CIGAR string

    Contains start(int), end(int), strand(bool True=forward), query position (int)
    """
    def __init__(self,start,end,strand,query_position):
        """Create Alignment instance

        Genomic interval is defined by start, end, and query_pos integers
        Strand is bool (True for forward)
        """
        self.start = start
        self.end = end
        self.strand = strand
        self.query_position = query_position
    def __str__(self):
        return ','.join([str(x) for x in [self.start, 
                                          self.end, 
                                          self.strand, 
                                          self.query_position]])

    def __repr__(self):
        return ('Alignment(%s,%s,%s,%s)' % \
                (self.start, self.end, self.strand,self.query_position))

class LongRead:
    """container of LongRead info

    Contains start(int), end(int), list of Alignments
    """

    def __init__(self,start,end,alignments):
        """Create LongRead instance

        Genomic interval is defined by start, end integers
        List of Alignments set by parameter
        """
        self.start = start
        self.end = end
        self.alignments = alignments
    def __str__(self):
        return ','.join([str(x) for x in [self.start, 
                                          self.end, 
                                          len(self.alignments)]])
    def __repr__(self):
        return ('LongRead(%s,%s,%s)' % \
                (self.start, self.end, len(self.alignments)))
    
def get_alignments_from_cigar(curr_pos, strand, cigartuples, reverse=False):
    """Breaks CIGAR string into individual Aignments

    Starting point within genome given by curr_pos and strand
    Set of CIGAR operations and lengths as pairs passed in as cigartuples
    Direction of alignment set to reverse with reverse boolean

    Return list of Alignments
    """
    alignments = []
    q_pos = 0
    if reverse:
        cigartuples = cigartuples[::-1]
    for op,length in cigartuples:
        if op in [CIGAR_MAP['M'],CIGAR_MAP['='],CIGAR_MAP['X']]:
            alignments.append(Alignment(curr_pos,
                                        curr_pos+length,
                                        strand,
                                        q_pos))
            curr_pos += length
            q_pos += length
        elif op == CIGAR_MAP['I']: 
            q_pos += length
        elif op == CIGAR_MAP['D']:
            curr_pos += length
        elif op == CIGAR_MAP['N']:
            curr_pos += length
        elif op == CIGAR_MAP['S']:
            q_pos += length
    return alignments

def get_cigartuples_from_string(cigarstring):
    """Extracts operations,lengths as tuples from cigar string"

    Returns list of tuples of [operation,length]
    """
    cigartuples = []
    for match in re.findall(r'(\d+)([A-Z]{1})', cigarstring):
        length = int(match[0])
        op = match[1]
        cigartuples.append((CIGAR_MAP[op], length))
    
    return cigartuples

def merge_alignments(min_gap, alignments):
    """Combines previously identified alignments if close together

    Alignments are combined if within min_gap distance
    
    Returns list of Alignments
    """

    merged_alignments = []

    for alignment in alignments:
        if len(merged_alignments) == 0:
            merged_alignments.append(alignment)
        else:
            if alignment.start < merged_alignments[-1].end + min_gap:
                merged_alignments[-1].end = alignment.end
            else:
                merged_alignments.append(alignment)
    return merged_alignments

def add_long_reads(read, long_reads, range_min, range_max, min_event_size):
    """Adds a (primary, non-supplementary, long) read to the long_reads list

    Read added to long_reads if within the inteval defined by range_min, range_max
    Alignments belonging to the LongRead instance combined if within the min_event_size distance apart
    """
    READ_TYPES_USED["Aligned long read"] = True

    if read.is_supplementary or read.is_secondary: return

    hp = 0 

    if read.has_tag('HP'):
        hp = int(read.get_tag('HP'))

    alignments = get_alignments_from_cigar(read.pos,
                                           not read.is_reverse,
                                           read.cigartuples)
    min_gap = min_event_size
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

def get_long_read_plan(read_name, long_reads, range_min, range_max):
    """Create a plan to render a long read

    Plan consists of the largest event within the read 
        (used to determine the y-axis position of read)
        and the alignment types for plotting each Alignment within 
        LongRead.alignments ALIGN, DUP, DEL, INVIN, INVOUT

    Returns plan
    """

    alignments = []

    for long_read in long_reads[read_name]:
        for alignment in long_read.alignments:
            alignments.append(alignment)

    if len(alignments) <= 0:
        return None
    alignments.sort(key=lambda x: x.query_position)

    gaps = []

    curr = alignments[0]

    # we set the primary strand to be the one with the longest alignment
    # this will affect which alignment is inverted. There are clearly edge
    # cases here that we will need to address as we get more examples 
    # of inversions

    longest_alignment = 0
    longest_alignment_i = -1
    for i in range(len(alignments)):
        if longest_alignment < alignments[i].end - alignments[i].start:
            longest_alignment = alignments[i].end - alignments[i].start
            longest_alignment_i = i
    primary_strand = alignments[longest_alignment_i].strand

    steps = []
    steps.append( [ [curr.start,curr.end], 'ALIGN' ] )

    for i in range(1,len(alignments)):
        last = alignments[i-1]
        curr = alignments[i]

        gap = 0

        # figure out what the event is

        # INV
        if (curr.strand != last.strand):
            gap = 0

            # it is possible that we have a complex even that 
            # is an inverted DUP
            if (curr.start  < last.end):
                steps.append( [ [last.end,curr.start], 'DUP' ] )

            if (curr.strand != primary_strand):
                # last (primary) | curr 
                # +++++++++++++++|-------
                #               ^.......^
                #             end           end

                # last (primary) | curr 
                # ---------------|+++++++
                #               ^.......^
                #             end           end
                steps.append( [ [last.end, curr.end], 'INVIN' ] )
                gap = abs(curr.end - last.end)
            else:
                if (curr.start  < last.end):
                    steps.append( [ [last.end,curr.start], 'DUP' ] )

                # last   | curr (primary)
                # +++++++|-------------
                # ^.......^
                # start   start

                # last   | curr (primary)
                # -------|+++++++++++++++
                # ^.......^
                # start   start
                steps.append( [ [last.start, curr.start], 'INVOUT' ] )
                gap = abs(curr.start - last.start)
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

def get_long_read_max_gap(read_name, long_reads):
    """Finds the largest gap between alignments in LongRead alignments, plus lengths of the  alignments 

    Returns the integer max gap
    """
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

def plot_variant(start, end, sv_type, ax, range_min, range_max):
    """Plots the variant bar at the top of the image

    """
    r=[float(int(start) - range_min)/float(range_max - range_min), \
        float(int(end) - range_min)/float(range_max - range_min)]
    ax.plot(r,[0,0],'-',color='black',lw=8,solid_capstyle="butt",alpha=0.5)
    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
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

def plot_confidence_interval(breakpoint,ci, ax, range_min, range_max):
    """Plots a confidence interval on the variant bar
    """
    r=[float(int(breakpoint)-int(ci[0]) - range_min)/float(range_max - range_min), \
        float(int(breakpoint)+int(ci[1]) - range_min)/float(range_max - range_min)]
    
    
    ax.plot(r,[0,0],'-',color='black',lw=.5, alpha=1)
    ax.axvline(r[0], color='black', lw=0.5,alpha=1, ymin=0.40, ymax=0.60)
    ax.axvline(r[1], color='black', lw=0.5,alpha=1, ymin=0.40, ymax=0.60)

    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x',length=0)
    ax.tick_params(axis='y',length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

def get_pair_event_type(pe_read):
    """Decide what type of event the read supports (del/normal, dup, inv)
    """
    event_by_strand = {
        (True, False): 'Deletion/Normal',
        (False, True): 'Duplication',
        (False, False): 'Inversion',
        (True, True): 'Inversion'
    }
    event_type = event_by_strand[pe_read[0].strand,pe_read[1].strand]
    return event_type

def points_in_window(points):
    """Checks whether these points lie within the window of interest

    Points is a list of one start, one end coordinate (ints)
    """
    if points[0] < -5 or points[1] < -5 or points[0] > 5 or points[1] > 5:
        return False
    return True


def plot_pair(pair, y, ax, range_min, range_max):
    """Plots a PairedEnd read at the y-position corresponding to insert size

    If read lies outside the range-min or range_max, it is not plotted
    """
    
    if pair[0].end < range_min or pair[1].start > range_max:
        return

    p = [float(pair[0].start - range_min)/float(range_max - range_min), \
         float(pair[1].end - range_min)/float(range_max - range_min)]

    # some points are far outside of the printable area, so we ignore them 
    if not points_in_window(p):
        return

    event_type = get_pair_event_type(pair)
    READ_TYPES_USED[event_type] = True
    color = COLORS[event_type]

    # plot the individual pair
    ax.plot(p,\
            [y,y],\
            '-',color=color, \
            alpha=0.25, \
            lw=0.5, \
            marker='s', \
            markersize=marker_size, zorder=10)

def plot_pairs(pairs,
               ax,
               range_min,
               range_max,
               curr_min_insert_size,
               curr_max_insert_size):
    """Plots all PairedEnd reads for the region
    """

    ## set the color of the pair based on the two strands
    for read_name in pairs:
        pair = pairs[read_name]

        if len(pair) != 2: continue

        if pair[0].MI or pair[1].MI: continue

        # y value is the insert size
        insert_size = pair[1].end - pair[0].start
        # use this to scale the y-axis
        points = [float(pair[0].start - range_min)/float(range_max - range_min), \
         float(pair[1].end - range_min)/float(range_max - range_min)]
        if points_in_window(points):
            if not curr_min_insert_size or curr_min_insert_size > insert_size:
                curr_min_insert_size = insert_size
            if not curr_max_insert_size or curr_max_insert_size < insert_size:
                curr_max_insert_size = insert_size

        plot_pair(pair, insert_size, ax, range_min, range_max)

    return [curr_min_insert_size, curr_max_insert_size]

def plot_linked_reads(pairs,
                      splits,
                      linked_reads,
                      ax,
                      range_min,
                      range_max,
                      curr_min_insert_size,
                      curr_max_insert_size):
    """Plots all LinkedReads for the region
    """

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

            if linked_pair[1].end > range_max or linked_pair[0].start < range_min:
                continue
            points = [float(linked_pair[0].start - range_min)/float(range_max - range_min), 
                    float(linked_pair[1].end - range_min)/float(range_max - range_min)]
            if points_in_window(points):
                gap_sizes.append(abs(linked_pair[1].end - linked_pair[0].start))

        for linked_split in linked_splits:
            alignments.append([linked_split[0].end,linked_split[1].start])
            if linked_split[1].start > range_max or linked_split[0].end < range_min:
                continue
            points = [float(linked_split[0].start - range_min)/float(range_max - range_min), 
                    float(linked_split[1].end - range_min)/float(range_max - range_min)]
            if points_in_window(points):
                gap_sizes.append(abs(linked_split[1].start - linked_split[0].end))

        if len(gap_sizes) == 0 : continue

        insert_size = max(gap_sizes)
        

        if not curr_min_insert_size or curr_min_insert_size > insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size < insert_size:
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
        
        #ignore points outside window
        if not points_in_window(p):
            continue

        ax.plot(p,
                [insert_size,insert_size],\
                '-',color='green', \
                alpha=0.75, \
                lw=0.25)

        for name in linked_reads[linked_read][0]:
            if name in pairs and len(pairs[name]) == 2:
                pair = pairs[name]
                plot_pair(pair, insert_size, ax, range_min, range_max)

        for name in linked_reads[linked_read][1]:
            if name in splits:
                split = splits[name]
                plot_split(split, insert_size, ax, range_min, range_max)

    return [curr_min_insert_size, curr_max_insert_size]

def get_split_event_type(split):
    """Decide what type of event the read supports (del/normal, dup, inv)
    """

    first = split[0]
    second = split[1]
    if first.start > second.end:
        second = split[0]
        first = split[1]

    # first.strand, second.strand, first.query<second.query,first.start<second.start
    event_type_by_strand_and_order = {
        (True, False)               : 'Inversion',       #mixed strands
        (False, True)               : 'Inversion',       #mixed strands 
        (True, True, True)          : 'Deletion/Normal', #forward strand
        (True, True, False)         : 'Duplication',     #forward strand
        (False, False, False, False): 'Deletion/Normal', #reverse strand
        (False, False, False, True) : 'Duplication',     #reverse strand
        (False, False, True, True)  : 'Deletion/Normal', #reverse strand
        (False, False, True, False)  : 'Duplication'     #reverse strand
    }
    orientations = [first.strand, second.strand]
    
    #if same strand, need query position info
    if orientations[0] == orientations[1]:
        #first query position smaller than second query position, normal for forward strand
        orientations.append(first.query_pos < second.query_pos)
        
        #reverse strand requires start position info
        if False in orientations[:2]:
            #first start smaller than second start, normal for forward strand
            orientations.append(first.start < second.start)
    event_type = event_type_by_strand_and_order[tuple(orientations)]
    return event_type
            

def plot_split(split, y, ax, range_min, range_max):
    """Plots a SplitRead at the y-position corresponding to insert size

    If read lies outside the range-min or range_max, it is not plotted
    """
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

    if not points_in_window(p):
        return
    event_type = get_split_event_type(split)
    color = COLORS[event_type]
    
    ax.plot(p,\
            [y,y],\
            ':',\
            color=color, \
            alpha=0.25, \
            lw=1, \
            marker='o',
            markersize=marker_size)

def plot_splits(splits,
                ax,
                range_min,
                range_max,
                curr_min_insert_size,
                curr_max_insert_size):
    """Plots all SplitReads for the region
    """

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
        points = [float(start.end - range_min)/float(range_max - range_min), \
                float(end.start - range_min)/float(range_max - range_min)]
        if points_in_window(points):
            if not curr_min_insert_size or curr_min_insert_size > insert_size:
                curr_min_insert_size = insert_size
            if not curr_max_insert_size or curr_max_insert_size < insert_size:
                curr_max_insert_size = insert_size

        plot_split(split, insert_size, ax, range_min, range_max)
    return [curr_min_insert_size, curr_max_insert_size]

def plot_long_reads(long_reads,
                    ax,
                    range_min,
                    range_max,
                    curr_min_insert_size,
                    curr_max_insert_size):
    """Plots all LongReads for the region
    """

    Path = mpath.Path

    colors = { 'ALIGN' : 'orange',
               'DEL' : 'black',
               'INVIN' : 'blue',
               'INVOUT' : 'blue',
               'DUP' : 'red' }

    for read_name in long_reads:
        long_read_plan = get_long_read_plan(read_name,
                                            long_reads,
                                            range_min,
                                            range_max)
        if long_read_plan is None:
            continue
        max_gap = long_read_plan[0]
        steps = long_read_plan[1]
        for step in steps:
            step_cords = step[0]
            step_type = step[1]

            p = [float(step_cords[0]-range_min)/float(range_max - range_min),
                 float(step_cords[1]-range_min)/float(range_max - range_min)]

            # some points are far outside of the printable area, so we ignore them 
            if not points_in_window(p):
                continue

            if step_type == 'ALIGN':
                ax.plot(p,
                        [max_gap,max_gap],
                        '-', 
                        color=colors[step_type],
                        alpha=0.25,
                        lw=1)

                if max_gap > curr_max_insert_size:
                    curr_max_insert_size = max_gap
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

                # all some room for the bend line
                if max_gap*1.1 > curr_max_insert_size:
                    curr_max_insert_size = max_gap*1.1

    return [curr_min_insert_size, curr_max_insert_size]

def plot_coverage(coverage, 
                ax, 
                range_min, 
                range_max, 
                hp_count,
                max_coverage, 
                tracktype,
                yaxis_label_fontsize):
    """Plots high and low quality coverage for the region

    User may specify a preference between stacked and superimposed 
        superimposed may cause unexpected behavior if low-quality depth is greater than high 
    """
   
    cover_x = []
    cover_y_lowqual = []
    cover_y_highqual = []
    cover_y_all = []

    for pos in range(range_min,range_max+1):
        if pos in coverage:
            cover_x.append(\
                    float(pos-range_min)/float(range_max - range_min))
            cover_y_all.append(coverage[pos][0] + coverage[pos][1])
            cover_y_highqual.append(coverage[pos][0])
            cover_y_lowqual.append(coverage[pos][1])
        else:
            cover_x.append(\
                    float(pos-range_min)/float(range_max - range_min))
            cover_y_lowqual.append(0)
            cover_y_highqual.append(0)
            cover_y_all.append(0)
    cover_y_lowqual = np.array(cover_y_lowqual)
    cover_y_highqual = np.array(cover_y_highqual)
    cover_y_all = np.array(cover_y_all)

    if max_coverage > 0:
        max_plot_depth = max_coverage
    elif cover_y_all.max() > 3 * cover_y_all.mean():
        max_plot_depth = max(np.percentile(cover_y_all, 99.5), np.percentile(cover_y_all, 99.5))
    else:
        max_plot_depth = max(cover_y_all.max(), cover_y_all.max())
    ax2 = ax.twinx()
    ax2.set_xlim([0,1])
    
    if 0 == max_plot_depth:
        max_plot_depth = 0.01

    ax2.set_ylim([0,max_plot_depth])
    bottom_fill = np.zeros(len(cover_y_all))
    if tracktype == "stack":
        ax2.fill_between(cover_x, \
                         cover_y_highqual, \
                         bottom_fill,\
                         color='darkgrey',
                         step="pre",
                         alpha=.4)

        ax2.fill_between(cover_x, \
                         cover_y_all, \
                         cover_y_highqual,
                         color='grey',
                         step="pre",
                         alpha=0.15)

        
    elif tracktype == "superimpose": 
        ax2.fill_between(cover_x, \
                         cover_y_lowqual, \
                         bottom_fill,\
                         color='grey',
                         step="pre",
                         alpha=.15)


        ax2.fill_between(cover_x, \
                         cover_y_highqual, \
                         cover_y_lowqual,\
                         color='darkgrey',
                         step="pre",
                         alpha=.4)

        ax2.fill_between(cover_x, \
                         cover_y_lowqual, \
                         bottom_fill,
                         color='grey',
                         step="pre",
                         alpha=0.15)
       
    #number of ticks should be 6 if there's one hp, 3 otherwise
    tick_count = 5 if hp_count==1 else 2
    tick_count = max(int(max_plot_depth/tick_count), 1)

    # set axis parameters
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(tick_count))
    ax2.tick_params(axis='y', colors='grey', labelsize=yaxis_label_fontsize)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis='x',length=0)
    ax2.tick_params(axis='y',length=0)

    return ax2

def get_pair_insert_sizes(pairs):
    """Extracts the integer insert sizes for all pairs

    Return list of integer insert sizes
    """
    pair_insert_sizes = []

    for hp in pairs:
        for read_name in pairs[hp]:
            if len(pairs[hp][read_name]) == 2:
                first = pairs[hp][read_name][0]
                second = pairs[hp][read_name][1]
                pair_insert_sizes.append(second.end - first.start)
    return pair_insert_sizes

def get_split_insert_size(splits):
    """Extracts the integer gap sizes for all split reads

    Return list of integer gap sizes
    """
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

def get_long_read_gap_sizes(long_reads):
    """Extracts the integer gap sizes for all long reads

    Return list of integer gap sizes
    """
    long_read_gap_sizes = [] 

    for hp in long_reads:
        for read_name in long_reads[hp]:
            long_read_gap_sizes.append(\
                    get_long_read_max_gap(read_name, long_reads[hp]))
    return long_read_gap_sizes

def pair(arg):
    """Defines behavior for ArgParse pairs 

    Pairs must be comma-separated list of two items
    """
    try:
        parsed_arg = [int(x) for x in arg.split(',')]
        if len(parsed_arg) == 2:
            return parsed_arg
        else:
            parser.error('Invalid number of pair values')
    except:
        parser.error('Invalid pair values')

def print_arguments(options):
    """Prints out the arguments to samplot as a json object

    Used as metadata for PlotCritic
    """
    if options.print_args or options.json_only:
        import json
        args_filename = os.path.splitext(options.output_file)[0] + ".json"
        args_info = {
            'titles': options.titles if options.titles else 'None',
            'reference': options.reference if options.reference else 'None',
            'bams': options.bams,
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
 

def setup_arguments():
    """Defines the allowed arguments for samplot
    """
    parser = argparse.ArgumentParser(
            description="SAMPLOT creates images of genome regions from CRAM/SAM alignments, "+\
                    "optimized for structural variant call review")
    parser.add_argument("--marker_size",
                      type=int,
                      default=3,
                      help="Size of marks on pairs and splits (default 3) ",
                      required=False)

    parser.add_argument("-n",
                      "--titles",
                      help="Space-delimited list of plot titles. "+\
                              "Use quote marks to include spaces (i.e. \"plot 1\" \"plot 2\")",
                      type=str,
                      nargs="+",
                      required=False);

    parser.add_argument("-r",
                      "--reference",
                      help="Reference file for CRAM, required if CRAM files used",
                      type=str,
                      required=False);

    parser.add_argument("-z",
                      "--z",
                      type=int,
                      default=4,
                      help="Number of stdevs from the mean (default 4)",
                      required=False)

    parser.add_argument("-b",
                      "--bams",
                      type=str,
                      nargs="+",
                      help="Space-delimited list of BAM/CRAM file names",
                      required=True)

    parser.add_argument("-o",
                      "--output_file",
                      type=str,
                      help="Output file name",
                      required=True)

    parser.add_argument("-s",
                      "--start",
                      type=int,
                      help="Start position of region/variant",
                      required=True)

    parser.add_argument("-e",
                      "--end",
                      type=int,
                      help="End position of region/variant",
                      required=True)

    parser.add_argument("-c",
                      "--chrom",
                      type=str,
                      help="Chromosome",
                      required=True)

    parser.add_argument("-w",
                      "--window",
                      type=int,
                      help="Window size (count of bases to include in view), default(0.5 * len)",
                      required=False)

    parser.add_argument("-d",
                      "--max_depth",
                      type=int,
                      help="Max number of normal pairs to plot",
                      default=100,
                      required=False)

    parser.add_argument("--minq",
                      type=int,
                      help="coverage from reads with MAPQ <= minq plotted in lighter grey."+\
                              " To disable, pass in negative value",
                      default=0,
                      required=False)

    parser.add_argument("-t",
                      "--sv_type",
                      type=str,
                      help="SV type. If omitted, plot is created without variant bar",
                      required=False)

    parser.add_argument("-T",
                      "--transcript_file",
                      help="GFF of transcripts",
                      required=False)

    parser.add_argument("-A",
                      "--annotation_files",
                      type=str,
                      nargs="+",
                      help="Space-delimited list of bed.gz tabixed files of annotations "+\
                              "(such as repeats, mappability, etc.)",
                      required=False)

    parser.add_argument("--coverage_tracktype",
                      type=str,
                      help="type of track to use for low MAPQ coverage plot.",
                      choices=['stack','superimpose'],
                      default="stack",
                      required=False)

    parser.add_argument("-a",
                      "--print_args",
                      action="store_true",
                      default=False,
                      help="Print commandline arguments",
                      required=False)

    parser.add_argument("-H",
                      "--plot_height",
                      type=int,
                      help="Plot height",
                      required=False)

    parser.add_argument("-W",
                      "--plot_width",
                      type=int,
                      help="Plot width",
                      required=False)

    parser.add_argument("-q",
                      "--min_mqual",
                      type=int,
                      help="Min mapping quality of reads to be included in plot",
                      required=False)

    parser.add_argument("-j",
                      "--json_only",
                      action="store_true",
                      default=False,
                      help="Create only the json file, not the image plot",
                      required=False)

    parser.add_argument("--start_ci",
                      help="confidence intervals of SV first breakpoint "+\
                              "(distance from the breakpoint). Must be a comma-separated pair of ints (i.e. 20,40)",
                      type=pair,
                      required=False)

    parser.add_argument("--end_ci",
                      help="confidence intervals of SV end breakpoint "+\
                              "(distance from the breakpoint). Must be a comma-separated pair of ints (i.e. 20,40)",
                      type=pair,
                      required=False)

    parser.add_argument("--long_read",
                      type=int,
                      default=1000,
                      help="Min length of a read to be treated as a long-read (default 1000)",
                      required=False)

    parser.add_argument("--min_event_size",
                      type=int,
                      default=100,
                      help="Min size of an event in long-read CIGAR to include (default 100)",
                      required=False)

    parser.add_argument("--xaxis_label_fontsize",
                      type=int,
                      default=6,
                      help="Font size for X-axis labels (default 6)",
                      required=False)

    parser.add_argument("--yaxis_label_fontsize",
                      type=int,
                      default=6,
                      help="Font size for Y-axis labels (default 6)",
                      required=False)

    parser.add_argument("--legend_fontsize",
                      type=int,
                      default=6,
                      help="Font size for legend labels (default 6)",
                      required=False)

    parser.add_argument("--annotation_fontsize",
                      type=int,
                      default=6,
                      help="Font size for annotation labels (default 6)",
                      required=False)

    parser.add_argument("--common_insert_size",
                      action="store_true",
                      default=False,
                      help="Set common insert size for all plots",
                      required=False)

    parser.add_argument("--hide_annotation_labels",
                      action="store_true",
                      default=False,
                      help="Hide the label (fourth column text) "+\
                              "from annotation files, useful for region with many annotations",
                      required=False)

    parser.add_argument("--coverage_only",
                      action="store_true",
                      default=False,
                      help="Hide all reads and show only coverage",
                      required=False)

    parser.add_argument("--same_yaxis_scales",
                      action="store_true",
                      default=False,
                      help="Set the scales of the Y axes to the max of all",
                      required=False)
    options = parser.parse_args()
    
    if options.print_args or options.json_only:
        print_arguments(options)
        if options.json_only:
            sys.exit(0)

    for bam in options.bams:
        if ".cram" in bam:
            if not options.reference:
                parser.print_help(sys.stderr)
                sys.exit("Error: Missing reference for CRAM")
    return options

def set_plot_dimensions(start, end, sv_type, arg_plot_height, arg_plot_width, 
        bams, annotation_files, transcript_file, arg_window):
    """Chooses appropriate dimensions for the plot

    Includes the number of samples, whether a variant type is included, and any annotations in height
    Includes the start, end, and window argument in width
    If height and width are chosen by used, these are used instead

    Return plot height, width, and window as integers
    """

    plot_height = 5
    plot_width = 8
    if arg_plot_height:
        plot_height = arg_plot_height
    else:
        num_subplots = len(bams)
        if annotation_files:
            num_subplots += .3*len(annotation_files)
        if transcript_file:
            num_subplots += .3
        plot_height = 2 + num_subplots

    if arg_plot_width:
        plot_width = arg_plot_width

    # if an SV type is given, then expand the window around its bounds
    window = 0
    if arg_window:
        window = arg_window
    elif sv_type:
        window = int((int(end) - int(start))/2)

    return plot_height,plot_width,window

def get_read_data(chrom, start, end, bams, reference, min_mqual, coverage_only, 
        long_read_length, same_yaxis_scales, max_depth, z_score):
    """Reads alignment files to extract reads for the region

    Region and alignment files given with chrom, start, end, bams
    If CRAM files are used, reference must be provided
    Reads with mapping quality below min_mqual will not be retrieved
    If coverage_only, reads are not kept and used only for checking coverage
    Reads longer than long_read_length will be treated as long reads
    Max coverages values will be set to same value for all samples if same_yaxis_scales
    If max_depth, only max_depth reads will be retrieved, although all will be included in coverage
    If PairedEnd read insert size is greater than z_score standard deviations from mean, read will be treated as discordant
    """
    
    all_pairs = []
    all_splits = []
    all_coverages = []
    all_long_reads = []
    all_linked_reads = []

    range_min = max(0,int(start) - window)
    range_max = int(end) + window
    min_insert_size = None
    max_insert_size = None
    max_coverage = 0

    for bam_file_name in bams:
        bam_file = None
        if not reference:
            bam_file = pysam.AlignmentFile(bam_file_name, "rb")
        else:
            bam_file = pysam.AlignmentFile(bam_file_name, \
                                           "rc", \
                                           reference_filename=reference)

        pairs = {}
        splits = {}
        long_reads = {}
        coverage = {}
        linked_reads = {}
        for read in bam_file.fetch(chrom,
                                   max(0,range_min-1000), 
                                   range_max+1000):
            if min_mqual and int(read.mapping_quality) < min_mqual:
                continue
            
            if not coverage_only:
                if read.query_length >= long_read_length:
                    add_long_reads(read, 
                        long_reads, 
                        range_min, 
                        range_max, 
                        options.min_event_size)
                else:
                    add_pair_end(read, pairs, linked_reads)
                    add_split(read, splits, bam_file, linked_reads)
            add_coverage(read, coverage, options.minq)

        pair_insert_sizes = get_pair_insert_sizes(pairs)
        split_insert_sizes = get_split_insert_size(splits)
        long_read_gap_sizes = get_long_read_gap_sizes(long_reads)

        insert_sizes = pair_insert_sizes + \
                       split_insert_sizes + \
                       long_read_gap_sizes
        if not insert_sizes or len(insert_sizes) == 0:
            if not coverage_only:
                print('Warning: No data returned from fetch in region  ' + \
                        chrom + ':' + str(start) + '-' + \
                        str(end) + \
                        ' from ' + bam_file_name, file=sys.stderr)
            insert_sizes.append(0)

        if not min_insert_size:
            min_insert_size = min(insert_sizes)
        else: 
            min_insert_size = min(min(insert_sizes), min_insert_size)

        if not max_insert_size:
            max_insert_size = max(insert_sizes)
        else: 
            max_insert_size = max(max(insert_sizes), min_insert_size)

        if same_yaxis_scales:
            for i in coverage:
                curr_max = max(max(coverage[i].values()))
                if curr_max > max_coverage:
                    max_coverage = curr_max
        all_coverages.append(coverage)
        all_pairs.append(pairs)
        all_splits.append(splits)
        all_long_reads.append(long_reads)
        all_linked_reads.append(linked_reads)
    
    read_data = {
            "all_pairs" :       all_pairs,
            "all_splits":       all_splits,
            "all_coverages":    all_coverages,
            "all_long_reads":   all_long_reads,
            "all_linked_reads": all_linked_reads
        }

    # Sample +/- pairs in the normal insert size range
    if max_depth:
        read_data['all_pairs'] = downsample_pairs(max_depth, z_score, read_data['all_pairs'])
    
    return read_data,max_coverage

def downsample_pairs(max_depth, z_score, all_pairs):
    """Downsamples to keep only max_depth normal pairs from all PairedEnd reads 
    """
    for i in range(len(all_pairs)):
        for hp in all_pairs[i]:
            all_pairs[i][hp] = sample_normal(max_depth,all_pairs[i][hp],z_score)
    return all_pairs

def set_haplotypes(curr_coverage):
    """Creates a list to manage counting haplotypes for subplots
    """
    hps = sorted(curr_coverage.keys(), reverse=True)
    #if there are multiple haplotypes, must have 0,1,2
    if len(hps) > 1 or (len(hps) == 1 and hps[0] != 0):
        if 0 not in hps:
            hps.append(0)
        if 1 not in hps:
            hps.append(1)
        if 2 not in hps:
            hps.append(2)
    elif 0 not in hps:
        hps.append(0)
    hps.sort(reverse=True)
    return hps


def plot_samples(read_data, 
        grid, 
        ax_i, 
        number_of_axes, 
        bams, 
        chrom, 
        coverage_tracktype, 
        titles, 
        same_yaxis_scales,
        xaxis_label_fontsize, 
        yaxis_label_fontsize, 
        annotation_files, 
        transcript_file,
        max_coverage):
    """Plots all samples
    """
    max_insert_size = 0
    for i in range(len(bams)):
        ax =  matplotlib.pyplot.subplot(grid[ax_i])
        hps = set_haplotypes(read_data['all_coverages'][i])
        inner_axs = gridspec.GridSpecFromSubplotSpec(len(hps), 
                                                     1,
                                                     subplot_spec=grid[ax_i],
                                                     wspace=0.0,
                                                     hspace=0.5)
        axs = {}
        for j in range(len(hps)):
            axs[j] = matplotlib.pyplot.subplot(inner_axs[hps[j]])
            
        curr_min_insert_size = None
        curr_max_insert_size = None

        cover_axs = {}
        curr_axs = ''
        overall_insert_size_range = []
        overall_coverage_range = []
        for hp in hps:
            curr_ax = axs[hp]

            curr_pairs = []
            if hp in read_data['all_pairs'][i]:
                curr_pairs = read_data['all_pairs'][i][hp]

            curr_splits = []
            if hp in read_data['all_splits'][i]:
                curr_splits = read_data['all_splits'][i][hp]

            curr_linked_reads = []
            if hp in read_data['all_linked_reads'][i]:
                curr_linked_reads = read_data['all_linked_reads'][i][hp]

            curr_long_reads = []
            if hp in read_data['all_long_reads'][i]:
                curr_long_reads = read_data['all_long_reads'][i][hp]

            curr_coverage = []
            if hp in read_data['all_coverages'][i]:
                curr_coverage = read_data['all_coverages'][i][hp]

            cover_ax = plot_coverage(curr_coverage, 
                    curr_ax, 
                    range_min, 
                    range_max,
                    len(hps), 
                    max_coverage, 
                    coverage_tracktype, 
                    yaxis_label_fontsize)
            
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

            cover_axs[hp] = cover_ax
            if curr_max_insert_size > max_insert_size:
                max_insert_size = curr_max_insert_size

        #{{{ set axis parameters
        #set the axis title to be either one passed in or filename
        curr_ax = axs[hps[0]]

        if titles and \
                len(titles) == len(bams):
            curr_ax.set_title(titles[ax_i-1], \
                         fontsize=8, loc='left')
        else:
            curr_ax.set_title(os.path.basename(bams[ax_i-1]), \
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
            if same_yaxis_scales:
                curr_ax.set_ylim([0,max_insert_size])
            else:
                curr_ax.set_ylim([0,curr_max_insert_size])
            curr_ax.spines['top'].set_visible(False)
            curr_ax.spines['bottom'].set_visible(False)
            curr_ax.spines['left'].set_visible(False)
            curr_ax.spines['right'].set_visible(False)
            curr_ax.tick_params(axis='y', labelsize=yaxis_label_fontsize)
            #if there's one hp, 6 ticks fit. Otherwise, do 3
            tick_count = 6 if len(hps)==1 else 3
            curr_ax.yaxis.set_major_locator(ticker.LinearLocator(tick_count))
            curr_ax.tick_params(axis='both', length=0)
            curr_ax.set_xticklabels([])
        last_sample_num = number_of_axes - 1
        if annotation_files:
            last_sample_num -= len(annotation_files)
        if transcript_file:
            last_sample_num -= 1
        
        if (ax_i == last_sample_num):
            curr_ax = axs[ hps[-1] ]
            labels = [int(range_min + l*(range_max-range_min)) \
                    for l in curr_ax.xaxis.get_majorticklocs()]
            curr_ax.set_xticklabels(labels, fontsize=xaxis_label_fontsize)
            curr_ax.set_xlabel('Chromosomal position on ' + chrom, fontsize=8)

        curr_ax = axs[hps[ int(len(hps)/2)    ]]
        curr_ax.set_ylabel('Insert size', fontsize=8)
        cover_ax = cover_axs[hps[ int(len(hps)/2)    ]]
        cover_ax.set_ylabel('Coverage', fontsize=8)
        #}}}

        ax_i += 1
    return ax_i


def plot_legend(fig, legend_fontsize):
    """Plots the figure legend
    """
    marker_colors = []
    marker_labels = []
    read_colors = {
        "Deletion/Normal":'black', 
        "Duplication":'red', 
        "Inversion":'blue', 
        "Aligned long read":'orange', 
        "Linked read":'green'
    }

    for read_type in READ_TYPES_USED:
        if read_type in read_colors:
            color = read_colors[read_type]
            flag = READ_TYPES_USED[read_type]
            if flag:
                marker_colors.append(color)
                marker_labels.append(read_type)
    legend_elements = []

    for color in marker_colors:
        legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1], \
                color=color,
                linestyle='-',
                lw=1)]
    if READ_TYPES_USED["Split-read"]:
        marker_labels.append("Split read")
        legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1], \
                    markerfacecolor="None",
                    markeredgecolor='grey',
                    color='grey',
                    marker='o', \
                    markersize=marker_size,
                    linestyle=':',
                    lw=1)]
    
    if READ_TYPES_USED["Paired-end read"] \
            or READ_TYPES_USED["Deletion/Normal"] \
            or READ_TYPES_USED["Inversion"]:
        marker_labels.append("Paired-end read")
        legend_elements += [matplotlib.pyplot.Line2D([0,0],[0,1], \
                    markerfacecolor="None",
                    markeredgecolor='grey',
                    color='grey',
                    marker='s', \
                    markersize=marker_size,
                    linestyle='-',
                    lw=1)]

    fig.legend( legend_elements ,
                marker_labels, 
                loc=1,
                fontsize = legend_fontsize,
                frameon=False)

def get_tabix_iter(chrom, pos, end, datafile):
    """Gets an iterator from a tabix BED/GFF/GFF3 file

    Used to avoid chrX vs. X notation issues when extracting data from  annotation files
    """
    tbx = pysam.TabixFile(datafile)
    itr = None
    try:
        itr = tbx.fetch(chrom, max(0,range_min-1000), range_max+1000)
    except ValueError:
        # try and account for chr/no chr prefix
        if chrom[:3] == 'chr':
            chrom = chrom[3:]
        else:
            chrom = 'chr' + chrom

        try:
            itr = tbx.fetch(chrom,
                            max(0,range_min-1000), 
                            range_max+1000)
        except ValueError:
            sys.exit('Warning: Could not fetch ' + \
                    chrom + ':' + pos + '-' + end + \
                    ' from ' + datafile)
    return itr
        


def plot_annotations(annotation_files, chrom, start, end, 
        hide_annotation_labels, annotation_fontsize, grid, ax_i):
    """Plots annotation information from region 
    """
    
    for annotation_file in annotation_files:
        itr = get_tabix_iter(chrom,start, end, annotation_file)
        ax =  matplotlib.pyplot.subplot(grid[ax_i])
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
                    ax.plot(r,[0,0],'-',color=str(v),lw=15)
                except ValueError:
                    ax.plot(r,[0,0],'-',color='black',lw=15)
                    if not hide_annotation_labels:
                        ax.text(r[0],
                                0 + 0.1,
                                A[3],
                                color='black', 
                                fontsize=annotation_fontsize)
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

        ax.tick_params(axis='x',length=0)
        ax.tick_params(axis='y',length=0)
        ax.set_xticklabels([])
        ax.set_yticklabels([])

def plot_transcript(transcript_file, chrom, start, end, 
        grid, annotation_fontsize, xaxis_label_fontsize):
    """Plots a transcript file annotation
    """
    genes = {}
    transcripts = {}
    cdss = {}

    # fetch and parse data from the tabixed gff file
    itr = get_tabix_iter(chrom, start, end, transcript_file)
    
    for row in itr:
        gene_annotation = row.split()

        if gene_annotation[2] == 'gene':
            info =  dict([list(val.split('=')) for val in gene_annotation[8].split(';')])
            info['strand'] = gene_annotation[6] == "+"

            genes[info['Name']] = [
                gene_annotation[0],
                int(gene_annotation[3]), 
                int(gene_annotation[4]),
                info
            ]

        elif gene_annotation[2] == 'transcript':
            info =  dict([list(val.split('=')) for val in gene_annotation[8].split(';')])
            info['strand'] = gene_annotation[6] == "+"
            
            if info['Parent'] not in transcripts:
                transcripts[info['Parent']] = {}
            transcripts[info['Parent']][info['ID']] = [
                gene_annotation[0],
                int(gene_annotation[3]), 
                int(gene_annotation[4]),
                info
            ]

        elif gene_annotation[2] == 'CDS':
            info =  dict([list(val.split('=')) for val in gene_annotation[8].split(';')])
            info['strand'] = gene_annotation[6] == "+"
            
            if info['Parent'] not in cdss:
                cdss[info['Parent']] = {}

            if info['ID'] not in cdss[info['Parent']]:
                cdss[info['Parent']][info['ID']] = []

            cdss[info['Parent']][info['ID']].append([gene_annotation[0], \
                                                     int(gene_annotation[3]), \
                                                     int(gene_annotation[4]), \
                                                     info])
    ax =  matplotlib.pyplot.subplot(grid[-1])

    transcript_idx = 0
    arrow_loc = 0.02
    for gene in genes:
        gene_id = genes[gene][3]['ID']
        if gene_id not in transcripts: continue
        for transcript in transcripts[gene_id]:
            t_start = max(range_min, transcripts[gene_id][transcript][1])
            t_end = min(range_max, transcripts[gene_id][transcript][2])
            r=[float(t_start - range_min)/float(range_max - range_min), \
               float(t_end - range_min)/float(range_max - range_min)]
            
            ax.plot(r,[transcript_idx,transcript_idx],'-',color='cornflowerblue',lw=0.5)

            ax.text(r[0],
                    transcript_idx + 0.02,
                    gene,
                    color='blue', 
                    fontsize=annotation_fontsize)

            if genes[gene][3]['strand']:
                ax.annotate("",
                    xy=(1,transcript_idx), 
                    xytext=(1-arrow_loc, transcript_idx),
                    arrowprops=dict(arrowstyle="->",color="cornflowerblue", lw=1), 
                    annotation_clip=True
                )
            else:
                ax.annotate("",
                    xy=(1-arrow_loc,transcript_idx), 
                    xytext=(1, transcript_idx),
                    arrowprops=dict(arrowstyle="->",color="cornflowerblue", lw=1), 
                    annotation_clip=True
                )
        
            if transcript in cdss:
                for cds in cdss[transcript]:
                    for exon in cdss[transcript][cds]:
                        e_start = max(range_min,exon[1])
                        e_end = min(range_max,exon[2])

                        r=[float(e_start - range_min)/float(range_max - range_min), \
                            float(e_end - range_min)/float(range_max - range_min)]
                        ax.plot(r,[transcript_idx,transcript_idx],'-',color='cornflowerblue',lw=4)
                transcript_idx += 1
        
    # set axis parameters
    ax.set_xlim([0,1])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='x',length=0)
    ax.tick_params(axis='y',length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(os.path.basename(transcript_file), fontsize=8, loc='left')

def create_variant_plot(grid, 
        ax_i, 
        start, 
        end, 
        sv_type, 
        range_min, 
        range_max, 
        start_ci, 
        end_ci):
    """Plots the pieces of the variant bar at the top, including bar and confidence intervals
    """
    ax = matplotlib.pyplot.subplot(grid[ax_i])
    plot_variant(start, end, sv_type, ax, range_min, range_max)
    ax_i += 1
    #plot confidence intervals if provided
    if start_ci and start_ci != None:
        plot_confidence_interval(start, start_ci, ax, range_min, range_max)
    if end_ci and end_ci != None:
        plot_confidence_interval(end, end_ci, ax, range_min, range_max)
    return ax_i

def create_gridspec(bams, transcript_file, annotation_files, sv_type ):
    """Helper function for creation of a correctly-sized GridSpec instance
    """
    # give one axis to display each sample
    num_ax = len(bams)

    # add another if we are displaying the SV
    if sv_type:
        num_ax+=1

    # add another if a annotation file is given
    if transcript_file:
        num_ax+=1

    if annotation_files:
        num_ax+=len(options.annotation_files)

    # set the relative sizes for each
    ratios = []
    if sv_type:
        ratios = [1] 
    
    for i in range(len(bams)):
        ratios.append( len(read_data['all_coverages'][i]) * 3 )
        if len(read_data['all_coverages']) > 0:
            ratios[-1] = 9
    
    if annotation_files:
        ratios += [.3]*len(annotation_files)
    if transcript_file:
        ratios.append(2)
    return gridspec.GridSpec(num_ax, 1, height_ratios = ratios), num_ax

########################################################################
# main block
########################################################################
if __name__ == '__main__':
    options = setup_arguments()
    
    # set up plot 
    plot_height,plot_width,window = set_plot_dimensions(options.start, 
            options.end, 
            options.sv_type, 
            options.plot_height, 
            options.plot_width, 
            options.bams, 
            options.annotation_files, 
            options.transcript_file, 
            options.window)

    marker_size = options.marker_size
    range_min = max(0,int(options.start) - window)
    range_max = int(options.end) + window

    # set up sub plots
    matplotlib.rcParams.update({'font.size': 12})
    fig = matplotlib.pyplot.figure(figsize=(plot_width, plot_height), dpi=300)

    # read alignment data
    read_data,max_coverage = get_read_data(options.chrom, 
            options.start, 
            options.end, 
            options.bams, 
            options.reference, 
            options.min_mqual, 
            options.coverage_only, 
            options.long_read, 
            options.same_yaxis_scales, 
            options.max_depth, 
            options.z)
    
    # set up grid organizer
    grid,num_ax = create_gridspec(options.bams, 
            options.transcript_file, 
            options.annotation_files, 
            options.sv_type )
    current_axis_idx = 0
    
    # plot variant on top
    if options.sv_type:
        current_axis_idx = create_variant_plot(grid, 
            current_axis_idx, 
            options.start, 
            options.end, 
            options.sv_type, 
            range_min, 
            range_max, 
            options.start_ci, 
            options.end_ci)
    
    # Plot each sample
    current_axis_idx = plot_samples(read_data, 
        grid, 
        current_axis_idx, 
        num_ax, 
        options.bams, 
        options.chrom,  
        options.coverage_tracktype, 
        options.titles, 
        options.same_yaxis_scales, 
        options.xaxis_label_fontsize, 
        options.yaxis_label_fontsize, 
        options.annotation_files, 
        options.transcript_file,
        max_coverage)   

    # plot legend
    plot_legend(fig, options.legend_fontsize)

    # Plot annotation files
    if options.annotation_files:
        plot_annotations(options.annotation_files, 
            options.chrom, 
            options.start, 
            options.end, 
            options.hide_annotation_labels, 
            options.annotation_fontsize, 
            grid, 
            current_axis_idx)

    # Plot sorted/bgziped/tabixed transcript file
    if options.transcript_file:
        plot_transcript(options.transcript_file, 
            options.chrom, 
            options.start, 
            options.end, 
            grid, 
            options.annotation_fontsize, 
            options.xaxis_label_fontsize)
    
    # save
    matplotlib.rcParams['agg.path.chunksize'] = 100000
    matplotlib.pyplot.tight_layout(pad=0.8,h_pad=.1, w_pad=.1)
    try:
        matplotlib.pyplot.savefig(options.output_file)
    except:
        print ("Failed to save figure " + options.output_file + ". Region may be too large")
    matplotlib.pyplot.savefig(options.output_file)
