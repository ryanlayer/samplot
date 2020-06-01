#!/usr/bin/env python
from __future__ import print_function

import os
import random
import re
import sys

import matplotlib
matplotlib.use("Agg") #must be before imports of submodules in matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pysam
from matplotlib.offsetbox import AnchoredText



INTERCHROM_YAXIS = 5000

COLORS = {
    "Deletion/Normal": "black",
    "Deletion": "black",
    "Duplication": "red",
    "Inversion": "blue",
    "InterChrmInversion": "blue",
    "InterChrm": "black",
}

READ_TYPES_USED = {
    "Deletion/Normal": False,
    "Duplication": False,
    "Inversion": False,
    "Aligned long read": False,
    "Linked read": False,
    "Split-read": False,
    "Paired-end read": False,
}

# pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
CIGAR_MAP = {
    "M": 0,
    "I": 1,
    "D": 2,
    "N": 3,
    "S": 4,
    "H": 5,
    "P": 6,
    "=": 7,
    "X": 8,
    "B": 9,
}

# {{{class plan_step:
class plan_step:
    step_events = ["Align", "ANNOTATION"]

    def __init__(self, start_pos, end_pos, event, info=None):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.event = event
        self.info = info

    def __str__(self):
        if self.info:
            return (
                "Step("
                + str(self.start_pos)
                + ", "
                + str(self.end_pos)
                + ", "
                + self.event
                + ", "
                + str(self.info)
                + ")"
            )
        else:
            return (
                "Step("
                + str(self.start_pos)
                + ", "
                + str(self.end_pos)
                + ", "
                + self.event
                + ")"
            )

    def __repr__(self):
        return str(self)


# }}}

# {{{class genome_interval:
class genome_interval:
    def __init__(self, chrm, start, end):
        self.chrm = chrm
        self.start = start
        self.end = end

    def __str__(self):
        return "(" + self.chrm + "," + str(self.start) + "," + str(self.end) + ")"

    def __repr__(self):
        return str(self)

    def __eq__(self, gi2):
        return self.chrm == gi2.chrm and self.start == gi2.start and self.end == gi2.end

    """ return -1 if before, 0 if in, 1 if after """

    def intersect(self, gi):
        if gi.chrm.strip("chr") < self.chrm.strip("chr") or gi.end < self.start:
            return -1
        elif gi.chrm.strip("chr") > self.chrm.strip("chr") or gi.start > self.end:
            return 1
        else:
            return 0


# }}}

# {{{def get_range_hit(ranges, chrm, point):
def get_range_hit(ranges, chrm, point):
    for j in range(len(ranges)):
        r = ranges[j]
        if (
            r.chrm.strip("chr") == chrm.strip("chr")
            and r.start <= point
            and r.end >= point
        ):
            return j
    return None


# }}}

# {{{def map_genome_point_to_range_points(ranges, chrm, point):
def map_genome_point_to_range_points(ranges, chrm, point):
    range_hit = get_range_hit(ranges, chrm, point)

    if range_hit == None:
        return None
    p = 1.0 / len(ranges) * range_hit + (1.0 / len(ranges)) * (
        float(point - ranges[range_hit].start)
        / float(ranges[range_hit].end - ranges[range_hit].start)
    )

    return p


# }}}

# {{{def points_in_window(points):
def points_in_window(points):
    """Checks whether these points lie within the window of interest

    Points is a list of one start, one end coordinate (ints)
    """
    if (
        None in points
        or points[0] < -5
        or points[1] < -5
        or points[0] > 5
        or points[1] > 5
    ):
        return False
    return True


# }}}

# {{{ def get_tabix_iter(chrm, start, end, datafile):
def get_tabix_iter(chrm, start, end, datafile):
    """Gets an iterator from a tabix BED/GFF/GFF3 file

    Used to avoid chrX vs. X notation issues when extracting data from
    annotation files
    """
    tbx = pysam.TabixFile(datafile)
    itr = None
    try:
        itr = tbx.fetch(chrm, max(0, start - 1000), end + 1000)
    except ValueError:
        # try and account for chr/no chr prefix
        if chrm[:3] == "chr":
            chrm = chrm[3:]
        else:
            chrm = "chr" + chrm

        try:
            itr = tbx.fetch(chrm, max(0, start - 1000), end + 1000)
        except ValueError as e:
            print(
                "Warning: Could not fetch "
                + chrm
                + ":"
                + str(start)
                + "-"
                + str(end)
                + " from "
                + datafile,
                file=sys.stderr,
            )
            print(e)
    return itr


# }}}

##Coverage methods
# {{{def add_coverage(bam_file, read, coverage, separate_mqual):
def add_coverage(bam_file, read, coverage, separate_mqual, ignore_hp):
    """Adds a read to the known coverage 
    
    Coverage from Pysam read is added to the coverage list
    Coverage list is a pair of high- and low-quality lists
    Quality is determined by separate_mqual, which is min quality
    """

    chrm = bam_file.get_reference_name(read.reference_id).strip("chr")

    hp = 0

    if not ignore_hp and read.has_tag("HP"):
        hp = int(read.get_tag("HP"))

    if hp not in coverage:
        coverage[hp] = {}
    if chrm not in coverage[hp]:
        coverage[hp][chrm] = {}

    curr_pos = read.reference_start
    if not read.cigartuples:
        return

    for op, length in read.cigartuples:
        if op in [CIGAR_MAP["M"], CIGAR_MAP["="], CIGAR_MAP["X"]]:
            for pos in range(curr_pos, curr_pos + length):
                if pos not in coverage[hp][chrm]:
                    coverage[hp][chrm][pos] = [0, 0]

                # the two coverage tracks are [0] high-quality and [1]
                # low-quality
                if separate_mqual and (read.mapping_quality > separate_mqual):
                    coverage[hp][chrm][pos][0] += 1
                else:
                    coverage[hp][chrm][pos][1] += 1
            curr_pos += length
        elif op == CIGAR_MAP["I"]:
            curr_pos = curr_pos
        elif op == CIGAR_MAP["D"]:
            curr_pos += length
        elif op == CIGAR_MAP["N"]:
            curr_pos = length
        elif op == CIGAR_MAP["S"]:
            curr_pos = curr_pos
        elif op == CIGAR_MAP["H"]:
            curr_pos = curr_pos
        else:
            curr_pos += length


# }}}

# {{{def plot_coverage(coverage,
def plot_coverage(
    coverage,
    ax,
    ranges,
    hp_count,
    max_coverage,
    tracktype,
    yaxis_label_fontsize,
    same_yaxis_labels=False,
):
    """Plots high and low quality coverage for the region

    User may specify a preference between stacked and superimposed 
    superimposed may cause unexpected behavior if low-quality depth is
    greater than high 
    """

    cover_x = []
    cover_y_lowqual = []
    cover_y_highqual = []
    cover_y_all = []

    for i in range(len(ranges)):
        r = ranges[i]
        for pos in range(r.start, r.end + 1):
            cover_x.append(map_genome_point_to_range_points(ranges, r.chrm, pos))
            if r.chrm in coverage and pos in coverage[r.chrm]:
                cover_y_all.append(coverage[r.chrm][pos][0] + coverage[r.chrm][pos][1])
                cover_y_highqual.append(coverage[r.chrm][pos][0])
                cover_y_lowqual.append(coverage[r.chrm][pos][1])
            else:
                cover_y_lowqual.append(0)
                cover_y_highqual.append(0)
                cover_y_all.append(0)
    cover_y_lowqual = np.array(cover_y_lowqual)
    cover_y_highqual = np.array(cover_y_highqual)
    cover_y_all = np.array(cover_y_all)

    if max_coverage > 0 and same_yaxis_labels:
        max_plot_depth = max_coverage
    elif cover_y_all.max() > 3 * cover_y_all.mean():
        max_plot_depth = max(
            np.percentile(cover_y_all, 99.5), np.percentile(cover_y_all, 99.5)
        )
    else:
        max_plot_depth = np.percentile(cover_y_all.max(), 99.5)
    ax2 = ax.twinx()
    ax2.set_xlim([0, 1])

    if 0 == max_plot_depth:
        max_plot_depth = 0.01

    ax2.set_ylim([0, max(1, max_plot_depth)])
    bottom_fill = np.zeros(len(cover_y_all))
    if tracktype == "stack":
        ax2.fill_between(
            cover_x,
            cover_y_highqual,
            bottom_fill,
            color="darkgrey",
            step="pre",
            alpha=0.4,
        )

        ax2.fill_between(
            cover_x, cover_y_all, cover_y_highqual, color="grey", step="pre", alpha=0.15
        )

    elif tracktype == "superimpose":
        ax2.fill_between(
            cover_x, cover_y_lowqual, bottom_fill, color="grey", step="pre", alpha=0.15
        )

        ax2.fill_between(
            cover_x,
            cover_y_highqual,
            cover_y_lowqual,
            color="darkgrey",
            step="pre",
            alpha=0.4,
        )

        ax2.fill_between(
            cover_x, cover_y_lowqual, bottom_fill, color="grey", step="pre", alpha=0.15
        )

    # number of ticks should be 6 if there's one hp, 3 otherwise
    tick_count = 5 if hp_count == 1 else 2
    tick_count = max(int(max_plot_depth / tick_count), 1)

    # set axis parameters
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(tick_count))
    ax2.tick_params(axis="y", colors="grey", labelsize=yaxis_label_fontsize)
    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.tick_params(axis="x", length=0)
    ax2.tick_params(axis="y", length=0)

    # break the variant plot when we have multiple ranges
    for i in range(1, len(ranges)):
        ax2.axvline(x=1.0 / len(ranges), color="white", linewidth=5)

    return ax2


# }}}

##Pair End methods
# {{{class PairedEnd:
class PairedEnd:
    """container of paired-end read info

    Contains start(int), end(int), strand(bool True=forward), MI (int
    molecular identifier), HP (int haplotype)
    """

    def __init__(self, chrm, start, end, is_reverse, MI_tag, HP_tag):
        """Create PairedEnd instance

        Genomic interval is defined by start and end integers
        Strand is opposite of is_reverse
        Molecular identifier and Haplotype are integers if present, else
        False
        """
        self.pos = genome_interval(chrm, start, end)
        self.strand = not (is_reverse)
        # molecular identifier - linked reads only
        self.MI = None
        # haplotype - phased reads only
        self.HP = 0

        if MI_tag:
            self.MI = MI_tag
        if HP_tag:
            self.HP = HP_tag

    def __repr__(self):
        return "PairedEnd(%s,%s,%s,%s,%s,%s)" % (
            self.pos.chrm,
            self.pos.start,
            self.pos.end,
            self.strand,
            self.MI,
            self.HP,
        )


# }}}

# {{{ def add_pair_end(bam_file, read, pairs, linked_reads):
def add_pair_end(bam_file, read, pairs, linked_reads, ignore_hp):
    """adds a (mapped, primary, non-supplementary, and paired) read to the
    pairs list

    Pysam read is added as simpified PairedEnd instance to pairs
    Also added to linked_reads list if there is an associated MI tag
    """

    if read.is_unmapped:
        return
    if not (read.is_paired):
        return
    if read.is_secondary:
        return
    if read.is_supplementary:
        return

    MI_tag = False
    HP_tag = False

    if read.has_tag("MI"):
        MI_tag = int(read.get_tag("MI"))
    if not ignore_hp and read.has_tag("HP"):
        HP_tag = int(read.get_tag("HP"))

    READ_TYPES_USED["Paired-end read"] = True

    pe = PairedEnd(
        bam_file.get_reference_name(read.reference_id),
        read.reference_start,
        read.reference_end,
        read.is_reverse,
        MI_tag,
        HP_tag,
    )

    if pe.HP not in pairs:
        pairs[pe.HP] = {}

    if read.query_name not in pairs[pe.HP]:
        pairs[pe.HP][read.query_name] = []

    if pe.MI:
        READ_TYPES_USED["Linked read"] = True
        if pe.HP not in linked_reads:
            linked_reads[pe.HP] = {}

        if pe.MI not in linked_reads[pe.HP]:
            linked_reads[pe.HP][pe.MI] = [[], []]
        linked_reads[pe.HP][pe.MI][0].append(read.query_name)

    pairs[pe.HP][read.query_name].append(pe)
    pairs[pe.HP][read.query_name].sort(key=lambda x: x.pos.start)


# }}}

# {{{def sample_normal(max_depth, pairs, z):
def sample_normal(max_depth, pairs, z):
    """Downsamples paired-end reads 
    
    Selects max_depth reads
    Does not remove discordant pairs, those with insert distance greater
    than z stdevs from mean

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
        lens = np.array(
            [pair[1].pos.end - pair[0].pos.start for pair in plus_minus_pairs.values()]
        )
        mean = np.mean(lens)
        stdev = np.std(lens)

        inside_norm = {}

        for read_name in pairs:
            pair = pairs[read_name]
            if len(pair) != 2:
                continue
            if pair[1].pos.end - pair[0].pos.start >= mean + z * stdev:
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


# }}}

# {{{def get_pairs_insert_sizes(pairs):
def get_pairs_insert_sizes(ranges, pairs):
    """Extracts the integer insert sizes for all pairs

    Return list of integer insert sizes
    """
    pair_insert_sizes = []

    for hp in pairs:
        for read_name in pairs[hp]:
            if len(pairs[hp][read_name]) == 2:
                size = get_pair_insert_size(ranges, pairs[hp][read_name])

                if size:
                    pair_insert_sizes.append(size)

    return pair_insert_sizes


# }}}

# {{{def get_pair_insert_size(ranges, pair):
def get_pair_insert_size(ranges, pair):
    """ Gives the outer distance
    """
    first = pair[0]
    second = pair[1]

    # make sure both sides are in range
    if (
        get_range_hit(ranges, first.pos.chrm, first.pos.start) != None
        or get_range_hit(ranges, first.pos.chrm, first.pos.end) != None
    ) and (
        get_range_hit(ranges, second.pos.chrm, second.pos.start) != None
        or get_range_hit(ranges, second.pos.chrm, second.pos.end) != None
    ):

        if first.pos.chrm == second.pos.chrm:
            return abs(second.pos.end - first.pos.start)
        else:
            return INTERCHROM_YAXIS
    else:
        return None


# }}}

# {{{ def get_pairs_plan(ranges, pairs, linked_plan=False):
def get_pairs_plan(ranges, pairs, linked_plan=False):
    steps = []
    max_event = 0

    insert_sizes = []

    for read_name in pairs:
        pair = pairs[read_name]

        plan = get_pair_plan(ranges, pair)

        if plan:
            insert_size, step = plan
            insert_sizes.append(insert_size)
            steps.append(step)

    if len(insert_sizes) > 0:
        max_event = max(insert_sizes)

    plan = [max_event, steps]

    return plan


# }}}

# {{{def get_pair_plan(ranges, pair, linked_plan=False):
def get_pair_plan(ranges, pair, linked_plan=False):
    if pair == None or len(pair) != 2:
        return None

    first = pair[0]
    second = pair[1]

    # see if they are part of a linked read
    if not linked_plan and (first.MI or second.MI):
        return None

    # make sure both ends are in the plotted region
    first_s_hit = get_range_hit(ranges, first.pos.chrm, first.pos.start)
    first_e_hit = get_range_hit(ranges, first.pos.chrm, first.pos.end)
    second_s_hit = get_range_hit(ranges, second.pos.chrm, second.pos.start)
    second_e_hit = get_range_hit(ranges, second.pos.chrm, second.pos.end)

    if (first_s_hit == None and first_e_hit == None) or (
        second_s_hit == None and second_e_hit == None
    ):
        return None

    insert_size = get_pair_insert_size(ranges, pair)

    first_hit = first_s_hit if first_s_hit != None else first_e_hit
    second_hit = second_e_hit if second_e_hit != None else second_s_hit

    start = genome_interval(
        first.pos.chrm,
        max(first.pos.start, ranges[first_hit].start),
        max(first.pos.start, ranges[first_hit].start),
    )

    end = genome_interval(
        second.pos.chrm,
        min(second.pos.end, ranges[second_hit].end),
        min(second.pos.end, ranges[second_hit].end),
    )

    step = plan_step(start, end, "PAIREND")

    event_type = get_pair_event_type(pair)
    step.info = {"TYPE": event_type, "INSERTSIZE": insert_size}

    return insert_size, step


# }}}

# {{{def get_pair_event_type(pe_read):
def get_pair_event_type(pe_read):
    """Decide what type of event the read supports (del/normal, dup, inv)
    """
    event_by_strand = {
        (True, False): "Deletion/Normal",
        (False, True): "Duplication",
        (False, False): "Inversion",
        (True, True): "Inversion",
    }
    event_type = event_by_strand[pe_read[0].strand, pe_read[1].strand]
    return event_type


# }}}

# {{{def plot_pair_plan(ranges, step, ax):
def plot_pair_plan(ranges, step, ax, marker_size):
    p = [
        map_genome_point_to_range_points(
            ranges, step.start_pos.chrm, step.start_pos.start
        ),
        map_genome_point_to_range_points(ranges, step.end_pos.chrm, step.end_pos.end),
    ]

    if None in p:
        return False

    # some points are far outside of the printable area, so we ignore them
    if not points_in_window(p):
        return False

    y = step.info["INSERTSIZE"]
    event_type = step.info["TYPE"]
    READ_TYPES_USED[event_type] = True
    color = COLORS[event_type]

    # plot the individual pair
    ax.plot(
        p,
        [y, y],
        "-",
        color=color,
        alpha=0.25,
        lw=0.5,
        marker="s",
        markersize=marker_size,
        zorder=10,
    )

    return True


# }}}

# {{{def plot_pairs(pairs,
def plot_pairs(
    pairs, ax, ranges, curr_min_insert_size, curr_max_insert_size, marker_size
):
    """Plots all PairedEnd reads for the region
    """

    plan = get_pairs_plan(ranges, pairs)

    if not plan:
        [curr_min_insert_size, curr_max_insert_size]

    max_event, steps = plan

    for step in steps:
        plot_pair_plan(ranges, step, ax, marker_size)

    if not curr_min_insert_size or curr_min_insert_size > max_event:
        curr_min_insert_size = max_event
    if not curr_max_insert_size or curr_max_insert_size < max_event:
        curr_max_insert_size = max_event

    return [curr_min_insert_size, curr_max_insert_size]


# }}}

##Split Read methods
# {{{class SplitRead:
class SplitRead:
    """container of split read info

    Contains start(int), end(int), strand(bool True=forward), query
    position (int), MI (int molecular identifier), HP (int haplotype)
    """

    def __init__(self, chrm, start, end, strand, query_pos, MI_tag=None, HP_tag=None):
        """Create SplitRead instance

        Genomic interval is defined by start, end, and query_pos integers
        Strand is opposite of is_reverse
        Molecular identifier and Haplotype are integers if present, else
        False
        """
        self.pos = genome_interval(chrm, start, end)
        self.strand = strand
        self.query_pos = query_pos
        # molecular identifier - linked reads only
        self.MI = None
        # haplotype - phased reads only
        self.HP = 0

        if MI_tag:
            self.MI = MI_tag
        if HP_tag:
            self.HP = HP_tag

    def __repr__(self):
        return "SplitRead(%s,%s,%s,%s,%s,%s,%s)" % (
            self.pos.chrm,
            self.pos.start,
            self.pos.end,
            self.strand,
            self.query_pos,
            self.MI,
            self.HP,
        )


# }}}

# {{{def calc_query_pos_from_cigar(cigar, strand):
def calc_query_pos_from_cigar(cigar, strand):
    """Uses the CIGAR string to determine the query position of a read

    The cigar arg is a string like the following: 86M65S
    The strand arg is a boolean, True for forward strand and False for
    reverse

    Returns pair of ints for query start, end positions
    """

    cigar_ops = [[int(op[0]), op[1]] for op in re.findall("(\d+)([A-Za-z])", cigar)]

    order_ops = cigar_ops
    if not strand:  # - strand
        order_ops = order_ops[::-1]

    qs_pos = 0
    qe_pos = 0
    q_len = 0

    for op_position in range(len(cigar_ops)):
        op_len = cigar_ops[op_position][0]
        op_type = cigar_ops[op_position][1]

        if op_position == 0 and (op_type == "H" or op_type == "S"):
            qs_pos += op_len
            qe_pos += op_len
            q_len += op_len
        elif op_type == "H" or op_type == "S":
            q_len += op_len
        elif op_type == "M" or op_type == "I" or op_type == "X":
            qe_pos += op_len
            q_len += op_len

    return qs_pos, qe_pos


# }}}

# {{{def add_split(read, splits, bam_file, linked_reads):
def add_split(read, splits, bam_file, linked_reads, ignore_hp):
    """adds a (primary, non-supplementary) read to the splits list

    Pysam read is added as simpified SplitRead instance to splits
    Also added to linked_reads list if there is an associated MI tag
    """
    if read.is_secondary:
        return
    if read.is_supplementary:
        return
    if not read.has_tag("SA"):
        return

    READ_TYPES_USED["Split-read"] = True
    qs_pos, qe_pos = calc_query_pos_from_cigar(read.cigarstring, (not read.is_reverse))

    HP_tag = False
    MI_tag = False
    if read.has_tag("MI"):
        MI_tag = int(read.get_tag("MI"))

    if not ignore_hp and read.has_tag("HP"):
        HP_tag = int(read.get_tag("HP"))
    sr = SplitRead(
        bam_file.get_reference_name(read.reference_id),
        read.reference_start,
        read.reference_end,
        not (read.is_reverse),
        qs_pos,
        MI_tag,
        HP_tag,
    )

    if sr.MI:
        if sr.HP not in linked_reads:
            linked_reads[sr.HP] = {}
        if sr.MI not in linked_reads[sr.HP]:
            linked_reads[sr.HP][sr.MI] = [[], []]
        linked_reads[sr.HP][sr.MI][1].append(read.query_name)

    if sr.HP not in splits:
        splits[sr.HP] = {}

    splits[sr.HP][read.query_name] = [sr]

    for sa in read.get_tag("SA").split(";"):
        if len(sa) == 0:
            continue
        A = sa.split(",")
        chrm = A[0]
        pos = int(A[1])
        strand = A[2] == "+"
        cigar = A[3]
        #mapq and nm are never used, annotating this for code readability 
        mapq = int(A[4])
        nm = int(A[5])
        qs_pos, qe_pos = calc_query_pos_from_cigar(cigar, strand)
        splits[sr.HP][read.query_name].append(
            SplitRead(chrm, pos, pos + qe_pos, strand, qs_pos)
        )

    if len(splits[sr.HP][read.query_name]) == 1:
        del splits[sr.HP][read.query_name]
    else:
        splits[sr.HP][read.query_name].sort(key=lambda x: x.pos.start)


# }}}

# {{{def get_split_event_type(split):
def get_split_event_type(split):
    """Decide what type of event the read supports (del/normal, dup, inv)
    """

    first = split[0]
    second = split[1]
    if first.pos.start > second.pos.end or first.pos.chrm > second.pos.chrm:
        second = split[0]
        first = split[1]

    # first.strand, second.strand,
    # first.query<second.query,first.start<second.start
    event_type_by_strand_and_order = {
        (True, False): "Inversion",  # mixed strands
        (False, True): "Inversion",  # mixed strands
        (True, True, True): "Deletion/Normal",  # forward strand
        (True, True, False): "Duplication",  # forward strand
        (False, False, False, False): "Deletion/Normal",  # reverse strand
        (False, False, False, True): "Duplication",  # reverse strand
        (False, False, True, True): "Deletion/Normal",  # reverse strand
        (False, False, True, False): "Duplication",  # reverse strand
    }
    orientations = [first.strand, second.strand]

    # if same strand, need query position info
    if orientations[0] == orientations[1]:
        # first query position smaller than second query position,
        # normal for forward strand
        orientations.append(first.query_pos < second.query_pos)

        # reverse strand requires start position info
        if False in orientations[:2]:
            # first start smaller than second start, normal for forward strand
            orientations.append(first.pos.start < second.pos.start)
    event_type = event_type_by_strand_and_order[tuple(orientations)]

    if first.pos.chrm != second.pos.chrm:
        if event_type == "Inversion":
            event_type = "InterChrmInversion"
        else:
            event_type = "InterChrm"

    return event_type


# }}}

# {{{def get_split_insert_sizes(splits):
def get_splits_insert_sizes(ranges, splits):
    """Extracts the integer gap sizes for all split reads
    Return list of integer gap sizes
    """
    split_insert_sizes = []

    for hp in splits:
        for read_name in splits[hp]:
            sizes = get_split_insert_sizes(ranges, splits[hp][read_name])
            if sizes:
                split_insert_sizes += sizes
    return split_insert_sizes


# }}}

# {{{def get_split_insert_sizes(ranges, srs):
def get_split_insert_sizes(ranges, srs):
    split_insert_sizes = []

    srs = sorted(srs, key=lambda x: x.query_pos)
    # only keep the alignments that are in a range
    srs = [
        sr
        for sr in srs
        if get_range_hit(ranges, sr.pos.chrm, sr.pos.start) != None
        or get_range_hit(ranges, sr.pos.chrm, sr.pos.end) != None
    ]

    if len(srs) < 2:
        return None

    last = srs[0]

    for i in range(1, len(srs)):
        curr = srs[i]

        # INTERCHROM
        if last.pos.chrm != curr.pos.chrm:
            split_insert_sizes.append(INTERCHROM_YAXIS)
        # Inversion
        elif last.strand != curr.strand:
            split_insert_sizes.append(abs(curr.pos.end - last.pos.end))
        # Duplication
        elif curr.pos.start < last.pos.end:
            split_insert_sizes.append(abs(last.pos.end - curr.pos.start))
        # Deletion
        elif curr.pos.start > last.pos.end:
            split_insert_sizes.append(abs(curr.pos.start - last.pos.end))
        else:
            sys.stderr.write(
                "WARNING: Could not classify event:" + str(last) + str(curr)
            )
        last = curr

    return split_insert_sizes


# }}}

# {{{def get_split_plan(ranges, split):
def get_split_plan(ranges, split, linked_plan=False):
    """
    There can be 2 or more alignments in a split. Plot only those that are in a
    range, and set the insert size to be the largest gap

    A split read acts like a long read, so we will covert the split read
    to a long read, then convert the long read plan back to a split read plan
    """

    alignments = []
    for s in split:
        # see if they are part of a linked read
        if not linked_plan and (s.MI):
            return None
        alignment = Alignment(s.pos.chrm, s.pos.start, s.pos.end, s.strand, s.query_pos)
        alignments.append(alignment)

    long_read = LongRead(alignments)
    long_reads = {}
    long_reads["convert"] = [long_read]
    plan = get_long_read_plan("convert", long_reads, ranges)

    if not plan:
        return None

    max_gap, lr_steps = plan

    if len(lr_steps) < 3:
        return None

    sr_steps = []

    # a split read will include 3 long read steps, align, event, align
    for i in range(0, len(lr_steps), 2):
        if i + 2 > len(lr_steps):
            break
        if (
            lr_steps[i].info["TYPE"] == "Align"
            and lr_steps[i + 1].info["TYPE"] != "Align"
            and lr_steps[i + 2].info["TYPE"] == "Align"
        ):
            start = genome_interval(
                lr_steps[i].end_pos.chrm,
                lr_steps[i].end_pos.end,
                lr_steps[i].end_pos.end,
            )
            end = genome_interval(
                lr_steps[i + 2].start_pos.chrm,
                lr_steps[i + 2].start_pos.start,
                lr_steps[i + 2].start_pos.start,
            )
            sr_steps.append(
                plan_step(
                    start,
                    end,
                    "SPLITREAD",
                    info={"TYPE": lr_steps[i + 1].info["TYPE"], "INSERTSIZE": max_gap},
                )
            )
    return max_gap, sr_steps


# }}}

# {{{def get_splits_plan(ranges, splits, linked_plan=False):
def get_splits_plan(ranges, splits, linked_plan=False):
    steps = []
    max_event = 0

    insert_sizes = []

    for read_name in splits:
        split = splits[read_name]

        plan = get_split_plan(ranges, split)

        if plan:
            insert_size, step = plan
            insert_sizes.append(insert_size)
            steps += step

    if len(insert_sizes) > 0:
        max_event = max(insert_sizes)

    plan = [max_event, steps]

    return plan


# }}}

# {{{def plot_split(split, y, ax, ranges):
def plot_split(split, y, ax, ranges, marker_size):
    """Plots a SplitRead at the y-position corresponding to insert size

    If read lies outside the range-min or range_max, it is not plotted
    """
    start = split[0]
    end = split[1]
    if start.pos.chrm > end.pos.chrm or start.pos.start > end.pos.end:
        end = split[0]
        start = split[1]

    p = [
        map_genome_point_to_range_points(ranges, start.pos.chrm, start.pos.end),
        map_genome_point_to_range_points(ranges, end.pos.chrm, end.pos.start),
    ]

    if not points_in_window(p):
        return
    event_type = get_split_event_type(split)
    color = COLORS[event_type]

    ax.plot(
        p,
        [y, y],
        ":",
        color=color,
        alpha=0.25,
        lw=1,
        marker="o",
        markersize=marker_size,
    )


# }}}

# {{{def plot_split(split, y, ax, ranges):
def plot_split_plan(ranges, step, ax, marker_size):
    p = [
        map_genome_point_to_range_points(
            ranges, step.start_pos.chrm, step.start_pos.start
        ),
        map_genome_point_to_range_points(ranges, step.end_pos.chrm, step.end_pos.end),
    ]

    if None in p:
        return False

    # some points are far outside of the printable area, so we ignore them
    if not points_in_window(p):
        return False

    y = step.info["INSERTSIZE"]
    event_type = step.info["TYPE"]
    READ_TYPES_USED[event_type] = True
    color = COLORS[event_type]

    ax.plot(
        p,
        [y, y],
        ":",
        color=color,
        alpha=0.25,
        lw=1,
        marker="o",
        markersize=marker_size,
    )


# }}}

# {{{def plot_splits(splits,
def plot_splits(
    splits, ax, ranges, curr_min_insert_size, curr_max_insert_size, marker_size
):
    """Plots all SplitReads for the region
    """
    plan = get_splits_plan(ranges, splits)

    if not plan:
        [curr_min_insert_size, curr_max_insert_size]

    max_event, steps = plan

    for step in steps:
        plot_split_plan(ranges, step, ax, marker_size)

    if not curr_min_insert_size or curr_min_insert_size > max_event:
        curr_min_insert_size = max_event
    if not curr_max_insert_size or curr_max_insert_size < max_event:
        curr_max_insert_size = max_event

    return [curr_min_insert_size, curr_max_insert_size]


# }}}

##Long Read methods
# {{{class Alignment:
class Alignment:
    """container of alignment info, from CIGAR string

    Contains start(int), end(int), strand(bool True=forward), query
    position (int)
    """

    def __init__(self, chrm, start, end, strand, query_position):
        """Create Alignment instance

        Genomic interval is defined by start, end, and query_pos integers
        Strand is bool (True for forward)
        """
        self.pos = genome_interval(chrm, start, end)
        self.strand = strand
        self.query_position = query_position

    def __str__(self):
        return ",".join(
            [
                str(x)
                for x in [
                    self.pos.chrm,
                    self.pos.start,
                    self.pos.end,
                    self.strand,
                    self.query_position,
                ]
            ]
        )

    def __repr__(self):
        return "Alignment(%s,%s,%s,%s,%s)" % (
            self.pos.chrm,
            self.pos.start,
            self.pos.end,
            self.strand,
            self.query_position,
        )


# }}}

# {{{class LongRead:
class LongRead:
    """container of LongRead info

    Contains start(int), end(int), list of Alignments
    """

    def __init__(self, alignments):
        """Create LongRead instance

        Genomic interval is defined by start, end integers
        List of Alignments set by parameter
        """
        self.alignments = alignments

    def __str__(self):
        return ",".join([str(x) for x in self.alignments])

    def __repr__(self):
        return "LongRead(" + str(self) + ")"


# }}}

# {{{def get_alignments_from_cigar(chrm,
def get_alignments_from_cigar(chrm, curr_pos, strand, cigartuples, reverse=False):
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

    for op, length in cigartuples:
        if op in [CIGAR_MAP["M"], CIGAR_MAP["="], CIGAR_MAP["X"]]:
            alignments.append(
                Alignment(chrm, curr_pos, curr_pos + length, strand, q_pos)
            )
            curr_pos += length
            q_pos += length
        elif op == CIGAR_MAP["I"]:
            q_pos += length
        elif op == CIGAR_MAP["D"]:
            curr_pos += length
        elif op == CIGAR_MAP["N"]:
            curr_pos += length
        elif op == CIGAR_MAP["S"]:
            q_pos += length
    return alignments


# }}}

# {{{def get_cigartuples_from_string(cigarstring):
def get_cigartuples_from_string(cigarstring):
    """Extracts operations,lengths as tuples from cigar string"

    Returns list of tuples of [operation,length]
    """
    cigartuples = []
    for match in re.findall(r"(\d+)([A-Z]{1})", cigarstring):
        length = int(match[0])
        op = match[1]
        cigartuples.append((CIGAR_MAP[op], length))

    return cigartuples


# }}}

# {{{def merge_alignments(min_gap, alignments):
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
            if (
                alignment.pos.chrm == merged_alignments[-1].pos.chrm
                and alignment.pos.start < merged_alignments[-1].pos.end + min_gap
            ):
                merged_alignments[-1].pos.end = alignment.pos.end
            else:
                merged_alignments.append(alignment)
    return merged_alignments


# }}}

# {{{def add_long_reads(bam_file, read, long_reads, min_event_size):
def add_long_reads(bam_file, read, long_reads, min_event_size, ignore_hp):
    """Adds a (primary, non-supplementary, long) read to the long_reads list

    Read added to long_reads if within the inteval defined by ranges
    Alignments belonging to the LongRead instance combined if within the
    min_event_size distance apart
    """
    READ_TYPES_USED["Aligned long read"] = True

    if read.is_supplementary or read.is_secondary:
        return

    hp = 0

    if not ignore_hp and read.has_tag("HP"):
        hp = int(read.get_tag("HP"))

    alignments = get_alignments_from_cigar(
        bam_file.get_reference_name(read.reference_id),
        read.pos,
        not read.is_reverse,
        read.cigartuples,
    )

    min_gap = min_event_size
    merged_alignments = merge_alignments(min_gap, alignments)

    read_strand = not read.is_reverse

    if read.has_tag("SA"):
        for sa in read.get_tag("SA").split(";"):
            if len(sa) == 0:
                continue

            rname, pos, strand, cigar, mapq, nm = sa.split(",")

            sa_pos = int(pos)
            sa_strand = strand == "+"
            strand_match = read_strand != sa_strand
            sa_cigartuples = get_cigartuples_from_string(cigar)
            sa_alignments = get_alignments_from_cigar(
                rname, sa_pos, sa_strand, sa_cigartuples, reverse=strand_match
            )

            sa_merged_alignments = merge_alignments(min_gap, sa_alignments)

            if len(sa_merged_alignments) > 0:
                merged_alignments += sa_merged_alignments

    if hp not in long_reads:
        long_reads[hp] = {}

    if read.query_name not in long_reads[hp]:
        long_reads[hp][read.query_name] = []

    long_reads[hp][read.query_name].append(LongRead(merged_alignments))


# }}}

# {{{def add_align_step(alignment, steps, ranges):
def add_align_step(alignment, steps, ranges):
    # alignment can span ranges
    start_range_hit_i = get_range_hit(ranges, alignment.pos.chrm, alignment.pos.start)
    end_range_hit_i = get_range_hit(ranges, alignment.pos.chrm, alignment.pos.end)

    # neither end is in range, add nothing
    if start_range_hit_i == None and end_range_hit_i == None:
        return

    # start is not in range, use end hit
    if start_range_hit_i == None:
        start = genome_interval(
            alignment.pos.chrm,
            max(alignment.pos.start, ranges[end_range_hit_i].start),
            max(alignment.pos.start, ranges[end_range_hit_i].start),
        )
        end = genome_interval(
            alignment.pos.chrm,
            min(alignment.pos.end, ranges[end_range_hit_i].end),
            min(alignment.pos.end, ranges[end_range_hit_i].end),
        )
        steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Align"}))
    # end is not in range, use start hit
    elif end_range_hit_i == None:
        start = genome_interval(
            alignment.pos.chrm,
            max(alignment.pos.start, ranges[start_range_hit_i].start),
            max(alignment.pos.start, ranges[start_range_hit_i].start),
        )
        end = genome_interval(
            alignment.pos.chrm,
            min(alignment.pos.end, ranges[start_range_hit_i].end),
            min(alignment.pos.end, ranges[start_range_hit_i].end),
        )
        steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Align"}))
    # both are in the same range
    elif start_range_hit_i == end_range_hit_i:
        start = genome_interval(
            alignment.pos.chrm,
            max(alignment.pos.start, ranges[start_range_hit_i].start),
            max(alignment.pos.start, ranges[start_range_hit_i].start),
        )
        end = genome_interval(
            alignment.pos.chrm,
            min(alignment.pos.end, ranges[end_range_hit_i].end),
            min(alignment.pos.end, ranges[end_range_hit_i].end),
        )
        steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Align"}))
    # in different ranges
    else:
        start_1 = genome_interval(
            alignment.pos.chrm,
            max(alignment.pos.start, ranges[start_range_hit_i].start),
            max(alignment.pos.start, ranges[start_range_hit_i].start),
        )
        end_1 = genome_interval(
            alignment.pos.chrm,
            ranges[start_range_hit_i].end,
            ranges[start_range_hit_i].end,
        )
        steps.append(plan_step(start_1, end_1, "LONGREAD", info={"TYPE": "Align"}))

        start_2 = genome_interval(
            alignment.pos.chrm,
            ranges[end_range_hit_i].start,
            ranges[end_range_hit_i].start,
        )
        end_2 = genome_interval(
            alignment.pos.chrm,
            min(alignment.pos.end, ranges[end_range_hit_i].end),
            min(alignment.pos.end, ranges[end_range_hit_i].end),
        )
        steps.append(plan_step(start_2, end_2, "LONGREAD", info={"TYPE": "Align"}))


# }}}

# {{{def get_long_read_plan(read_name, long_reads, ranges):
def get_long_read_plan(read_name, long_reads, ranges):
    """Create a plan to render a long read

    Plan consists of the largest event within the read 
        (used to determine the y-axis position of read)
        and the alignment types for plotting each Alignment within 
        LongRead.alignments Align, Duplication, Deletion, Inversion,
        Inversion,
        InterChrmInversion, InterChrm

    Returns plan
    """

    alignments = []

    # only keep alignments that intersect a range
    seen = {}

    if read_name not in long_reads:
        sys.stderr.write("ERROR: Read name " + read_name + " not in list of long reads")
        sys.exit(1)

    for long_read in long_reads[read_name]:
        for alignment in long_read.alignments:
            if alignment.query_position in seen:
                continue
            seen[alignment.query_position] = 1
            # check to see if any part of this alignment overlaps a plot
            # range
            in_range = False
            for r in ranges:
                if r.intersect(alignment.pos) == 0:
                    in_range = True
            if in_range:
                alignments.append(alignment)

    if len(alignments) <= 0:
        return None
    alignments.sort(key=lambda x: x.query_position)

    # we set the primary strand to be the one with the longest alignment
    # this will affect which alignment is inverted. There are clearly edge
    # cases here that we will need to address as we get more examples
    # of inversions
    longest_alignment = 0
    longest_alignment_i = -1
    for i in range(len(alignments)):
        l = alignments[i].pos.end - alignments[i].pos.start
        if longest_alignment < l:
            longest_alignment = l
            longest_alignment_i = i
    primary_strand = alignments[longest_alignment_i].strand

    steps = []
    # long aglinments may spill over the edges, so we will clip that starts
    curr = alignments[0]

    add_align_step(curr, steps, ranges)

    for i in range(1, len(alignments)):
        last = alignments[i - 1]
        curr = alignments[i]

        # figure out what the event is

        # INTER CHROM
        if curr.pos.chrm != last.pos.chrm:
            if curr.strand != last.strand:
                start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)

                end = genome_interval(curr.pos.chrm, curr.pos.end, curr.pos.end)

                info = {"TYPE": "InterChrmInversion"}
                steps.append(plan_step(start, end, "LONGREAD", info=info))
            else:
                start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
                end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
                info = {"TYPE": "InterChrm"}
                steps.append(plan_step(start, end, "LONGREAD", info=info))

            add_align_step(curr, steps, ranges)
        # Inversion
        elif curr.strand != last.strand:
            # it is possible that we have a complex even that
            # is an inverted Duplication
            if curr.pos.start < last.pos.end:
                start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
                end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
                info = {"TYPE": "Deletion"}
                steps.append(plan_step(start, end, "LONGREAD", info=info))
            if curr.strand != primary_strand:
                # last (primary) | curr
                # +++++++++++++++|-------
                #               ^.......^
                #             end           end

                # last (primary) | curr
                # ---------------|+++++++
                #               ^.......^
                #             end           end

                start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
                end = genome_interval(curr.pos.chrm, curr.pos.end, curr.pos.end)
                info = {"TYPE": "Inversion"}
                steps.append(plan_step(start, end, "LONGREAD", info=info))
            else:
                if curr.pos.start < last.pos.end:
                    start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
                    end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
                    info = {"TYPE": "Duplication"}
                    steps.append(plan_step(start, end, "LONGREAD", info=info))

                # last   | curr (primary)
                # +++++++|-------------
                # ^.......^
                # start   start

                # last   | curr (primary)
                # -------|+++++++++++++++
                # ^.......^
                # start   start

                start = genome_interval(last.pos.chrm, last.pos.start, last.pos.start)
                end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
                info = {"TYPE": "Inversion"}
                steps.append(plan_step(start, end, "LONGREAD", info=info))

            add_align_step(curr, steps, ranges)
        # Duplication
        elif curr.pos.start < last.pos.end:
            start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
            end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
            info = {"TYPE": "Duplication"}
            steps.append(plan_step(start, end, "LONGREAD", info=info))
            add_align_step(curr, steps, ranges)
        # Deletion
        else:
            start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
            end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
            info = {"TYPE": "Deletion"}
            # steps.append(plan_step(start, end, 'LONGREAD', info=info))
            steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Deletion"}))
            add_align_step(curr, steps, ranges)

        # if either end is in a range, then add its gap to the list

    max_gap = None

    chrms = set([s.start_pos.chrm for s in steps] + [s.end_pos.chrm for s in steps])

    # set interchrm dist to 5000
    if len(chrms) > 1:
        max_gap = INTERCHROM_YAXIS
    else:
        step_sizes = [
            abs(step.end_pos.end - step.start_pos.start)
            for step in steps
            if step.info["TYPE"] != "Align"
            and get_range_hit(ranges, step.start_pos.chrm, step.start_pos.start) != None
            and get_range_hit(ranges, step.end_pos.chrm, step.end_pos.end) != None
        ]

        max_gap = max(step_sizes) if len(step_sizes) > 0 else 0

    plan = [max_gap, steps]

    return plan


# }}}

# {{{def get_long_read_max_gap(read_name, long_reads):
def get_long_read_max_gap(read_name, long_reads):
    """Finds the largest gap between alignments in LongRead alignments,
    plus lengths of the  alignments 

    Returns the integer max gap
    """
    alignments = []
    for long_read in long_reads[read_name]:
        for alignment in long_read.alignments:
            alignments.append(alignment)
    alignments.sort(key=lambda x: x.query_position)

    # find biggest gap
    max_gap = 0
    for i in range(1, len(alignments)):
        curr_gap = abs(alignments[i].start - alignments[i - 1].end)
        max_gap = max(max_gap, curr_gap)
    return max_gap


# }}}

##Variant methods
# {{{def plot_variant(sv, sv_type, ax, ranges):
def plot_variant(sv, sv_type, ax, ranges):
    """Plots the variant bar at the top of the image

    """

    r = [
        map_genome_point_to_range_points(ranges, sv[0].chrm, sv[0].start),
        map_genome_point_to_range_points(ranges, sv[-1].chrm, sv[-1].end),
    ]

    ax.plot(r, [0, 0], "-", color="black", lw=8, solid_capstyle="butt", alpha=0.5)

    ax.set_xlim([0, 1])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ## make SV title
    sv_title = ""
    if sv[0].chrm == sv[-1].chrm:
        sv_size = float(sv[0].end) - float(sv[0].start)
        sv_size_unit = "bp"

        if sv_size > 1000000:
            sv_size = "{0:0.2f}".format(sv_size / 1000000.0)
            sv_size_unit = "mb"
        elif sv_size > 1000:
            sv_size = "{0:0.2f}".format(sv_size / 1000.0)
            sv_size_unit = "kb"

        sv_title = str(sv_size) + " " + sv_size_unit + " " + sv_type
    else:
        sv_title = sv_type

    ax.set_title(sv_title, fontsize=8)


# }}}

# {{{def plot_confidence_interval(chrm, breakpoint,ci, ax, ranges):
def plot_confidence_interval(chrm, breakpoint, ci, ax, ranges):
    """Plots a confidence interval on the variant bar
    """

    r = [
        map_genome_point_to_range_points(ranges, chrm, breakpoint - int(ci[0])),
        map_genome_point_to_range_points(ranges, chrm, breakpoint + int(ci[1])),
    ]

    ax.plot(r, [0, 0], "-", color="black", lw=0.5, alpha=1)
    ax.axvline(r[0], color="black", lw=0.5, alpha=1, ymin=0.40, ymax=0.60)
    ax.axvline(r[1], color="black", lw=0.5, alpha=1, ymin=0.40, ymax=0.60)

    ax.set_xlim([0, 1])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])


# }}}

# {{{def create_variant_plot(grid,
def create_variant_plot(grid, ax_i, sv, sv_type, ranges, start_ci, end_ci):
    """Plots the pieces of the variant bar at the top, including bar and
    confidence intervals 
    """
    ax = plt.subplot(grid[ax_i])
    plot_variant(sv, sv_type, ax, ranges)
    ax_i += 1
    # plot confidence intervals if provided
    if start_ci and start_ci != None:
        plot_confidence_interval(sv[0].chrm, sv[0].start, start_ci, ax, ranges)
    if end_ci and end_ci != None:
        plot_confidence_interval(sv[-1].chrm, sv[-1].end, end_ci, ax, ranges)

    # break the variant plot when we have multiple ranges
    for i in range(1, len(ranges)):
        ax.axvline(x=1.0 / len(ranges), color="white", linewidth=5)
        ax.text(
            1.0 / len(ranges),
            0,
            "...",
            fontsize=6,
            fontdict=None,
            horizontalalignment="center",
        )

    return ax_i


# }}}

# Linked Reads methods
# {{{ def get_linked_plan(ranges, pairs, splits, linked_reads, gem_name):
def get_linked_plan(ranges, pairs, splits, linked_reads, gem_name):
    insert_sizes = []

    gem_poss = [[] for i in range(len(ranges))]

    linked_pair_steps = []
    # collect all the pairs in a gem
    for name in linked_reads[gem_name][0]:
        if name in pairs and len(pairs[name]) == 2:
            pair = pairs[name]
            plan = get_pair_plan(ranges, pair, linked_plan=True)
            if plan:
                insert_size, step = plan
                insert_sizes.append(insert_size)
                linked_pair_steps.append(step)

    # collect all the splits in a gem
    linked_split_steps = []
    for name in linked_reads[gem_name][1]:
        if name in splits:
            split = splits[name]
            plan = get_split_plan(ranges, split, linked_plan=True)
            if plan:
                insert_size, steps = plan
                insert_sizes.append(insert_size)
                linked_split_steps += steps

    if len(linked_split_steps) == 0 and len(linked_pair_steps) == 0:
        return None

    for step in linked_split_steps + linked_pair_steps:
        poss = [
            (step.start_pos.chrm, step.start_pos.start),
            (step.start_pos.chrm, step.start_pos.end),
            (step.end_pos.chrm, step.end_pos.start),
            (step.end_pos.chrm, step.end_pos.end),
        ]
        for pos in poss:
            hit = get_range_hit(ranges, pos[0], pos[1])
            if hit > -1:
                gem_poss[hit].append(pos[1])

    max_event_size = max(insert_sizes)

    gem_steps = []

    for i in range(len(ranges)):
        if len(gem_poss[i]) == 0:
            continue
        start = genome_interval(ranges[i].chrm, min(gem_poss[i]), min(gem_poss[i]))
        end = genome_interval(ranges[i].chrm, max(gem_poss[i]), max(gem_poss[i]))
        gem_steps.append(plan_step(start, end, "LINKED"))

    # if the gem extends beyond the range, then push the end pos to the
    # end/begining of the range
    if len(gem_steps) > 1:
        gem_steps[0].end_pos.start = ranges[0].end
        gem_steps[0].end_pos.end = ranges[0].end

        gem_steps[1].start_pos.start = ranges[1].start
        gem_steps[1].start_pos.end = ranges[1].start

    info = {
        "INSERTSIZE": max_event_size,
        "PAIR_STEPS": linked_pair_steps,
        "SPLIT_STEPS": linked_split_steps,
    }

    gem_steps[0].info = info

    return max(insert_sizes), gem_steps


# }}}

# {{{ def plot_linked_reads(pairs,
def plot_linked_reads(
    pairs,
    splits,
    linked_reads,
    ax,
    ranges,
    curr_min_insert_size,
    curr_max_insert_size,
    marker_size,
):
    """Plots all LinkedReads for the region
    """
    for linked_read in linked_reads:
        plan = get_linked_plan(ranges, pairs, splits, linked_reads, linked_read)

        if not plan:
            continue

        insert_size, steps = plan

        if not curr_min_insert_size or curr_min_insert_size > insert_size:
            curr_min_insert_size = insert_size
        if not curr_max_insert_size or curr_max_insert_size < insert_size:
            curr_max_insert_size = insert_size

        for step in steps:
            p = [
                map_genome_point_to_range_points(
                    ranges, step.start_pos.chrm, step.start_pos.start
                ),
                map_genome_point_to_range_points(
                    ranges, step.end_pos.chrm, step.end_pos.end
                ),
            ]
            # ignore points outside window
            if not points_in_window(p):
                continue

            ax.plot(
                p, [insert_size, insert_size], "-", color="green", alpha=0.75, lw=0.25
            )

        for pair_step in steps[0].info["PAIR_STEPS"]:
            pair_step.info["INSERTSIZE"] = insert_size
            plot_pair_plan(ranges, pair_step, ax, marker_size)

        for split_step in steps[0].info["SPLIT_STEPS"]:
            split_step.info["INSERTSIZE"] = insert_size
            plot_split_plan(ranges, split_step, ax, marker_size)

    return [curr_min_insert_size, curr_max_insert_size]


# }}}

# {{{def plot_long_reads(long_reads,
def plot_long_reads(long_reads, ax, ranges, curr_min_insert_size, curr_max_insert_size):
    """Plots all LongReads for the region
    """

    Path = mpath.Path

    colors = {
        "Align": "orange",
        "Deletion": "black",
        "Inversion": "blue",
        "Duplication": "red",
        "InterChrm": "black",
        "InterChrmInversion": "blue",
    }

    for read_name in long_reads:
        long_read_plan = get_long_read_plan(read_name, long_reads, ranges)

        if long_read_plan is None:
            continue
        max_gap = long_read_plan[0]
        steps = long_read_plan[1]
        for step in steps:

            p = [
                map_genome_point_to_range_points(
                    ranges, step.start_pos.chrm, step.start_pos.start
                ),
                map_genome_point_to_range_points(
                    ranges, step.end_pos.chrm, step.end_pos.end
                ),
            ]

            # some points are far outside of the printable area, so we
            # ignore them
            if not points_in_window(p):
                continue

            if step.info["TYPE"] == "Align":
                ax.plot(
                    p,
                    [max_gap, max_gap],
                    "-",
                    color=colors[step.info["TYPE"]],
                    alpha=0.25,
                    lw=1,
                )

                if curr_max_insert_size and (max_gap > curr_max_insert_size):
                    curr_max_insert_size = max_gap
            else:
                x1 = p[0]
                x2 = p[1]

                pp = mpatches.PathPatch(
                    Path(
                        [
                            (x1, max_gap),
                            (x1, max_gap * 1.1),
                            (x2, max_gap * 1.1),
                            (x2, max_gap),
                        ],
                        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
                    ),
                    fc="none",
                    color=colors[step.info["TYPE"]],
                    alpha=0.25,
                    lw=1,
                    ls=":",
                )
                ax.add_patch(pp)

                # add some room for the bend line
                if (
                    curr_max_insert_size is None
                ) or max_gap * 1.1 > curr_max_insert_size:
                    curr_max_insert_size = max_gap * 1.1

    return [curr_min_insert_size, curr_max_insert_size]


# }}}

##Setup
# {{{def pair(arg):
def pair(arg):
    """Defines behavior for ArgParse pairs 

    Pairs must be comma-separated list of two items
    """
    try:
        parsed_arg = [int(x) for x in arg.split(",")]
        if len(parsed_arg) == 2:
            return parsed_arg
        else:
            sys.exit("Invalid number of pair values")
    except:
        sys.exit("Invalid pair values")


# }}}

# {{{def print_arguments(options):
def print_arguments(options):
    """Prints out the arguments to samplot as a json object

    Used as metadata for PlotCritic
    """
    if options.print_args or options.json_only:
        import json

        args_filename = os.path.splitext(options.output_file)[0] + ".json"
        args_info = {
            "titles": options.titles if options.titles else "None",
            "reference": options.reference if options.reference else "None",
            "bams": options.bams,
            "output_file": options.output_file,
            "start": options.start,
            "end": options.end,
            "chrom": options.chrom,
            "window": options.window,
            "max_depth": options.max_depth if options.max_depth else "None",
            "sv_type": options.sv_type,
            "transcript_file": options.transcript_file
            if options.transcript_file
            else "None",
        }
        with open(args_filename, "w") as outfile:
            json.dump(args_info, outfile)


# }}}


# {{{def setup_arguments():
def add_plot(parent_parser):
    """Defines the allowed arguments for plot function
    """
    parser = parent_parser.add_parser(
        "plot",
        help="Plot an image of a genome region from "
        + "CRAM/SAM alignments, "
        + "optimized for structural variant call review",
    )

    parser.add_argument(
        "-n",
        "--titles",
        help="Space-delimited list of plot titles. "
        + "Use quote marks to include spaces "
        + '(i.e. "plot 1" "plot 2")',
        type=str,
        nargs="+",
        required=False,
    )

    parser.add_argument(
        "-r",
        "--reference",
        help="Reference file for CRAM, required if " + "CRAM files used",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-z",
        "--z",
        type=int,
        default=4,
        help="Number of stdevs from the mean (default 4)",
        required=False,
    )

    def bam_file(bam):
        if not os.path.isfile(bam):
            parser.error("alignment file {} does not exist or is not a valid file".format(bam))
        options = ["sam", "bam", "cram"]
        idx_options = ["sai", "bai", "crai"]
        fields = os.path.splitext(bam)
        ext = fields[1][1:].lower()
        if ext not in options:
            parser.error("alignment file {} is not in SAM/BAM/CRAM format".format(bam))
        idx_type = idx_options[options.index(ext)]
        if not os.path.isfile(bam + "." + idx_type):
            parser.error("alignment file {} has no index".format(bam))
        return bam


    parser.add_argument(
        "-b",
        "--bams",
        type=bam_file,
        nargs="+",
        help="Space-delimited list of BAM/CRAM file names",
        required=True,
    )

    parser.add_argument(
        "-o", 
        "--output_file", 
        type=str, 
        help="Output file name/type. "
        +"Defaults to {type}_{chrom}_{start}_{end}.png",
        required=False,
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        default=".",
        help="Output directory name. Defaults to working dir. "
        +"Ignored if --output_file is set",
        required=False,
    )

    parser.add_argument(
        "-s",
        "--start",
        type=int,
        help="Start position of region/variant",
        action="append",
        required=True,
    )

    parser.add_argument(
        "-e",
        "--end",
        type=int,
        help="End position of region/variant",
        action="append",
        required=True,
    )

    parser.add_argument(
        "-c",
        "--chrom", type=str,
        help="Chromosome",
        action="append",
        required=True
    )

    parser.add_argument(
        "-w",
        "--window",
        type=int,
        help="Window size (count of bases to include " + "in view), default(0.5 * len)",
        required=False,
    )

    parser.add_argument(
        "-d",
        "--max_depth",
        type=int,
        help="Max number of normal pairs to plot",
        default=1,
        required=False,
    )

    parser.add_argument(
        "-t",
        "--sv_type",
        type=str,
        help="SV type. If omitted, plot is created " + "without variant bar",
        required=False,
    )
    
    def gff_file(transcript_file):
        if not os.path.isfile(transcript_file):
            parser.error("transcript file {} does not exist or is not a valid file".format(transcript_file))
        options = ["gff", "gff3"]
        fields = os.path.splitext(transcript_file)
        ext = fields[1][1:]
        if ext == "gz":
            ext = os.path.splitext(fields[0])[1][1:]
        ext = ext.lower()
        if ext not in options:
            parser.error("transcript file {} is not in GFF/GFF3 format".format(transcript_file))

        idx_file = transcript_file + ".tbi"
        if not os.path.isfile(idx_file):
            parser.error("transcript file {} is missing .tbi index file".format(transcript_file))
        return transcript_file

    parser.add_argument(
        "-T", "--transcript_file",
        help="GFF3 of transcripts",
        required=False,
        type=gff_file,
    )

    parser.add_argument(
        "--transcript_filename",
        help="Name for transcript track",
        required=False,
        type=str,
    )

    def bed_file(annotation_file):
        if not os.path.isfile(annotation_file):
            parser.error("annotation file {} does not exist or is not a valid file".format(annotation_file))
        fields = os.path.splitext(annotation_file)
        ext = fields[1][1:]
        if ext == "gz":
            ext = os.path.splitext(fields[0])[1][1:]
        ext = ext.lower()
        if ext != "bed":
            parser.error("annotation file {} is not in BED format".format(annotation_file))

        idx_file = annotation_file + ".tbi"
        if not os.path.isfile(idx_file):
            parser.error("annotation file {} is missing .tbi index file".format(annotation_file))
        return annotation_file

    parser.add_argument(
        "-A",
        "--annotation_files",
        type=bed_file,
        nargs="+",
        help="Space-delimited list of bed.gz tabixed "
        + "files of annotations (such as repeats, "
        + "mappability, etc.)",
        required=False,
    )
    
    parser.add_argument(
        "--annotation_filenames",
        type=str,
        nargs="+",
        help="Space-delimited list of names for the tracks in --annotation_files",
        required=False,
    )

    parser.add_argument(
        "--coverage_tracktype",
        type=str,
        help="type of track to use for low MAPQ " + "coverage plot.",
        choices=["stack", "superimpose"],
        default="stack",
        required=False,
    )

    parser.add_argument(
        "-a",
        "--print_args",
        action="store_true",
        default=False,
        help="Print commandline arguments",
        required=False,
    )

    parser.add_argument(
        "-H", "--plot_height", type=int, help="Plot height", required=False
    )

    parser.add_argument(
        "-W", "--plot_width", type=int, help="Plot width", required=False
    )

    parser.add_argument(
        "-q",
        "--include_mqual",
        type=int,
        help="Min mapping quality of reads to be included in plot",
        default=1,
        required=False,
    )

    parser.add_argument(
        "--separate_mqual",
        type=int,
        help="coverage from reads with MAPQ <= separate_mqual "
        + "plotted in lighter grey. To disable, "
        + "pass in negative value",
        default=0,
        required=False,
    )

    parser.add_argument(
        "-j",
        "--json_only",
        action="store_true",
        default=False,
        help="Create only the json file, not the " + "image plot",
        required=False,
    )

    parser.add_argument(
        "--start_ci",
        help="confidence intervals of SV first "
        + "breakpoint (distance from the "
        + "breakpoint). Must be a "
        + "comma-separated pair of ints (i.e. 20,40)",
        type=pair,
        required=False,
    )

    parser.add_argument(
        "--end_ci",
        help="confidence intervals of SV end "
        + "breakpoint (distance from the "
        + "breakpoint). Must be a "
        + "comma-separated pair of ints (i.e. 20,40)",
        type=pair,
        required=False,
    )

    parser.add_argument(
        "--long_read",
        type=int,
        default=1000,
        help="Min length of a read to be treated as a " + "long-read (default 1000)",
        required=False,
    )

    parser.add_argument(
        "--ignore_hp",
        action="store_true",
        help="Choose to ignore HP tag in alignment files",
        required=False,
    )
    parser.add_argument(
        "--min_event_size",
        type=int,
        default=20,
        help="Min size of an event in long-read " + "CIGAR to include (default 20)",
        required=False,
    )

    parser.add_argument(
        "--xaxis_label_fontsize",
        type=int,
        default=6,
        help="Font size for X-axis labels (default 6)",
        required=False,
    )

    parser.add_argument(
        "--yaxis_label_fontsize",
        type=int,
        default=6,
        help="Font size for Y-axis labels (default 6)",
        required=False,
    )

    parser.add_argument(
        "--legend_fontsize",
        type=int,
        default=6,
        help="Font size for legend labels (default 6)",
        required=False,
    )

    parser.add_argument(
        "--annotation_fontsize",
        type=int,
        default=6,
        help="Font size for annotation labels (default 6)",
        required=False,
    )

    parser.add_argument(
        "--common_insert_size",
        action="store_true",
        default=False,
        help="Set common insert size for all plots",
        required=False,
    )

    parser.add_argument(
        "--hide_annotation_labels",
        action="store_true",
        default=False,
        help="Hide the label (fourth column text) "
        + "from annotation files, useful for regions "
        + "with many annotations",
        required=False,
    )

    parser.add_argument(
        "--coverage_only",
        action="store_true",
        default=False,
        help="Hide all reads and show only coverage",
        required=False,
    )

    parser.add_argument(
        "--same_yaxis_scales",
        action="store_true",
        default=False,
        help="Set the scales of the Y axes to the " + "max of all",
        required=False,
    )

    parser.add_argument(
        "--marker_size",
        type=int,
        default=3,
        help="Size of marks on pairs and splits (default 3)",
        required=False,
    )
    
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Dots per inches (pixel count, default 300)",
        required=False,
    )
    parser.add_argument(
        "--zoom",
        type=int,
        default=500000,
        help="Only show +- zoom amount around breakpoints, "
            +"much faster for large regions. "
            +"Ignored if region smaller than --zoom (default 500000)",
        required=False,
    )

    parser.add_argument(
        "--debug",
        type=str,
        help="Print debug statements",
        required=False
    )

    parser.set_defaults(func=plot)


# }}}

# {{{def estimate_fragment_len(bam)
def estimate_fragment_len(bam, reference):
    try:
        if not reference:
            bam_file = pysam.AlignmentFile(bam, "rb")
        else:
            bam_file = pysam.AlignmentFile(bam, "rc", reference_filename=reference)
    except Exception as err:
        print("Error:", err, file=sys.stderr)
        sys.exit(1)

    frag_lens = []

    for i, read in enumerate(bam_file):
        if i >= 10000:
            break
        frag_lens.append(abs(read.tlen))
    if len(frag_lens) >= 5000:
        return np.median(frag_lens)
    else:
        print(
            "Insufficient reads for fragment length estimate.\nContinuing with unmodified window size",
            file=sys.stderr,
        )
        return 0


# {{{def set_plot_dimensions(sv,
def set_plot_dimensions(
    sv,
    sv_type,
    arg_plot_height,
    arg_plot_width,
    bams,
    reference,
    annotation_files,
    transcript_file,
    arg_window,
    zoom,
):
    """Chooses appropriate dimensions for the plot

    Includes the number of samples, whether a variant type is included, and
    any annotations in height Includes the start, end, and window argument
    in width If height and width are chosen by used, these are used instead

    Return plot height, width, and window as integers
    """

    plot_height = 5
    plot_width = 8
    if arg_plot_height:
        plot_height = arg_plot_height
    else:
        num_subplots = len(bams)
        if annotation_files:
            num_subplots += 0.3 * len(annotation_files)
        if transcript_file:
            num_subplots += 0.3
        plot_height = 2 + num_subplots

    if arg_plot_width:
        plot_width = arg_plot_width

    window = 0
    ranges = []
    if arg_window:
        window = arg_window

    """ 
    Several things determine the window size. 
    1) SV is not given, window = 0
    1) SV is given
        1) it is directly set
        2) it is not directly set
           2.1) single interval SV
           2.2) zoom set
           2.3) 2-interval SV
    """
    # if an SV type is given, then expand the window around its bounds
    if sv_type:
        # if the sv has one interval then set the window proportional
        # to sv size and set one range
        if len(sv) == 1:
            if arg_window:
                window = arg_window
            else:
                window = int((sv[0].end - sv[0].start) / 2)
                frag_len = estimate_fragment_len(bams[0], reference)

                if (0 < frag_len) and (window < 1.5 * frag_len):
                    old_window = window
                    window = int(1.5 * frag_len)
                    print(
                        "Window size is under 1.5x the estimated fragment length "
                        + "and will be resized to {}. Rerun with -w {} to override".format(
                            window, old_window
                        ),
                        file=sys.stderr,
                    )

            ranges = [
                genome_interval(
                    sv[0].chrm, max(0, sv[0].start - window), sv[0].end + window
                )
            ]

            # if region is larger than zoom, set window to zoom and set two ranges
            if window >= zoom:
                window = zoom
                ranges = [
                    genome_interval(
                        sv[0].chrm,
                        max(0, sv[0].start - window),
                        sv[0].start + window,
                    ),
                    genome_interval(
                        sv[0].chrm, max(0, sv[0].end - window), sv[0].end + window
                    ),
                ]
        elif len(sv) == 2:
            if arg_window:
                window = arg_window
            elif zoom:
                window = zoom
            else:
                window = 1000

            ranges = [
                genome_interval(
                    sv[0].chrm, max(0, sv[0].start - window), sv[0].start + window
                ),
                genome_interval(
                    sv[1].chrm, max(0, sv[1].end - window), sv[1].end + window
                ),
            ]
        else:
            sys.stderr.write(str(len(sv)) + " genome splits are not supported")
            sys.exit(1)
    else:
        ranges = [genome_interval(sv[0].chrm, sv[0].start, sv[0].end)]

    return plot_height, plot_width, window, ranges


# }}}

# {{{def get_read_data(ranges,
def get_read_data(
    ranges,
    bams,
    reference,
    separate_mqual,
    include_mqual,
    coverage_only,
    long_read_length,
    min_event_size,
    same_yaxis_scales,
    max_depth,
    z_score,
    ignore_hp,
):
    """Reads alignment files to extract reads for the region

    Region and alignment files given with chrom, start, end, bams
    If CRAM files are used, reference must be provided
    Reads with mapping quality below include_mqual will not be retrieved
    If coverage_only, reads are not kept and used only for checking
    coverage Reads longer than long_read_length will be treated as long
    reads Max coverages values will be set to same value for all samples if
    same_yaxis_scales If max_depth, only max_depth reads will be retrieved,
    although all will be included in coverage If PairedEnd read insert size
    is greater than z_score standard deviations from mean, read will be
    treated as discordant
    """

    all_pairs = []
    all_splits = []
    all_coverages = []
    all_long_reads = []
    all_linked_reads = []

    max_coverage = 0

    for bam_file_name in bams:
        bam_file = None
        try:
            if not reference:
                bam_file = pysam.AlignmentFile(bam_file_name, "rb")
            else:
                bam_file = pysam.AlignmentFile(
                    bam_file_name, "rc", reference_filename=reference
                )
        except Exception as err:
            print("Error:", err, file=sys.stderr)
            sys.exit(1)

        pairs = {}
        splits = {}
        long_reads = {}
        coverage = {}
        linked_reads = {}

        for r in ranges:
            try:
                bam_iter = bam_file.fetch(r.chrm, max(0, r.start - 1000), r.end + 1000)
            except ValueError:
                chrm = r.chrm
                if chrm[:3] == "chr":
                    chrm = chrm[3:]
                else:
                    chrm = "chr" + chrm
                bam_iter = bam_file.fetch(chrm, max(0, r.start - 1000), r.end + 1000)

            for read in bam_iter:
                if (
                    read.is_qcfail
                    or read.is_unmapped
                    or read.is_duplicate
                    or int(read.mapping_quality) < include_mqual
                ):
                    continue

                if not coverage_only:
                    if read.query_length >= long_read_length:
                        add_long_reads(bam_file, read, long_reads, min_event_size, ignore_hp)
                    else:
                        add_pair_end(bam_file, read, pairs, linked_reads, ignore_hp)
                        add_split(read, splits, bam_file, linked_reads, ignore_hp)
                add_coverage(bam_file, read, coverage, separate_mqual, ignore_hp)

        if (
            len(pairs) == 0
            and len(splits) == 0
            and len(long_reads) == 0
            and len(linked_reads) == 0
        ):
            if not coverage_only:
                print(
                    "Warning: No data returned from fetch in "
                    + "regions  "
                    + " ".join([str(r) for r in ranges])
                    + " from "
                    + bam_file_name,
                    file=sys.stderr,
                )

        for chrm in coverage:
            for pos in coverage[chrm]:
                sn_coverages = [
                    v for values in coverage[chrm][pos].values() for v in values
                ]
                curr_max = 0
                if len(sn_coverages) > 0:
                    curr_max = np.percentile(sn_coverages, 99.5)
                if curr_max > max_coverage:
                    max_coverage = curr_max
        all_coverages.append(coverage)
        all_pairs.append(pairs)
        all_splits.append(splits)
        all_long_reads.append(long_reads)
        all_linked_reads.append(linked_reads)

    read_data = {
        "all_pairs": all_pairs,
        "all_splits": all_splits,
        "all_coverages": all_coverages,
        "all_long_reads": all_long_reads,
        "all_linked_reads": all_linked_reads,
    }

    # Sample +/- pairs in the normal insert size range
    if max_depth:
        read_data["all_pairs"] = downsample_pairs(
            max_depth, z_score, read_data["all_pairs"]
        )

    return read_data, max_coverage


# }}}

# {{{def downsample_pairs(max_depth, z_score, all_pairs):
def downsample_pairs(max_depth, z_score, all_pairs):
    """Downsamples to keep only max_depth normal pairs from all PairedEnd
    reads 
    """
    for bam_i in range(len(all_pairs)):
        for hp_i in all_pairs[bam_i]:
            all_pairs[bam_i][hp_i] = sample_normal(
                max_depth, all_pairs[bam_i][hp_i], z_score
            )
    return all_pairs


# }}}

# {{{def set_haplotypes(curr_coverage):
def set_haplotypes(curr_coverage):
    """Creates a list to manage counting haplotypes for subplots
    """
    hps = sorted(curr_coverage.keys(), reverse=True)
    # if there are multiple haplotypes, must have 0,1,2
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


# }}}

# {{{def plot_samples(ranges,
def plot_samples(
    ranges,
    read_data,
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
    max_coverage,
    marker_size,
):

    """Plots all samples
    """
    max_insert_size = 0
    for i in range(len(bams)):
        #ax is never used, annotating this for readability
        ax = plt.subplot(grid[ax_i])
        hps = set_haplotypes(read_data["all_coverages"][i])
        inner_axs = gridspec.GridSpecFromSubplotSpec(
            len(hps), 1, subplot_spec=grid[ax_i], wspace=0.0, hspace=0.5
        )
        axs = {}
        for j in range(len(hps)):
            axs[j] = plt.subplot(inner_axs[hps[j]])

        curr_min_insert_size = None
        curr_max_insert_size = 0

        cover_axs = {}
        for hp in hps:
            curr_ax = axs[hp]

            curr_splits = []
            if hp in read_data["all_splits"][i]:
                curr_splits = read_data["all_splits"][i][hp]

            curr_linked_reads = []
            if hp in read_data["all_linked_reads"][i]:
                curr_linked_reads = read_data["all_linked_reads"][i][hp]

            curr_long_reads = []
            if hp in read_data["all_long_reads"][i]:
                curr_long_reads = read_data["all_long_reads"][i][hp]

            curr_pairs = []
            if hp in read_data["all_pairs"][i]:
                curr_pairs = read_data["all_pairs"][i][hp]

            curr_coverage = {}
            if hp in read_data["all_coverages"][i]:
                curr_coverage = read_data["all_coverages"][i][hp]

            cover_ax = plot_coverage(
                curr_coverage,
                curr_ax,
                ranges,
                len(hps),
                max_coverage,
                coverage_tracktype,
                yaxis_label_fontsize,
                same_yaxis_scales,
            )

            if len(curr_linked_reads) > 0:
                curr_min_insert_size, curr_max_insert_size = plot_linked_reads(
                    curr_pairs,
                    curr_splits,
                    curr_linked_reads,
                    curr_ax,
                    ranges,
                    curr_min_insert_size,
                    curr_max_insert_size,
                    marker_size,
                )
            elif len(curr_long_reads) > 0:
                curr_min_insert_size, curr_max_insert_size = plot_long_reads(
                    curr_long_reads,
                    curr_ax,
                    ranges,
                    curr_min_insert_size,
                    curr_max_insert_size,
                )
            else:
                curr_min_insert_size, curr_max_insert_size = plot_pairs(
                    curr_pairs,
                    curr_ax,
                    ranges,
                    curr_min_insert_size,
                    curr_max_insert_size,
                    marker_size,
                )

                curr_min_insert_size, curr_max_insert_size = plot_splits(
                    curr_splits,
                    curr_ax,
                    ranges,
                    curr_min_insert_size,
                    curr_max_insert_size,
                    marker_size,
                )

            cover_axs[hp] = cover_ax
            if curr_max_insert_size and (curr_max_insert_size > max_insert_size):
                max_insert_size = curr_max_insert_size

        # {{{ set axis parameters
        # set the axis title to be either one passed in or filename
        curr_ax = axs[hps[0]]

        if titles and len(titles) == len(bams):
            curr_ax.set_title(titles[ax_i - 1], fontsize=8, loc="left")
        else:
            curr_ax.set_title(os.path.basename(bams[ax_i - 1]), fontsize=8, loc="left")

        if len(axs) > 1:
            for j in axs:
                curr_ax = axs[j]
                fp = dict(size=8, backgroundcolor="white")
                text = "HP: "
                if j == 0:
                    text += "Undef"
                else:
                    text += str(j)
                at = AnchoredText(
                    text, loc=2, prop=fp, borderpad=0, pad=0, frameon=False
                )
                curr_ax.add_artist(at)

        for j in hps:
            curr_ax = axs[j]
            curr_ax.set_xlim([0, 1])
            if same_yaxis_scales:
                curr_ax.set_ylim([0, max(1, max_insert_size * 1.10)])
            else:
                curr_ax.set_ylim([0, max(1, curr_max_insert_size * 1.10)])
            curr_ax.spines["top"].set_visible(False)
            curr_ax.spines["bottom"].set_visible(False)
            curr_ax.spines["left"].set_visible(False)
            curr_ax.spines["right"].set_visible(False)
            curr_ax.tick_params(axis="y", labelsize=yaxis_label_fontsize)
            # if there's one hp, 6 ticks fit. Otherwise, do 3
            tick_count = 6 if len(hps) == 1 else 3
            curr_ax.yaxis.set_major_locator(ticker.LinearLocator(tick_count))
            curr_ax.tick_params(axis="both", length=0)
            curr_ax.set_xticklabels([])
        last_sample_num = number_of_axes - 1
        if annotation_files:
            last_sample_num -= len(annotation_files)
        if transcript_file:
            last_sample_num -= 1

        if ax_i == last_sample_num:
            curr_ax = axs[hps[-1]]

            labels = []
            if len(ranges) == 1:
                labels = [
                    int(ranges[0].start + l * (ranges[0].end - ranges[0].start))
                    for l in curr_ax.xaxis.get_majorticklocs()
                ]
            elif len(ranges) == 2:
                x_ticks = curr_ax.xaxis.get_majorticklocs()
                labels_per_range = int(
                    len(curr_ax.xaxis.get_majorticklocs()) / len(ranges)
                )
                labels = [
                    int(ranges[0].start + l * (ranges[0].end - ranges[0].start))
                    for l in x_ticks[:labels_per_range]
                ]
                try:
                    labels += [
                        int(ranges[-1].start + l * (ranges[-1].end - ranges[-1].start))
                        for l in x_ticks[labels_per_range:]
                    ]
                except:
                    sys.exit(labels_per_range)
            else:
                sys.stderr.write("Ranges greater than 2 are not supported\n")
                sys.exit(1)

            curr_ax.set_xticklabels(labels, fontsize=xaxis_label_fontsize)
            chrms = [x.chrm for x in ranges]
            curr_ax.set_xlabel("Chromosomal position on " + "/".join(chrms), fontsize=8)

        curr_ax = axs[hps[int(len(hps) / 2)]]
        curr_ax.set_ylabel("Insert size", fontsize=8)
        cover_ax = cover_axs[hps[int(len(hps) / 2)]]
        cover_ax.set_ylabel("Coverage", fontsize=8)
        # }}}

        ax_i += 1
    return ax_i


# }}}

# {{{def plot_legend(fig, legend_fontsize):
def plot_legend(fig, legend_fontsize, marker_size):
    """Plots the figure legend
    """
    marker_colors = []
    marker_labels = []
    read_colors = {
        "Deletion/Normal": "black",
        "Duplication": "red",
        "Inversion": "blue",
        "Aligned long read": "orange",
        "Linked read": "green",
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
        legend_elements += [
            plt.Line2D([0, 0], [0, 1], color=color, linestyle="-", lw=1)
        ]
    if READ_TYPES_USED["Split-read"]:
        marker_labels.append("Split read")
        legend_elements += [
            plt.Line2D(
                [0, 0],
                [0, 1],
                markerfacecolor="None",
                markeredgecolor="grey",
                color="grey",
                marker="o",
                markersize=marker_size,
                linestyle=":",
                lw=1,
            )
        ]

    if (
        READ_TYPES_USED["Paired-end read"]
        or READ_TYPES_USED["Deletion/Normal"]
        or READ_TYPES_USED["Inversion"]
    ):
        marker_labels.append("Paired-end read")
        legend_elements += [
            plt.Line2D(
                [0, 0],
                [0, 1],
                markerfacecolor="None",
                markeredgecolor="grey",
                color="grey",
                marker="s",
                markersize=marker_size,
                linestyle="-",
                lw=1,
            )
        ]

    fig.legend(
        legend_elements, marker_labels, loc=1, fontsize=legend_fontsize, frameon=False
    )


# }}}

# {{{def create_gridspec(bams, transcript_file, annotation_files, sv_type ):
def create_gridspec(bams, transcript_file, annotation_files, sv_type, read_data):
    """Helper function for creation of a correctly-sized GridSpec instance
    """
    # give one axis to display each sample
    num_ax = len(bams)

    # add another if we are displaying the SV
    if sv_type:
        num_ax += 1

    # add another if a annotation file is given
    if transcript_file:
        num_ax += 1

    if annotation_files:
        num_ax += len(annotation_files)

    # set the relative sizes for each
    ratios = []
    if sv_type:
        ratios = [1]

    for i in range(len(bams)):
        ratios.append(len(read_data["all_coverages"][i]) * 3)
        if len(read_data["all_coverages"]) > 0:
            ratios[-1] = 9

    if annotation_files:
        ratios += [0.3] * len(annotation_files)
    if transcript_file:
        ratios.append(2)

    return gridspec.GridSpec(num_ax, 1, height_ratios=ratios), num_ax


# }}}

##Annotations/Transcript methods
# {{{def get_plot_annotation_plan(ranges, annotation_file):
def get_plot_annotation_plan(ranges, annotation_file):
    annotation_plan = []
    for r in ranges:
        itr = get_tabix_iter(r.chrm, r.start, r.end, annotation_file)
        if not (itr):
            continue
        for row in itr:
            A = row.rstrip().split()
            A[0] = A[0].replace("chr", "")
            chrm = A[0]
            start = int(A[1])
            end = int(A[2])

            interval = genome_interval(chrm, start, end)

            # check to see if any part of this alignment overlaps a plot
            # range
            in_range = False
            for r in ranges:
                if r.intersect(interval) == 0:
                    in_range = True
            if in_range:
                step = plan_step(
                    genome_interval(chrm, start, start),
                    genome_interval(chrm, end, end),
                    "ANNOTATION",
                )
                if len(A) > 3:
                    try:
                        v = float(A[3])
                        step.event = "FLOAT_ANNOTATION"
                        step.info = v
                    except ValueError:
                        step.event = "STRING_ANNOTATION"
                        step.info = A[3]

                annotation_plan.append(step)
    return annotation_plan


# }}}

# {{{def plot_annotations(annotation_files, chrom, start, end,
def plot_annotations(
    annotation_files, annotation_filenames, ranges, hide_annotation_labels, annotation_fontsize, grid, ax_i
):
    """Plots annotation information from region 
    """
    if not annotation_filenames:
        annotation_filenames = []
        for annotation_file in annotation_files:
            annotation_filenames.append(os.path.basename(annotation_file))

    for i,annotation_file in enumerate(annotation_files):
        annotation_plan = get_plot_annotation_plan(ranges, annotation_file)
        annotation_filename = annotation_filenames[i]

        if len(annotation_plan) == 0:
            continue
        ax = plt.subplot(grid[ax_i])
        ax_i += 1

        for step in annotation_plan:
            p = [
                map_genome_point_to_range_points(
                    ranges, step.start_pos.chrm, step.start_pos.start
                ),
                map_genome_point_to_range_points(
                    ranges, step.end_pos.chrm, step.end_pos.end
                ),
            ]
            # if an annotation lies outside the window, its coordinate will be None, so we trim to the window
            if p[0] is None:
                p[0] = 0
            if p[1] is None:
                p[1] = 1

            if step.event == "ANNOTATION":
                ax.plot(p, [0, 0], "-", color="black", lw=5)
            elif step.event == "FLOAT_ANNOTATION":
                ax.plot(p, [0, 0], "-", color=str(step.info), lw=15)
            elif step.event == "STRING_ANNOTATION":
                ax.plot(p, [0, 0], "-", color="black", lw=15)
                if step.info and not hide_annotation_labels:
                    ax.text(
                        p[0],
                        0 + 0.1,
                        step.info,
                        color="black",
                        fontsize=annotation_fontsize,
                    )
            else:
                sys.stderr.write("Unsupported annotation type: " + step.event)
                sys.exit(1)

            # set axis parameters
            ax.set_xlim([0, 1])
            ax.spines["top"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_title(annotation_filename, fontsize=8, loc="left")
            ax.tick_params(axis="x", length=0)
            ax.tick_params(axis="y", length=0)
            ax.set_xticklabels([])
            ax.set_yticklabels([])


# }}}

# {{{def get_interval_range_plan_start_end(ranges, interval):
def get_interval_range_plan_start_end(ranges, interval):

    # transcript can span ranges
    start_range_hit_i = get_range_hit(ranges, interval.chrm, interval.start)
    end_range_hit_i = get_range_hit(ranges, interval.chrm, interval.end)

    if start_range_hit_i is None and end_range_hit_i is None:
        for i, range_item in enumerate(ranges):
            if (
                (range_item.chrm.strip("chr") == interval.chrm.strip("chr"))
                and (interval.start <= range_item.start <= interval.end)
                and (interval.start <= range_item.end <= interval.end)
            ):
                start_range_hit_i = i
                end_range_hit_i = i

    start = None
    end = None
    # neither end is in range, add nothing
    if start_range_hit_i == None and end_range_hit_i == None:
        return None, None
    # start is in, end is not
    elif end_range_hit_i == None:
        start = genome_interval(
            interval.chrm,
            max(interval.start, ranges[start_range_hit_i].start),
            max(interval.start, ranges[start_range_hit_i].start),
        )
        end = genome_interval(
            interval.chrm, ranges[start_range_hit_i].end, ranges[start_range_hit_i].end
        )
    # end is in, start is not
    elif start_range_hit_i == None:
        start = genome_interval(
            interval.chrm, ranges[end_range_hit_i].start, ranges[end_range_hit_i].start
        )
        end = genome_interval(
            interval.chrm,
            min(interval.end, ranges[end_range_hit_i].end),
            min(interval.end, ranges[end_range_hit_i].end),
        )
    # in same range or in different ranges
    else:
        start = genome_interval(
            interval.chrm,
            max(interval.start, ranges[start_range_hit_i].start),
            max(interval.start, ranges[start_range_hit_i].start),
        )
        end = genome_interval(
            interval.chrm,
            min(interval.end, ranges[end_range_hit_i].end),
            min(interval.end, ranges[end_range_hit_i].end),
        )
    return start, end


# }}}

# {{{def get_transcript_plan(ranges, transcript_file):
def get_transcript_plan(ranges, transcript_file):
    transcript_plan = []

    genes = {}
    transcripts = {}
    cdss = {}

    for r in ranges:
        itr = get_tabix_iter(r.chrm, r.start, r.end, transcript_file)
        for row in itr:
            gene_annotation = row.rstrip().split()

            if gene_annotation[2] == "gene":
                info = dict(
                    [list(val.split("=")) for val in gene_annotation[8].split(";")]
                )

                info["strand"] = gene_annotation[6] == "+"

                genes[info["Name"]] = [
                    genome_interval(
                        gene_annotation[0],
                        int(gene_annotation[3]),
                        int(gene_annotation[4]),
                    ),
                    info,
                ]
            elif gene_annotation[2] in ["transcript", "mRNA"]:
                info = dict(
                    [list(val.split("=")) for val in gene_annotation[8].split(";")]
                )
                info["strand"] = gene_annotation[6] == "+"

                if info["Parent"] not in transcripts:
                    transcripts[info["Parent"]] = {}
                transcripts[info["Parent"]][info["ID"]] = [
                    genome_interval(
                        gene_annotation[0],
                        int(gene_annotation[3]),
                        int(gene_annotation[4]),
                    ),
                    info,
                ]
            elif gene_annotation[2] == "CDS":
                info = dict(
                    [list(val.split("=")) for val in gene_annotation[8].split(";")]
                )
                info["strand"] = gene_annotation[6] == "+"

                if info["Parent"] not in cdss:
                    cdss[info["Parent"]] = {}

                if info["ID"] not in cdss[info["Parent"]]:
                    cdss[info["Parent"]][info["ID"]] = []

                cdss[info["Parent"]][info["ID"]].append(
                    genome_interval(
                        gene_annotation[0],
                        int(gene_annotation[3]),
                        int(gene_annotation[4]),
                    )
                )
    transcript_plan = []
    for gene in genes:
        gene_id = genes[gene][1]["ID"]
        if gene_id not in transcripts:
            continue
        for transcript in transcripts[gene_id]:
            interval, info = transcripts[gene_id][transcript]
            start, end = get_interval_range_plan_start_end(ranges, interval)

            if not start or not end:
                continue

            step = plan_step(start, end, "TRANSCRIPT")
            step.info = {"Name": None, "Strand": None, "Exons": None}
            step.info["Name"] = info["Name"]
            step.info["Strand"] = info["strand"]

            exons = []
            if transcript in cdss:
                for cds in cdss[transcript]:
                    for exon in cdss[transcript][cds]:
                        start, end = get_interval_range_plan_start_end(ranges, exon)
                        if start and end:
                            exons.append(plan_step(start, end, "EXON"))
            if len(exons) > 0:
                step.info["Exons"] = exons

            transcript_plan.append(step)
    return transcript_plan


# }}}

# {{{ def plot_transcript(transcript_file, chrom, start, end,
def plot_transcript(
    transcript_file, transcript_filename, ranges, grid, annotation_fontsize, xaxis_label_fontsize
):
    """Plots a transcript file annotation
    """
    if not transcript_filename:
        transcript_filename = os.path.basename(transcript_file)
    transcript_idx = 0
    arrow_loc = 0.02
    ax = plt.subplot(grid[-1])

    transcript_plan = get_transcript_plan(ranges, transcript_file)

    for step in transcript_plan:
        p = [
            map_genome_point_to_range_points(
                ranges, step.start_pos.chrm, step.start_pos.start
            ),
            map_genome_point_to_range_points(
                ranges, step.end_pos.chrm, step.end_pos.end
            ),
        ]
        # if an annotation lies outside the window, its coordinate will be None, so we trim to the window
        if p[0] is None:
            p[0] = 0
        if p[1] is None:
            p[1] = 0

        ax.plot(
            p, [transcript_idx, transcript_idx], "-", color="cornflowerblue", lw=0.5
        )

        ax.text(
            p[0],
            transcript_idx + 0.02,
            step.info["Name"],
            color="blue",
            fontsize=annotation_fontsize,
        )

        if step.info["Strand"]:
            ax.annotate(
                "",
                xy=(1, transcript_idx),
                xytext=(1 - arrow_loc, transcript_idx),
                arrowprops=dict(arrowstyle="->", color="cornflowerblue", lw=1),
                annotation_clip=True,
            )
        else:
            ax.annotate(
                "",
                xy=(1 - arrow_loc, transcript_idx),
                xytext=(1, transcript_idx),
                arrowprops=dict(arrowstyle="->", color="cornflowerblue", lw=1),
                annotation_clip=True,
            )

        if step.info["Exons"]:
            for exon in step.info["Exons"]:
                p = [
                    map_genome_point_to_range_points(
                        ranges, exon.start_pos.chrm, exon.start_pos.start
                    ),
                    map_genome_point_to_range_points(
                        ranges, exon.end_pos.chrm, exon.end_pos.end
                    ),
                ]
                if not points_in_window(p):
                    continue

                ax.plot(
                    p,
                    [transcript_idx, transcript_idx],
                    "-",
                    color="cornflowerblue",
                    lw=4,
                )

        transcript_idx += 1

    # set axis parameters
    ax.set_xlim([0, 1])
    ax.set_ylim([transcript_idx * -0.1, 0.01+(transcript_idx * 1.01)])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(transcript_filename, fontsize=8, loc="left")


# }}}


########################################################################
# main block
########################################################################
def plot(parser):
    """
    To support translocations, the SVs are specified as an array of 
    genome_interval. For now we let that arry be size 1 or 2.
    """
    options = parser.parse_args()

    if options.print_args or options.json_only:
        print_arguments(options)
        if options.json_only:
            sys.exit(0)

    if options.output_file:
        output_file = options.output_file
    else:
        if not os.path.isdir(options.output_dir):
            os.mkdir(options.output_dir)
        name_fields = [
            options.sv_type,
            "-".join(options.chrom),
            "-".join([str(s) for s in options.start]),
            "-".join([str(e) for e in options.end]),
        ]
        if options.sv_type:
            output_file = os.path.join(options.output_dir, "_".join(name_fields))
        else:
            output_file = os.path.join(options.output_dir, "_".join(name_fields[1:]))
    if (options.annotation_files 
            and options.annotation_filenames 
            and len(options.annotation_files) != len(options.annotation_filenames)):
        print("annotation filenames do not match annotation files", file=sys.stderr)
        sys.exit()

    for bam in options.bams:
        if ".cram" in bam:
            if not options.reference:
                parser.print_help(sys.stderr)
                sys.exit("Error: Missing reference for CRAM")

    if len(options.chrom) != len(options.start) != len(options.end):
        print("The number of chromosomes, starts, and ends do not match.", file=sys.stderr)
        sys.exit()

    sv = []
    for i in range(len(options.chrom)):
        options.chrom[i] = options.chrom[i].replace("chr", "")
        sv.append(genome_interval(options.chrom[i], options.start[i], options.end[i]))
    # set up plot
    plot_height, plot_width, window, ranges = set_plot_dimensions(
        sv,
        options.sv_type,
        options.plot_height,
        options.plot_width,
        options.bams,
        options.reference,
        options.annotation_files,
        options.transcript_file,
        options.window,
        options.zoom,
    )

    marker_size = options.marker_size

    # set up sub plots
    matplotlib.rcParams.update({"font.size": 12})
    fig = plt.figure(figsize=(plot_width, plot_height), dpi=options.dpi)

    # read alignment data
    read_data, max_coverage = get_read_data(
        ranges,
        options.bams,
        options.reference,
        options.separate_mqual,
        options.include_mqual,
        options.coverage_only,
        options.long_read,
        options.min_event_size,
        options.same_yaxis_scales,
        options.max_depth,
        options.z,
        options.ignore_hp,
    )

    # set up grid organizer
    grid, num_ax = create_gridspec(
        options.bams,
        options.transcript_file,
        options.annotation_files,
        options.sv_type,
        read_data,
    )
    current_axis_idx = 0

    # plot variant on top
    if options.sv_type:
        current_axis_idx = create_variant_plot(
            grid,
            current_axis_idx,
            sv,
            options.sv_type,
            ranges,
            options.start_ci,
            options.end_ci,
        )

    # Plot each sample
    current_axis_idx = plot_samples(
        ranges,
        read_data,
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
        max_coverage,
        marker_size,
    )
    # plot legend
    plot_legend(fig, options.legend_fontsize, marker_size)

    # Plot annotation files
    if options.annotation_files:
        plot_annotations(
            options.annotation_files,
            options.annotation_filenames,
            ranges,
            options.hide_annotation_labels,
            options.annotation_fontsize,
            grid,
            current_axis_idx,
        )

    # Plot sorted/bgziped/tabixed transcript file
    if options.transcript_file:
        plot_transcript(
            options.transcript_file,
            options.transcript_filename,
            ranges,
            grid,
            options.annotation_fontsize,
            options.xaxis_label_fontsize,
        )

    # save
    matplotlib.rcParams["agg.path.chunksize"] = 100000
    plt.tight_layout(pad=0.8, h_pad=0.1, w_pad=0.1)
    try:
        plt.savefig(output_file)
    except:
        print(
            "Failed to save figure " + output_file
        )
# }}}
