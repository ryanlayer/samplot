"""
Module containing types and functions for working
with genomic coordinates/intervals in genome and
plot space.

TODO
At the moment, there will be some redundancy between
some of the utilities here and what can be found in
samplot.py.  Eventually, those will be replaced if
applicable with these.
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass
class GenomeCoord:
    """
    Simple struct to store a genomic coordinate.
    """

    chrm: str
    at: int


# TODO make dataclass
class GenomeInterval:
    def __init__(self, chrm, start, end):
        self.chrm = chrm
        self.start = start
        self.end = end

    @property
    def startCoord(self):
        return GenomeCoord(self.chrm, self.start)

    @property
    def endCoord(self):
        return GenomeCoord(self.chrm, self.end)

    def __str__(self):
        return f"{self.chrm}:{self.start}-{self.end}"
        # return "(" + self.chrm + "," + str(self.start) + "," + str(self.end) + ")"

    def __repr__(self):
        return str(self)

    def __eq__(self, gi2):
        return self.chrm == gi2.chrm and self.start == gi2.start and self.end == gi2.end

    def intersect(self, gi):
        """ return -1 if before, 0 if in, 1 if after """
        if (
            gi.chrm.replace("chr", "") < self.chrm.replace("chr", "")
            or gi.end < self.start
        ):
            return -1
        elif (
            gi.chrm.replace("chr", "") > self.chrm.replace("chr", "")
            or gi.start > self.end
        ):
            return 1
        else:
            return 0

    def __contains__(self, query: GenomeCoord | "GenomeInterval") -> bool:
        """A simple boolean alternative to the intersect function
        using the 'in' operator. Called with `query in self`. The
        operation is commutative.
        """
        if self.chrm != query.chrm:
            return False
        # TODO when python 3.10 gets better support
        # in conda/mamba, switch to structural pattern match
        if isinstance(query, GenomeCoord):
            return self.start <= query.at <= self.end
        elif isinstance(query, GenomeInterval):
            return (
                self.start <= query.start <= self.end
                or self.start <= query.end <= self.end
            )
        else:
            raise TypeError(
                "GenomeIterval.__contains__ only accepts GenomeCoord or GenomeInterval"
            )


def coord_in_range(ranges: list[GenomeInterval], coord: GenomeCoord) -> int | None:
    """
    Check if a coordinate is in a list of genomic intervals.
    Returns the index of the interval in ranges that contains
    the coord, and None if it was not contained in any.

    TODO I'm a bit iffy on this function as well as get_range_hit
    I'm mimicing its functionality, but I think that there is a better
    overall way.
    """
    for i, r in enumerate(ranges):
        if coord in r:
            return i
    return


def interval_in_range(
    ranges: list[GenomeInterval], interval: GenomeInterval
) -> int | None:
    """
    Check if an interval is in a list of genomic intervals.
    Returns the index of the interval in ranges that contains
    the interval, and None if it was not contained in any.
    """
    for i, r in enumerate(ranges):
        if interval in r:
            return i
    return


def genomic_to_axes(ranges: list[GenomeInterval], coord: GenomeCoord,) -> float | None:
    """
    Map a genomic coordinate to plot space (x-axis).
    """
    # TODO there's redundant range checking here when its also done
    # during plan creation. keep as is for now, but later figure
    # out a way to add this information to the plan
    range_hit = coord_in_range(ranges, coord)

    # TODO there are ways to handle optional values using a monad
    # like pattern, but leave unchanged for now
    if range_hit is None:
        return
    r = ranges[range_hit]
    margin = float(r.end - r.start)
    r_scale = (1.0/len(ranges)) # ex: 1/3 for 3 ranges
    r_offset = r_scale * range_hit
    p = r_scale * float(coord.at - r.start) / margin + r_offset

    return p


# def points_in_window(points): # TODO add param for tolerance
#     """
#     Checks whether these points lie within the window of interest
#     Points is a list of one start, one end coordinate (ints)
#     """
#     if (
#         None in points
#         or points[0] < -5
#         or points[1] < -5
#         or points[0] > 5
#         or points[1] > 5
#     ):
#         return False
#     return True
