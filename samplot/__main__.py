#!/usr/bin/env python
import argparse
import sys

from .__init__ import __version__
from .samplot import add_plot
from .samplot_vcf import add_vcf


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="samplot", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Installed version",
        action="version",
        version="%(prog)s " + str(__version__),
    )
    sub = parser.add_subparsers(title="[sub-commands]", dest="command")
    sub.required = True

    # TODO instead of doing this, just check the subcommand 
    # and call the appropriate function.
    # These functions just add args to the parser. 
    add_plot(sub)
    add_vcf(sub)

    # TODO this is such a opaque way to call the function,
    # change this to be as explicit as possible
    args,extra_args = parser.parse_known_args(args)
    args.func(parser, args, extra_args)


if __name__ == "__main__":
    sys.exit(main() or 0)
