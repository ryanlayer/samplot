#!/usr/bin/env python
import argparse
import logging
import sys

from .__init__ import __version__
from .samplot import add_plot
from .samplot_vcf import add_vcf


def main(args=None):
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(module)s - %(levelname)s: %(message)s")
    
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

    add_plot(sub)
    add_vcf(sub)

    args,extra_args = parser.parse_known_args(args)
    args.func(parser, args, extra_args)


if __name__ == "__main__":
    sys.exit(main() or 0)
