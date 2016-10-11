#!/usr/bin/env python2
import sys, os.path, io, argparse, re

import qform
import libm
import matrix
import mpmat
import conmax

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", metavar="test.q",
                        help="write to file with Windows newlines")
    parser.add_argument("-m", metavar="libm", nargs="+",
                        help="only include tests from specified modules")
    args = parser.parse_args()

    if args.o:
        qform.output_file = io.open(args.o, "w", newline="\r\n")

    template = open(os.path.join(os.path.dirname(__file__), "test.q.template"))
    for line in template:
        line = line.rstrip()
        match = re.match(r"^---(\w+)---$", line)
        if match:
            module, = match.groups()
            if not args.m or module in args.m:
                print >>sys.stderr, "generating", module
                globals()[module].tests()
        else:
            qform.output(line)
