#!/usr/bin/env python2

import os.path, io

import qform
import libm
import matrix
import mpmat
import conmax

if __name__ == '__main__':
    with io.open("test.q", "w", newline="\r\n") as qform.output_file:
        with open(os.path.join(os.path.dirname(__file__), "test.q.template"),
                  "rU") as read:
            for line in read:
                line = line.rstrip()
                if line == "---LIBM---":
                    libm.tests()
                elif line == "---MATRIX---":
                    matrix.tests()
                elif line == "---MPMAT---":
                    mpmat.tests()
                elif line == "---CONMAX---":
                    conmax.tests()
                else:
                    qform.output(line)
