#!/usr/bin/env python

import os.path

import matrix
import conmax

if __name__ == '__main__':
    with open("test.q", "w") as write:
        with open(os.path.join(os.path.dirname(__file__), "test.q.template"),
                  "rU") as read:
            def output(s):
                write.write("\r\n".join(s.splitlines()) + "\r\n")

            for line in read:
                line = line.rstrip()
                if line == "---MATRIX---":
                    matrix.output = output
                    matrix.tests()
                elif line == "---CONMAX---":
                    conmax.output = output
                    conmax.tests()
                else:
                    output(line)
