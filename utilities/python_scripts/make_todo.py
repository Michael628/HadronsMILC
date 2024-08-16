#! /usr/bin/env python

# Python 3 version

import sys

# Script for creating todo file entries with precessed source times
# C. DeTar 

######################################################################
def main():

    if len(sys.argv) < 4:
        print("Usage", sys.argv[0], "<series>, <cfglow>, <cfghi>, <cfgstep>")
        sys.exit(1)

    ( series, cfglow, cfghi, cfgstep ) = sys.argv[1:5]

    cfgnoRange = range(int(cfglow),int(cfghi),int(cfgstep))

    for cfgno in cfgnoRange:
        print("{0:s}.{1:d} S 0 E 0".format(series, cfgno))

######################################################################
main()
