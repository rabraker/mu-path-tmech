#!/usr/bin/python

import matlab.engine as me
import sys

script = sys.argv[1]

eng = me.connect_matlab()

eng.run(script, nargout=0)

# close all the leftover figures.
eng.close('all')
