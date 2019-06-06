#!/usr/bin/python

import matlab.engine as me
import sys
import os

script = sys.argv[1]


cur_dir = os.getcwd()

eng = me.connect_matlab()

eng.cd(cur_dir, nargout=0)
eng.run(script, nargout=0)

# close all the leftover figures.
eng.close('all')
