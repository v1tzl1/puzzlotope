#!/usr/bin/python3

""" import standard python libraries """
import os.path
import sys


""" add puzzlotope folder to sys.path if we are in a folder called demo and there is a puzzlotope folder next to it """
# allows to run the demo without moving puzzlotope to system folders
if os.path.basename(os.path.dirname(os.path.abspath(__file__)))=='demo' and os.path.isdir('../puzzlotope'):
    sys.path.append(os.path.dirname('../puzzlotope'))
import puzzlotope

""" Set parameters """

# maximum mass difference between measured spectrum and theoretical value to still consider a match
tolerance = 0.3

# neglect isotopes combinations if they have a probability of less than this value. Set to 0 to not discard at all
cutoff_percent=1e-5

# consider the mass range from until (both boundaries are included)
limits=[414, 425]

# filename of spectrum measurement file
spectrum_file = 'spectrum.xy'

# project name (temporary files will contain this in their file name)
proj_name = 'demo'

# skip computation steps
compute_from=0 # compute everything again every run
#compute_from=3 # do not recompute combinations and spectra, rather load them from disk. Do not use this when values above were changes


""" Run isotope puzzle """

puzzlotope.run(proj_name, spectrum_file, tolerance, cutoff_percent/100.0, limits, recompute_from=compute_from, do_plots=False)
