#!/usr/bin/python3

""" import standard python libraries """
import os.path
import sys


""" add puzzlotope folder to sys.path if we are in a folder called demo and there is a puzzlotope folder next to it """
# allows to run the demo without moving puzzlotope to system folders
if os.path.basename(os.path.dirname(os.path.abspath(__file__)))=='demo' and os.path.isdir('../puzzlotope'):
    sys.path.append(os.path.dirname('../puzzlotope'))
import puzzlotope
from puzzlotope.Chem import Block as Block
from puzzlotope.Chem import BlockTag as Tag

""" Set parameters """

# maximum mass difference between measured spectrum and theoretical value to still consider a match
tolerance = 0.3

# neglect isotopes combinations if they have a probability of less than this value. Set to 0 to not discard at all
cutoff_percent=1e-4

# consider the mass range from until (both boundaries are included)
limits=[414, 425]

# filename of spectrum measurement file
spectrum_file = 'spectrum.xy'

# project name (temporary files will contain this in their file name)
# use the name of this file without the .py extension
proj_name = os.path.splitext(os.path.basename(__file__))[0]

# Set to True if combinations should be recmputed
recompute=False

# Neglect Mg and Cl to speed up computation
Blocks = (
            Block('Ni'     , [0,1,2,3,4]   ),
#            Block('Mg'     , [2,]          ),
            Block('C6H5'   , [-1,]         , tags=[Tag.optional, Tag.canOxydize]),
            Block('C5H8'   , [0,]          , tags=[Tag.optional]),
#            Block('Cl'     , [-1,]         , tags=[Tag.optional]),
            Block('OH'     , [-1,]         , tags=[Tag.optional]),
            Block('O'      , [-2,]        , tags=[Tag.optional])
        )

puzzlotope.setBlocks(Blocks)

""" Run isotope puzzle """

puzzlotope.run(proj_name, spectrum_file, tolerance, cutoff_percent/100.0, limits, recompute=recompute, do_plots=False)
