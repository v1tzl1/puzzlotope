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
# the file should contain two numeric values per line. The values should
# be separated by a space. The first value is a mass, the second one is
# a measured intensity (relative to other intensities in the same file)
spectrum_file = 'spectrum.xy'

# project name (temporary files will contain this in their file name)
# as this is used in filenames, spaces and umlauts should be avoided
#
# use the name of this file without the .py extension
proj_name = os.path.splitext(os.path.basename(__file__))[0]
# alternatively set a name yourself
# proj_name = 'myawesomeproject'

# Set to True if combinations should be recmputed
# otherwise they are loaded from a stored file if they have been
# computed before
#
# ! Important !
# When the Blocks, spectrum, or cidpeaks variables are changed a recomputation
# had to be performed (this is not detected automatically) otherwise the
# output is still for the old values.
recompute=True

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

# Measured CID peaks
cidpeaks = [357.12, 346.99, 338.01, 299.18, 289.05, 280.08, 269.95, 260.97, 231.12, 222.14, 212.01, 203.04, 192.91, 183.93, 154.08, 145.10, 134.97, 126.00, 115.87,  77.04,  68.06,  57.94]

# How many mass units can a measured and predicted peak be apart to
# still count as the same peak?
cidtolerance = 0.3

puzzlotope.setBlocks(Blocks)

""" Run isotope puzzle """

puzzlotope.run(proj_name, spectrum_file, tolerance, cutoff_percent/100.0, limits, cidpeaks, cidtolerance, recompute=recompute, do_plots=False, num_threads=1)
