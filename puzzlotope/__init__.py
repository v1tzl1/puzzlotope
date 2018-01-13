#/usr/bin/python3
import sys
import pickle
import multiprocessing
import os.path
from timeit import default_timer as timer

try:
    import numpy as np
except ImportError:
    print('Numpy is required')
    sys.exit(1)

try:
    from scipy.stats import stats
except ImportError:
    print('Scipy is required')
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print('Matplotlib is required')
    sys.exit(1)

import puzzlotope.Solver
import puzzlotope.Chem
import puzzlotope.Measurement
import puzzlotope.gmm
#import puzzlotope.probability

def setBlocks(blocks):
	for block in blocks:
		if not isinstance(block, puzzlotope.Chem.Block):
			raise TypeError('Element in blocks list is not of type Chem.Block')
	
	puzzlotope.Chem.Blocks = blocks

def run(proj_name, spectrum_fname, weight_tolerance, prob_cutoff, xlims, recompute=False, do_plots=False, num_threads=multiprocessing.cpu_count()):
	Chem.update()
	sep = '-'*15
	
	print(sep)
	print('\n  Calculating project \'%s\'\n' % proj_name)
	print(sep)
	
	### Step 0: Read measured Spectrum
	x, y = Measurement.parseFile(spectrum_fname, xmin=xlims[0]-0.5, xmax=xlims[1]+0.5)
	sigma, diracs = gmm.estimateSigma(x, y)
	masses_meas=[p.mass for p in diracs]
	weights_meas=[p.prob for p in diracs]
	mass_main_peak=masses_meas[weights_meas.index(max(weights_meas))]
	print('Step 0: Read spectrum from "%s"' % spectrum_fname)
	print('    using %d mass peaks' % len(diracs))
	print('    main mass was detected at %f' % mass_main_peak)
	print(sep)
	
	### Step 1: Computing possible combinations
	print('Step 1: Computing main isotope combinations')
	s1_fname='data/%s_combinations.p' % proj_name
	s1_txtname='data/%s_combinations.txt' % proj_name
	print('    storing results in %s' % s1_fname)
	print('    a detailed list of combinations is in %s' % s1_txtname)
	print('    target mass is %f' % mass_main_peak)
	
	if recompute or not os.path.isfile(s1_fname):
		allCombinations = puzzlotope.Solver.solveCombination(mass_main_peak)
		with open(s1_txtname,'w') as txtf:
			filteredCombinations = puzzlotope.Solver.CombinationResult.filter(allCombinations, weight_tolerance, file_out=txtf)
		with open(s1_fname, 'wb') as _f:
			pickle.dump(filteredCombinations, _f)
	else:
		with open(s1_fname, 'rb') as _f:
			filteredCombinations = pickle.load(_f)
			print('    loaded %d combinations' % len(filteredCombinations))
	print('    Using %d combinations:' % len(filteredCombinations))
	for comb in filteredCombinations:
		print('      ' + comb.getCombString())
	print(sep)
	
	### Step 2: Computing spectra for found combinations
	print('Step 2: Computing spectrum for each combination')
	s2_fname='data/%s_spectra.p' % proj_name
	print('    storing results in %s' % s2_fname)
	print('    number of combinations is %d' % len(filteredCombinations))
	if recompute or not os.path.isfile(s2_fname):
		start = timer()
		spectra = puzzlotope.Solver.buildSpectrums(filteredCombinations, weight_tolerance, prob_cutoff, num_threads, verbose=True, prefix='    ')
		end = timer()
		ellapsed = (end-start)
		print('    Time for spectrum computation: %f s' % ellapsed)
		print('      thats %f spectra per minute' % (len(filteredCombinations)/ellapsed*60.0))
		with open(s2_fname, 'wb') as _f:
			pickle.dump(spectra, _f)
	else:
		with open(s2_fname, 'rb') as _f:
			spectra = pickle.load(_f)
			print('    loaded %d spectra' % len(spectra))
	print(sep)
	
	
	### Step 3: Analyse spectrum mass peaks
	print('Step 3: Analyse spectrum mass peaks')
	s3_txtname = 'data/%s_spectra.txt' % proj_name
	s3_csvname = 'data/%s_spectra.csv' % proj_name
	with open(s3_txtname,'w') as txtf:
		with open(s3_csvname, 'w') as csvf:
			puzzlotope.Solver.Spectrum.printSpectra(spectra, xlims[0], xlims[1], diracs, prefix='    ', f_txt=txtf, f_csv=csvf)
	print(sep)
	
	"""
	### Step 4: Analyse predicted spectral measurements
	print('Step 4: Analyse predicted spectral measurements')
	
	pdf_meas=y/sum(y)
	pdfs=[pdf_meas,]
	labels=['Measurement',]
	print('Comparing spectra to measurement')
	for spectrum in spectra:
		_masses = [p.mass for p in spectrum.pms]
		ofst = gmm.estimateOffset(masses_meas, _masses, weight_tolerance)
		_pdf = gmm.buildPDF(x, spectrum.pms, -ofst, sigma)
		pdfs.append(_pdf)
		labels.append(spectrum.comb.getCombString(space=False))
		print('  %s: RMS=%12.5f' % (spectrum.comb.getCombString(), Measurement.getRMS(pdf_meas, _pdf)))
	_reproj_pdf = gmm.buildPDF(x, diracs, 0.0, sigma)
	print('Reprojection RMS')
	print('  %-*s: RMS=%12.5f' % (len(spectra[0].comb.getCombString()), 'Measurement', Measurement.getRMS(pdf_meas, _reproj_pdf)))
	print(sep)
	"""
	
	if do_plots:
		gmm.plotPDFs(x, pdfs, labels, 'Spectrum comparison')
	
