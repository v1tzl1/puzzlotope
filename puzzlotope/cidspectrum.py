#!/usr/bin/env python3

import copy
import time
import multiprocessing
import puzzlotope.Solver
import puzzlotope.Chem as Chem

# For test routine
from puzzlotope.Chem import Block as Block
from puzzlotope.Chem import BlockTag as Tag

def _buildSpectrumWorker(comb_queue, return_queue, peaks, tolerance):
	while True:
		try:
			comb = comb_queue.get(block=True, timeout=None)
			
			if comb == 'QUIT':
				break
			
			cidspectrum = CIDSpectrum(comb, peaks, tolerance)
			return_queue.put(cidspectrum, block=True, timeout=None)
		except KeyboardInterrupt:
			print('Keyboard interrupt, quitting')
			break

def buildSpectrums(combinations, cid_peaks, cid_tolerance, num_threads, verbose=False, prefix=''):
	spectra=[]
	
	if num_threads == 1:
		""" Single Thread version """
		for comb in combinations:
			if verbose:
				print(prefix+'Building CID spectrum for ' + str(comb))
			spectra.append(CIDSpectrum(comb, cid_peaks, cid_tolerance))
		
		if verbose:
			print(prefix+'  done building spectra')
		
	else:
		""" Multi Thread Option """
		comb_queue = multiprocessing.Queue(len(combinations)+num_threads)
		return_queue = multiprocessing.Queue(len(combinations))
		
		# queue combinations
		for comb in combinations:
			comb_queue.put(comb)
		
		# queue QUIT for every Thread
		for _ in range(num_threads):
			comb_queue.put( 'QUIT' )
		
		if verbose:
			print(prefix+'Starting %d parallel threads to compute CID spectra' % num_threads)
			
		# start processes
		processes = []
		for i in range(num_threads):
			p = multiprocessing.Process(target=_buildSpectrumWorker, args=(comb_queue, return_queue, cid_peaks, cid_tolerance))
			p.start()
			processes.append(p)
		
		# collect computed spectras
		while len(spectra) < len(combinations):
			_spectrum = return_queue.get(block=True, timeout=None)
			spectra.append(_spectrum)
			if verbose:
				time_str = time.strftime("%H:%M:%S")
				num_threads_active = sum([1 for p in processes if p.is_alive()])
				progress_total = 100.0*len(spectra)/len(combinations)
				print(prefix+'	%s CID spectrum for %s computed, %5.1f%% total, %2d/%d threads active' % (time_str, _spectrum.original_comb.getCombString(), progress_total, num_threads_active, num_threads))
		
		for p in processes:
			p.join()
		
		if verbose:
			print(prefix+'  done building spectra')
	
	return spectra

def printSpectra(cidspectra, peaks, tolerance, prefix='', f_txt=None, f_csv=None):
	TOPN=10 if len(cidspectra)>=10 else len(cidspectra)
	cidspectrasorted = sorted(cidspectra, key=lambda x: x.metric)
	
	if f_txt:
		print(cidspectra[0].getSummaryHeader(prefix), file=f_txt)
		for _s in cidspectrasorted:
			print(_s.getSummary(prefix), file=f_txt)
	
	if f_csv:
		print('"Combination","Metric","Measured but not predicted", "Predicted but not measured"', file=f_csv)
		for _s in cidspectrasorted:
			
			#print('"%s","%f","%s","%s"' % (_s.original_comb.getCombString(space=False), _s.metric, ' and '.join(_s.meas_only), ' and '.join(_s.pred_only)), file=f_csv)
			print('"%s","%f","%s","%s"' % (_s.original_comb.getCombString(space=False), _s.metric, _s.getMeasOnlyStr(), _s.getPredOnlyStr()), file=f_csv)
	
	print('%sTop %d of %d CID spectra' % (prefix, TOPN, len(cidspectra)))
	if f_txt or f_csv:
		print('%s  detailed lists are in:' % prefix)
		if f_txt:
			print('%s	%s' % (prefix, f_txt.name))
		if f_csv:
			print('%s	%s' % (prefix, f_csv.name))
		print('')
	print(cidspectra[0].getSummaryHeader(prefix))
	for i in range(TOPN):
		print(cidspectrasorted[i].getSummary(prefix))

class CIDSpectrum:
	original_comb=None
	combs=None
	peaks=None
	breakoffs=None
	remainings=None
	measured_peaks=None
	tolerance=None
	metric=None
	pred_only=None
	meas_only=None
	pred_only_name=None
	meas_only_name=None
	
	def __init__(self, combination, measured_peaks, tolerance):
		cstr = combination.getCombString()
		self.original_comb=combination
		self.measured_peaks = measured_peaks
		self.tolerance = tolerance
		
		#print("Breaking combination %s" % cstr)
		
		self.combs = CIDSpectrum.getCIDCombinations(combination)
		self.remainings = []
		self.breakoffs = []
		self.peaks = []
		
		for c in self.combs:
			self.remainings.append(c.getCombString(space=False))
			
			_breakoff = puzzlotope.Solver.CombinationResult.getDifference(combination, c)
			self.breakoffs.append(_breakoff.getCombString(space=False))
			
			self.peaks.append(c.weight)
		
		self.computeMetric();
	
	def iterate(self):
		
		yield self.original_comb.getCombString(space=False), '', self.original_comb.weight
		
		for rem, cut, w in zip(self.remainings, self.breakoffs, self.peaks):
			yield rem, cut, w
		
		return
	
	def computeMetric(self):
		# penalty if a peak is predicted, but not measured
		_w_pred_only=0.1 # can be fine, maybe the bond is so strong, that in practive this breakoff never happens and is thus not measured
		
		# penalty if a peak was measured, but not predicted
		_w_meas_only=1.0 # this is bad, as it means there is a split happening, but it cannot be explained by out building blocks
		
		
		if len(self.peaks) == 0:
			raise ValueError('No predicted peaks at all')
		
		_rems = copy.copy(self.remainings)
		_pred = copy.copy(self.peaks)
		_meas = copy.copy(self.measured_peaks)
		
		# Keep in mind that _meas could be empty
		self.pred_only	  = [_pred[i] for i in range(len(_pred)) if min([abs(_pred[i]-_m) for _m in _meas])>self.tolerance] if _meas else _pred
		self.pred_only_name = [_rems[i] for i in range(len(_pred)) if min([abs(_pred[i]-_m) for _m in _meas])>self.tolerance] 
		
		self.meas_only	  = [_meas[i] for i in range(len(_meas)) if min([abs(_p-_meas[i]) for _p in _pred])>self.tolerance] if _meas else []
		self.meas_only_name = [_rems[i] for i in range(len(_meas)) if min([abs(_p-_meas[i]) for _p in _pred])>self.tolerance] if _meas else []
		
		self.metric = _w_pred_only*len(self.pred_only) + _w_meas_only*len(self.meas_only)
	
	def getSummaryHeader(self, prefix):
		return '%s%-*s | %7s | M only | P only' % (prefix, self.original_comb.getStringLen(), 'Combination', 'Metric')
	
	def getSummary(self, prefix=''):
		if self.metric is None:
			self.computeMetric()
		
		#_ponly_out = [(_p, ('P:%6.2f' % _p)) for _p in _ponly]
		#_monly_out = [(_m, ('M:%6.2f' % _m)) for _m in _monly]
		#_out = sorted(_ponly_out+_monly_out, key=lambda x: x[0], reverse=True)
		#_out_str = ', '.join([_o[1] for _o in _out])
		return '%s%s | %7.2f |  % 3d   |  % 3d' % (prefix, self.original_comb.getCombString(), self.metric, len(self.meas_only), len(self.pred_only))
	
	@staticmethod
	def _getPeaksToString(peaks, names):
		if not len(peaks) == len(names):
			ValueError('Length of peaks and names does not match')
		
		if peaks:
			return ('%d peaks missing: '%len(peaks)) + ' and '.join(['%s with mass %f' % (_name, _mass) for _name, _mass in zip(names, peaks)])
		else:
			return 'None'
	
	def getPredOnlyStr(self):
		return CIDSpectrum._getPeaksToString(self.pred_only, self.pred_only_name)
	
	def getMeasOnlyStr(self):
		return CIDSpectrum._getPeaksToString(self.meas_only, self.meas_only_name)
	
	@staticmethod
	def getCIDCombinations(combination):
		_nodes = [combination.N,]
		_leafs = []
		
		while _nodes:
			_newnodes = []
			for n in _nodes:
				_newnodes += CIDSpectrum._createChildren(n)
				
			_leafs = CIDSpectrum._getuniques(_newnodes + _leafs)
			_nodes = CIDSpectrum._getuniques(_newnodes)
		
		combs = [puzzlotope.Solver.CombinationResult(combination.total, N=_N) for _N in sorted(_leafs, reverse=True)]
		return combs
	
	@staticmethod
	def _createChildren(combN):
		ret = []
		
		if sum(combN) == 1:
			return ret
		
		for i in range(Chem.getNumBlocks()):
			if combN[i] == 0:
				continue
			
			ret.append(CIDSpectrum._splitoff(combN, i))
		
		return ret
	
	@staticmethod
	def _splitoff(N, indx):
		if indx >= len(N):
			raise ValueError("Index out of bounds")
		if N[indx] == 0:
			raise ValueError("Attempting to split of a non existant block")

		newN = copy.copy(N)
		newN[indx] -= 1

		return newN
	
	@staticmethod
	def _getuniques(Ns):
		_d = dict()
		
		for N in Ns:
			_Nstr = '['+','.join([str(n) for n in N])+']'
			_d[_Nstr] = N
		
		return list(_d.values())
	
	def printSpectrum(self):
		l = self.original_comb.getStringLen()
		for r, c, w in sorted(list(self.iterate()), key=lambda x: x[2], reverse=True):
			print('%-*s has weight %12.8f (without %-*s)' % (l, r, w, l, c))
	
	
	"""
	
	
	children=None
	comb=None
	def __init__(self, combination):
		self.children=[]
		self.comb=combination

	def buildTree(self):
		if len(self.children) > 0:
			raise ValueError("Attempt to build existing tree")
		
		_nodes = [self,]
		while True:
			_newnodes = []
			for n in _nodes:
				_newnodes += n.createChildren()
			
			if not _newnodes:
				return set(_nodes)
			
			print('New nodes:   %s' % _newnodes)
			print('New uniques: %s' % set(_newnodes))
			_nodes = set(_newnodes)

	def createChildren(self):
		if sum(self.comb.N) <= 1:
			return []
		
		for i in range(len(self.comb.N)):
			if self.comb.N[i] == 0:
				continue
			
			self.children.append(CIDNode(self.comb.getSplitoff(i)))
		print('created %d new children: %s' % (len(self.children), str(self.children)))
		return self.children
		
	def __repr__(self):
		return 'TreeNode (%s)' % (self.comb.getCombString(space=False))
	"""	
	
if __name__ == '__main__':
	Blocks = (
			Block('Ni'	 , [0,1,2,3,4]   ),
			Block('OH'	 , [-1,]		 , tags=[Tag.optional]),
		)
	puzzlotope.setBlocks(Blocks)
	
	N = [2,3]
	mpeaks = [167, 150, 132.9, 116, 109, 92, 75, 34, 17, 42] # One peak missing at 57.93, one additional at 42
	
	# compute weight
	comb = puzzlotope.Solver.CombinationResult(0.0, N=N)
	
	# use computed weight as target
	comb = puzzlotope.Solver.CombinationResult(total=comb.weight, weight=comb.weight, N=N)
	
	s = CIDSpectrum(comb, mpeaks, 0.3)
	
	s.printResult()
	
