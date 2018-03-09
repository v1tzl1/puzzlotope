import math
import copy
import numpy as np
import multiprocessing
import time
import sys

import puzzlotope.Chem as Chem
import puzzlotope.Measurement
import puzzlotope.probability
from puzzlotope.Element import ProbMass

def solveCombination(total, verbose=False):
	
	tmp = [CombinationResult(total), ]
	
	for level in range(Chem.getNumBlocks()):
		if verbose:
			print('Solving level %d of %d' % (level, Chem.getNumBlocks()))
			print('  inputs=%s' % str(tmp))
		
		results = []
		for e in tmp:
			results += e.iterate(verbose=verbose)
		tmp=results
	return results

def _buildSpectrumWorker(comb_queue, return_queue, tolerance, cutoff):
	while True:
		try:
			comb = comb_queue.get(block=True, timeout=None)
			
			if comb == 'QUIT':
				break
			
			spectrum = Spectrum.computeSpectrum(comb, tolerance, cutoff)
			return_queue.put(spectrum, block=True, timeout=None)
		except KeyboardInterrupt:
			print('Keyboard interrupt, quitting')
			break

def buildSpectrums(combinations, tolerance, cutoff, num_threads, prefix='', verbose=False):
	spectra=[]
	
	if num_threads == 1:
		""" Single Thread version """
		for comb in combinations:
			if verbose:
				print(prefix+'Building Isotope spectrum for ' + str(comb))
			spectra.append(Spectrum.computeSpectrum(comb, tolerance, cutoff))
		
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
			print(prefix+'Starting %d parallel threads to compute spectra' % num_threads)
			
		# start processes
		processes = []
		for i in range(num_threads):
			p = multiprocessing.Process(target=_buildSpectrumWorker, args=(comb_queue, return_queue, tolerance, cutoff))
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
				print(prefix+'    %s Spectrum for %s computed, %5.1f%% total, %2d/%d threads active' % (time_str, _spectrum.comb.getCombString(), progress_total, num_threads_active, num_threads))
		
		for p in processes:
			p.join()
		
		if verbose:
			print(prefix+'  done building spectra')
	
	return spectra

class Spectrum:
	comb=None
	pms=[]
	scale=None
	metrics=None
	metrics_name=None
	
	@staticmethod
	def computeSpectrum(comb, tolerance, cutoff):
		tree=comb.buildIsotopeTree(cutoff)
		
		masses=[]
		probs=[]
		
		for leaf in tree.getLeafs():
			mass=0.0
			for isotope in leaf.getBranchString():
				mass+=Chem.getIsotopeMass(isotope)
			masses.append(mass)
			probs.append(leaf.prob)
		
		return Spectrum(comb, probs, masses, tolerance, cutoff)
		
	"""
	def __init__(self, combination, tolerance, cutoff):
		tree=combination.buildIsotopeTree(cutoff)
		
		self.comb = combination
		self.pms=[]
		
		for leaf in tree.getLeafs():
			mass=0.0
			for isotope in leaf.getBranchString():
				mass+=Chem.getIsotopeMass(isotope)
			self.pms.append(ProbMass(leaf.prob, mass))
		self.__normalize(tolerance, cutoff)
	"""
	
	def __init__(self, combination, probs, masses, tolerance, cutoff):
		#print('Creating spectrum with %d masses and %d probabilities' % (len(masses), len(probs)))
		self.comb = combination
		self.pms = [ProbMass(_p, _m) for _p, _m in zip(probs, masses)]
		self.__normalize(tolerance, cutoff)
		
	"""
	def __init__(self, combination, tolerance, cutoff):
		tmp = [ProbMass(),]
		for i in range(Chem.getNumBlocks()):
			if combination.N[i]>0:
				tmp = ProbMass.merge(tmp, Chem.Blocks[i].getMasses(cutoff), combination.N[i], cutoff)
		
		self.comb = combination
		self.pms=tmp
		#print('finished spectrum construction')
		self.__normalize(tolerance, cutoff)
	 """ 
	
	def __normalize(self, tolerance, cutoff):
		__pms=sorted(self.pms, reverse=True)
		__new_pms=[]
		while len(__pms)>0:
			pivot = __pms.pop(0)
			flag=True
			while flag:
				flag=False
				for i in range(len(__pms)):
					if pivot.getMassDiff(__pms[i])<tolerance:
						rem = __pms.pop(i)
						pivot.prob += rem.prob
						flag=True
						break
			__new_pms.append(pivot)
		#print('normalization reduced spectrum from %d elements to %d' % (len(self.pms), len(__new_pms)))
		self.pms=__new_pms
		self.scale=self.getScale()
        
	
	@staticmethod
	def _extractSpectrum(weights, spectrum):
		_tmp = {w:0.0 for w in weights}
		_other_count=0
		_other_prob = 0.0
		for pm in spectrum.pms:
			_k=round(pm.mass)
			if _k in _tmp:
				_tmp[_k] += pm.prob
			else:
				_other_count += 1
				_other_prob += pm.prob
		
		_vec=[]
		for w in weights:
			_vec.append(spectrum.scale*_tmp[w])
		_vec=np.array(_vec)

		return _vec, _other_count, _other_prob
		
	
	@staticmethod
	def printSpectra(spectra, w_min, w_max, spectrum_meas=None, prefix='', f_txt=None, f_csv=None):
		TOPN=10 if len(spectra)>=10 else len(spectra)
		
		weights=range(w_min, w_max+1)
		norms_title= ('RMS', 'L_inf', 'L_1')
		sort_by_indx=0 # sort by RMS
		
		print_normst = lambda x: '  '.join(['%8s' % _t for _t in x])
		print_norms  = lambda x: '  '.join(['%8.4f' % _w for _w in x])
		print_weights = lambda x: '  '.join(['%5.1f' %  _w for _w in x])
		print_normst_csv = lambda x: ','.join(['"%s"' % _t for _t in x])
		print_norms_csv  = lambda x: ','.join(['"%f"' % _w for _w in x])
		print_weightst_csv = lambda x: ','.join(['"%d"' %  _w for _w in x])
		print_weights_csv = lambda x: ','.join(['"%f"' %  _w for _w in x])
		name_len=len(spectra[0].comb.getCombString())
		
		if f_txt:		
			print('%-*s:  %s' % (name_len, 'Combination', print_normst(norms_title)), end='', file=f_txt)
			for w in weights:
				print('  %5d' % w, end='', file=f_txt)
			print('', file=f_txt)
		if f_csv:		
			print('"%s",%s,%s,"remarks"' % ('Combination', print_normst_csv(norms_title), print_weightst_csv(weights)), file=f_csv)
		print(prefix + 'Top %d of %d combinations.' % (TOPN, len(spectra)))
		if f_txt or f_csv:
			print(prefix + '  detailed lists are in:')
			if f_txt:			
				print(prefix + '    ' + f_txt .name)
			if f_csv:
				print(prefix + '    ' + f_csv.name)
			print('')
		print(prefix + '%-*s:  %s' % (name_len, 'Combination', print_normst(norms_title)))
		
		if spectrum_meas:
			if f_txt:
				print('%-*s:  ' %(name_len, 'Measurements'), end='', file=f_txt)
			_tmp = {w:0.0 for w in weights}
			for pm in spectrum_meas:
				_k=round(pm.mass)
				if _k in _tmp:
					_tmp[_k] += pm.prob
			
			vec_meas=[]
			for w in weights:
				vec_meas.append(_tmp[w])
			vec_meas=np.array(vec_meas)
			vec_meas = vec_meas*100.0/np.max(vec_meas)
			
			_empty_norm_str = print_normst(['',]*len(norms_title))
			if f_txt:
				print('%s  ' % _empty_norm_str, end='', file=f_txt)
				print(print_weights(vec_meas), file=f_txt)
			print(prefix+'%-*s:  %s' %(name_len, 'Measurements', _empty_norm_str))
			

		results = []
		for spectrum in spectra:
			_vec, _other_count, _other_prob = Spectrum._extractSpectrum(weights, spectrum)
			_norm_rms = np.asscalar(np.sqrt(np.sum((_vec-vec_meas)**2.0)))
			_norm_inf = np.asscalar(np.max(np.abs(_vec-vec_meas)))
			_norm_max = np.asscalar(np.sqrt(np.sum(np.abs(_vec-vec_meas))))
			
			_name = spectrum.comb.getCombString()
			_norms = [_norm_rms, _norm_inf, _norm_max]
			
			if not len(_norms) == len(norms_title):
				raise ValueError('Mismatch in norm titles and implementation')
			
			spectrum.metrics_name = norms_title
			spectrum.metrics = _norms
			
			_peaks = [_m for _m in np.nditer(_vec)]
			 
			if _other_count > 0:
				_other_str = '  %d more with %8.3f' % (_other_count, spectrum.scale*_other_prob)
			else:
				_other_str = ''
			
			results.append( (_name, _norms, _peaks, _other_str) )
		
		if f_txt:
			for _name, _norms, _peaks, _other_str in sorted(results, key=lambda x:x[1][sort_by_indx]):
				print('%-*s:  %s  %s  %s' % (name_len, _name, print_norms(_norms), print_weights(_peaks), _other_str), file=f_txt)
		if f_csv:
			for _name, _norms, _peaks, _other_str in sorted(results, key=lambda x:x[1][sort_by_indx]):
				print('"%s",%s,%s,"%s"' % (_name, print_norms_csv(_norms), print_weights_csv(_peaks-vec_meas), _other_str), file=f_csv)
		for _name, _norms, _peaks, _other_str in sorted(results, key=lambda x:x[1][sort_by_indx])[:TOPN]:
				print(prefix+'%-*s:  %s' % (name_len, _name, print_norms(_norms)))
	
	def getTopN(self,N=4):
		ret = sorted(self.pms, reverse=True)
		return ret[0:N]
	
	def getScale(self):
		return 100.0/max([p.prob for p in self.pms])
	
	def __str__(self):
		scale=self.getScale()
		spec_scaled=[s.getScaledString(scale) for s in sorted(self.getTopN(9), key=lambda x: x.mass)]
		return '%s: %s' % (self.comb.getCombString(space=True), str(spec_scaled))

class CombinationResult:
	def __init__(self, total, weight=None, N=[]):
		self.N = copy.copy(N)
		self.total=total
		if weight:
			self.weight=copy.copy(weight)
		else:
			self.weight = 0.0
			for i in range(len(self.N)):
				self.weight += self.N[i] * Chem.getBlockMasses(i)
	
	@staticmethod
	def getDifference(a, b):
		newN = [na-nb for na, nb in zip(a.N, b.N)]
		return CombinationResult(a.total, N=newN)


	def buildIsotopeTree(self, cutoff):
		elements=self.getElements()
		tree=puzzlotope.probability.TreeNode()
		for elem in elements:
			probs, labels = Chem.getElement(elem).getProbsAndLabels()
			tree.addLevel(probs, labels)
			tree.cutoff(cutoff)
		return tree
	
	def spawn(self, newN, verbose):
		level = len(self.N)
		ret = CombinationResult(self.total, self.weight, self.N)
		ret.weight += Chem.getBlockMasses(level)*newN
		ret.N.append(newN)
		if verbose:
			print('spawning %s => %s' % (self, ret))
		
		return ret
	
	def iterate(self, verbose):
		if self.isComplete():
			return
		else:
			level = len(self.N)
			maxN = math.ceil((self.total-self.weight)/Chem.getBlockMasses(level))+1
			
			if verbose:
				print('Iterating level %d, N=%s, w=%f' % (level, str(self.N), self.weight))
			for i in range(maxN+1):
				if verbose:
					print('  trying %d for level %d' % (i, level))
				yield self.spawn(i, verbose)
		return
	
	def getElements(self):
		elements=[]
		for i in range(len(self.N)):
			if self.N[i]>0:
				elements += Chem.getBlock(i).getElements()*self.N[i]
		return sorted(elements)
	
	def isComplete(self):
		return len(self.N) >= Chem.getNumBlocks()

	""" Test if combination result is a valid result """
	def isValid(self, tolerance):
		reasons = []
		
		""" Combination has to be complete """
		if not self.isComplete():
			return False, reasons
		
		""" Combined weight has to match total weight up to tolerance """
		if self.getRes() > tolerance:
			return False, reasons
		
		""" At lease one required block has to be present in combination """
		flag=False
		for _N, _block in zip(self.N, Chem.Blocks):
			if _N > 0 and (not _block.hasTag(Chem.BlockTag.optional)):
				flag=True
				break
		if not flag:
			reasons.append('Only contains optional blocks')
		
		""" Charge numbers has to be -1, to simplify we check whether -1 is within the possible interval """
		charges_min, charges_max = self.getCharges()
		if charges_min>-1 or charges_max<-1:
			reasons.append('Charge value has to be -1, but is between %d and %d' % (charges_min, charges_max))
		
		""" if there is no neglection reason so far, assume combination is valid """
		return (len(reasons)==0), reasons
	
	def getCharges(self):
		_c_min=0
		_c_max=0
		_num_oxydized=0
		_indx_O=None
		
		for _N, _block, i in zip(self.N, Chem.Blocks, range(Chem.getNumBlocks())):
			# treat oxygen separately
			if _block.name == 'O':
				_indx_O=i
				continue
			_c_min += _N*_block.getMinCharge()
			_c_max += _N*_block.getMaxCharge()
			_num_oxydized += _N if _block.hasTag(Chem.BlockTag.canOxydize) else 0
		
		if _indx_O is not None:
			# there is oxygen, compute its charge
			_N=self.N[_indx_O]
			_block=Chem.Blocks[_indx_O]
			rem=_N-_num_oxydized
			_c_min += rem*_block.getMinCharge() # + 0*_num_oxydized
			_c_max += _*_block.getMaxCharge()
		
		return _c_min, _c_max
		
	def __lt__(self, obj):
		if self.getLevel() < obj.getLevel():
			return True
		elif self.getLevel() > obj.getLevel():
			return False
		else:
			return (self.getRes() < obj.getRes())

	def getLevel(self):
		return len(self.N)

	def getRes(self):
		return abs(self.total-self.weight)

	def __repr__(self):
		return '[result (%s): res=%8.4f, w=%8.4f, N=%s]' % (('complete' if self.isComplete() else 'incomplete'), self.total-self.weight, self.weight, str(self.N))
	
	def getStringLen(self, space=True):
		if space:
			return len(self.N)*5 - 1 + sum([len(b.toSymbol()) for b in Chem.Blocks])
		else:
			return len(self.getCombString(spcae))
	
	def getCombString(self, space=True):
		string=[]
		for i in range(Chem.getNumBlocks()):
			_N = self.N[i]
			if space:
				if _N > 0:
					string.append('%3d %s' % (_N, Chem.Blocks[i].toSymbol()))
				else:
					string.append('%3s %s' % ('', ' '*len(Chem.Blocks[i].toSymbol())))
			else:
				if _N > 0:
					string.append('%d %s' % (_N, Chem.Blocks[i].toSymbol()))
		return ' '.join(string)
	
	@staticmethod
	def filter(combinations, tolerance, file_out=sys.stdout):
		print('Filtering %d combination resuls' % len(combinations), file=file_out)
		_valid=[]
		_invalid_count=0
		_neglected=[]
		_neglected_reasons=[]
		for comb in combinations:
			valid, reasons = comb.isValid(tolerance)
			if valid:
				_valid.append(comb)
			elif reasons:
				_neglected.append(comb)
				_neglected_reasons.append(' AND '.join(reasons))
			else:
				_invalid_count += 1
		print('  %d invalid combinations' % _invalid_count, file=file_out)
		print('  %d neglected:' % len(_neglected), file=file_out)
		for _n, _reason  in zip(_neglected, _neglected_reasons):
			print('	%s: %s' % (_n.getCombString(), _reason), file=file_out)
		print('  %d combinations remaining:' % len(_valid), file=file_out)
		for _c in _valid:
			print('	' + _c.getCombString(), file=file_out)
		return _valid
	
if __name__ == '__main__':
	Chem.update()
	__N=([0,]*Chem.getNumBlocks())
	__N[-1]=2
	cr = CombinationResult(415, weight=415, N=__N)
	print(Spectrum.computeSpectrum(cr,0.3,1e-7))
