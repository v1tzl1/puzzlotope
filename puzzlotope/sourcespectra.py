import puzzlotope.Solver
import puzzlotope.Chem
import numpy as np
import sys
import puzzlotope.multinomial

def buildSpectra(combinations, tolerance, cutoff, num_threads, prefix='', verbose=False):
	num_elements = dict()
	out = []
	count_dicts = []
	
	for comb in combinations:
		_tmp_dic = dict()
		for _e in comb.getElements():
			if not _e in _tmp_dic:
				_tmp_dic[_e] = 0
			_tmp_dic[_e] = _tmp_dic[_e]+1
		count_dicts.append(_tmp_dic)
		
		for _e, _c in _tmp_dic.items():
			if _e not in num_elements:
				num_elements[_e] = set()
			num_elements[_e].add(_c)
	
	isotope_tree = dict()
	for _e, _v in num_elements.items():
		isotope_tree[_e] = buildElementTree(_e, _v, cutoff)
	
	for i in range(len(combinations)):
		_tmp_tree = getEmptyTree()
		for _e, _c in count_dicts[i].items():
			_tmp_tree = mergeTrees(_tmp_tree, isotope_tree[_e][_c], cutoff)
		
		print(combinations[i])
		out.append(puzzlotope.Solver.Spectrum(combinations[i], _tmp_tree[0], _tmp_tree[1], tolerance, cutoff))
		#print(_tmp_tree)
		#	mergeTrees(_tmp_tree, isotope_tree[]
	#print(mergeTrees(isotope_tree['C'][5], isotope_tree['O'][3]))
	#sys.exit(0)
	return out

def mergeTrees(a, b, cutoff):
	_a_probs, _a_weights, _a_labels = a
	_b_probs, _b_weights, _b_labels = b
	
	_out_probs   = []
	_out_weights = []
	_out_labels  = []
	for i in range(len(_a_probs)):
		for j in range(len(_b_probs)):
			if _a_probs[i]*_b_probs[j] >= cutoff:
				_out_probs.append(_a_probs[i]*_b_probs[j])
				_out_weights.append(_a_weights[i]+_b_weights[j])
				_out_labels.append(_a_labels[i]+' '+_b_labels[j])
	
	return _out_probs, _out_weights, _out_labels

def getEmptyTree():
	return [1.0,], [0.0,], ['',]

"""
Element tree:
tree = tuple( probabilities, weights,  labels)
probabilities = tuple(p_1, p_2, ..., p_N)    0 <= p <= 1,    sum(p) = 1.0
weights tuple with N weights
labels tuple with N labels. e.g. "5C12 3C13"
"""

def buildElementTree(element_symb, nums, cutoff):
	element = puzzlotope.Chem.getElement(element_symb)
	probs, labels, masses = element.getProbsAndLabels()
	
	out = dict()
	
	#for num in range(max_num+1):
	for num in sorted(nums):
		coeff_dict = puzzlotope.multinomial.multinomial_coefficients(len(labels), num)
		
		_prbs = []
		_wgts = []
		_lbls = []
		
		for _factors, _num in coeff_dict.items():
			_l = []
			_p = 1.0;
			_w = 0.0;
			for i, _factor in zip(range(len(_factors)), _factors):
				if _factor > 0.0:
					_l.append(str(_factor)+labels[i])
					_p *= (probs[i])**_factor
					_w += _factor*masses[i]
			
			if _num*_p >= cutoff:
				_lbls.append(' '.join(_l))
				_prbs.append(_num*_p)
				_wgts.append(_w)
		
		if len(_prbs) > 0:
			out[num] = tuple(_prbs), tuple(_wgts), tuple(_lbls)
		else:
			out[num] = getEmptyTree()
		#print(str(out))
		
	return out
	
