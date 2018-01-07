import puzzlotope.gmm
import re
import numpy as np
#import matplotlib.pyplot as plt

def parseFile(filename, xmin=None, xmax=None, cutoff=0.0):
	linere = re.compile('([0-9]*.[0-9]*)[ 	]*([0-9.]+)')
	x=[]
	y=[]
	
	with open(filename, 'r') as f:
		for line in f:
			m=linere.match(line)
			if m:
				x.append(float(m.group(1)))
				y.append(float(m.group(2)))
			else:
				print('Cannot parse line %s' % line)
	nx = np.array(x)
	ny = np.array(y)
	
	if xmin is None and xmax is None:
		mask = np.ones(nx.shape, dtype=bool)
	elif xmax is None:
		mask = (nx >= xmin)
	elif xmin is None:
		mask = (nx <= xmax)
	else:
		mask = ( (nx <= xmax) & (nx >= xmin) )
	 
	return nx[mask], ny[mask]

def getRMS(a, b):
	return np.asscalar(np.sqrt(np.sum((a-b)**2)))
	 	 	 
if __name__=='__main__':
	filename='spectrum.xy'
	cutoff=0.01
	
	x, y = parseFile(filename, xmin=414.5, xmax=423.5)
	yscaled=100.0/max(y) * y
	mask = yscaled >= cutoff
	
	xmasked = x[mask]
	ymasked = yscaled[mask]
	
	sigma, diracs = gmm.estimateSigma(xmasked, ymasked, verbose=True)
	print()
	print(sigma)
	print(diracs)
	reconstr_pdf = gmm.buildPDF(xmasked, diracs, 0.0, sigma)
	gmm.plotPDFs(xmasked, [ymasked, reconstr_pdf], ['Measurement','Reconstruction'], '')
