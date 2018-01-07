import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import puzzlotope.Input
import puzzlotope.Chem as Chem

try:
	from sklearn.mixture import GaussianMixture
except ImportError:
	try:
		from sklearn.mixture import GMM as GaussianMixture
		print('Your version of sklearn is too old.')
	except:
		raise
	raise


def estimateOffset(ma, mb, maxdiff):
	na=np.array(ma)
	nb=np.array(mb)
	ncomb=np.hstack((na, nb))
	masses=np.round(ncomb)
	#masses=np.round(np.concatenate(np.array(ma), np.array(mb)))
	#masses=(np.concatenate(np.round(ma), np.round(mb)))
	masses=np.unique(masses)
	
	diff=[]
	for _m in np.nditer(masses):
		idxa = (np.abs(ma-_m)).argmin()
		idxb = (np.abs(mb-_m)).argmin()
		
		if (np.abs(ma[idxa]-_m) < 0.5) and (np.abs(mb[idxb]-_m) < 0.5):
			diff.append(ma[idxa]-mb[idxb])
	ofst=np.mean(np.array(diff))
	if ofst > maxdiff:
		ofst = maxdiff
	if ofst < -maxdiff:
		ofst = -maxdiff
	#print('estimated offset: %f'%ofst)
	return ofst
	
def estimateSigma(x, y, num_samples=10000, verbose=False):
	xmin=np.ceil(min(x))
	xmax=np.floor(max(x))
	num=int(xmax-xmin+1)
	
	yscaled=np.round(num_samples*(y/sum(y)))
	samples=[]
	for i in range(yscaled.size):
		if yscaled[i]>0:
			samples += [x[i],]*int(yscaled[i])
	ns=np.array(samples)
	
	means_init=np.linspace(xmin, xmax, num).reshape(-1,1)
	gmm = GaussianMixture(n_components=num, covariance_type='tied', means_init=means_init)
	gmm.fit(ns.reshape(-1,1))
	
	sigma=np.asscalar(np.sqrt(gmm.covariances_))
	means=gmm.means_.reshape(-1)
	weights=gmm.weights_.reshape(-1)
	
	if verbose:
		print('GMM fit with %d components:' % num)
		print('  converged: %s'%gmm.converged_)
		print('  num iterations: %d'%gmm.n_iter_)
		print('  means: %s' % means)
		print('  weights: %s' % weights)
		print('  sigma: %f' % sigma)
	
	pms=[Chem.ProbMass(_p, _m) for _m, _p in zip(means, weights)]
		
	return sigma, pms
	
def buildPDF(xmeas, diracs, offset, sigma):
	x = np.array(xmeas)+offset
	
	mus=[d.mass for d in diracs]
	w=np.array([d.prob for d in diracs])
	w=w/sum(w)
	
	pdf=np.zeros(x.shape)
	for i in range(len(mus)):
		pdf += w[i]*norm.pdf(x,loc=mus[i], scale=sigma)
		
	return pdf

def plotPDFs(x, pdfs, labels, title):
	fig, ax = plt.subplots(1, 1)
	for pdf, label in zip(pdfs, labels):
		scale=100.0/max(pdf)
		ax.plot(x, scale*pdf, label=label)
	
	ax.legend(loc='best')
	ax.set_xlim(min(x), max(x))
	plt.title(title)
	plt.show()

if __name__=='__main__':
	pdfs=[]
	labels=[]
	
	x=np.linspace(414.0, 425.0,1000)
	
	diracs=(
		Chem.ProbMass(100, 415.0),
		Chem.ProbMass(24.64, 416.0),
		Chem.ProbMass(79.43, 417.0),
		Chem.ProbMass(22.33, 418.0),
		Chem.ProbMass(28.18, 419.0),
		Chem.ProbMass(7.58, 420.0),
		Chem.ProbMass(4.79, 421.0),
		Chem.ProbMass(1.01, 422.0),
		Chem.ProbMass(0.29, 423.0)
				)
	
	params=(
				#[0.0, 0.2],
				#[0.0, 0.225],
				#[0.0, 0.25],
				#[0.0, 0.275],
				[0.0, 0.3],
				[0.0, 0.325],
				[0.0, 0.35],
				[0.0, 0.375],
				[0.0, 0.4],
				[0.0, 0.425],
				[0.0, 0.45],
				[0.0, 0.475],
				[0.0, 0.5],
				[0.0, 0.525],
				[0.0, 0.55],
				[0.0, 0.575],
				#[0.0, 0.5],
				#[0.0, 0.6],
				#[0.0, 0.8]
					)
	
	for ofst, sigma in params:
		pdfs.append(buildPDF(x, diracs, ofst, sigma))
		labels.append('Delta=%f, sigma=%f' % (ofst, sigma))
	
	plotPDFs(x, pdfs, labels, '2Ni 3C6H5 1C5H8')
