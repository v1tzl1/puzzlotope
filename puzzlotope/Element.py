class Element:
    symbol=None
    isotopes=[]
    
    def __init__(self, symbol, isotopes):
        self.symbol=symbol
        
        masses = [x[0] for x in isotopes]
        probs  = [x[1] for x in isotopes]
        sort_indx = sorted(range(len(probs)), key=probs.__getitem__, reverse=True)
        
        self.isotopes=[]
        for i in range(len(isotopes)):
            self.isotopes.append(ProbMass(probs[sort_indx[i]], masses[sort_indx[i]]))
        
        
    def getMass(self):
        return self.isotopes[0].mass
    
    def getMasses(self):
        return self.isotopes
    
    def getProbsAndLabels(self):
    	probs=[]
    	labels=[]
    	for pm in sorted(self.isotopes):
    		probs.append(pm.prob)
    		labels.append(self.symbol + str(round(pm.mass)))
    	return probs, labels

class ProbMass:
    prob=None
    mass=None
    
    def __init__(self, prob=1.0, mass=0.0):
        self.prob=prob
        self.mass=mass
    
    def addProbMass(self, num, add, cutoff):
        combs = combinations.getNofM(num, len(add))
        for comb in combs:
            __mass = self.mass
            __probs = self.prob
            for i in range(len(comb)):
        	    for j in range(comb[i]):
        		    __mass += add[i].mass
        		    __probs *= add[i].prob
            if __probs > cutoff:
                yield ProbMass(__probs, __mass)
        return
        
    def getScaledString(self, scale):
        return 'm=%8.4f (%6.2f)' % (self.mass, self.prob*scale)
    	
    def __repr__(self):
        return 'm=%8.4f (%6.2f%%)' % (self.mass, 100.0*self.prob)
    
    @staticmethod
    def merge(a,b, num_b, cutoff):
        merged = []
        for pma in a:
            merged += list(pma.addProbMass(num_b, b, cutoff))
        
        return merged
    
    def getMassDiff(self, b):
        return abs(self.mass-b.mass)
        
    def __lt__(self, b):
        return (self.prob < b.prob)
        #return (self.mass < b.mass)

