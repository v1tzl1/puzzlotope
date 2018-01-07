import re
import puzzlotope.Input
import puzzlotope.combinations
import enum

Elements=None
Blocks=None

regexp_elem = re.compile(r"([A-Z][a-z]*)([0-9]*)")

def update():
    global Elements
    global Blocks
    
    if Elements is None:
        Elements = {e.symbol: e for e in puzzlotope.Input.Elements}
    
    if Blocks is None:
        Blocks = puzzlotope.Input.Blocks

def getBlockMasses(i=None):
    global Blocks
    if i is None:
        return [b.getMass() for b in Blocks]
    else:
        return Blocks[i].getMass()

def getNumBlocks():
    global Blocks
    return len(Blocks)

def getBlock(i):
	global Blocks
	return Blocks[i]

def getIsotopeMass(iso_str):
	global regexp_elem
	m = regexp_elem.match(iso_str)
	
	if not m:
		raise ValueError('Cannot find isotope %s' % iso_str)
	
	elem = getElement(m.group(1))
	if m.group(2) == '':
		raise ValueError('Cannot determine Isotope Mass from %s' % iso_str)
	num=int(m.group(2))
	for m in elem.getMasses():
		if round(m.mass)==num:
			return m.mass
	raise ValueError('Isotope %s is not defined' % iso_str)

def getElement(sym):
    global Elements
    update()
    
    return Elements[sym]

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


class BlockTag(enum.Enum):
	""" indicate that this block is a byproduct. Combinations of only optional blocks are discarded """
	optional=1
	
	""" blocks with this tag can insert with oxygen such that the oxygen has a charge value of 0 instead of -2 """
	canOxydize=2
	
	
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
        
class Block:
    name=None
    elements=[]
    elements_num=[]
    charges=[]
    tags=set()
    
    def __init__(self, name, charges, tags=[]):
        global regexp_elem

        """ Split name into elements """
        match = regexp_elem.findall(name)
    
        if not match:
            raise 'Cannot parse \'%s\'' % name
        
        #print('Parsing %s' % name)
        #print(match)
        
        self.elements = []
        self.elements_num = []
        
        for element, number in match:
            number = int(number) if not number == '' else 1
            
            #print('Found %2d %s' % (number, element))
            
            self.elements.append(element)
            self.elements_num.append(number)

        """ Set other options """
        self.tags=set(tags)
        self.charges=charges
    
    def hasTag(self, t):
    	return (t in self.tags)
    
    def getMass(self):
        mass = 0
        for i in range(len(self.elements)):
            mass += self.elements_num[i]*getElement(self.elements[i]).getMass()
        return mass
    
    def getMasses(self, cutoff):
        pms=[ProbMass(1.0,0.0), ]
        for i in range(len(self.elements)):
            pms_tmp=[]
            for pm in pms:
                pms_tmp+= list(pm.addProbMass(self.elements_num[i], getElement(self.elements[i]).getMasses(), cutoff))
            pms=pms_tmp
        
        return pms
     
    def getElements(self):
    	elements=[]
    	for i in range(len(self.elements)):
    		for j in range(self.elements_num[i]):
    			elements.append(self.elements[i])
    	return elements 
    
    def getMinCharge(self):
        return min(self.charges)
    
    def getMaxCharge(self):
        return max(self.charges)
    
    def toSymbol(self):
        return ''.join([self.elements[i]+(str(self.elements_num[i]) if self.elements_num[i] > 1 else '') for i in range(len(self.elements))])
    
    def __repr__(self):
        return self.toSymbol() + (' (required)' if self.required else '')
