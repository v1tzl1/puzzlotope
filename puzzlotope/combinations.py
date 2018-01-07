import copy

def getNofM(N, M):
	prefixes=[[],]
	for i in range(M):
		tmp=[]
		for prefix in prefixes:
			tmp+=list(__getOptions(N, M-i, prefix))
		prefixes=tmp
	return prefixes

def __getOptions(N, remDim, prefix):
	remValues=N-sum(prefix)
	if remDim==1:
		yield __append(prefix, remValues)
		return
	for num in range(remValues+1):
		yield __append(prefix, num)
		
def __append(a,b):
	ret=copy.copy(a)
	ret.append(b)
	return ret

if __name__ == "__main__":
	
	res=getNofM(3, 5)
	print('result:\n%s' % '\n'.join([str(r) for r in res]))