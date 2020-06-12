import os

# self-explanatory convenience functions

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def makedir(d):
	if not os.path.isdir(d):
		os.system('mkdir '+d)

def getInputFileExtension(fn):
	fn = fn.lower()
	isPICKLE = (fn[-2:] == '.p')
	isFASTQ  = (fn[-3:] == '.fq' or fn[-6:] == '.fq.gz' or fn[-6:] == '.fastq' or fn[-9:] == '.fastq.gz')
	isFASTA  = (fn[-3:] == '.fa' or fn[-6:] == '.fa.gz' or fn[-6:] == '.fasta' or fn[-9:] == '.fasta.gz')
	isGZIP   = (fn[-3:] == '.gz')
	return (isPICKLE, isFASTQ, isFASTA, isGZIP)

RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join([RC_DICT[n] for n in s[::-1]])

def read_fq_entry(fq_file):
	myName = fq_file.readline().strip()[1:]
	if not myName:
		return ('','','')
	myRDat = fq_file.readline().strip()
	skip   = fq_file.readline().strip()
	myQDat = fq_file.readline().strip()
	return (myRDat, myQDat, myName)

def read_fa_entry(fa_file, dummyVal=120):
	myName = fa_file.readline().strip()[1:]
	if not myName:
		return ('','','')
	myRDat = fa_file.readline().strip()
	myQDat = ''.join([chr(dummyVal) for n in myRDat])	# dummy value
	return (myRDat, myQDat, myName)

def print_nicely(d,order=None,spacing=2,end_boof=0,custom_tot=None):
	b1 = max([len(str(n)) for n in d.keys()])+spacing
	b2 = max([len(str(n)) for n in d.values()])+spacing
	if custom_tot == None:
		tot = float(sum(d.values()))
	else:
		tot = float(custom_tot)
	if order == None: order = sorted(d.keys())
	for k in order:
		outstr = str(k)+':'+' '*(b1-len(str(k)))+str(d[k])+' '*(b2-len(str(d[k])))+'({0:0.2f}%)'.format(100.*(d[k]/tot))
		print(outstr+' '*(max([0,end_boof-len(outstr)])))

def l_2_cd(l):
	"""
	count how many times each item in a list occurs

	l - list of stuff
	returns: dictionary where keys are list elements and values their count
	"""
	cd = {}
	for n in l:
		if n not in cd:
			cd[n] = 0
		cd[n] += 1
	sorted_res = sorted([(cd[k],k) for k in cd.keys()],reverse=True)
	return [(n[1],n[0]) for n in sorted_res]

# longest common subsequence
def lcs(s1,s2):
	m = len(s1)
	n = len(s2)
	L = [[None]*(n + 1) for i in range(m + 1)] 
	for i in range(m + 1): 
		for j in range(n + 1): 
			if i == 0 or j == 0 : 
				L[i][j] = 0
			elif s1[i-1] == s2[j-1]: 
				L[i][j] = L[i-1][j-1]+1
			else: 
				L[i][j] = max(L[i-1][j], L[i][j-1]) 
	return L[m][n] 

if __name__ == '__main__':
	pass
	