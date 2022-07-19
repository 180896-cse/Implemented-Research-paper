#! /bin/python
import random
import sys

class CSP:
	def __init__(self, n, cst, prc, cor, mem, stor, b = 1, perf=1, ra=1, rel=1):
		self.name = n
		self.cost = cst
		self.price = prc
		self.cores = cor
		self.memory = mem
		self.storage = stor
		self.brandValue = b
		self.Qp = perf
		self.Qra = ra
		self.Qrel = rel
		self.Qscal_c = 0
		self.Qscal_m = 0
		self.Qscal_s = 0
		self.x = [0,0,0]
		self.satc = 0
		self.satm = 0
		self.sats = 0
	
	def getName(self):
		return self.name	
		
	def getCost(self):
		return self.cost
		
	def getPrice(self):
		return self.price
		
	def getCores(self):
		return self.cores
		
	def getMemory(self):
		return self.memory
		
	def getStorage(self):
		return self.storage
		
	def getBrandValue(self):
		return self.brandValue
		
	def getQp(self):
		return self.Qp
		
	def getQra(self):
		return self.Qra
		
	def getQrel(self):
		return self.Qrel
		
	def getQ(self):
		return (self.Qp + self.Qra + self.Qrel)/3.0
		
	def getQscal(self):
		return (self.Qscal_c + self.Qscal_m + self.Qscal_s)/3.0
			
	def getX(self):
		return self.x
		
	def getSat(self):
		return (self.satc + self.satm + self.sats)/3.0
	
	def getTrust(self):
		return self.getQrel()
	
		

SP = []
SP.append(CSP("US East (N. Virginia)", [0.0665, 0.133, 0.266], [0.133, 0.266, 0.532], 32, 120, 640, 0.678, 0.689, 0.687, 0.692))
SP.append(CSP("US West (Northern California)", [0.077, 0.154, 0.308], [0.154, 0.308, 0.616], 32, 120, 640, 0.786, 0.736, 0.732, 0.739))
SP.append(CSP("EU (Ireland)", [0.073, 0.146, 0.293], [0.146, 0.293, 0.585], 32, 120, 640, 0.745, 0.717, 0.713, 0.725))
SP.append(CSP("Asia Pacific (Singapore)", [0.098, 0.196, 0.392], [0.196, 0.392, 0.784], 32, 120, 640, 1.000, 0.931, 0.924, 0.927))
SP.append(CSP("Asia Pacific (Tokyo)", [0.0865, 0.193, 0.385], [0.193, 0.385, 0.77], 32, 120, 640, 0.985, 0.889, 0.882, 0.885))
SP.append(CSP("Asia Pacific (Sydney)", [0.093, 0.186, 0.372], [0.186, 0.372, 0.745], 32, 120, 640, 0.949, 0.824, 0.817, 0.813))
SP.append(CSP("South America (Sao Paulo)", [0.095, 0.19, 0.381], [0.19, 0.381, 0.761], 32, 120, 640, 0.969, 0.858, 0.854, 0.851))
SP.append(CSP("AWS GovCloud (US)", [0.084, 0.168, 0.336], [0.168, 0.336, 0.672], 32, 120, 640, 0.857, 0.793, 0.785, 0.787))
#Untrusted CSPs 
for i in range(8):
	SP.append(CSP("Untrusted CSP "+str(i+1), [0.155*random.random(), 0.309*random.random(), 0.618*random.random()], [0.155, 0.309, 0.618], 32, 120, 640, 0.39+0.1*random.random(), 0.39+0.1*random.random(), 0.39+0.1*random.random(), 0.39+0.1*random.random()))
	
GrandFederation = SP
	
def powerset(seq):
    """
    Returns all the subsets of this set. This is a generator.
    """
    if len(seq) <= 1:
        yield seq
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item
            

def Banzhaf(spi, GrandFederation, request):
	m = len(GrandFederation)
	newset = [t for t in GrandFederation if t != spi]
	Fspace = [x for x in powerset(newset)]
	total = sum([v(F+[spi], request) - v(F, request) for F in Fspace])
	return total/2**(m-1)
	
def NormalizedBanzhaf(spi, GrandFederation, request):
	vspi = Banzhaf(spi, GrandFederation, request)
	vall = sum([Banzhaf(tsp, GrandFederation, request) for tsp in GrandFederation])
	if vall>0:
		return vspi/vall
	else:
		return 0

from lpsolve55 import *

def v(F, request):
	if len(F)==0:
		return 0
		
	lp = lpsolve('make_lp', 0, 3*len(F))
	lpsolve('set_verbose', lp, IMPORTANT)
	ret = lpsolve('set_obj_fn', lp, [(spi.getPrice()[0] - spi.getCost()[0]) for spi in F]+[(spi.getPrice()[1] - spi.getCost()[1]) for spi in F]+[(spi.getPrice()[2] - spi.getCost()[2]) for spi in F])
	vector = [0]*len(F)*3
	for i in range(len(F)):
		vector[i] = 2
		vector[len(F)+i] = 4
		vector[2*len(F)+i] = 8
		ret = lpsolve('add_constraint', lp, vector, LE, F[i].getCores())
		vector[i] = 0
		vector[len(F)+i] = 0
		vector[2*len(F)+i] = 0
	for i in range(len(F)):
		vector[i] = 7.5
		vector[len(F)+i] = 15
		vector[2*len(F)+i] = 30
		ret = lpsolve('add_constraint', lp, vector, LE, F[i].getMemory())
		vector[i] = 0
		vector[len(F)+i] = 0
		vector[2*len(F)+i] = 0 	
	for i in range(len(F)):
		vector[i] = 32
		vector[len(F)+i] = 80
		vector[2*len(F)+i] = 160
		ret = lpsolve('add_constraint', lp, vector, LE, F[i].getStorage())
		vector[i] = 0
		vector[len(F)+i] = 0
		vector[2*len(F)+i] = 0
		
	for j in range(3):
		for i in range(len(F)):
			vector[j*len(F)+i] = 1
				
		ret = lpsolve('add_constraint', lp, vector, EQ, request[j])
			
		for i in range(len(F)):
			vector[j*len(F)+i] = 0
				
	for k in range(3*len(F)):
		vector[k] = 1
		ret = lpsolve('add_constraint', lp, vector, GE, 0)
		vector[k] = 0
		
	for i in range(len(F)):
		for j in range(3):
			vector[j*len(F)+i] = 1
		ret = lpsolve('add_constraint', lp, vector, GE, 1)
		for j in range(3):
			vector[j*len(F)+i] = 0
			
	for i in range(3*len(F)):
		ret = lpsolve('set_int', lp, i+1, 1)
	
	ret = lpsolve('set_maxim', lp)	
	ret = lpsolve('solve', lp)
	
	if(ret != 0):
		return 0
		
	obj = lpsolve('get_objective', lp)
	allocation = lpsolve('get_variables', lp)[0]
	for i in range(len(F)):
		for j in range(3):
			F[i].x[j] = allocation[j*len(F)+i]
		F[i].satc= (2*allocation[i] + 4*allocation[len(F)+i] + 8*allocation[2*len(F)+i])/F[i].getCores()
		F[i].satm= (7.5*allocation[i] + 15*allocation[len(F)+i] + 30*allocation[2*len(F)+i])/F[i].getMemory()
		F[i].sats= (32*allocation[i] + 80*allocation[len(F)+i] + 160*allocation[2*len(F)+i])/F[i].getStorage()
		F[i].Qscal_c = 1-F[i].satc
		F[i].Qscal_m = 1-F[i].satm
		F[i].Qscal_s = 1-F[i].sats
		
					
	lpsolve('delete_lp', lp)
			
	return obj	 
		


def MergeCompare(union, part1, part2, request):
	return (v(union, request)>v(part1, request) and v(union, request)>v(part2, request))

def SplitCompare(union, part1, part2, request):
	return (v(part1, request)>=v(union, request) or v(part2)>=v(union, request))


def removeDuplicates(L):
	temp = []
	for item in L:
		if item not in temp:
			temp.append(item)
	return temp

def MergeFederations(FS, checked, request):
	#checked = []
	while True:
		feasible = [x for x in FS if x not in checked]
		#print "Feasible list for MergeFederations:"
		#print [sp.getName() for x in feasible for sp in x]
		if len(feasible) == 0:
			return
		mergeFound = False
		
		for (fi, fj) in [(a,b) for a in feasible for b in feasible if a!=b]:
			checked.append(tuple(set(fi+fj)))
			checked = list(set(checked))
			#for f in checked:
				#print [sp.getName() for sp in f]
			if MergeCompare(list(set(list(fi)+list(fj))),fi,fj, request):
				#fi = list(fi)+list(fj)
				#fi = set(fi)
				#FS.remove(fi)
				FS.remove(fj)
				FS.append(list(set(list(fi)+list(fj))))
				mergeFound = True
				#print "Merge found"
				break
		if not mergeFound:
			break

def SplitFederations(FS, checked, request):
	splitFound = False
	for fi in [x for x in FS if len(x)>1]:
		for (fj, fk) in [(x,y) for x in powerset(fi) for y in powerset(fi) if x+y==fi]:
			checked.append(tuple(fj))
			checked.append(tuple(fk))
			#checked += [fj, fk]
			checked = list(set(checked))
			if SplitCompare(fi,fj,fk, request):
				#fi = fj
				FS.append(fj)
				FS.append(fk)
				FS.remove(fi)
				#FS=set(FS)
				FS = removeDuplicates(FS)
				splitFound = True
				#print "Split found"
				break
	return splitFound
	
def getQp(F):
	if sum([(ci.getX()[0]+2*ci.getX()[1]+4*ci.getX()[2]) for ci in F]) == 0:
		return 0
	return sum([(ci.getX()[0]+2*ci.getX()[1]+4*ci.getX()[2])*ci.getQp() for ci in F]) / sum([(ci.getX()[0]+2*ci.getX()[1]+4*ci.getX()[2]) for ci in F])
	
def getQra(F):
	product = 1
	for ci in F:
		product *= 1 - ci.getQra()
	return 1 - product
	
def getQrel(F):
	return sum([ci.getQrel() for ci in F])/len(F)
	
def getQscal(F):
	scal_c = sum([ci.getCores() - 2*ci.getX()[0]-4*ci.getX()[1]-8*ci.getX()[2] for ci in F])/sum([ci.getCores() for ci in F])
	scal_m = sum([ci.getMemory() - 7.5*ci.getX()[0]-15*ci.getX()[1]-30*ci.getX()[2] for ci in F])/sum([ci.getMemory() for ci in F])
	scal_s = sum([ci.getStorage() - 32*ci.getX()[0]-80*ci.getX()[1]-160*ci.getX()[2] for ci in F])/sum([ci.getStorage() for ci in F])
	return (scal_c + scal_m + scal_s)/3.0
	
def getQ(F):
	return (getQp(F) + getQra(F) + getQrel(F))/3.0
	
def sat(v, vmax, vmin):
	if vmin<=v and v <= vmax:
		return 0.01 + 0.99*abs((v-vmin)/(vmax-vmin))
	elif v > vmax:
		return 1
	else:
		return 0

def getFSPsat(fd, spi, request):
	if len(fd) == 0:
		return 0
	norm = NormalizedBanzhaf(spi, fd, request)
	profit = norm*totalprofit
	return 0.5*sat(profit, 2*v(fd, request)/len(fd), 0.01*v(fd, request)/len(fd)) + 0.5*sat(getQra(fd), 1, 0.2)
	
def getFsat(fd, request):
	if len(fd) == 0:
		return 0
	return 0.5*sat(v(fd, request), sum([2*sum([spi.getX()[j]*spi.getCost()[j] for j in range(3)]) for spi in fd])/len(fd), 0.01*sum([sum([spi.getX()[j]*spi.getCost()[j] for j in range(3)]) for spi in fd])/len(fd)) + 0.5*sat(getQra(fd), 1, 0.2)
	
		
checked = []
request = [int(arg) for arg in sys.argv[1:]]
if(len(request))!= 3:
	print "Enter a 3-component vector for the request"
	quit()
	
import time

time1=time.clock()
print "Mechanism\tAvg F.sat\tAvg F.SP.sat\tAvg F.Q\tAvg F.A\tAvg F.Profit\tAvg F.SP.Profit\tAvg F.Size\tExecution Time(s)"
fspprofit = 0
fspsat = 0
fsize = 0
fsat = 0
fq = 0
fa = 0
fprofit = 0
count = 0


while True:

	FS = [[x] for x in GrandFederation]

	for i in range(1):
		MergeFederations(FS, checked, request)
		if not SplitFederations(FS, checked, request):
			break

	maxval = 0
	finalfed = []
	#print "Final Federation Structure:"
	for fed in FS:
		#print [spi.getName() for spi in fed]
		if v(fed, request) > maxval:
			finalfed = fed
			maxval = v(fed, request)

	#allocate resources using finalfed

	totalprofit = v(finalfed, request)
	if totalprofit == 0:
		break
	
	#print "---FINAL SELECTED FEDERATION---"
	#print [spi.getName() for spi in finalfed]
	#print "Number of requests for VM instances: ", request
	#print "Total profit of the final federation: ", totalprofit

	#print "request[0]\trequest[1]\trequest[2]\tCSP\tCSP.Q\tCSP.Qs\tCSP.x[0]\tCSP.x[1]\tCSP.x[2]\tProfit of CSP\tCSP.sat\tF.Q\tF.Qs\tProfit of Federation\tFD.SP.sat\tFD.sat\tFD.A\tRunning Time (s)"
	
	for ci in finalfed:
		X = list(ci.getX())
		norm = NormalizedBanzhaf(ci, finalfed, request)
		#print "Banzhaf value of CSP %s: " % spi.getName(), Banzhaf(spi, finalfed, request)
		#print "Normalized Banzhaf value of CSP %s: " % spi.getName(), norm
		#print "Share of CSP %s in the profit of the final federation: " % spi.getName(), norm*totalprofit
		profit = norm*totalprofit
		fspprofit = fspprofit + profit
		fspsat = fspsat + getFSPsat(finalfed, ci, request)
		#print "%d\t%d\t%d\t%s\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (request[0], request[1], request[2],  ci.getName(), ci.getQ(), ci.getQscal(), X[0], X[1], X[2], profit, ci.getSat(), getQ(finalfed), getQscal(finalfed), totalprofit, getFSPsat(finalfed, ci, request), getFsat(finalfed, request), getQra(finalfed), runtime)
		GrandFederation.remove(ci)
		
	count = count+1
	fsize = fsize + len(finalfed)
	fsat = fsat + getFsat(finalfed, request)
	fq = fq + getQ(finalfed)
	fa = fa + getQra(finalfed)
	fprofit = fprofit + totalprofit
	
time2 = time.clock()
runtime = time2-time1

fspprofit = fspprofit / fsize
fspsat = fspsat / fsize
fsat = fsat / count
fq = fq / count
fa = fa / count
fprofit = fprofit / count
fsize = fsize / count

print "Lena (50%% Untrusted)\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" %(fsat, fspsat, fq, fa, fprofit, fspprofit, fsize, runtime)



