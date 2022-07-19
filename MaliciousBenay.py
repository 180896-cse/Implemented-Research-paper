#! /bin/python
import random
import sys
import math

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

	
from lpsolve55 import *

wquality = 0.5

wprofit = 0.5

ThresholdTrust = 0.5

def getBV(F):
	if len(F) == 0:
		return 0
	return max([ci.getBrandValue() for ci in F])
	
def getQp(F):
	#if sum([(ci.getX()[0]+2*ci.getX()[1]+4*ci.getX()[2]) for ci in F]) == 0:
		#return 0
	#return sum([(ci.getX()[0]+2*ci.getX()[1]+4*ci.getX()[2])*ci.getQp() for ci in F]) / sum([(ci.getX()[0]+2*ci.getX()[1]+4*ci.getX()[2]) for ci in F])
	return sum([ci.getQp() for ci in F])/len(F)
		
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
		
#print [C.getName() for C in sorted(SP, key = lambda C : C.cost[0])]
		
def getRCost(fd, j):
	N = len(fd)
	temp =	sorted(fd, key = lambda C : C.cost[j])
	s = temp[0].getQra()*temp[0].cost[j]
	for i in range(1, N):
		addend = 1
		for k in range(i):
			addend = addend*(1-temp[k].getQra())*temp[i].getQra()*temp[i].cost[j]
		s = s+addend
	return s
	
def getRPrice(fd, j):
	return getQra(fd)*sum([C.price[j] for C in fd])/len(fd)
	 		
def getRProfit(fd, j):
	return  (4/len(fd))*(getRPrice(fd, j) - getRCost(fd, j))

def v(fd):
	if fd == []:
		return 0
	else:
		return sum([getRProfit(fd, j) for j in range(3)])

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
            	
#O(2^N) complexity
def getShapleyPayoff(Eta, C):
	s = 0
	if C in Eta:
		Eta.remove(C)
		newEta = list(Eta)
		Eta.append(C)
	else:
		newEta = list(Eta)
	l1=len(Eta)
	for fd in powerset(newEta):
		l2 = len(fd)
		addend = math.factorial(l2)*math.factorial(l1-l2-1)*(v(fd+[C])-v(fd))/math.factorial(l1)
		s = s + addend
	return s
	
def getEgalitarianPayoff(fd, C):
	index = 0
	while fd[index].getName() != C.getName():
		index = index + 1
	return v([C]) + (v(fd) - sum([v([spi]) for spi in fd]))/len(fd)
	
def sat(v, vmax, vmin):
	if vmin<=v and v <= vmax:
		return 0.01 + 0.99*abs((v-vmin)/(vmax-vmin))
	elif v > vmax:
		return 1
	else:
		return 0

def getSPsat(spi):
	return 0.5*sat(v([spi]), 2*4*sum(spi.getCost()), 0.01*4*sum(spi.getCost())) + 0.5*sat(spi.getQra(), 1, 0.2)
	
def getFSPsat(Eta, fd, spi):
	if len(fd) == 0:
		return 0
	return 0.5*sat(getShapleyPayoff(Eta, spi), 2*v(fd)/len(fd), 0.01*v(fd)/len(fd)) + 0.5*sat(getQra(fd), 1, 0.2)
	#return 0.5*sat(getEgalitarianPayoff(fd, spi), 2*v(fd)/len(fd), 0.01*v(fd)/len(fd)) + 0.5*sat(getQra(fd), 1, 0.2)
	
def getFsat(fd):
	if len(fd) == 0:
		return 0
	return 0.5*sat(v(fd), 4*sum([2*sum(spi.getCost()) for spi in fd])/len(fd), 0.01*4*sum([sum(spi.getCost()) for spi in fd])/len(fd)) + 0.5*sat(getQra(fd), 1, 0.2)
	


def CGCFF(Eta):
	partition = []
	for spi in Eta:
		join(spi, partition)
	PartitionUnchanged = False # needed only to enter the loop (REPEAT-UNTIL would have handled this messy situation easily)
	while not PartitionUnchanged:
		PartitionUnchanged = True 
		for fd in partition:
						
			for sp in fd:
				singleton = False
				fd.remove(sp)
				if len(fd) == 0:
					singleton = True
					partition.remove(fd)
				
				join(sp, partition)
				
				ProviderStabilized = False
				
				if singleton:
					if [sp] in partition:
						ProviderStabilized = True
					
				else:
					if sp in fd:
						ProviderStabilized = True
					
				PartitionUnchanged = PartitionUnchanged and ProviderStabilized

	return partition
		
		
def join(sp, partition):
	maxSat = getSPsat(sp)
	bestFD = []
	for fd in partition:
		dummyFD = fd+[sp]
		Eta = [x for fed in partition for x in fed]+[sp]
		if getFSPsat(Eta, dummyFD, sp) >= maxSat and sp.getTrust()>=ThresholdTrust and sp.getQra() <= getQra(dummyFD) and getFsat(dummyFD) >= getFsat(fd) and getQra(dummyFD) >= getQra(fd):
			maxSat = getFSPsat(Eta, dummyFD, sp)
			bestFD = fd
	if bestFD!=[]:
		bestFD.append(sp)
		
	else:
		newFD = [sp]
		partition.append(newFD)
		
def beliefT(Si, Sj):
	#return (Si.getBrandValue()+Sj.getBrandValue())/2 - 0.3
	return random.random()
	
def beliefM(Si, Sj):
	return 1 - beliefT(Si, Sj)
	
		
import time		

time1 = time.clock()
partition = CGCFF(SP)
time2 = time.clock()
runtime = time2 - time1
print "CSP Name\tShapley Payoff\tFD.SP.sat\tSP.sat\tFD.sat\tFD.Profit\tSP.Profit\tFD.A\tFD.Q\tSP.a\tFD Size\tExecution Time(s)"
for fd in partition:
	print "\n"
	N = len(fd)
	if N>1:
		correctionFactor = sum([beliefT(Si, Sj) for Si in fd for Sj in fd if Si.getName()!= Sj.getName()])/(N*(N-1))
	else:
		correctionFactor = 1
	for sp in fd:
		print "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f" % (sp.getName(), getShapleyPayoff(SP, sp)*correctionFactor, getFSPsat(SP, fd, sp)*correctionFactor, getSPsat(sp), getFsat(fd)*correctionFactor, v(fd)*correctionFactor, v([sp]), getQra(fd)*correctionFactor, getQ(fd)*correctionFactor, sp.getQra(), len(fd), runtime) 
	

		
	


		
		 
