'''Sequential Piecewise Linear Discriminant Plane Seperation
December 2014 
Author: Joshua G. Albert (albert@strw.leidenuniv.nl)

Basic usage: see example breast_cancer_test()
	Given an N-dimensional dataspace suppose you have two distinct sets, A and B, of datapoints. Suppose you have K datapoints in A, and M datapoints in B. Organize them as 1xN rows in two arrays, i.e.
	
	A = [[a11,a12,...,a1N],[a21,a22,...,a2N],...,[aK1,aK2,...,aKN]]
	B = [[b11,b12,...,b1N],[b21,b22,...,b2N],...,[bK1,bK2,...,bKN]]

	Then the following steps are typical:
	1) Make a training set, and test set from A, and B, with say 75% of the data in the training set.
	
	A_train,B_train,A_test,B_test = make_train_set(A,B,0.75)

	2) Calculate, and save a solution.
	
	sols = msm_solve(A_train,B_train)
	print "Solution:",sols
	save_sol("solution.npy",sols)

	3) Test how well your solution works, by measuring false positive rate, and plot this.

	plot_msm_test_data(sols,A_test,B_test)
	
Extended Usage:
	Suppose you have L datasets instead of just two. You can extend this method in the fairly straightforward way by treating the problem as a "log(L)/log(2)-bit" problem. In the basic usage you have a "1bit" problem.

Note: this requires the accompanying C program "msm" to compile and work on your system. See "Makefile".'''

import time
import numpy as np
import pylab as plt
from os import system

def msm_solve (A_,B_):
	'''Takes data sets A_ and B_ and creates sequences of piecewise linear discriminant planes which seperates the data.'''
	def sep(x,w,alpha,beta):
		if np.dot(x,w) > beta:
			return 1
		elif np.dot(x,w) < alpha:
			return -1
		else:
			return 0
	A = []
	B = []
	print "A: %d x %d"%(len(A_),len(A_[0]))
	print "B: %d x %d"%(len(B_),len(B_[0]))
	print "Reducing to make dis-joint"
	def vecinarray(v,A):
		i = 0
		B = A
		while i < len(v):
			B = B[(v[i] == B[:,i]),:]
			if np.size(B) == 0:
				return False
			i += 1
		return True

	count = 0
	for a in A_:
		if vecinarray(a,B_):
			#	print "removing",a
			count += 1
		else:
			A.append(a)
	for b in B_:
		if vecinarray(b,A_):
			pass
		else:
			B.append(b)
	if count > 0:
		print "Dis-joint: Removed %d duplicates"%(count)
		print "*A: %d x %d"%(len(A),len(A[0]))
		print "*B: %d x %d"%(len(B),len(B[0]))
	else:
		print "Already dis-joint"
	
	def make_mod(A,B,n=0):
		M1 = len(A)
		M2 = len(B)
		N = len(A[0])
		nrows = M1+M2
		ncols = N+2
		n1 = 1
		n2 = 0
		mn = []
		while n1 <= N:
			f = open("msm_mod","w+")
			f.write("N %d %d\n"%(nrows,ncols)) #set minimization and nrows, ncols
			col = 1
			while col <= N:
				if col == n1:
					if n == 0:
						f.write("d %d %f %f\n"%(col,-1.,1.)) #set bounds col -1 to 1
					else:
						f.write("f %d %f\n"%(col,(-1.)**n2)) # set col n1 = +-1.
					n1 += n2
					n2 += 1
					n2 = (n2%2)
				else:
					f.write("d %d %f %f\n"%(col,-1.,1.)) #set bounds col -1 to 1
				col += 1
			f.write("r %d\n"%(N+1)) #set alpha and beta free range
			f.write("r %d\n"%(N+2))
			f.write("C %d %f\n"%(N+1,-1.))#set C[alpha]=-1
			f.write("C %d %f\n"%(N+2,1.))#set C[beta]=1
			row = 1
			while row <= M1:
				col = 1
				while col <= N:
					f.write("A %d %d %f\n"%(row,col,A[row-1][col-1]))# Set A
					col += 1
				f.write("A %d %d %f\n"%(row,N+1,-1.))
				f.write("L %d %f\n"%(row,0.))# A w - alpha >= eps
				row += 1
			while row <= M1 + M2:
				col = 1
				while col <= N:
					f.write("A %d %d %f\n"%(row,col,-B[row-1-M1][col-1]))# set B
					col += 1
				f.write("A %d %d %f\n"%(row,N+2,1.))
				f.write("L %d %f\n"%(row,0.))# -B w + beta >= eps
				row += 1
			f.close()
			# Run msm.c -> create msm_mod.mps and msm_mod_results
		#	print "Running linear program msm"
			system("./msm msm_mod")
			D = np.genfromtxt("msm_mod_results" , dtype=('double'))
			w = np.array(D[0:N,1]).reshape([N,1])
			alpha = D[N,1]
			beta = D[N+1,1]
			obj = -alpha + beta
		#	if obj < 0.:
		#		print "found convex solution"
			s = {}
			s['obj'] = obj
			s['alpha'] = alpha
			s['beta'] = beta
			s['w'] = w
			miss = 0
			for a in A:
				if np.dot(a,w) <= beta:
					miss += 1
			for b in B:
				if np.dot(b,w) >= alpha:
					miss += 1
			s['miss'] = miss
			mn.append(miss)
			if miss == min(mn):
				minsol = s
			if n == 0:
				break
		return minsol
	sols = {}
	sol = make_mod(A,B,n=0)
	if sol['obj'] < 0. and sol['miss'] == 0:
		print "Convex Solution"
		sols[0] = sol
		sols['jmax'] = 0
		return sols
	else:
		nonconvex = True
	Total = len(A) + len(B)
	j = 0
	while nonconvex:
		print "Working on plane: %d"%(j)
		print "*A: %d x %d"%(len(A),len(A[0]))
		print "*B: %d x %d"%(len(B),len(B[0]))
		M1 = len(A)
		M2 = len(B)
		N = len(A[0])
		sol = make_mod(A,B,n=1)
		obj = sol['obj']
		miss = sol['miss']
		if miss == M1+M2:
			if M1 > M2: #remove larger first
				print "-->Degeneracy: Resolving via A"
				for a in A:
					sol = make_mod([a],B,n=0)
					if sol['obj'] < 0. and sol['miss'] == 0:
						sol['alpha'] = float('-inf')
						break

			else:
				print "-->Degeneracy: Resolving via B"
				for b in B:
					sol = make_mod(A,[b],n=0)
					if sol['obj'] < 0. and sol['miss'] == 0:
						sol['beta'] = float("inf")
						break
		w = sol['w']
		alpha = sol['alpha']
		beta = sol['beta']
		miss = 0.
		A__ = []
		B__ = []
		for a in A:
			if np.dot(a,w) <= beta:
				A__.append(a)
				miss += 1
		for b in B:
			if np.dot(b,w) >= alpha:
				B__.append(b)
				miss += 1
		sol['miss'] = miss
		print "Seperation: %f | Unclassified: %d | Classified by this plane: %d"%(1. - float(miss)/(Total),miss, M1+M2-miss)
	#	if obj < 0.:
	#		print "Found convex solution"
	#		nonconvex = False
		if len(A__) == 0:
			if len(B__) == 0:
				nonconvex = False
			else:
				print "A empty, B: %d x %d -> convex solution exists"%(len(B__),len(B__[0]))
				sol=make_mod(A,B,n=0)
#				print sol['miss']
				if sol['obj'] < 0. and sol['miss'] == 0:
					nonconvex=False
				else:
					print "Failed to find convex"
					exit()
		if len(B__) == 0:
			if len(A__) == 0:
				nonconvex=False
			else:
				print "B empty, A: %d x %d -> convex solution exists"%(len(A__),len(A__[0]))
				sol=make_mod(A,B,n=0)
#				print sol['miss']
				if sol['obj'] < 0. and sol['miss'] == 0:
					nonconvex=False
				else:
					print "Failed to find convex"
					exit()
		A = A__
		B = B__
		sols[j] = sol
		j += 1
	sols['jmax'] = j - 1
	print "Finished at plane: %d"%(j-1)
	return sols

def sep(sols,v):
	'''Use solution sols to seperate v. If decides v in A, returns positive plane of seperation (1 indexed). If decides v in B, returns negative plane of seperation (1 indexed).'''
	jmax = sols['jmax']#zero indexed
	j = 0
	while j < jmax:
		w = sols[j]['w']
		alpha = sols[j]['alpha']
		beta = sols[j]['beta']
		if np.dot(v,w) > beta:
		#	print "seperated in plane",j #(zero indexed)
			return j + 1
		elif np.dot(v,w) < alpha:
		#	print "seperated in plane",j #(zero indexed)
			return -1 - j
		j += 1
	w = sols[jmax]['w']
	alpha = sols[jmax]['alpha']
	beta = sols[jmax]['beta']
	if np.dot(v,w) > (alpha+beta)/2.:
#		print "seperated in plane",jmax
		return jmax + 1
	else:
#		print "seperated in plane",jmax
		return -1 - jmax

def plot_msm_test_data(sols,A,B,plot=True):
	'''Use solution sols to verify on test sets A,B. Will verify the number of false positives, and plot data.'''
	global plot_num_i
	global plot_num_j
	global Nw
	global ax2
	Nw = len(sols[0]['w'])
	plot_data=True
	if plot_data:
		if plot:
			plt.figure()
			ax3 = plt.axes((0.01,0.01,0.05,0.04))
			ax1 = plt.axes((0.06,0.01,0.05,0.04))
			ax2 = plt.axes((0.09,0.08,0.9,0.9))
			next_btn = plt.Button(ax1,"Next")
			prev_btn = plt.Button(ax3,"Prev")
		plot_num_i = 0
		plot_num_j = 1
		Am = []
		Af = []
		Bm = []
		Bf = []
		fail = 0
		for a in A:
			p = sep(sols,a)
			if p < 0:
				fail += 1
				Am.append(a)
			else:
				Af.append(a)
		print "A: false positives: %d from: %d total test data"%(fail,len(A))
		fail = 0
		for b in B:
			p = sep(sols,b)
			if p > 0:
				fail += 1
				Bm.append(b)
			else:
				Bf.append(b)		
		print "B: false positives: %d from: %d total test data"%(fail,len(B))

		def plot_next(event):
			global plot_num_i
			global plot_num_j
			global ax2
			global Nw
		#	print event
			if plot_num_i >= (Nw-1.):
				print "No more iterations left"
				return
			print "plotting %d vs %d"%(plot_num_i,plot_num_j)
			plt.axes(ax2)
			plt.cla()
			if len(Am) > 0:
				plt.scatter(np.array(Am)[:,plot_num_i],np.array(Am)[:,plot_num_j],c='yellow',label='missed a')
			if len(Af) > 0:
				plt.scatter(np.array(Af)[:,plot_num_i],np.array(Af)[:,plot_num_j],c='blue',label='found a')
			if len(Bm) > 0:
				plt.scatter(np.array(Bm)[:,plot_num_i],np.array(Bm)[:,plot_num_j],c='red',label='missed b')
			if len(Bf) > 0:
				plt.scatter(np.array(Bf)[:,plot_num_i],np.array(Bf)[:,plot_num_j],c='pink',label='found b')			
			plt.legend()
			plt.xlabel("Data %d"%plot_num_i,fontsize=12)
			plt.ylabel("Data %d"%plot_num_j,fontsize=12)
			plot_num_j += 1
			if plot_num_j == Nw:
				plot_num_i += 1
				plot_num_j = plot_num_i + 1
			plt.draw()
		def plot_prev(event):
			global plot_num_i
			global plot_num_j
			global Nw
			if plot_num_i == 0 and plot_num_j == 1:
				return
			plot_num_j += -2
			if plot_num_j <= plot_num_i:
				plot_num_i += -1
				plot_num_j = Nw - 1
			plot_next(event)
		if plot:
			next_btn.on_clicked(plot_next)
			prev_btn.on_clicked(plot_prev)
			plot_next("First plot up")
			plt.show()

def make_train_set(A,B,p):
	'''Takes full A,B sets and returns p*100% of the set in A_train,B_train and the rest in A_test,B_test''' 
	A_test = []
	B_test = []
	A_train = []
	B_train = []
	for a in A:
		if np.random.uniform() >= p:
			A_test.append(a)
		else:
			A_train.append(a)
	for b in B:
		if np.random.uniform() >= p:
			B_test.append(b)
		else:
			B_train.append(b)
	return np.array(A_train),np.array(B_train),np.array(A_test),np.array(B_test)

def save_sol(fname,sol):
	'''Saves your solution to fname. Appends .npy if not already there.'''
	print "Saving solution to:",fname
	np.save(fname,sol)

def load_sol(fname):
	'''Loads the solution in fname'''
	k = np.load(fname)
	return k.item(0)

def breast_cancer_test():
	D = np.genfromtxt("wdbc.data",delimiter=",")[:,2:]
	#Uncomment to see what happens when there is no structure in the data. Hint try running twice and comparing solutions.
#	D = np.random.uniform(low=0.,high=1.,size=np.shape(D))
	M = np.genfromtxt("wdbc.data",delimiter=",",usecols=[1],dtype=('S'))
	m = 0
	A = []
	B = []
	while m < len(M):
		if M[m] == 'M':
			A.append(D[m,:])	
		elif M[m] == 'B':
			B.append(D[m,:])
		m += 1
	A = np.array(A)
	B = np.array(B)
	A_train,B_train,A_test,B_test = make_train_set(A,B,0.5)
	sols = msm_solve(A_train,B_train)
	print "Breast Cancer Solution:",sols
	save_sol("breastcancersolution.npy",sols)
	plot_msm_test_data(sols,A_test,B_test)
	return sols,A,B

class KbitPattern():
	def __init__(self,fname):
		self.fname = fname
		self.load()
	def save(self):
		fname = "Kbit-%s.npy"%(self.fname)
		print "Saving problem:",fname
		np.save(fname,{"splits":self.splits,"sols":self.sols})
	def load(self):
		fname = "Kbit-%s.npy"%(self.fname)
		d = np.DataSource()
		if d.exists(fname):
			print "Loading:",fname
			k = np.load(fname)
			K = k.item(0)
			self.splits=K['splits']
			self.sols=K['sols']
		else:
			print "No file:",fname
	def add_data(self,D,decisioncol,trainpercent=0.7):
		A = np.arange(np.size(D[0]))
		self.data = D[:,A!=decisioncol]
		self.decisions = D[:,decisioncol]
		self.train_data = []
		self.test_data = []
		self.train_decisions = []
		self.test_decisions = [] 
		i = 0
		while i < len(self.decisions):
			if np.random.uniform() < trainpercent:
				self.train_data.append(self.data[i])
				self.train_decisions.append(self.decisions[i])
			else:
				self.test_data.append(self.data[i])
				self.test_decisions.append(self.decisions[i])
			i += 1
		self.train_decisions = np.array(self.train_decisions)
		self.test_decisions = np.array(self.test_decisions)
	def recursive_split(self,p,minp):
		split = np.median(p)
		plow = p[p<split]
		phigh = p[p>=split]
		if np.size(plow) >= minp:
			#sucess split is a divider
			self.add_split(split)
			self.recursive_split(plow,minp)
			self.recursive_split(phigh,minp)
		else:
			return
	def partition_data(self,minpercent=0.1):
		'''Calculates an even set of partitions which have at least minpercent sample size in deciesion variable x.'''
		self.splits=[]
		self.recursive_split(self.train_decisions,np.size(self.train_decisions)*minpercent)
		self.datasets = {'train':{},'test':{}}
		b = 0
		while b < len(self.splits):
			i = 0
			while i < len(self.train_decisions):
				pat = self.get_pattern(self.train_decisions[i])
				if self.splits[b] not in self.datasets['train'].keys():
					self.datasets['train'][self.splits[b]] = {'A':[],'B':[]}
				if pat[b]:
					self.datasets['train'][self.splits[b]]['A'].append(self.train_data[i])

				else:
					self.datasets['train'][self.splits[b]]['B'].append(self.train_data[i])
				i += 1
			i = 0
			while i < len(self.test_decisions):
				pat = self.get_pattern(self.test_decisions[i])

				if self.splits[b] not in self.datasets['test'].keys():
					self.datasets['test'][self.splits[b]] = {'A':[],'B':[]}
				if pat[b]:
					self.datasets['test'][self.splits[b]]['A'].append(self.test_data[i])
					
				else:
					self.datasets['test'][self.splits[b]]['B'].append(self.test_data[i])	
				i += 1
			b += 1
	def sep(self,x):
		b = 0
		pat = []
		while b < len(self.splits):
			p = sep(self.sols[self.splits[b]],x)
			if p > 0:
				pat.append(1)
			else:
				pat.append(0)
			b += 1
		lb,ub=self.get_interval(pat)
		return lb,ub
	def test_data_sep(self):
		i = 0
		plt.figure()
		lb=[]
		ub = []
		while i < len(self.test_decisions):
			lb_,ub_ = self.sep(self.test_data[i])
			lb.append(lb_)
			ub.append(ub_)
			i += 1
		plt.plot(np.arange(len(self.test_decisions)),lb,'g',label='Model lower bound')
		plt.plot(np.arange(len(self.test_decisions)),ub,'b',label="Model upper bound")
		plt.scatter(np.arange(len(self.test_decisions)),self.test_decisions,c='r',label="Measurment")
		plt.legend()
		plt.show()
		mean = []
		i = 0
		while i < len(lb):
			if np.isinf(lb[i]):
				mean.append(ub[i])
			elif np.isinf(ub[i]):
				mean.append(lb[i])
			else:
				mean.append((np.array(lb[i])+np.array(ub[i]))/2.)
			i += 1
		dev = (self.test_decisions - np.array(mean))
		rms = np.sqrt(np.sum(dev**2)/len(lb))
		#plt.figure()
		#plt.scatter(self.test_decisions,np.abs(dev))
		#plt.show()
		print rms
	def solve(self):
		b = 0 
		self.sols = {}
		while b < len(self.splits):
			A_train = np.array(self.datasets['train'][self.splits[b]]['A'])
			B_train = np.array(self.datasets['train'][self.splits[b]]['B'])
		#	A_test = np.array(self.datasets['test'][self.splits[b]]['A'])
		#	B_test = np.array(self.datasets['test'][self.splits[b]]['B'])
			print "Solving split:",self.splits[b]
			sols = msm_solve(A_train,B_train)
#			save_sol("Kbit-%s-%.1f.npy"%(self.fname,self.splits[b]),sols)
			plot_msm_test_data(sols,A_train,B_train,plot=False)
			self.sols[self.splits[b]] = sols
			b += 1
		self.save()
	def add_split(self,split):
		self.splits.append(split)
	def get_pattern(self,x):
		pattern = []
		lb = -np.inf
		ub = np.inf
		for b in self.splits:
			if x >= b:
				lb = max(b,lb)
				pattern.append(1)
			else:
				ub = min(b,ub)
				pattern.append(0)
		return pattern
	def num2pat(self,num):
		b = 0
		bits = []
		k = num
		while b < len(self.splits):
			bits.append(k % 2)
			k = (k - bits[-1])/2
			b += 1
		return bits
	def pat2num(self,pattern):
		b = 0
		num=0
		while b < len(self.splits):
			num += pattern[b]*2**b
			b += 1
		return num
	def get_interval(self,pattern):
		#fix this
		lb = -np.inf
		ub = np.inf
		#lb = self.splits[0]
		#ub = self.splits[0]
		b = 0
		while b < len(self.splits):
			if pattern[b]:
				lb_ = max(self.splits[b],lb)
				if lb_ <= ub:
					lb = lb_
			else:
				ub_ = min(self.splits[b],ub)
				if ub_ >= lb:
					ub = ub_
			b += 1
		return lb,ub

def bodyfat_test1():
	D = np.genfromtxt("bodyfat.csv",delimiter=",",usecols=(1,2,3,4,5,6),comments="#")
	K = KbitPattern("bodyfat-women")
	K.add_data(D,1,1.)
	K.partition_data(0.1)
	K.solve()
	K.add_data(D,1,0.)
	K.test_data_sep()

def bodyfat_test2():  
	'''0 Case Number
1 Percent body fat using Brozek's equation, 457/Density - 414.2
2 Percent body fat using Siri's equation, 495/Density - 450
3 Density (gm/cm^3)
4 Age (yrs)
5 Weight (lbs)
6 Height (inches)
7 Adiposity index = Weight/Height^2 (kg/m^2)
8 Fat Free Weight = (1 - fraction of body fat) * Weight, using Brozek's formula (lbs)
9 Neck circumference (cm)
10 Chest circumference (cm)
11 Abdomen circumference (cm) "at the umbilicus and level with the iliac crest"
12 Hip circumference (cm)
13 Thigh circumference (cm)
14 Knee circumference (cm)
15 Ankle circumference (cm)
16 Extended biceps circumference (cm)
17 Forearm circumference (cm)
18 Wrist circumference (cm) "distal to the styloid processes"'''
	#age(years), weight(lbs), height(in), neckcirc (cm), chest circ (cm), abdomen at belly button circ (cm), hip circ (cm), thigh circ (cm), knee circ (cm), ankle circ (cm), extended biceps circ (cm), forarm circ (cm), wrist circ (cm)
	D = np.genfromtxt("fat.dat.txt",usecols=(1,4,5,6,9,10,11,12,13,14,15,16,17,18),comments="#")
	K = KbitPattern("bodyfat-men")
	K.add_data(D,0,1.)
	K.partition_data(0.05)
	K.solve()
	K.add_data(D,0,0.)
	K.test_data_sep()
	print K.decisions, K.splits

if __name__ == '__main__':
	#bodyfat_test1()

	#bodyfat_test2()
	breast_cancer_test()
