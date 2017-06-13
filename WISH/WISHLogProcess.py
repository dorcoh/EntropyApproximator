#--------------------------------------------------------------------------------------
# Copyright, 2013:
#
# Stefano Ermon 		- Cornell University 			, ermonste@cs.cornell.edu
# Ashish Sabharwal 		- IBM Watson Research Center 	, ashish.sabharwal@us.ibm.com
#---------------------------------------------------------------------------------------

import sys
import math
import random
import os
import re
import numpy
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import pylab
import os



def median(alist):
	
	if len(alist)==1:
		return alist[0]
	else:
		srtd = sorted(alist) # returns a sorted copy
		mid = len(alist)/2   # remember that integer division truncates

		if len(alist) % 2 == 0:  # take the avg of middle two
			return (srtd[mid-1] + srtd[mid]) / 2.0
		else:
			return srtd[mid]
			
def logsumexp(M,maxdepth):
	curmax = M[0]*math.log(10)
	for i in range(1,maxdepth+1):
		if (M[i]*math.log(10)+(i-1)*math.log(2)) > curmax:
			curmax = M[i]*math.log(10)+(i-1)*math.log(2)
			
	W = 0.0
	for i in range(0,maxdepth+1):
		W = W +math.exp(M[i]*math.log(10)+(i-1)*math.log(2)-curmax)
	return math.log(W)+curmax
	
def process_logs(folder):

	## w_i per level. format is (i,t,w_i^t, true|false)
	## true when certificate is found (solved to optimality)
	w=[]
	Samples=[] 
	maxdepth =-1

	os.chdir(folder)
	for files in os.listdir("."):
		if files.endswith(".LOG"):
			print "processing " + files
			m = re.search('.xor(\d+).loglen(\d+).(\d+).',files)
			if m is not None:
				#print m.group(1)+","+m.group(2)+","+m.group(3)
				# read in the optimal val
				i = int(m.group(1))							# level
				if i>maxdepth:
					maxdepth = i
				loglen = int(m.group(2))
				t = int(m.group(3))							# sample number
				with open(files, 'r') as f:
					 #check if solved to optimality
					 lines = f.read()
					 entries = re.search("Optimum: (\d+) log10like: ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)", lines)
					 if entries is not None:
						w.append([i,t,float(entries.group(2)),True])
						entriesSol = re.search("Optimal solution: (\d+)", lines)
						if entriesSol is not None:
							Samples.append([i,t,float(entries.group(2)),True,entriesSol.group(1)])
					 else:
						entries2 = re.findall("New solution: (\d+) log10like: ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)", lines)
						if entries2:
							# find last entry
							#print entries2[-1]
							w.append([i,t,float(entries2[-1][1]),False])
							entriesSol = re.findall("Current solution: (\d+)", lines)
							if entriesSol:
								#print entriesSol[-1]
								Samples.append([i,t,float(entries2[-1][1]),False,entriesSol[-1]])
						else:
							# Inconsistency detected!
							entries3 = re.search("No solution in (\d+) backtracks and", lines)
							entries4 = re.search("Inconsistency detected!", lines)
							if (entries3 is not None or entries4 is not None):
								#nothing = 34
								w.append([i,t,float("-inf"),True])			# problem is unsat
							else:
								print "---> Could not read solution value"
			else:
				print "Error " + files

	x=[]
	y=[]

	xOpt=[]
	xSubOpt=[]
	yOpt=[]
	ySubOpt=[]
	# compute medians
	M = [0] * (maxdepth+1)
	for i in range(maxdepth+1):
		results = [t[2] for t in w if t[0] == i]
		resultsOpt = [t[2] for t in w if ((t[0] == i) & (t[3]==True))]
		resultsSubOpt = [t[2] for t in w if ((t[0] == i) & (t[3]==False))]
		if len(results)>0:
			x = x+ results
			xOpt = xOpt + resultsOpt 
			xSubOpt = xSubOpt + resultsSubOpt 
			y = y + [-i]*len(results)
			yOpt = yOpt + [-i]*len(resultsOpt)
			ySubOpt = ySubOpt + [-i]*len(resultsSubOpt)
			M[i] = median(results)
		else:
			M[i] =float("-inf")
			print "Level " +str(i) + " is empty"
		print str(i) + "," + str(M[i])+";"+str(len(results)-results.count(float("-inf")))


	print "Final log-estimate: " + str(logsumexp(M,maxdepth))


	pylab.plot(xOpt,yOpt,'yo',xSubOpt,ySubOpt,'go')
	pylab.xlabel('Log Score')
	pylab.ylabel('-Depth')
	pylab.title('Distribution of w_i^t (CDF)')
	pylab.grid(True)
	pylab.savefig('cdf.plot.pdf')
	pylab.show()

def process_logs_cplex_LB(folder):
	w=[]
	Samples=[] 
	maxdepth =-1
	import os
	os.chdir(folder)
	for files in os.listdir("."):
		if files.endswith(".LOG"):
			print "processing " + files
			m = re.search('.xor(\d+).loglen(\d+).(\d+).',files)
			if m is not None:
				#print m.group(1)+","+m.group(2)+","+m.group(3)
				# read in the optimal val
				i = int(m.group(1))							# level
				if i>maxdepth:
					maxdepth = i
				loglen = int(m.group(2))
				t = int(m.group(3))							# sample number
				with open(files, 'r') as f:
					 #check if solved to optimality
					 lines = f.read()
					 #entries = re.search("Optimum: (\d+) log10like: ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) prob: ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) in (\d+) backtracks and (\d+) nodes and ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?) seconds", lines)
					 entriesoo = re.search("Solution status = Optimal", lines)
					 entries = re.search("Solution value log10lik = ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)", lines)
					 if entries is not None and entriesoo is not None:
						print entries.group(1)
						print entries.group(2)
						w.append([i,t,float(entries.group(1)),True])
						entriesSol = re.search("Optimal solution: (\d+)", lines)
						if entriesSol is not None:
							Samples.append([i,t,float(entries.group(1)),True,entriesSol.group(1)])
					 else:
						entries2 = re.findall("\n(\s+)(\d+)(\s+)(\d+)(\s+)([-+]?\d+\.\d+)(\s+)(\d+)(\s+)(\s+|[-+]?\d+\.\d+)(\s+)([-+]?\d+\.\d+)", lines)
						if entries2:
							# find last entry
							#print entries2
							if is_number(entries2[-1][-3]):
								#print "found something " + str(float(entries2[-1][-3]))
								w.append([i,t,float(entries2[-1][-3]),False])
							else:
								w.append([i,t,float("-inf"),True])			# problem is unsat
							entriesSol = re.findall("Current solution: (\d+)", lines)
							if entriesSol:
								#print entriesSol[-1]
								Samples.append([i,t,float(entries2[-1][1]),False,entriesSol[-1]])
						else:
							# Inconsistency detected!
							entries3 = re.search("No solution in (\d+) backtracks and", lines)
							entries4 = re.search("Inconsistency detected!", lines)
							entries5 =	re.search("Infeasibility row", lines)
							entries6 =	re.search("Failed to optimize LP", lines)							
							entries7 =	re.search("Infeasible column", lines)
							
							if (entries3 is not None or entries4 is not None or entries5 is not None or entries6 is not None or entries7 is not None):
								#nothing = 34
								w.append([i,t,float("-inf"),True])			# problem is unsat
							else:
								print "ERROR reading logs", files
					#New solution: 0 log10like: 0 prob: 1 (0 backtracks, 24 nodes, depth 25)
					#Optimum: 91499957 log10like: -3.97379 prob: 0.00010622 in 299 backtracks and 600 nodes and 0.08 seconds.
			else:
				print "Error " + files



	x=[]
	y=[]

	xOpt=[]
	xSubOpt=[]
	yOpt=[]
	ySubOpt=[]
	# compute medians
	M = [0] * (maxdepth+1)
	for i in range(maxdepth+1):
		results = [t[2] for t in w if t[0] == i]
		resultsOpt = [t[2] for t in w if ((t[0] == i) & (t[3]==True))]
		resultsSubOpt = [t[2] for t in w if ((t[0] == i) & (t[3]==False))]
		if len(results)>0:
			x = x+ results
			xOpt = xOpt + resultsOpt 
			xSubOpt = xSubOpt + resultsSubOpt 
			y = y + [-i]*len(results)
			yOpt = yOpt + [-i]*len(resultsOpt)
			ySubOpt = ySubOpt + [-i]*len(resultsSubOpt)
			M[i] = median(results)
			print str(i) + "," + str(M[i])
		else:
			M[i] =float("-inf")
			#print "Level " +str(i) + " is empty"
		

	# 
	# W = pow(10,M[0])
	# for i in range(1,maxdepth+1):
		# print str(i) + "," + str(M[i])
		# W = W +pow(10,M[i])*pow(2,i-1)

	print "Final log-estimate: " + str(logsumexp(M,maxdepth))

	#W = math.exp(logsumexp(M))

def process_logs_cplex_UB(folder):	
	w=[]
	Samples=[] 
	maxdepth =-1
	#import os
	#os.chdir(folder)
	for files in os.listdir("."):
		if files.endswith(".LOG"):
			print "processing " + files
			m = re.search('.xor(\d+).loglen(\d+).(\d+).',files)
			if m is not None:
				#print m.group(1)+","+m.group(2)+","+m.group(3)
				# read in the optimal val
				i = int(m.group(1))							# level
				if i>maxdepth:
					maxdepth = i
				loglen = int(m.group(2))
				t = int(m.group(3))							# sample number
				with open(files, 'r') as f:
					 #check if solved to optimality
					 lines = f.read()
					 entries = re.search("Solution value log10lik = ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)", lines)
					 entriesoo = re.search("Solution status = Optimal", lines)
					 if entries is not None and entriesoo is not None:
						w.append([i,t,float(entries.group(1)),True])
						entriesSol = re.search("Optimal solution: (\d+)", lines)
						if entriesSol is not None:
							Samples.append([i,t,float(entries.group(1)),True,entriesSol.group(1)])
					 else:

						entries2 = re.findall("\n(\s+)(\d+)(\s+)(\d+)(\s+)([-+]?\d+\.\d+)(\s+)(\d+)(\s+)(\s+|[-+]?\d+\.\d+)(\s+)([-+]?\d+\.\d+)", lines)
						if entries2:
							w.append([i,t,float(entries2[-1][-1]),False])
							entriesSol = re.findall("Current solution: (\d+)", lines)
							if entriesSol:
								Samples.append([i,t,float(entries2[-1][1]),False,entriesSol[-1]])
						else:
							entries3 = re.search("No solution in (\d+) backtracks and", lines)
							entries4 = re.search("Inconsistency detected!", lines)
							entries5 =	re.search("Infeasibility row", lines)
							entries6 =	re.search("Failed to optimize LP", lines)							
							entries7 =	re.search("Infeasible column", lines)							
							
							if (entries3 is not None or entries4 is not None or entries5 is not None or entries6 is not None or entries7 is not None):
								w.append([i,t,float("-inf"),True])			# problem is unsat
							else:
								print "Error processing log file"

			else:
				print "Error " + files



	x=[]
	y=[]

	xOpt=[]
	xSubOpt=[]
	yOpt=[]
	ySubOpt=[]
	# compute medians
	M = [0] * (maxdepth+1)
	for i in range(maxdepth+1):
		results = [t[2] for t in w if t[0] == i]
		resultsOpt = [t[2] for t in w if ((t[0] == i) & (t[3]==True))]
		resultsSubOpt = [t[2] for t in w if ((t[0] == i) & (t[3]==False))]
		if len(results)>0:
			x = x+ results
			xOpt = xOpt + resultsOpt 
			xSubOpt = xSubOpt + resultsSubOpt 
			y = y + [-i]*len(results)
			yOpt = yOpt + [-i]*len(resultsOpt)
			ySubOpt = ySubOpt + [-i]*len(resultsSubOpt)
			M[i] = median(results)
			print str(i) + "," + str(M[i])
		else:
			M[i] =float("-inf")
			#print "Level " +str(i) + " is empty"
	
	
	print "Upper bound log-estimate: " + str(logsumexp(M,maxdepth))
	
	
if __name__ == "__main__":
	print "Processing LOGs.."
	process_logs_cplex_LB(sys.argv[1])	