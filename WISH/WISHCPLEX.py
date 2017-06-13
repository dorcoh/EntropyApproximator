#----------------------------------------------------------------------------------------
# Copyright, 2013:
#
# Stefano Ermon 		- Cornell University 			, ermonste@cs.cornell.edu
# Ashish Sabharwal 		- IBM Watson Research Center 	, ashish.sabharwal@us.ibm.com
#----------------------------------------------------------------------------------------

import sys
import math
import random
import os
import argparse
from WISHLogProcess import process_logs 
from WISHLogProcess import process_logs_cplex_LB
from WISHLogProcess import process_logs_cplex_UB

# version number
__version__ = '1.0'

#########################################
# Usage Information:
# run "python WISH.py -h" for help
#########################################

parser = argparse.ArgumentParser(description='Estimate the partition function using the WISH algorithm and CPLEX for the optimization.')

parser.add_argument('-v', '--version',      action='version', version='%(prog)s ' + __version__)

parser.add_argument("infile", help="Graphical model (in UAI format)")
parser.add_argument("outfolder", help="Folder where logs are stored")

parser.add_argument('-alpha', '--alpha', type=float, help="Accuracy alpha", default=1.0)

parser.add_argument('-delta', '--delta', type=float, help="Failure probability delta", default=0.1)

parser.add_argument('-timeout', '--timeout', type=int, help="Timeout for each optimization instance (seconds)", default=10)

args = parser.parse_args()


print "Reading factor graph from " + args.infile	
inputfile = open(args.infile, "r")

fileName, fileExtension = os.path.splitext(args.infile)

ind = 0
origNbrFactor = 0
origNbrVar = 0
for l in inputfile:
	if not l.strip()=='':
		ind = ind +1
		if ind==2:
			origNbrVar=int(l)
		elif ind==3:
			l = l.rstrip("\n")
		elif ind==4:			## add xor cpt tabe
			origNbrFactor = int(l)
		elif ind>5:
			break
print "Model with " + str(origNbrVar) + "variables and "+str(origNbrFactor) +" factors"

depth = origNbrVar

T = int(math.ceil(math.log(origNbrVar)*math.log(1.0/args.delta)/args.alpha))

print "Using " + str(T) +" samples per level"

os.system("mkdir "+args.outfolder)
 
for i in range(0,depth+1):			## main for loop
	if i==0:
		sampnum=1
	else:
		sampnum=T
	for t in range(1,sampnum+1):			## main for loop
		outfilenamelog = "%s.xor%d.loglen%d.%d.ILOGLUE.uai.LOG" % (os.path.basename(fileName) , i , 0 , t)
		cmdline = ("timeout %d ./WH_cplex -paritylevel 1 -number %d -seed 10 %s > %s") % (args.timeout , i , args.infile , args.outfolder +"/"+ outfilenamelog)
		os.system(cmdline)
		## Parallel execution:
		##
		## assign this job to a separate core (a system dependent script is needed here)
		## we provide an example based on Torque/PBS:
		##
		## os.system("qsub -v basedir="+basedir+",file="+infile+",level="+str(i)+",len="+str(0)+",outdir="+outdir+",sample="+str(t)+",timeout=900s"+" LaunchIloglue.sh")

		
process_logs_cplex_LB(args.outfolder)
process_logs_cplex_UB(args.outfolder)