from mpi4py import MPI
from subprocess import Popen, PIPE, STDOUT, call, check_output, CalledProcessError
import random
import time

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

sim_path = '/sciclone/home00/geogdan/geoMatch_testing/SimTests/geoMatch_simTests.R'

iterations = 1000

for c in range(0,size):
	if rank == c:
		for i in range((iterations*rank/size + 1), (iterations*(rank+1)/size+1)):
			print "Worker - rank %d on %s."%(rank, name) 
			print "i:%d"%i
			
			try:
				R_ret = check_output("Rscript " +
				sim_path,
				stderr=STDOUT,
				shell=True)

				print i
				print R_ret
				print "========================"
			
			except CalledProcessError as sts_err:
				print sts_err
