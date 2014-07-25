import os, string, threading, sys, shutil, math, pickle, random, re, copy
import numpy
import filetypes, lammps, utils

def bend():
	by_type = []
	for t in range(3,18):
		run_name = 'smd19_'+str(t)
		energies = []
		for line in open('lammps/'+run_name+'.log'):
			try:
				a = line.split()
				if len(a)==3:
					energies.append( tuple([float(s) for s in a]) )
			except ValueError: pass
		energies = [(abs(e[0]-energies[0][0]),e[1],e[2]) for e in energies]
		#energies = [e[1] for e in energies]
		
		spring_e = [e[1] for e in energies]
		spring_x = [e[0] for e in energies]
		
		'''
		spring_data = [e for e in energies if e[0]<1.0]
		
		#coeffs = numpy.polyfit(spring_x, spring_e, 2)
		#p = numpy.poly1d(coeffs)
		
		#spring_x2 = [x**2 for x in spring_x]
		#A = numpy.vstack([spring_x2, numpy.ones(len(spring_x2))]).T
		#m, c = numpy.linalg.lstsq(A, spring_e)[0]
		
		k = 1.0
		best_error = 1e10
		for step in range(0):
			k_new = random.gauss(k,k*0.2)
			error = sum([ (k_new*spring_x[i]**2 - spring_e[i])**2 for i in range(len(spring_x)) ])
			if error<best_error:
				best_error = error
				k = k_new
		
		#print best_error**0.5/len(spring_x), k
		
		p = numpy.poly1d([k,0.0,0.0])
		
		#print [p(x/10.0) for x in range(10)]
		
		section = [e for e in energies if abs(e[0]-1.0)<0.1]
		average = sum([e[1] for e in section])/len(section)
		print average
		'''
		#exponential_moving_average = [sum([e[1]*math.exp(-(e[0]-energies[i][0])**2) for e in energies]) for i in range(len(energies))]
		#energies = [ (energies[i][0], exponential_moving_average[i]) for i in range(len(energies))]
		
		N = 10
		medians = []
		for i in range(N/2, len(energies)-N/2):
			median_e = sorted(spring_e[i-N/2:i+N/2])[N/2-1]
			median_x = sorted(spring_x[i-N/2:i+N/2])[N/2-1]
			medians.append( (median_x, median_e) )
		
		by_type.append(energies[::25])
		print run_name, len(energies)
	f = open('out.txt', 'w')
	for row in range(len(by_type[0])):
		#if by_type[0][row][0]>1.0:
		#	continue
		for col in range(len(by_type)):
			f.write( ('%f\t') % by_type[col][row] )
		f.write('\n')

def smd():
	by_type = []
	for t in range(3,13):
		run_name = 'smd8_'+str(t)
		energies = []
		for line in open('lammps/'+run_name+'.log'):
			try:
				energies.append( float(line) )
			except ValueError: pass
		by_type.append(energies)
		print run_name, len(energies)
	f = open('out.txt', 'w')
	for row in range(len(by_type[0])):
		for col in range(len(by_type)):
			f.write('%f\t' % by_type[col][row])
		f.write('\n')

bend()

'''
for step in range(5):
	run_name = 'azoto_flat__'+str(step)
	energies = []
	for line in open('lammps/'+run_name+'.log'):
		try:
			a = line.split()
			if len(a)==3:
				energies.append( [float(s) for s in a] )
		except ValueError: pass
	energies = energies[1000:]
	avg_epair = sum([e[0] for e in energies])/len(energies)
	avg_pe = sum([e[1] for e in energies])/len(energies)
	avg_e = sum([e[2] for e in energies])/len(energies)
	print len(energies), step, avg_epair, avg_pe, avg_e
'''

'''
data = []
for ri in range(6):
	r = 20.0 + 10*ri
	run_name = 'azo_hex_r'+str(r)
	energies = []
	for line in open('lammps/'+run_name+'.log'):
		if 'atoms in group test_molecules' in line:
			n_molecules = int(line.split()[0])
		try:
			for s in line.split():
				energies.append( float(s) )
		except ValueError: pass
	energies = energies[100:]
	print len(energies), r, n_molecules, sum(energies)/len(energies), sum(energies)/len(energies)/n_molecules
	data.append( (n_molecules, sum(energies)/len(energies)) )

for i in range(len(data)-1):
	d_mol = data[i+1][0]-data[i][0]
	d_e = data[i+1][1]-data[i][1]
	
	print d_e/d_mol
'''
