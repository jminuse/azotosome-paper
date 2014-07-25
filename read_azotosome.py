import os, string, threading, sys, shutil, math, pickle, random, re, copy, cPickle
import numpy
import filetypes, lammps, utils

def bend():
	by_type = []
	for t in range(3,18):
		run_name = 'smd19_'+str(t)
		energies = []
		for line in open('lammps/'+run_name+'.log').read().split('push[6] push[7] zrest')[1].splitlines():
			try:
				a = line.split()
				if len(a)==3:
					energies.append( tuple([float(s) for s in a]) )
			except ValueError: pass
		#energies = [(abs(e[0]-energies[0][0]),e[1],e[2]) for e in energies]
		energies = [e[2] for e in energies]
		
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

def ti():
	by_type = []
	for t in range(3,18):
		run_name = 'ti0_'+str(t)
		average_energies = []
		average_forces = []
		for section in open('lammps/'+run_name+'.log').read().split('PotEng hold[4] hold[6] zrest \n')[1:]:
			energies = []
			forces = []
			for line in section.splitlines():
				try:
					a = line.split()
					if len(a)==4:
						energies.append( float(a[0]) )
						forces.append( float(a[1]) )
				except ValueError: pass
			average_energies.append( sum(energies)/len(energies) )
			average_forces.append( sum(forces)/len(forces) )
		
		by_type.append(average_forces)
		print run_name, len(energies)
	f = open('out.txt', 'w')
	for row in range(len(by_type[0])):
		#if by_type[0][row][0]>1.0:
		#	continue
		for col in range(len(by_type)):
			f.write( ('%f\t') % by_type[col][row] )
		f.write('\n')

def restart_ti():
	by_type = []
	for t in [3,4,5,6,7,8,9,10,11,14,16]:
		run_name = 'ti16_'+str(t)
		average_energies = []
		average_forces = []
		for step in range(-10,51):
			energies = []
			forces = []
			#for line in open('lammps/'+run_name+'__'+str(step)+'.log').read().split('average[ average[ T/CPU ')[1].splitlines():
			for line in open('lammps/'+run_name+'__'+str(step)+'.log').read().split('average[ average[ ')[1].splitlines():
				try:
					a = line.split()
					if len(a)==2:
						energies.append( float(a[0]) )
						forces.append( float(a[1]) )
				except ValueError: pass
			#energies = energies[2:-1]
			#forces = forces[2:-1]
			energies = energies[-50:-1]
			forces = forces[-50:-1]
			try:
				average_energies.append( sum(energies)/len(energies) )
				average_forces.append( sum(forces)/len(forces) )
			except:
				print 'lammps/'+run_name+'__'+str(step)+'.log'
			
			#average_energies.append( sorted(energies)[len(energies)/2] ) #median
		
		by_type.append( average_forces )
		# k from force
		#by_type.append([max(average_forces), (average_forces.index(max(average_forces))-2)*0.5, (average_energies.index(min(average_energies))-2)*0.5 ])
		# k from energy
		#by_type.append([average_energies[ average_forces.index(max(average_forces)) ], (average_forces.index(max(average_forces))-2)*0.5, (average_energies.index(min(average_energies))-2)*0.5 ])
		
		print run_name, len(energies)
	f = open('out.txt', 'w')
	for row in range(len(by_type[0])):
		for col in range(len(by_type)):
			f.write( ('%f\t') % by_type[col][row] )
		f.write('\n')

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

def umbrella():
	by_type = []
	for t in [3,4,5,6,7,8,9,10,11,14,16]:
		run_name = 'ti19_'+str(t)
		data = []
		for step in range(-10,51):
			print step
			count_lines = 0
			for line in open('lammps/'+run_name+'__'+str(step)+'.log'):
				try:
					a = line.split()
					if len(a)==4: # thermo_p hold[4] fz_pushe z_pushed
						count_lines+=1
						if count_lines > 50000: data.append( [float(s) for s in a]  )
				except ValueError: pass
		data.sort(key=lambda a:a[3])
		
		bins = []
		z = data[0][3]
		z_increment = 12.0/200
		count = 0
		sums = [0.0, 0.0, 0.0, 0.0]
		for i,d in enumerate(data):
			sums[0] += d[0]
			sums[1] += d[1]
			sums[2] += d[2]
			sums[3] += d[3]
			count += 1
			if d[3] > z+z_increment:
				sums = [s/count for s in sums]
				bins.append(sums + [count])
				#print bins[-1]
				#print count
				count = 0
				z += z_increment
				sums = [0.0, 0.0, 0.0, 0.0]
		
		print 'lammps/'+run_name+'__'+str(step)+'.log', len(data)
		
		by_type.append( bins )
		f = open(str(t)+'_data.txt', 'w')
		for b in bins:
			f.write( ('%f\t'*5+'\n') % tuple(b) )
		f.close()
		#print bins
		
	cPickle.dump(by_type, open('by_type.dat', 'w'))
	
	f = open('out.txt', 'w')
	f.write( str(by_type) )
	return
	
	for row in range(len(by_type[0])):
		for col in range(len(by_type)):
			#f.write( ('%f\t')*5 % tuple(by_type[col][row]) )
			f.write( str(by_type[col][row]) )
		f.write('\n')

def write_by_type():
	by_type = cPickle.load( open('by_type.dat') )
	f = open('out.txt', 'w')
	for t in by_type:
		#f.write( str(t)+'\n' )
		f.write( '\n'.join( [' '.join([str(x) for x in s]) for s in t] )+'\n' )
		f.write('\n\n')
	
import numpy
import matplotlib.pyplot as plt
by_type = cPickle.load( open('by_type.dat') )
types_by_index = [3,4,5,6,7,8,9,10,11,14,16]
type_names = ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane', 'aminobutane', 'aminopentane', 'aminohexane', 'hexane', 'acrylonitrile', 'acetonitrile']
z_offsets = [-4.5, 4.5, 5.0, 6.5, 2.5, 3.5, 5.0, 6.5, 9.0, -4.0, 2.0]

def potential_energy_():
	csv = open('data.csv', 'w')

	for type_index, t in enumerate(by_type):
		#if type_names[type_index] not in ('hexane'):
		#	continue
		min_PE = min([data[0] for data in t])
	
		print type_names[type_index], max([data[0] for data in t[3:-3]])-min([data[0] for data in t[3:-3]])
	
		PE = numpy.array([data[0]-min_PE for data in t])
		forces = numpy.array([data[2]-data[1] for data in t])
		zs = numpy.array([data[3] for data in t])
		
		csv.write(type_names[type_index].capitalize()+',,,\n')
		csv.write( '"z (Angstroms)","PE (kcal/mol)","F (kcal/mol-Angstrom)","Number of samples"\n' )
		for i in range(len(zs)):
			csv.write( '%g,%g,%g,%d\n' % (zs[i], PE[i], forces[i], t[i][-1]) )
		
		#plt.plot(zs, PE, label=type_names[type_index])
		#plt.legend()
		#plt.xlabel('Distance (Angstroms)')
		#plt.ylabel('Energy (kcal/mol)')
		#plt.show()

	sys.exit()

potential_energy_()

'''
for type_index, t in enumerate(by_type):
	forces = numpy.array([data[2]-data[1] for data in t])
	zs = numpy.array([data[3] for data in t])
	pmf = [0.0]
	for i in range(len(zs)-1):
		pmf.append( pmf[-1] + (forces[i]+forces[i+1])*0.5*(zs[i+1]-zs[i]) )
	plt.plot(zs, pmf, label=type_names[type_index])
	plt.legend()
	plt.show()
	print type_names[type_index], max(pmf)-min(pmf)

sys.exit()
'''

for type_index, t in enumerate(by_type):
	potential_energy = [data[0] for data in t]
	potential_energy_minimum = potential_energy.index( min(potential_energy) )
	direction_of_motion = abs(max([data[3] for data in t])) - abs(min([data[3] for data in t]))
	direction_of_motion = direction_of_motion/abs(direction_of_motion)
	#print direction_of_motion, potential_energy_minimum
	#continue
	
	#print t[0][3], t[200][3]
	#continue
	
	best_k = 0.0
	best_c = 0.0
	best_r2 = -1e10
	for z_index in range( len(t) ):
		end_point = z_index
		length_required = 2.0
		while end_point < len(t) and abs(t[end_point][3] - t[z_index][3]) < length_required:
			end_point += 1
		if end_point==len(t) or abs(t[end_point][3] - t[z_index][3]) < length_required:
			break
		
		#find out whether it's on the forward side of PE min?
		#if ( (z_index+end_point)/2 - potential_energy_minimum ) * direction_of_motion < 0.0:
		#	continue
		
		#if ( (z_index+end_point)/2 - len(t)*0.5 ) * direction_of_motion < 0.0:
		
		if False: #pulling out
			mid = (z_index+end_point)/2
			if direction_of_motion > 0.0 and ( mid < len(t)*0.3 or mid > len(t)*0.5):
				continue
			if direction_of_motion < 0.0 and ( mid > len(t)*0.7 or mid < len(t)*0.5):
				continue
		else: #pushing in
			mid = (z_index+end_point)/2
			if direction_of_motion > 0.0 and ( mid > len(t)*0.1 or mid > len(t)*0.4): #0.2 for propanenitrile, 0.1 for aminohexane
				continue
			if direction_of_motion < 0.0 and ( mid > len(t)*0.9 or mid < len(t)*0.6):
				continue
		
		data_points = t[z_index:end_point]
		forces = numpy.array([data[2]-data[1] for data in data_points])
		zs = numpy.array([data[3] for data in data_points])
		#zs_matrix = zs[:,numpy.newaxis] #numpy.newaxis forces y-intercept of zero - unwanted here
		zs_matrix = numpy.array([ (z, 1.0) for z in zs])
		result = numpy.linalg.lstsq(zs_matrix, forces)
		k, y_intercept = result[0]
		mean_force = sum(forces)/len(forces)
		r2 = 1 - result[1][0]/ sum( [(force-mean_force)**2 for force in forces] )
		
		if r2>best_r2 and k < 0.0: # k > 0.0 for pulling out, < 0.0 for pushing in
			best_k = k
			best_r2 = r2
			best_c = y_intercept
		
		'''
		forces = numpy.array([data[2]-data[1] for data in t])
		zs = numpy.array([data[3] for data in t])
		zs_matrix = zs[:,numpy.newaxis] #numpy.array([ (z, 1.0) for z in zs]) #0.0 forces y-intercept of zero?
		result = numpy.linalg.lstsq(zs_matrix, forces)
		k = result[0]
		mean_force = sum(forces)/len(forces)
		r2 = 1 - result[1][0]/ sum( [(force-mean_force)**2 for force in forces] )
		'''
	k = best_k
	r2 = best_r2
	c = best_c
	
	print type_names[type_index].capitalize(), k, r2
	#sys.exit()
	
	if True:
		forces = numpy.array([data[2]-data[1] for data in t])
		zs = numpy.array([data[3] for data in t])
		plt.title(type_names[type_index])
		plt.plot(zs, forces, 'o', label='Original data', markersize=10)
		plt.plot(zs, k*zs + c, 'r', label='Fitted line')
		plt.legend()
		plt.xlabel('z (Angstroms)')
		plt.ylabel('F (kcal/mol-Angstrom)')
		plt.ylim([min(forces)-1,max(forces)+1])
		plt.show()
		#plt.savefig('test.png')
		#sys.exit()
	

