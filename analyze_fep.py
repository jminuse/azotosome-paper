import os, string, threading, sys, shutil, math, cPickle, random, re
import filetypes, lammps, utils

kT = 0.199
n_steps = 100

os.chdir('lammps')

def run(molecule_name, solvent_ratio):
	#run_prefix = molecule_name+'_'+('_'.join([str(x) for x in solvent_ratio]))+'__'
	run_prefix = molecule_name+'_solid2__'
	files = [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith(run_prefix) and f.endswith('.log')]
	files.sort(key=lambda f: float(f[len(run_prefix):-4]))

	free_energies_forward = []
	free_energies_backward = []
	average_dHs = []

	for log_file in files:
		step = float(log_file[len(run_prefix):-4])
		if step==0: step = 1e-4
		strength = 1.0*step/n_steps
		sum_exp_forward = 0.0
		sum_exp_backward = 0.0
		sum_dH = 0.0
		count = 0
		skip = 10000
		for line in open(log_file):
			try:
				eng = float(line)
				if skip>0: skip-=1; continue
				E0 = eng/strength
				new_energy = E0/n_steps
				sum_exp_forward += math.exp( -new_energy/kT )
				sum_exp_backward += math.exp( new_energy/kT )
				sum_dH += E0 #if H = lambda * E0, then dH/d(lambda) = E0
				count += 1
			#except ValueError: pass
			except: pass
		#print count, 'samples'
		b = -kT * math.log(sum_exp_backward/count)
		f = -kT * math.log(sum_exp_forward/count)
		if abs(b)>kT*0.5 or abs(f)>kT*0.5 or abs(abs(b)-abs(f))>kT*0.1:
			pass
			#print 'Bad step %d: %.2f vs %.2f' % (step, b, f)
			#intermediate_steps = int(2*abs(abs(b)-abs(f))/kT)
			#print intermediate_steps
			#b = b/abs(b) * abs(b+f)/2
			#f = -b #CHEESY
		if int(step)>0-1 and int(step)<n_steps: free_energies_forward.append( f )
		if int(step)>1-1: free_energies_backward.append( b ) #skip first
		#if int(step)==n_steps: free_energies_backward.append( -kT * math.log(sum_exp_backward/count) )
		if int(step)>0-1: average_dHs.append(sum_dH/count)

	print molecule_name
	print '%.3f,' % sum(free_energies_forward),
	print '%.3f,' % sum([(average_dHs[i]+average_dHs[i+1])*0.5/n_steps for i in range(len(average_dHs)-1)]),
	print '%.3f' % sum(free_energies_backward)

	#print 'FEP: %.3f,' % sum(free_energies_forward),
	#print 'TI: %.3f,' % sum([(average_dHs[i]+average_dHs[i+1])*0.5/n_steps for i in range(len(average_dHs)-1)]),
	#print 'FEP backwards: %.3f' % sum(free_energies_backward)

#for molecule_name in ['propanenitrile','N2','acetylene','acrylonitrile','hexane','acetonitrile']:
for molecule_name in ['propanenitrile']:
	run(molecule_name, [100, 0, 0])


'''
for hydrocarbons_ch4_fraction_i in range(0,2+1):
	hydrocarbons_ch4_fraction = 1.0 - hydrocarbons_ch4_fraction_i*0.1
	for n2_percent in range(0,10+1,10):
		ch4_percent = int(round( hydrocarbons_ch4_fraction*(100-n2_percent) ))
		c2h6_percent = 100 - ch4_percent - n2_percent
		#print [ch4_percent, c2h6_percent, n2_percent]
		run([ch4_percent, c2h6_percent, n2_percent])
'''
