import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run(mode, timescale, step, step_size, makefile):
	run_name = 'ti10_'+str(mode)+'__'+str(step)
	z_separation = step * step_size
	solvent = utils.Molecule('methane.arc')
	#solvent = utils.Molecule('benzene2.arc')
	#solvent = utils.Molecule('hexane2.arc')
	box_size = [21.,21.,42.]
	S = 3.5

	if mode==3:
		tail = utils.Molecule('propanenitrile.arc')
		z_offset = -4.5
	elif mode==4:
		tail = utils.Molecule('butanenitrile.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 4.5
	elif mode==5:
		tail = utils.Molecule('pentanenitrile.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 5.0
	elif mode==6:
		tail = utils.Molecule('hexanenitrile.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 6.5
	elif mode==7:
		tail = utils.Molecule('aminopropane2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 2.5
	elif mode==8:
		tail = utils.Molecule('aminobutane2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, -a.x
		z_offset = 3.5
	elif mode==9:
		tail = utils.Molecule('aminopentane2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 5.0
	elif mode==10:
		tail = utils.Molecule('aminohexane2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 6.5
	elif mode==11:
		tail = utils.Molecule('hexane2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 9.0
	elif mode==12:
		tail = utils.Molecule('hcn2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		S = 3.0
		z_offset = -0.5
	elif mode==13:
		tail = utils.Molecule('cyanoacetylene2.arc') #unfixed dihedrals irrelevant, since angle=180
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = -2.0
	elif mode==14:
		tail = utils.Molecule('acrylonitrile2.arc') #unfixed dihedrals irrelevant, since angle=180
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = -4.0
	elif mode==15:
		tail = utils.Molecule('cyanoallene2.arc') #added one dihedral param (18   19   37  100), other two free, ( 36   37  100   37 ) added too, or modified? Basically, doesn't work in OPLS; allene dihedrals are wrong. 
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 4.0
	elif mode==16:
		tail = utils.Molecule('acetonitrile2.arc') #all dihedrals free
		z_offset = 2.0
	elif mode==17:
		tail = utils.Molecule('hc5n2.arc') #doesn't work well in OPLS
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = -4.0

	atoms = []
	bonds = []
	angles = []
	dihedrals = []

	theta = 0; z = 0.0
	molecules_not_added = 0
	
	held_atoms = []
	pushed_atoms = []
	
	if mode==11:
		for xi in range(box_size[0]/S):
			x = xi*S - box_size[0]/2
			for yi in range(box_size[1]/S):
				y = yi*S - box_size[1]/2

				#tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				if xi==3 and yi==3:
					tail.add_to(x, y, z+z_separation*(1 if z_offset<0 else -1), atoms, bonds, angles, dihedrals)
					pushed_atoms.append( atoms[-15] )
				else:
					tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
					if xi==0 or yi==0:
						held_atoms.append( atoms[-15] )

				old_atom_positions = [(a.x, a.y, a.z) for a in tail.atoms]
				for a in tail.atoms:
					a.z = -a.z+z_offset
					a.y = -a.y
				tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]
	else:
		for xi in range(box_size[0]/S):
			x = xi*S - box_size[0]/2
			for yi in range(box_size[1]/S):
				y = yi*S - box_size[1]/2
				if (xi%2)==(yi%2):
					if xi==3 and yi==3:
						tail.add_to(x, y, z+z_separation*(1 if z_offset<0 else -1), atoms, bonds, angles, dihedrals)
					else:
						tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				else:
					old_atom_positions = [(a.x, a.y, a.z) for a in tail.atoms]
					for a in tail.atoms:
						a.z = -a.z+z_offset
						a.y = -a.y
					tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
					for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]
				
				for a in atoms[::-1]:
					if a.element=='N':
						if xi==3 and yi==3:
							pushed_atoms.append( a )
						elif xi==0 or yi==0:
							held_atoms.append( a )
						break

	#filetypes.write_xyz('', atoms, xyz_file)
	print pushed_atoms[0].z
	return

	print len(atoms)/len(tail.atoms)
	atom_count = len(atoms)
	solvent_spacing = [(max(solvent.atoms, key=lambda a:a.x).x-min(solvent.atoms, key=lambda a:a.x).x), (max(solvent.atoms, key=lambda a:a.y).y-min(solvent.atoms, key=lambda a:a.y).y), (max(solvent.atoms, key=lambda a:a.z).z-min(solvent.atoms, key=lambda a:a.z).z)]
	solvent_spacing = [x+1.0 for x in solvent_spacing]
	max_vdw_r = max(solvent.atoms, key=lambda a:a.type.vdw_r).type.vdw_r
	solvent_spacing = [max(x,max_vdw_r) for x in solvent_spacing]
	for xi in range(box_size[0]/solvent_spacing[0]):
		for yi in range(box_size[1]/solvent_spacing[1]):
			for zi in range(box_size[2]/solvent_spacing[2]):
				x, y, z = xi*solvent_spacing[0]-box_size[0]/2,  yi*solvent_spacing[1]-box_size[1]/2,  zi*solvent_spacing[2]-box_size[2]/2
				try:
					for a in atoms[:atom_count]:
						if (a.x-x)**2 + (a.y-y)**2 + (a.z-z)**2 < max_vdw_r**2:
							raise Exception()
					solvent.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				except: pass

	T = 94.0
	P = 1.45

	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	is_charged = True

	directory = 'lammps'
	os.chdir(directory)

	lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types)
	os.system('cp ../'+sys.argv[0]+' '+run_name+'.py')

	f = open(run_name+'.in', 'w')
	f.write('units	real\natom_style	full #bonds, angles, dihedrals, impropers, charges\n')

	if is_charged:
		f.write('pair_style lj/cut/coul/long 8.0\n')
	else:
		f.write('pair_style lj/cut 8.0\n')
	if bonds: f.write('bond_style harmonic\n')
	if angles: f.write('angle_style harmonic\n')
	if dihedrals: f.write('dihedral_style opls\n')
	if is_charged: f.write('kspace_style pppm 1.0e-3\n')

	f.write('special_bonds lj/coul 0.0 0.0 0.5\nread_data	'+run_name+'.data\n')

	f.write('''
dump	1 all xyz '''+str(int(10000*timescale))+''' '''+run_name+'''.xyz

thermo_style custom etotal ke temp pe ebond eangle edihed epair press lx ly lz tpcpu
thermo_modify	line multi format float %14.6f

group held_atoms id '''+(' '.join([str(a.index) for a in held_atoms]))+'''
group pushed_atoms id '''+(' '.join([str(a.index) for a in pushed_atoms]))+'''

fix no_z held_atoms planeforce 0.0 0.0 1.0
fix hold pushed_atoms smd cvel 100.0 0.0 tether NULL NULL 100.0 0.0

fix no_z_temp pushed_atoms planeforce 0.0 0.0 1.0 #since smd is not enforced during minimization
minimize 0.0 1.0e-8 1000 100000
unfix no_z_temp

timestep  2.0
fix relax all npt temp 10.0 10.0 50 aniso '''+str(P)+' '+str(P)+''' 500
run 1000
unfix relax
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
velocity held_atoms set 0.0 0.0 0.0
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000

fix average all ave/time 1 1000 1000 c_thermo_pe f_hold[4] f_hold[6]
thermo 1000
thermo_style 	custom f_average[1] f_average[2] f_average[3]

run '''+str(int(50000*timescale))+'''
write_restart '''+run_name+'.restart\n')

	f.close()
	if True: #start multiple processes
		#os.system('nohup ~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
		#print 'Running', run_name
		makefile.write( '%s:\n\t~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null\n' % (run_name, run_name, run_name)  )
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (run_name, run_name) )
		sys.exit()
	os.chdir('..')


def restart(mode, timescale, step, step_size, makefile):
	os.chdir('lammps')
	
	energies = []
	zs = []
	for line in open('ti18_'+str(mode)+'__0.log').read().split('average[ average[ ')[1].splitlines():
		try:
			a = line.split()
			if len(a)==2:
				energies.append( float(a[0]) )
				zs.append( float(a[1]) )
		except ValueError: pass
	zs = zs[len(zs)/2:-1]
	average_z = sum(zs)/len(zs)
	old_run_name = 'ti18_'+str(mode)+'__0'
	z_sign = {3:-1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1, 12:-1, 13:-1, 14:-1, 16:1}
	z_separation = step * step_size * -z_sign[mode] + average_z
	T = 94.0
	P = 1.45
	new_run_name = 'ti19_'+str(mode)+'__'+str(step)
	f = open(new_run_name+'.in', 'w')
	f.write('''read_restart '''+old_run_name+'''.restart
kspace_style pppm 1.0e-3
fix no_z held_atoms planeforce 0.0 0.0 1.0
set group pushed_atoms z '''+str(z_separation)+'''
fix hold pushed_atoms smd cvel 100.0 0.0 tether NULL NULL 100.0 0.0
fix relax all nve/limit 0.01
fix drag all viscous 0.05

compute z_pushed_atoms pushed_atoms property/atom z
compute z_pushed all reduce sum c_z_pushed_atoms
compute fz_pushed_atoms pushed_atoms property/atom fz
compute fz_pushed all reduce sum c_fz_pushed_atoms
#fix average all ave/time 1 1000 1000 c_thermo_pe f_hold[4] c_fz_pushed c_z_pushed
thermo 1
#thermo_style 	custom f_average[1] f_average[2] f_average[3] f_average[4]
thermo_style custom c_thermo_pe f_hold[4] c_fz_pushed c_z_pushed
#dump	1 all xyz 1000 '''+new_run_name+'''.xyz

timestep 2.0
run 2500
unfix relax
fix relax all nve
run 25000
unfix relax
unfix drag
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
velocity held_atoms set 0.0 0.0 0.0
run 250000
write_restart '''+new_run_name+'.restart\n')
	f.close()
	if True: #start multiple processes
		makefile.write( '%s:\n\t~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null\n' % (new_run_name, new_run_name, new_run_name)  )
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (new_run_name, new_run_name) )
		sys.exit()
	
	os.chdir('..')

def restart_minimized(mode, timescale, step, step_size):
	os.chdir('lammps')
	old_run_name = 'ti12_'+str(mode)+'__'+str(step)
	T = 94.0
	P = 1.45
	new_run_name = 'ti18_'+str(mode)+'__'+str(step)
	f = open(new_run_name+'.in', 'w')
	f.write('''read_restart '''+old_run_name+'''.restart
kspace_style pppm 1.0e-3
fix no_z held_atoms planeforce 0.0 0.0 1.0
fix press all npt temp '''+str(T)+''' '''+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
compute z_pushed_atoms pushed_atoms property/atom z
compute z_pushed all reduce sum c_z_pushed_atoms
fix average all ave/time 1 1000 1000 c_thermo_pe c_z_pushed
thermo 1000
thermo_style 	custom f_average[1] f_average[2]
timestep 2.0
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
velocity held_atoms set 0.0 0.0 0.0
run 250000
write_restart '''+new_run_name+'.restart\n')
	f.close()
	if True: #start multiple processes
		os.system('nohup ~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null &' % (new_run_name, new_run_name) )
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (new_run_name, new_run_name) )
		sys.exit()
	os.chdir('..')

def view_xyz(mode, timescale, step, step_size):
	#os.chdir('lammps')
	old_run_name = 'ti15_'+str(mode)+'__'+str(step)
	z_sign = {3:-1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1, 12:-1, 13:-1, 14:-1, 16:1}
	z_separation = step * step_size * -z_sign[mode]
	
	print z_separation
	return
	
	T = 94.0
	P = 1.45
	new_run_name = 'ti14_'+str(mode)+'__'+str(step)
	f = open(new_run_name+'.in', 'w')
	f.write('''read_restart '''+old_run_name+'''.restart
kspace_style pppm 1.0e-3
dump	1 all xyz 1 '''+new_run_name+'''.xyz
run 0''')
	f.close()
	if False: #start multiple processes
		makefile.write( '%s:\n\t~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null\n' % (new_run_name, new_run_name, new_run_name)  )
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (new_run_name, new_run_name) )
		#sys.exit()
	
	os.chdir('..')

def view_dimensions(mode, step, step_size):
	os.chdir('lammps')
	old_run_name = 'ti19_'+str(mode)+'__'+str(step)
	z_sign = {3:-1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1, 12:-1, 13:-1, 14:-1, 16:1}
	z_separation = step * step_size * -z_sign[mode]
	
	T = 94.0
	P = 1.45
	new_run_name = 'ti20_'+str(mode)+'__'+str(step)
	f = open(new_run_name+'.in', 'w')
	f.write('''read_restart '''+old_run_name+'''.restart
kspace_style pppm 1.0e-3
fix no_z held_atoms planeforce 0.0 0.0 1.0
fix press all npt temp '''+str(T)+''' '''+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
compute z_pushed_atoms pushed_atoms property/atom z
compute z_pushed all reduce sum c_z_pushed_atoms
variable xbox equal lx
variable ybox equal ly
fix average all ave/time 1 1000 1000 v_xbox v_ybox
thermo 1000
thermo_style 	custom f_average[1] f_average[2]
timestep 2.0
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
velocity held_atoms set 0.0 0.0 0.0
run 2000''')
	f.close()
	if False: #start multiple processes
		makefile.write( '%s:\n\t~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null\n' % (new_run_name, new_run_name, new_run_name)  )
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (new_run_name, new_run_name) )
		#sys.exit()
		log = open(new_run_name+'.log').read()
		match = re.search('average\[ average\[\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)', log)
		lx, ly = match.group(1), match.group(2)
	
	os.chdir('..')
	return lx, ly

'''
makefile = open('lammps/Makefile19', 'w')
makefile.write('all :')
for mode in [3,4,5,6,7,8,9,10,11,14,16]:
	for step in range(-10,51):
		run_name = 'ti19_'+str(mode)+'__'+str(step)
		makefile.write(' %s' % run_name )
makefile.write('\n\n')
'''

'''
utils.Molecule.set_params('oplsaa4.prm')
os.chdir('lammps')
for mode in [3]:
	s = 'cat '
	for step in range(-10,51):
		#view_xyz(mode, 1.0, step, 0.2)
		old_run_name = 'ti14_'+str(mode)+'__'+str(step)+'.xyz'
		s = s + ' ' + old_run_name
	s = s + '> out.xyz'
	os.system(s)
'''

'''
utils.Molecule.set_params('oplsaa4.prm')
for mode in [3,4,5,6,7,8,9,10,11,14,16]:
	for step in range(-10,51):
		restart(mode, 1.0, step, 0.2, makefile)
'''

mode = 3
step = 0
step_size = 0.2

#print view_dimensions(mode, step, step_size)
#sys.exit()

dimensions = []
for step in range(-10,0):
	dimensions.append( view_dimensions(mode, step, step_size) )
for d in dimensions:
	print ','.join(d)

