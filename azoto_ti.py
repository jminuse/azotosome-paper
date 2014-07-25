import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run(mode, timescale):
	run_name = 'ti6_'+str(mode)
	utils.Molecule.set_params('oplsaa4.prm')
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
		z_offset = 8.0
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

	directory = 'lammps'
	os.chdir(directory)

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

				tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)

				old_atom_positions = [(a.x, a.y, a.z) for a in tail.atoms]
				for a in tail.atoms:
					a.z = -a.z+z_offset
					a.y = -a.y
				tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]
				
				if xi==3 and yi==3:
					pushed_atoms.append( atoms[-15] )
				elif xi==0 or yi==0:
					held_atoms.append( atoms[-15] )
	else:
		for xi in range(box_size[0]/S):
			x = xi*S - box_size[0]/2
			for yi in range(box_size[1]/S):
				y = yi*S - box_size[1]/2
				if (xi%2)==(yi%2):
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

	#filetypes.write_xyz(run_name, atoms)
	#os.chdir('..')
	#return

	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	is_charged = True

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

minimize 0.0 1.0e-8 1000 100000

timestep  2.0
fix relax all npt temp 10.0 10.0 50 aniso '''+str(P)+' '+str(P)+''' 500
run 1000
unfix relax
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
run '''+str(int(50000*timescale))+'''
group held_atoms id '''+(' '.join([str(a.index) for a in held_atoms]))+'''
group pushed_atoms id '''+(' '.join([str(a.index) for a in pushed_atoms]))+'''

fix zrest held_atoms spring/self 100.0 z

fix push pushed_atoms smd cvel 100.0 '''+str((1 if z_offset<0 else -1)*0.00001/timescale)+''' tether NULL NULL 100.0 0.0
run '''+str(int(100000*timescale))+''' #push in 2 Angstroms, then pull out 12
fix push pushed_atoms smd cvel 100.0 '''+str((-1 if z_offset<0 else 1)*0.00001/timescale)+''' tether NULL NULL 100.0 0.0
write_restart '''+run_name+'''.0.restart
''')
	for step in range(1,121):
		f.write('run '+str(int(5000*timescale))+'''
write_restart '''+run_name+'.'+str(step)+'''.restart
''')

	f.close()
	if True: #start multiple processes
		os.system('nohup ~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
		print 'Running', run_name
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (run_name, run_name) )
		sys.exit()
	os.chdir('..')

makefile = open('lammps/Makefile', 'w')
makefile.write('all :')
for mode in range(3,17):
	for step in range(0,121,2):
		makefile.write(' %d_%d' % (mode,step) )
makefile.write('\n\n')

def restart(mode):
	os.chdir('lammps')
	old_run_name = 'ti6_'+str(mode)
	T = 94.0
	P = 1.45
	for step in range(0,121,2):
		new_run_name = 'ti8_'+str(mode)+'__'+str(step)
		f = open(new_run_name+'.in', 'w')
		f.write('''read_restart '''+old_run_name+'.'+str(step)+'''.restart
kspace_style pppm 1.0e-3
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
fix zrest held_atoms spring/self 100.0 z
fix hold pushed_atoms smd cvel 100.0 0.0 tether NULL NULL 100.0 0.0
fix average all ave/time 1 1000 1000 c_thermo_pe f_hold[4] f_hold[6] f_zrest
thermo 1000
thermo_style 	custom f_average[1] f_average[2] f_average[3] f_average[4]
velocity all create '''+str(T)+''' 9999 rot yes dist gaussian
timestep 2.0
run 100000
''')
		f.close()
		if True: #start multiple processes
			#os.system('nohup ~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null &' % (new_run_name, new_run_name) )
			#print 'Running', new_run_name
			makefile.write( '%d_%d:\n\t~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null\n' % (mode, step, new_run_name, new_run_name)  )
		else:
			os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (new_run_name, new_run_name) )
			sys.exit()
	
	os.chdir('..')

for mode in range(3,17):
	restart(mode)

#for mode in range(3,17):
#	run(mode, 1.0)

