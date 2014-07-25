import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

utils.Molecule.set_params('oplsaa4.prm')
solvent = utils.Molecule('methane.arc')

def run(mode, timescale):
	run_name = 'bend11_'+str(mode)

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
		z_offset = -3.5
	elif mode==13:
		tail = utils.Molecule('cyanoacetylene2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 0.0
	elif mode==14:
		tail = utils.Molecule('acrylonitrile2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = -4.0
	elif mode==15:
		tail = utils.Molecule('cyanoallene2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = 4.0
	elif mode==16:
		tail = utils.Molecule('acetonitrile2.arc')
		z_offset = 0.0
	elif mode==17:
		tail = utils.Molecule('hc5n2.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_offset = -4.0

	directory = 'lammps'
	os.chdir(directory)

	atoms = []
	bonds = []
	angles = []
	dihedrals = []

	S_count = 12
	S = 3.5
	box_size = [ S_count*S ]*3
	theta = 0; z = 0.0
	molecules_not_added = 0
	
	pushable = []
	
	if mode==11:
		for xi in range(S_count):
			x = xi*S - box_size[0]/2
			for yi in range(S_count):
				y = yi*S - box_size[1]/2

				tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)

				push_r = S * ( (xi-S_count/2)**2 + (yi-S_count/2)**2 )**0.5
				if xi==0 or yi==0:
					push_r = 100.0
					atoms[-15].push_r = push_r
					pushable.append( atoms[-15] )

				old_atom_positions = [(a.x, a.y, a.z) for a in tail.atoms]
				for a in tail.atoms:
					a.z = -a.z+z_offset
					a.y = -a.y
				tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]
				
				atoms[-15].push_r = push_r
				pushable.append( atoms[-15] )
	else:
		for xi in range(S_count):
			x = xi*S - box_size[0]/2
			for yi in range(S_count):
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
						a.push_r = 100.0 if (xi==0 or yi==0) else S * ( (xi-S_count/2)**2 + (yi-S_count/2)**2 )**0.5
						pushable.append(a)
						break

	#for a in pushable:
	#	a.element = str(int(a.push_r))
	#for a in atoms:
	#	print a.element
	#filetypes.write_xyz('out', atoms)
	#sys.exit(0)

	print len(atoms)/len(tail.atoms)
	atom_count = len(atoms)
	solvent_spacing = [(max(solvent.atoms, key=lambda a:a.x).x-min(solvent.atoms, key=lambda a:a.x).x), (max(solvent.atoms, key=lambda a:a.y).y-min(solvent.atoms, key=lambda a:a.y).y), (max(solvent.atoms, key=lambda a:a.z).z-min(solvent.atoms, key=lambda a:a.z).z)]
	solvent_spacing = [x+1.0 for x in solvent_spacing]
	max_vdw_r = max(solvent.atoms, key=lambda a:a.type.vdw_r).type.vdw_r
	solvent_spacing = [max(x,max_vdw_r) for x in solvent_spacing]
	for xi in range(int(box_size[0]/solvent_spacing[0])):
		for yi in range(int(box_size[1]/solvent_spacing[1])):
			for zi in range(int(box_size[2]/solvent_spacing[2])):
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

	f.write('''thermo		0
dump	1 all xyz '''+str(int(10000*timescale))+''' '''+run_name+'''.xyz

thermo_style custom etotal ke temp pe ebond eangle edihed epair press lx ly lz tpcpu
thermo_modify	line multi format float %14.6f

minimize 0.0 1.0e-8 1000 100000

timestep  4.0
neigh_modify check yes every 1 delay 0
run_style respa 4 2 2 2 inner 2 4.5 6.0 middle 3 7.0 8.0 outer 4
thermo		1000
fix relax all npt temp 10.0 10.0 50 aniso '''+str(P)+' '+str(P)+''' 500
run    1000
unfix relax
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
run '''+str(int(25000*timescale))+'''
''')
	push_groups = [[a for a in pushable if a.push_r<=R and a.push_r>(R-0.01)] for R in [0.01*ri for ri in range(0,1000)]]
	push_groups = [p for p in push_groups if p]
	for p_i,p in enumerate(push_groups):
		f.write('group push_g'+str(p_i)+' id '+(' '.join([str(a.index) for a in p]))+'\n')
		r = p[0].push_r
		max_r = 10.0
		v = 0.5 * 1e-6*( max_r**2 - r**2 )**0.5 / timescale
		print r, v, 250000*timescale*4*v
		f.write('fix push_%d push_g%d smd cvel 20.0 %e tether NULL NULL 100.0 0.0\n' % (p_i, p_i, v) )
	
	f.write('group held_atoms id '+(' '.join([str(a.index) for a in pushable if a.push_r==100.0]))+'\nfix zrest held_atoms spring/self 20.0 z\n')
	f.write('thermo_style 	custom ' + ( ' '.join(['f_push_%d[6] f_push_%d[7]' % (p_i,p_i) for p_i,p in enumerate(push_groups)]) ) + '\n')

	f.write('''
thermo  '''+str(int(1000*timescale))+'''
run '''+str(int(250000*timescale))+'''
''')
	f.close()
	if True: #start multiple processes
		os.system('nohup ~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
		print 'Running', run_name
	else:
		os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (run_name, run_name) )
		sys.exit()
	os.chdir('..')

for timescale in [1.0]:
	#for mode in range(3,18):
	for mode in range(12,13):
		run(mode, timescale)

