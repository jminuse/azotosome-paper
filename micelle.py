import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run(run_name):
	utils.Molecule.set_params('oplsaa4.prm')

	tail = utils.Molecule('acrylonitrile2.arc')

	for a in tail.atoms:
		a.x, a.z = a.z, a.x-2.0
	
	atoms = []
	bonds = []
	angles = []
	dihedrals = []

	added_atoms = []
	S = 2.4
	R = 8.0
	for phi_i in range(math.pi*R/S):
		phi = S/R*0.5 + S/R*phi_i
		z = R*math.cos(phi)
		r = (R**2 - z**2)**0.5
		N_theta = int(math.pi*2*r/S)
		for theta_i in range(N_theta):
			theta = S/r + (math.pi*2)*theta_i/N_theta
			x, y = R*math.cos(theta)*math.sin(phi), R*math.sin(theta)*math.sin(phi)
		
			weighted_charge = 0.0
		
			for added in added_atoms:
				weight = 1/( (x-added.x)**2 + (y-added.y)**2 + (z-added.z)**2 )
				weighted_charge += weight*added.charge
		
			charge = 1 if weighted_charge<0 else -1
		
			added_atoms.append( utils.Struct( x=x, y=y, z=z, theta=theta, phi=phi, charge=charge ) )

	for added in added_atoms:
		x, y, z, theta, phi = added.x, added.y, added.z, added.theta, added.phi

		old_atom_positions = [(a.x, a.y, a.z) for a in tail.atoms]
		if added.charge>0:
			for a in tail.atoms:
				a.x, a.y, a.z = -a.z-7, a.y, -a.x
			c1,c2,s1,s2 = math.cos(theta),math.cos(phi+math.pi/2),math.sin(theta),math.sin(phi+math.pi/2)
			m = [ [c1*c2, -s1, c1*s2], [c2*s1, c1, s1*s2], [-s2, 0., c2] ]
			tail.rotate(m)
			tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
		else:
			for a in tail.atoms:
				a.x, a.y, a.z = -a.z-7, a.y, -a.x
			c1,c2,s1,s2 = math.cos(theta),math.cos(phi+math.pi/2),math.sin(theta),math.sin(phi+math.pi/2)
			m = [ [c1*c2, -s1, c1*s2], [c2*s1, c1, s1*s2], [-s2, 0., c2] ]
			tail.rotate(m)
			x*=0.7; y*-0.7; z*=0.7
			tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
		for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]

	box_size = [R*2+30]*3

	T = 94.0
	P = 1.45

	filetypes.write_xyz('out', atoms)
	sys.exit()

	directory = 'lammps'
	os.chdir(directory)

	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	is_charged = True

	lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types)
	os.system('cp ../'+sys.argv[0]+' '+run_name+'.py')

	f = open(run_name+'.in', 'w')
	f.write('units	real\natom_style	full #bonds, angles, dihedrals, impropers, charges\n')

	if is_charged:
	#	f.write('pair_style lj/cut/coul/long 8.0\n')
		f.write('pair_style lj/cut/coul/cut 8.0\n')
	else:
		f.write('pair_style lj/cut 8.0\n')
	if bonds: f.write('bond_style harmonic\n')
	if angles: f.write('angle_style harmonic\n')
	if dihedrals: f.write('dihedral_style opls\n')
	#if is_charged: f.write('kspace_style pppm 1.0e-3\n')

	f.write('special_bonds lj/coul 0.0 0.0 0.5\nread_data	'+run_name+'.data\n')

	f.write('''thermo		0
dump	1 all xyz 100 '''+run_name+'''.xyz
thermo_modify	line multi format float %14.6f
thermo 1000
minimize 0.0 1.0e-8 1000 100000
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix dynamics all nve
fix solvent all langevin '''+str(T)+' '+str(5*T)+''' 100 1337
timestep  2.0
neigh_modify check yes every 1 delay 0
run 100
	''')
	f.close()
	if False: #start multiple processes
		os.system('nohup ~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
		print 'Running', run_name
	else:
		os.system('~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log' % (run_name, run_name) )
		sys.exit()
	os.chdir('..')

run('micelle_demo')

