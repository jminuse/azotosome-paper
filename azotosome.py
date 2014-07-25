import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

utils.Molecule.set_params('oplsaa4.prm')
tail = utils.Molecule('hexanenitrile.arc')
solvent = utils.Molecule('methane.arc')

def run(R):
	atoms = []
	bonds = []
	angles = []
	dihedrals = []

	added_atoms = []

	S = 3.5
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
	
	added_atoms.append( utils.Struct( x=0., y=0., z=0., theta=0., phi=0., charge=1 ) )
	
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
				a.x, a.y, a.z = a.z, a.y, a.x
			c1,c2,s1,s2 = math.cos(theta),math.cos(phi+math.pi/2),math.sin(theta),math.sin(phi+math.pi/2)
			m = [ [c1*c2, -s1, c1*s2], [c2*s1, c1, s1*s2], [-s2, 0., c2] ]
			tail.rotate(m)
			tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
		for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]

	atom_count = len(atoms)

	box_size = [R*2+20]*3

	solvent_spacing = [(max(solvent.atoms, key=lambda a:a.x).x-min(solvent.atoms, key=lambda a:a.x).x), (max(solvent.atoms, key=lambda a:a.y).y-min(solvent.atoms, key=lambda a:a.y).y), (max(solvent.atoms, key=lambda a:a.z).z-min(solvent.atoms, key=lambda a:a.z).z)]
	solvent_spacing = [x+0.0 for x in solvent_spacing]
	max_vdw_r = max(solvent.atoms, key=lambda a:a.type.vdw_r).type.vdw_r
	solvent_spacing = [max(x,max_vdw_r) for x in solvent_spacing]
	for x in range(-box_size[0]/2+solvent_spacing[0], box_size[0]/2, solvent_spacing[0]):
		for y in range(-box_size[1]/2+solvent_spacing[1], box_size[1]/2, solvent_spacing[1]):
			for z in range(-box_size[2]/2+solvent_spacing[2], box_size[2]/2, solvent_spacing[2]):
				try:
					for a in atoms[:atom_count]:
						if (a.x-x)**2 + (a.y-y)**2 + (a.z-z)**2 < max_vdw_r**2:
							raise Exception()
					solvent.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				except: pass

	T = 94.0
	P = 1.5

	directory = 'lammps'
	if not os.path.isdir(directory):
		os.mkdir(directory)
	os.chdir(directory)

	run_name = 'azo_hex_r'+str(R)

	atom_types = dict( [(t.type,True) for t in atoms] ).keys()
	atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
	is_charged = any([a.charge!=0.0 for a in atoms])

	lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types)
	os.system('cp ../'+sys.argv[0]+' '+run_name+'.py')

	f = open(run_name+'.in', 'w')
	f.write('''units	real\natom_style	full #bonds, angles, dihedrals, impropers, charges\n''')

	f.write('pair_style lj/cut/coul/long 8.0\n')
	if bonds: f.write('bond_style harmonic\n')
	if angles: f.write('angle_style harmonic\n')
	if dihedrals: f.write('dihedral_style opls\n')
	f.write('kspace_style pppm 1.0e-3\n')

	f.write('''special_bonds lj/coul 0.0 0.0 0.5
read_data	'''+run_name+'''.data\n

neigh_modify exclude molecule all

group test_molecules id <= '''+str(atom_count)+'''

compute 1 all pair/local eng
compute 2 test_molecules reduce sum c_1

minimize 0.0 1.0e-8 1000 100000

velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix dynamics all nvt temp '''+str(T)+' '+str(T)+''' 100
timestep 4.0
run_style respa 4 2 2 2 inner 2 4.5 6.0 middle 3 7.0 8.0 outer 4
run 2500

thermo_style custom c_2
thermo 1

run 1000
''')
	f.close()
	if True: #start multiple processes
		os.system('nohup ~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
		print 'Running', run_name
	else:
		os.system('~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log' % (run_name, run_name) )
		sys.exit()
	os.chdir('..')

for ri in range(6):
	r = 20.0 + 10*ri
	run(r)

