import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run():
	utils.Molecule.set_params('oplsaa4.prm')
	solvent = utils.Molecule('methane.arc')

	if False:
		tail = utils.Molecule('propanenitrile.arc')
		z_off_1 = 3.0
		z_off_2 = -5.5
	else:
		tail = utils.Molecule('hexanenitrile.arc')
		for a in tail.atoms:
			a.x, a.z = a.z, a.x
		z_off_1 = 0.0
		z_off_2 = 6.5

	directory = 'lammps'
	os.chdir(directory)

	n_steps = 1
	for step in range(0,n_steps):

		atoms = []
		bonds = []
		angles = []
		dihedrals = []

		S = 3.5
		box_size = [21.,21.,35.]
		theta = 0; z = 0.0
		molecules_not_added = 0
		
		held_Ns = []
		pushed_Ns = []
		
		for xi in range(box_size[0]/S):
			x = xi*S - box_size[0]/2
			for yi in range(box_size[1]/S):
				y = yi*S - box_size[1]/2
				if (xi%2)==(yi%2):
					tail.add_to(x, y, z+z_off_1, atoms, bonds, angles, dihedrals)
				else:
					old_atom_positions = [(a.x, a.y, a.z) for a in tail.atoms]
					for a in tail.atoms:
						a.z = -a.z+z_off_2
						a.y = -a.y
					tail.add_to(x, y, z, atoms, bonds, angles, dihedrals)
					for i,a in enumerate(tail.atoms): a.x, a.y, a.z = old_atom_positions[i]
				if (xi,yi) in [(3,3)]:
					for a in atoms[::-1]:
						if a.element=='N':
							pushed_Ns.append(a)
							break
				elif (xi,yi) in [(0,0)]:
					for a in atoms[::-1]:
						if a.element=='N':
							held_Ns.append(a)
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
		T_hot = T*4.2
		P = 1.45

		#filetypes.write_xyz('out', atoms)
		#sys.exit(0)

		run_name = 'azoto_hot__0'
	
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
dump	1 all xyz 1000 '''+run_name+'''.xyz

thermo_style custom etotal ke temp pe ebond eangle edihed epair press lx ly lz tpcpu
thermo_modify	line multi format float %14.6f

minimize 0.0 1.0e-8 1000 100000

timestep  4.0
neigh_modify check yes every 1 delay 0
run_style respa 4 2 2 2 inner 2 4.5 6.0 middle 3 7.0 8.0 outer 4
thermo		1000
fix relax all npt temp 10.0 10.0 50 aniso '''+str(P)+' '+str(P)+''' 500
run 1000
unfix relax
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
run 1000
group held_Ns id '''+(' '.join([str(N.index) for N in held_Ns]))+'''
group pushed_Ns id '''+(' '.join([str(N.index) for N in pushed_Ns]))+'''

#fix push1 held_Ns smd cvel 20.0 -0.00001 tether NULL NULL 100.0 0.0
#fix push2 pushed_Ns smd cvel 20.0 0.00001 tether NULL NULL 100.0 0.0

unfix press
group membrane id <= '''+str(atom_count)+'''
group solvent subtract all membrane

fix press1 membrane npt temp '''+str(T_hot)+' '+str(T_hot)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
fix press2 solvent nvt temp '''+str(T_hot)+' '+str(T_hot)+''' 100

#thermo_style 	custom epair f_push1[7] f_push2[7]
thermo 100
run 10000
''')
		f.close()
		if False: #start multiple processes
			os.system('nohup ~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
			print 'Running', run_name
		else:
			os.system('~/lammps2/src/lmp_g++ -in %s.in -log %s.log' % (run_name, run_name) )
			sys.exit()
	os.chdir('..')

run()

