import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run(solvent_ratio):
	utils.Molecule.set_params('oplsaa4.prm')

	tail = utils.Molecule('hexanenitrile.arc')
	solvent = utils.Molecule('methane.arc')

	for a in tail.atoms:
		a.x, a.z = a.z, a.x

	directory = 'lammps'
	os.chdir(directory)

	n_steps = 1
	for step in range(0,n_steps):

		atoms = []
		bonds = []
		angles = []
		dihedrals = []

		S = 3.5
		box_size = [20,20,20]
		theta = 0; z = 0.0

		print len(atoms)/len(tail.atoms)
		atom_count = len(atoms)
		solvent_spacing = [(max(solvent.atoms, key=lambda a:a.x).x-min(solvent.atoms, key=lambda a:a.x).x), (max(solvent.atoms, key=lambda a:a.y).y-min(solvent.atoms, key=lambda a:a.y).y), (max(solvent.atoms, key=lambda a:a.z).z-min(solvent.atoms, key=lambda a:a.z).z)]
		solvent_spacing = [x+1.0 for x in solvent_spacing]
		max_vdw_r = max(solvent.atoms, key=lambda a:a.type.vdw_r).type.vdw_r
		solvent_spacing = [max(x,max_vdw_r) for x in solvent_spacing]
		for x in utils.frange(-box_size[0]/2, box_size[0]/2-solvent_spacing[0], solvent_spacing[0]):
			for y in utils.frange(-box_size[1]/2, box_size[1]/2-solvent_spacing[1], solvent_spacing[1]):
				for z in utils.frange(-box_size[2]/2, box_size[2]/2-solvent_spacing[2], solvent_spacing[2]):
					try:
						for a in atoms[:atom_count]:
							if (a.x-x)**2 + (a.y-y)**2 + (a.z-z)**2 < max_vdw_r**2:
								raise Exception()
						solvent.add_to(x, y, z, atoms, bonds, angles, dihedrals)
					except: pass

		T = 94.0
		P = 1.45

		print (len(atoms)-atom_count)/len(solvent.atoms)
		filetypes.write_xyz('out', atoms)

		run_name = 'azoto_solv__'+str(int(step))
	
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
dump	1 all xyz 10000 '''+run_name+'''.xyz

thermo_style custom etotal ke temp pe ebond eangle edihed evdwl ecoul elong press lx ly lz tpcpu
thermo_modify	line multi format float %14.6f

minimize 0.0 1.0e-8 1000 100000
velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 500
timestep  4.0
neigh_modify check yes every 1 delay 0
run_style respa 4 2 2 2 inner 2 4.5 6.0 middle 3 7.0 8.0 outer 4
thermo		100
run 25000
unfix press
fix dynamics all nvt temp '''+str(T)+' '+str(T)+''' 100
thermo_style 	custom epair pe etotal
thermo 1
run 25000
''')
		f.close()
		if True: #start multiple processes
			os.system('nohup ~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
			print 'Running', run_name
		else:
			os.system('~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log' % (run_name, run_name) )
			sys.exit()
	os.chdir('..')

run([100, 0, 0])

