import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run_dissolved(test_molecule_filename, solvent_ratio=[100,0,0]):
	utils.Molecule.set_params('oplsaa4.prm')
	test_molecule = utils.Molecule(test_molecule_filename)
	test_molecule_name = test_molecule_filename[ test_molecule_filename.rindex('/')+1 :-4]
	if test_molecule_name[-2:]=='e2': test_molecule_name = test_molecule_name[:-1]

	N2 = utils.Molecule('N2.arc')
	#N2 params: http://www.sciencedirect.com/science/article/pii/002240739290142Q
	N2.atoms[0].type.vdw_r = 3.35
	N2.atoms[0].type.vdw_e = 0.0721
	#http://cccbdb.nist.gov/exp2.asp?casno=7727379
	N2.bonds[0].type.r = 1.0977

	solvents = [utils.Molecule('methane.arc'), utils.Molecule('ethane.arc'), copy.deepcopy(N2)]

	water_types = []

	atoms = []
	bonds = []
	angles = []
	dihedrals = []

	test_molecule.add_to(0.0, 0.0, 0.0, atoms, bonds, angles, dihedrals)

	box_size = (20,20,20)
	atom_count = len(atoms)
	solvent_amount_to_add = [ratio for ratio in solvent_ratio]
	solvent = solvents[0] #cheesy
	solvent_spacing = [(max(solvent.atoms, key=lambda a:a.x).x-min(solvent.atoms, key=lambda a:a.x).x), (max(solvent.atoms, key=lambda a:a.y).y-min(solvent.atoms, key=lambda a:a.y).y), (max(solvent.atoms, key=lambda a:a.z).z-min(solvent.atoms, key=lambda a:a.z).z)]
	max_vdw_r = max(solvent.atoms, key=lambda a:a.type.vdw_r).type.vdw_r
	solvent_spacing = [max(x,max_vdw_r) for x in solvent_spacing]

	for x in range(-box_size[0]/2+solvent_spacing[0], box_size[0]/2, solvent_spacing[0]):
		for y in range(-box_size[1]/2+solvent_spacing[1], box_size[1]/2, solvent_spacing[1]):
			for z in range(-box_size[2]/2+solvent_spacing[2], box_size[2]/2, solvent_spacing[2]):
				try:
					for a in atoms[:atom_count]:
						if (a.x-x)**2 + (a.y-y)**2 + (a.z-z)**2 < max_vdw_r**2:
							raise AssertionError()
					while True:
						index_to_add = random.randrange(len(solvents))
						if solvent_amount_to_add[ index_to_add ] > 0:
							solvents[ index_to_add ].add_to(x, y, z, atoms, bonds, angles, dihedrals)
							solvent_amount_to_add[ index_to_add ] -= 1
							if sum(solvent_amount_to_add) == 0:
								solvent_amount_to_add = [ratio for ratio in solvent_ratio]
							break
				except AssertionError: pass

	T = 94.0
	P = 1.45

	filetypes.write_xyz('out', atoms)

	os.chdir('lammps')
	
	original_vdw_e = [a.type.vdw_e for a in test_molecule.atoms]
	original_charges = [a.type.charge for a in test_molecule.atoms]

	n_steps = 100
	for step in range(0,n_steps+1):
		for i,x in enumerate(test_molecule.atoms):
			a = atoms[i]
			a.type.vdw_e = original_vdw_e[i] * (1.0*step/n_steps)**2 #since LJ eps is geometric mean
			a.type.charge = original_charges[i] * step/n_steps
			a.charge = a.type.charge
	
		#run_name = test_molecule_name+'_'+('_'.join([str(x) for x in solvent_ratio]))+'__'+str(int(step))
		run_name = test_molecule_name+'_gas__'+str(int(step))
	
		atom_types = dict( [(t.type,True) for t in atoms] ).keys()
		atom_type_numbers = dict( [(t,i+1) for i,t in enumerate(atom_types)] )
		is_charged = any([a.charge!=0.0 for a in atoms])

		lammps.write_data_file_general(atoms, bonds, angles, dihedrals, box_size, run_name, atom_types=atom_types)
		pickle.dump((atoms, bonds, angles, dihedrals), open(run_name+'.pickle', 'w'))
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
	
		test_molecule_types = dict( [(t.type,True) for t in atoms[:len(test_molecule.atoms)]] ).keys()

		f.write('''thermo		0
dump	1 all xyz 10000 '''+run_name+'''.xyz

group test_molecule id <= '''+str(len(test_molecule.atoms))+'''
group others subtract all test_molecule

compute SolEng test_molecule group/group others
thermo_style 	custom c_SolEng pe
#thermo_style custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press vol c_SolEng tpcpu
#thermo_modify	line multi format float %14.6f

minimize 0.0 1.0e-8 1000 100000

velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 iso '''+str(P)+' '+str(P)+''' 1000
timestep 2.0
thermo 1
run 250000
''')
		f.close()
		if True: #start multiple processes
			os.system('nohup ~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
			print 'Running', run_name
		else:
			os.system('~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log' % (run_name, run_name) )
			sys.exit()
	os.chdir('..')


for filename in ['propanenitrile.arc', 'acrylonitrile2.arc']:
	run_dissolved(filename)

