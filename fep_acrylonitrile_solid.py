import os, string, threading, sys, shutil, math, pickle, random, re, copy
import filetypes, lammps, utils

def run(test_molecule_filename, makefile):
	utils.Molecule.set_params('oplsaa4.prm')
	test_molecule = utils.Molecule(test_molecule_filename)
	test_molecule_name = test_molecule_filename[ test_molecule_filename.rindex('/')+1 :-4]
	if test_molecule_name[-2:]=='e2': test_molecule_name = test_molecule_name[:-1]
	solvent = copy.deepcopy(test_molecule)

	solvent_spacing = [(max(solvent.atoms, key=lambda a:a.x).x-min(solvent.atoms, key=lambda a:a.x).x), (max(solvent.atoms, key=lambda a:a.y).y-min(solvent.atoms, key=lambda a:a.y).y), (max(solvent.atoms, key=lambda a:a.z).z-min(solvent.atoms, key=lambda a:a.z).z)]
	solvent_spacing = [x+2.0 for x in solvent_spacing]
	max_vdw_r = max(solvent.atoms, key=lambda a:a.type.vdw_r).type.vdw_r
	solvent_spacing = [max(x,max_vdw_r) for x in solvent_spacing]

	solvent_spacing[0] += 0.5
	solvent_spacing[1] -= 1.0
	solvent_spacing[2] -= 0.5

	L = 6
	box_size = [L*x for x in solvent_spacing]

	atoms = []
	bonds = []
	angles = []
	dihedrals = []
	
	#add staggering between layers (like shingles)
	
	for xi in range(L):
		for yi in range(L):
			for zi in range(L):
				x = (xi+zi*0.4)*solvent_spacing[0] - box_size[0]/2
				y = yi*solvent_spacing[1] - box_size[1]/2
				z = zi*solvent_spacing[2] - box_size[2]/2
				
				if xi==0 and yi==0 and zi==0:
					test_molecule.add_to(x, y, z, atoms, bonds, angles, dihedrals)
				else:
					if True: #(xi+yi+zi)%2==0:
						solvent.add_to(x, y, z, atoms, bonds, angles, dihedrals)
					else:
						for a in solvent.atoms:
							a.z = -a.z
						solvent.add_to(x, y, z, atoms, bonds, angles, dihedrals)
						for a in solvent.atoms:
							a.z = -a.z

	T = 94.0
	P = 1.45

	#filetypes.write_xyz('out', atoms)
	#sys.exit()

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
	
		run_name = test_molecule_name+'_solid2__'+str(int(step))
	
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
thermo_style 	custom c_SolEng
#thermo_style custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press vol c_SolEng tpcpu
#thermo_modify	line multi format float %14.6f

minimize 0.0 1.0e-8 1000 100000

fix relax all nve/limit 0.01
fix drag all viscous 0.05
timestep 2.0
run 2500
unfix relax
fix relax all nve
run 25000
unfix relax
unfix drag

velocity all create '''+str(T)+''' 1 rot yes dist gaussian
fix press all npt temp '''+str(T)+' '+str(T)+''' 100 aniso '''+str(P)+' '+str(P)+''' 1000
timestep 2.0
thermo 1
run 250000
''')
		f.close()
		if True: #start multiple processes
			#os.system('nohup ~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log &> /dev/null &' % (run_name, run_name) )
			makefile.write( '%s:\n\t~/lammps2/src/lmp_g++ -in %s.in -log %s.log &> /dev/null\n' % (run_name, run_name, run_name)  )
			print 'Running', run_name
		else:
			os.system('~/lammps/src/lmp_g++_no_cuda -in %s.in -log %s.log' % (run_name, run_name) )
			sys.exit()
	os.chdir('..')

makefile = open('lammps/Makefile_solid2', 'w')
makefile.write('all :')
for test_molecule_filename in ['acrylonitrile2.arc']:
	test_molecule_name = test_molecule_filename[ test_molecule_filename.rindex('/')+1 :-4]
	if test_molecule_name[-2:]=='e2': test_molecule_name = test_molecule_name[:-1]
	for step in range(0,101):
		run_name = test_molecule_name+'_solid__'+str(int(step))
		makefile.write(' %s' % run_name )
makefile.write('\n\n')

for filename in ['acrylonitrile2.arc']:
	run(filename, makefile)

