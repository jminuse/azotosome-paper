import os, string, threading, sys, shutil, math, pickle, random, re
import filetypes, lammps, utils, copy

utils.Molecule.set_params('oplsaa4.prm')

def run(filename):
	tail = utils.Molecule(filename+'.arc')
	
	atoms = []
	bonds = []
	angles = []
	dihedrals = []
	
	tail.add_to(0., 0., 0., atoms, bonds, angles, dihedrals)
	tail.add_to(0., 5., 0., atoms, bonds, angles, dihedrals)

	box_size = [40.]*3

	T = 100.0

	directory = 'lammps'
	if not os.path.isdir(directory):
		os.mkdir(directory)
	os.chdir(directory)

	run_name = 'bind3_'+filename

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

dump	1 all xyz 2500 '''+run_name+'''.xyz

thermo_modify	line multi format float %14.6f
minimize 0.0 1.0e-8 1000 100000

group mobile id > 1
velocity mobile create '''+str(T)+''' 1
run_style respa 4 2 2 2 inner 2 4.5 6.0 middle 3 7.0 8.0 outer 4
fix anneal mobile nvt temp '''+str(T)+''' 1.0 50.0
timestep 4.0
run 250000
minimize 0.0 1.0e-8 1000 100000

group aa id <= '''+str(len(atoms)/2)+'''
group bb subtract all aa
compute 1 aa group/group bb kspace yes
thermo_style custom c_1

run 0
''')
	f.close()
	os.system('lammps -in '+run_name+'.in -log '+run_name+'.log')
	os.chdir('..')
	
	
	for i,line in enumerate(open(directory+'/'+run_name+'.xyz').readlines()[ -len(atoms) : ]):
		columns = line.split()
		index, x, y, z = columns
		atoms[i].x, atoms[i].y, atoms[i].z = float(x), float(y), float(z)
	
	filetypes.write_xyz(run_name, atoms)

for n in ['cyanoallene2']: #['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane', 'hcn2', 'cyanoacetylene2', 'acrylonitrile2', 'cyanoallene2', 'acetonitrile2', 'hc5n2']:
	run(n)

