import re, os, math, copy, subprocess
import gaussian, filetypes, utils

def minimize():
	#for f in ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane']:
	for f in ['hcn34', 'cyanoacetylene34', 'acrylonitrile34', 'cyanoallene34', 'acetonitrile34', 'hc5n34']:
		atoms = utils.Molecule(f+'.arc').atoms
		gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', f+'_0', 'Opt', procs=1)
		gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', f+'_1', 'Opt SCRF(Solvent=n-Hexane)', procs=1)

def compare():
	utils.Molecule.set_params('oplsaa4.prm')
	for f in ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane', 'hcn2', 'cyanoacetylene2', 'acrylonitrile2', 'cyanoallene2', 'acetonitrile2', 'hc5n2']:
		atoms0 = utils.Molecule(f+'.arc').atoms
		atoms1 = gaussian.parse_atoms('gaussian/'+f+'_0.log')[1]
		for i in range(len(atoms0)):
			atoms0[i].index = i+1
			atoms1[i].index = i+1
			
			#atoms0[i].x = atoms1[i].xp
			#atoms0[i].y = atoms1[i].y
			#atoms0[i].z = atoms1[i].z
	
		bonds0 = filetypes.get_bonds(atoms0)
		bonds1 = filetypes.get_bonds(atoms1)
		
		#filetypes.write_arc(f[:-1]+'3', atoms0)
		#continue
		
		angles0, dihedrals0 = filetypes.get_angles_and_dihedrals(atoms0, bonds0)
		angles1, dihedrals1 = filetypes.get_angles_and_dihedrals(atoms1, bonds1)
		
		bond_error = 0.0
		for i in range(len(bonds0)):
			bond_error += abs(bonds0[i].d - bonds1[i].d)
		
		angles0.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		angles1.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		
		angle_error = 0.0
		for i in range(len(angles0)):
			angle_error += min( abs(angles0[i].theta - angles1[i].theta), abs(angles0[i].theta - angles1[i].theta+math.pi*2), abs(angles0[i].theta - angles1[i].theta-math.pi*2) )
			
		dihedrals0.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		dihedrals1.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		
		angles_by_atoms = {}
		for angle in angles1:
			angles_by_atoms[ tuple(angle.atoms) ] = angle
			angles_by_atoms[ tuple(reversed(angle.atoms)) ] = angle
		
		if dihedrals0:
			dihedral_error = 0.0
			for i in range(len(dihedrals0)):
				#if dihedral only controls a 180-degree angle, ignore comparison
				
				angle1 = angles_by_atoms[ tuple(dihedrals1[i].atoms[1:]) ]
				angle2 = angles_by_atoms[ tuple(dihedrals1[i].atoms[:-1]) ]
				
				#print angle1.theta, angle2.theta,
				
				if abs(angle1.theta-180)<2 or abs(angle1.theta+180)<2 or abs(angle2.theta-180)<2 or abs(angle2.theta+180)<2:
					#print angle1.theta, angle2.theta, 
					continue
				
				#if dihedral is single-bonded free rotor, ignore comparison
				if (dihedrals0[i].atoms[1].type.index,dihedrals0[i].atoms[1].type.index) == (78,78):
					continue
				#else make comparison
				try:
					a,b = utils.dihedral_angle(*dihedrals0[i].atoms)[0], utils.dihedral_angle(*dihedrals1[i].atoms)[0]
					dihedral_error += min( abs(a-b), abs(a-b+math.pi*2), abs(a-b-math.pi*2) )
					
					#print [aa.element for aa in dihedrals0[i].atoms], min( abs(a-b), abs(a-b+math.pi*2), abs(a-b-math.pi*2) )
				except ZeroDivisionError: pass
		
			print f, bond_error/len(bonds0), angle_error/len(angles0), 180/math.pi*dihedral_error/len(dihedrals0)
		else:
			print f, bond_error/len(bonds0), angle_error/len(angles0), 0.0

def compare_distances_only():
	utils.Molecule.set_params('oplsaa4.prm')
	for f in ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane', 'hcn2', 'cyanoacetylene2', 'acrylonitrile2', 'cyanoallene2', 'acetonitrile2', 'hc5n2']:
		atoms0 = utils.Molecule(f+'.arc').atoms
		atoms1 = gaussian.parse_atoms('gaussian/'+f+'_1.log')[1]
		for i in range(len(atoms0)):
			atoms0[i].index = i+1
			atoms1[i].index = i+1
			
			#atoms0[i].x = atoms1[i].xp
			#atoms0[i].y = atoms1[i].y
			#atoms0[i].z = atoms1[i].z
	
		bonds0 = filetypes.get_bonds(atoms0)
		bonds1 = filetypes.get_bonds(atoms1)
		
		#filetypes.write_arc(f[:-1]+'3', atoms0)
		#continue
		
		angles0, dihedrals0 = filetypes.get_angles_and_dihedrals(atoms0, bonds0)
		angles1, dihedrals1 = filetypes.get_angles_and_dihedrals(atoms1, bonds1)
		
		bond_error = 0.0
		for i in range(len(bonds0)):
			bond_error += abs(bonds0[i].d - bonds1[i].d)
		
		angles0.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		angles1.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		
		angle_error = 0.0
		for i in range(len(angles0)):
			d0 = utils.dist(angles0[i].atoms[0], angles0[i].atoms[2])
			d1 = utils.dist(angles1[i].atoms[0], angles1[i].atoms[2])
			angle_error += abs(d0-d1)
			
		dihedrals0.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		dihedrals1.sort(key=lambda a: tuple([n.index for n in a.atoms]) )
		
		if dihedrals0:
			dihedral_error = 0.0
			for i in range(len(dihedrals0)):
				#if dihedral is single-bonded free rotor, ignore comparison
				if (dihedrals0[i].atoms[1].type.index,dihedrals0[i].atoms[1].type.index) == (78,78):
					continue
				d0 = utils.dist(dihedrals0[i].atoms[0], dihedrals0[i].atoms[3])
				d1 = utils.dist(dihedrals1[i].atoms[0], dihedrals1[i].atoms[3])
				dihedral_error += abs(d0-d1)
		
			print f, bond_error/len(bonds0), angle_error/len(angles0), 180/math.pi*dihedral_error/len(dihedrals0)
		else:
			print f, bond_error/len(bonds0), angle_error/len(angles0), 0.0

def binding():
	#for n in ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane', 'hcn2', 'cyanoacetylene2', 'acrylonitrile2', 'cyanoallene2', 'acetonitrile2', 'hc5n2']:
	for n in ['aminobutane2']:
		
		contents = open('lammps/bind3_'+n+'.log').read()
		
		opls_binding_e = float(re.search('\n1\s+(\S+)\s+Loop time of', contents).group(1))
		
		atoms = filetypes.parse_xyz('bind3_'+n+'.xyz')
		
		#opt0: use B97D
		#gaussian.job(atoms, 'B97D/TZVP/Fit', 'batch', 'bind3_'+n+'_opt0', 'Opt=Cartesian', procs=1)
		
		#opt1: use M062X
		#os.system('cp gaussian/bind3_'+n+'_opt0.chk gaussian/bind3_'+n+'_opt1.chk')
		#gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1', 'Opt=Cartesian Geom=AllCheck Guess=Read', procs=1)
		#opt2: M062X with solvent
		#os.system('cp gaussian/bind3_'+n+'_opt0.chk gaussian/bind3_'+n+'_opt2.chk')
		#gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2', 'Opt=Cartesian Geom=AllCheck Guess=Read SCRF(Solvent=n-Hexane)', procs=1)
		
		'''
		jlist = subprocess.Popen('jlist', shell=True, stdout=subprocess.PIPE).communicate()[0]
		if ' bind3_'+n+'_opt1 ' not in jlist and ' bind3_'+n+'_opt1.1 ' not in jlist and gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.log') is None and (not os.path.isfile('gaussian/bind3_'+n+'_opt1.1.log') or gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.1.log') is None):
			os.system('cp gaussian/bind3_'+n+'_opt1.chk gaussian/bind3_'+n+'_opt1.1.chk')
			#gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1.1', 'Opt=Restart', procs=1)
			#gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1.1', 'Opt=Cartesian Geom=AllCheck Guess=Read', procs=1)
			print 'bind3_'+n+'_opt1.1'
		if ' bind3_'+n+'_opt2 ' not in jlist and ' bind3_'+n+'_opt2.1 ' not in jlist and gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2.log') is None and (not os.path.isfile('gaussian/bind3_'+n+'_opt2.1.log') or gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2.1.log') is None):
			os.system('cp gaussian/bind3_'+n+'_opt2.chk gaussian/bind3_'+n+'_opt2.1.chk')
			#gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2.1', 'Opt=Restart SCRF(Solvent=n-Hexane)', procs=1)
			#gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2.1', 'Opt=Cartesian Geom=AllCheck Guess=Read SCRF(Solvent=n-Hexane)', procs=1)
			print 'bind3_'+n+'_opt2.1'
		continue
		'''
		
		#Find BSSE without solvent (on opt1)
		#os.system('cp gaussian/'+n+'_0.chk gaussian/bind3_'+n+'_opt1a.chk')
		#os.system('cp gaussian/'+n+'_0.chk gaussian/bind3_'+n+'_opt1b.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt1.chk gaussian/bind3_'+n+'_opt1a_Bq.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt1.chk gaussian/bind3_'+n+'_opt1b_Bq.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt1.1.chk gaussian/bind3_'+n+'_opt1a_Bq.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt1.1.chk gaussian/bind3_'+n+'_opt1b_Bq.chk')
		try:
			atoms = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.1.log', check_convergence=False)[1]
		except:
			atoms = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.log', check_convergence=False)[1]
		#gaussian.job(atoms[ : len(atoms)/2 ], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1a', 'Guess=Read', procs=1)
		#gaussian.job(atoms[ len(atoms)/2 : ], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1b', 'Guess=Read', procs=1)
		atoms_a = copy.deepcopy(atoms)
		for a in atoms_a[ : len(atoms)/2 ]:
			a.element = a.element + '-Bq'
		#gaussian.job(atoms_a, 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1a_Bq', 'Guess=Read', procs=1)
		for a in atoms[ len(atoms)/2 : ]:
			a.element = a.element + '-Bq'
		#gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt1b_Bq', 'Guess=Read', procs=1)
		
		
		#Find BSSE with solvent (on opt2)
		#os.system('cp gaussian/'+n+'_1.chk gaussian/bind3_'+n+'_opt2a.chk')
		#os.system('cp gaussian/'+n+'_1.chk gaussian/bind3_'+n+'_opt2b.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt2.chk gaussian/bind3_'+n+'_opt2a_Bq.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt2.chk gaussian/bind3_'+n+'_opt2b_Bq.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt2.1.chk gaussian/bind3_'+n+'_opt2a_Bq.chk')
		#os.system('cp gaussian/bind3_'+n+'_opt2.1.chk gaussian/bind3_'+n+'_opt2b_Bq.chk')
		try:
			atoms = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2.1.log', check_convergence=False)[1]
		except:
			atoms = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2.log', check_convergence=False)[1]
		#gaussian.job(atoms[ : len(atoms)/2 ], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2a', 'Guess=Read SCRF(Solvent=n-Hexane)', procs=1)
		#gaussian.job(atoms[ len(atoms)/2 : ], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2b', 'Guess=Read SCRF(Solvent=n-Hexane)', procs=1)
		atoms_a = copy.deepcopy(atoms)
		for a in atoms_a[ : len(atoms)/2 ]:
			a.element = a.element + '-Bq'
		gaussian.job(atoms_a, 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2a_Bq', 'Guess=Read SCRF(Solvent=n-Hexane)', procs=1)
		for a in atoms[ len(atoms)/2 : ]:
			a.element = a.element + '-Bq'
		#gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'_opt2b_Bq', 'Guess=Read SCRF(Solvent=n-Hexane)', procs=1)
		
		continue
		
		#atoms = atoms[ : len(atoms)/2 ]
		#atoms = atoms[ len(atoms)/2 : ]
		
		#for a in atoms[ len(atoms)/2 : ]:
		#	a.element = a.element + '-Bq'
		
		#gaussian.job(atoms, 'M062X/cc-pVDZ', 'batch', n+'_sp4', 'SP', procs=1) #sp0 is pair, sp1,sp2 are single, sp3,sp4 with ghost atoms
		
		os.system('cp gaussian/'+n+'_sp0.chk gaussian/'+n+'_scrf0.chk')
		gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', n+'_scrf0', 'SP SCRF(Solvent=n-Hexane) Guess=Read', procs=1) #scrf0, scrf1, scrf2
		atoms1 = copy.deepcopy(atoms)
		for a in atoms1[ : len(atoms)/2 ]:
			a.element = a.element + '-Bq'
		os.system('cp gaussian/'+n+'_sp0.chk gaussian/'+n+'_scrf1.chk')
		gaussian.job(atoms1, 'M062X/aug-cc-pVDZ', 'batch', n+'_scrf1', 'SP SCRF(Solvent=n-Hexane) Guess=Read', procs=1)
		atoms2 = copy.deepcopy(atoms)
		for a in atoms2[ len(atoms)/2 : ]:
			a.element = a.element + '-Bq'
		os.system('cp gaussian/'+n+'_sp0.chk gaussian/'+n+'_scrf2.chk')
		gaussian.job(atoms2, 'M062X/aug-cc-pVDZ', 'batch', n+'_scrf2', 'SP SCRF(Solvent=n-Hexane) Guess=Read', procs=1)
		
		#pair = gaussian.parse_atoms('gaussian/'+n+'_sp0.log')[0]
		#single = gaussian.parse_atoms('gaussian/'+n+'_sp1.log')[0]
		
		#print n, opls_binding_e, 627.5*(pair - single*2)
		
		'''
		atoms = filetypes.parse_xyz(n+'_bound.xyz')
		#gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', n+'_opt', 'Opt', procs=1)
		
		for i,a in enumerate(atoms):
			if i<len(atoms)/2: a.fragment = 1
			else: a.fragment = 2
		
		gaussian.job(atoms, 'M062X/aug-cc-pVDZ', 'batch', n+'_cp', 'Counterpoise=2', procs=1, charge_and_multiplicity='0,1 0,1 0,1')
		'''
		
def binding2():
	for n in ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane', 'hcn2', 'cyanoacetylene2', 'acrylonitrile2', 'cyanoallene2', 'acetonitrile2', 'hc5n2']:
		'''
		try:
			contents = open('lammps/bind3_'+n+'.log').read()
			opls_binding_e = float(re.search('\n1\s+(\S+)\s+Loop time of', contents).group(1))
			
			pair = gaussian.parse_atoms('gaussian/'+n+'_opt.log')[0]
			lone1 = gaussian.parse_atoms('gaussian/'+n+'_0.log')[0]
			lone2 = gaussian.parse_atoms('gaussian/'+n+'_0.log')[0]
		
			print n, opls_binding_e, 627.5*(pair - lone1 - lone2)
		except: print n
		'''
		#(Optimized pair energy) - 2*(optimized single energy) - (BSSE of optimized pair)
		try:
			contents = open('lammps/bind3_'+n+'.log').read()
			opls_pair = float(re.findall('TotEng\s+=\s+(\S+)', contents)[-1])
			opls_single = float(re.findall('TotEng\s+=\s+(\S+)', contents)[-1])
			
			pair = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.log')[0]
			single = gaussian.parse_atoms('gaussian/'+n+'_0.log')[0]
		
			print n, opls_binding_e, 627.5*(pair - lone1 - lone2)
		except: print n
		
def restart():
	for n in ['aminopentane2_opt2']:
		os.system('cp gaussian/bind3_'+n+'.chk gaussian/bind3_'+n+'.1.chk')
		if n[-1]=='1':
			gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'.1', 'Opt=Restart', procs=1)
		elif n[-1]=='2':
			gaussian.job([], 'M062X/aug-cc-pVDZ', 'batch', 'bind3_'+n+'.1', 'Opt=Restart SCRF(Solvent=n-Hexane)', procs=1)

#minimize()
#compare()

#binding()

#restart()

def analyze_binding():
	for n in ['propanenitrile', 'butanenitrile', 'pentanenitrile', 'hexanenitrile', 'aminopropane2', 'aminobutane2', 'aminopentane2', 'aminohexane2', 'hexane2', 'methane', 'hcn2', 'cyanoacetylene2', 'acrylonitrile2', 'cyanoallene2', 'acetonitrile2', 'hc5n2']:
	
		contents = open('lammps/bind3_'+n+'.log').read()
		
		opls_binding_e = float(re.search('\n1\s+(\S+)\s+Loop time of', contents).group(1))
		
		try:
			e_pair = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.1.log', check_convergence=False)[0]
		except:
			e_pair = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1.log', check_convergence=False)[0]
			
		try:
			e_pair_scrf = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2.1.log', check_convergence=False)[0]
		except:
			e_pair_scrf = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2.log', check_convergence=False)[0]

		e_single = gaussian.parse_atoms('gaussian/'+n+'_0.log')[0]
		e_single_scrf = gaussian.parse_atoms('gaussian/'+n+'_1.log')[0]

		try:
			bsse = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1a.log')[0] + gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1b.log')[0] - gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1a_Bq.log')[0] - gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1b_Bq.log')[0]
		except:
			print 1, n
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1a.log') is None
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1b.log') is None
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1a_Bq.log') is None
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt1b_Bq.log') is None
			bsse = 0.0
		
		try:
			bsse_scrf = gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2a.log')[0] + gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2b.log')[0] - gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2a_Bq.log')[0] - gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2b_Bq.log')[0]
		except:
			print 2, n
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2a.log') is None
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2b.log') is None
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2a_Bq.log') is None
			#print gaussian.parse_atoms('gaussian/bind3_'+n+'_opt2b_Bq.log') is None
			bsse_scrf = 0.0
			
		#No such file or directory: 'gaussian/bind3_aminopentane2_opt2b_Bq.log
		#No such file or directory: 'gaussian/bind3_aminobutane2_opt2a_Bq.log'
		
		
		#print n, opls_binding_e, e_pair, e_single, bsse, e_pair_scrf, e_single_scrf, bsse_scrf
		
		print n, ',', opls_binding_e, ',', (e_pair-2*e_single+bsse)*627, ',',  (e_pair_scrf-2*e_single_scrf+bsse_scrf)*627

#binding()		

#analyze_binding()

compare()

#compare_distances_only()

