from iobasic import Base_Info


class Select_Extract(Base_Info):
	def __init__(self, _base_dir, _path_id, _replica_id, NPath, NReplicas, sysname, whichstate,  Single_path = False, Single_replica = False):
		Base_Info.__init__(self, _base_dir, _path_id, _replica_id, NPath, NReplicas, sysname, whichstate,  Single_path = False, Single_replica = False)
		atom_label, chrg = Base_Info.load_chrg(self)
		ene_pathway, barrier_energy, reverse_barrier = Base_Info.load_ene(self)

		self.u = Base_Info.load_unverise(self)
		self.NQMAtom = len(atom_label)
		self.chrg = chrg
		self.ene_pathway = ene_pathway           # the energy profile for the whole pathway
		self.barrier_energy = barrier_energy     # the barrier energy 
		self.reverse_barrier = reverse_barrier   # the reverse barrier energy
		self.qm_atom_format = self.get_qm_select_format() # get the qm atom format for selecting qm atoms by MDanalysis


	def get_qm_select_format(self):
		"""
		obtain the qm atoms format which can used for MDanalysis.
		The format is like (sigid Q{padthid} and ((resname *** and atom **) or () or ...)
		"""
		atom_label, _ = Base_Info.load_chrg(self)
		QM_atoms_selection = f"segid Q{self.replica_id} and ("
		for i in range(len(atom_label)):
			name_id_atom = atom_label[i].split('.')
			if i != (len(atom_label) - 1):             # if i is not the last one.
				id_and_atomtype = f"(resid {name_id_atom[1]} and name {name_id_atom[2]}) or "
			else:
				id_and_atomtype = f"(resid {name_id_atom[1]} and name {name_id_atom[2]}))"
#			print(id_and_atomtype)
			QM_atoms_selection = QM_atoms_selection + id_and_atomtype
#		print(QM_atoms_selection)
			
		return QM_atoms_selection
	
	def qm_atom_sel(self):
		"""
		selecting the qm atom using MDanalysis
		"""
		qmatomsel = (self.u).select_atoms(self.qm_atom_format)
		# ATOM NUMBER CHECK
		if (((qmatomsel.positions).shape)[0] != self.NQMAtom ):
			print("The QM atom number in selection is not consitant with the atom number we save in the provious step")
			exit()
		return qmatomsel

	def qmhvy_sel(self, return_str = False):
		"""
		selecting the heavy atoms of qm region
		for kpc_wt, res 72 is PHE, its qm region contains 50 heavy atoms
		for kpc_y72, res 72 is TYR irs qm region contains 51 heavy atoms
		"""
		if self.sysname == 'kpc_wt':
			qmhvysel =  f"segid Q{self.replica_id} and (" \
				          f"    (resid  69 and (name  CB or name  OG))" \
			              f" or (resid  71 and (name  CB or name  CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ))" \
                          f" or (resid  72 and (name  CB or name  CG or name  CD or name  CE or name  NZ))" \
                          f" or (resid 129 and (name  CB or name  OG))" \
                          f" or (resid 131 and (name  CB or name  CG or name OD1 or name ND2))" \
                          f" or (resid 165 and (name  CB or name  CG or name  CD or name OE1 or name OE2))" \
                          f" or (resid 169 and (name  CB or name  CG or name OD1 or name ND2 ))" \
                          f" or (resid 292 and (name OH2))" \
                          f" or (resid 291 and (name   C1 or name   C2 or name   C3 or name   N4 or name   C5 or "       \
                          f"                    name   C6 or name   C7 or name   O7 or name   S8 or name   C9 or "       \
                          f"                    name  C10 or name  N11 or name  C12 or name  N13 or name  C14 or "  \
                          f"                    name O14A or name O14B or name  C15 or name  O15 or name  C16 ))" \
                          f")"
			qmhvyatom = (self.u).select_atoms(qmhvysel)
			if qmhvyatom.n_atoms != 50 or qmhvyatom.n_residues != 9 :
				raise ValueError(f"qmhvy selected {qmhvyatom.n_residues}/{qmhvyatom.n_atoms} residues/atoms, should be 9/50.")

		if self.sysname == 'kpc_y72':
			qmhvysel =  f"segid Q{self.replica_id} and (" \
				          f"    (resid  69 and (name  CB or name  OG))" \
			              f" or (resid  71 and (name  CB or name  CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ or name OH))" \
                          f" or (resid  72 and (name  CB or name  CG or name  CD or name  CE or name  NZ))" \
                          f" or (resid 129 and (name  CB or name  OG))" \
                          f" or (resid 131 and (name  CB or name  CG or name OD1 or name ND2))" \
                          f" or (resid 165 and (name  CB or name  CG or name  CD or name OE1 or name OE2))" \
                          f" or (resid 169 and (name  CB or name  CG or name OD1 or name ND2 ))" \
                          f" or (resid 292 and (name OH2))" \
                          f" or (resid 291 and (name   C1 or name   C2 or name   C3 or name   N4 or name   C5 or "       \
                          f"                    name   C6 or name   C7 or name   O7 or name   S8 or name   C9 or "       \
                          f"                    name  C10 or name  N11 or name  C12 or name  N13 or name  C14 or "  \
                          f"                    name O14A or name O14B or name  C15 or name  O15 or name  C16 ))" \
                          f")"
			qmhvyatom = (self.u).select_atoms(qmhvysel)
			if qmhvyatom.n_atoms != 51 or qmhvyatom.n_residues != 9 :
				raise ValueError(f"qmhvy selected {qmhvyatom.n_residues}/{qmhvyatom.n_atoms} residues/atoms, should be 9/50.")

		if return_str == True:
			return qmhvysel # return the qmhvy_sel string, for the convience of selecting the hbonds, since HBA reuires string format
		else:
			return qmhvyatom

	def qmrxhdy_sel(self, return_str = False):
		"""
		Here we select the H atom involving the reaction process.
		resname LYS and resid 72 and name HZ1
		resname TIP3 and resid 292 and name H1 H2
		"""
		selcmd = f"segid Q{self.replica_id} and ("                          \
			f"    (resname  LYS and resid  72 and name HZ1)" \
            f" or (resname TIP3 and resid 292 and name  H1)" \
            f" or (resname TIP3 and resid 292 and name  H2)" \
            f")"

		if return_str == True:
			return selcmd
		
		qmrxhyd = (self.u).select_atoms(selcmd)
		if qmrxhyd.n_residues != 2 or qmrxhyd.n_atoms != 3:
			print([i for i in qmrxhyd])
			raise ValueError(f"qmhvy selected {qmrxhyd.n_residues}/{qmrxhyd.n_atoms} residues/atoms, should be 2/3.")
		return qmrxhyd
	
	def qmhdy_sel(self, return_str = False):
		"""
		sele the hydrogen atoms in the QM region in one replica
		"""
		selcmd = f"segid Q{self.replica_id} and name H* and (" \
                 f"    (resid  69 and resname SER)"  \
                 f" or (resid  71 and (resname PHE or resname TYR))"  \
                 f" or (resid  72 and resname LYS)"  \
                 f" or (resid 129 and resname SER)"  \
                 f" or (resid 131 and resname ASN)"  \
                 f" or (resid 165 and resname GLU)"  \
                 f" or (resid 169 and resname ASN)"  \
                 f" or (resid 291 and resname IMI)"  \
                 f" or (resid 292 and resname TIP3)" \
                 f")"
		if return_str == True:
			return selcmd
		
		qmhyd = (self.u).select_atoms(selcmd)
		return qmhyd

	def qmneighbor_sel(self, qmatomsel):
		"""
		select the neighboring atoms () for the specific atom
		"""
		qmatom = (self.u).select_atoms(f"segid Q{self.replica_id} and {qmatomsel}")
		qmhvy = self.qmhvy_sel()
		qmrxhyd = self.qmrxhdy_sel()
		qmneighbor = (self.u).select_atoms(f"(around 4 group qmatom) and (group qmhvy or group qmrxhyd)",
										   qmatom = qmatom, qmhvy = qmhvy, qmrxhyd = qmrxhyd)
		return qmneighbor
	   											   
		
if __name__ == "__main__":
	path_id = 1
	replica_id = 1
	sysname = 'kpc_wt'
	whichstate = 'd2'
	_basedir = "/users/chaoy/scratch/0.proj_kpc/1.sample"
	NPath = 200
	NReplicas = 36
	qmatomsel = "resid 72 and name NZ"
#	for i in range(1, NReplicas+1):
	X = Select_Extract(_basedir, path_id, replica_id, NPath, NReplicas, sysname, whichstate)
	#	ene_pathway, barrier_energy, reverse_barrier  = X.load_ene()
	#	U = X.load_unverise()
	c = X.qmneighbor_sel(qmatomsel)
	b = X.qmhdy_sel(True)
	a = X.qmhvy_sel(True)
	d = X.heavy_atoms_bonds()
	print(d,len(d))
	e = X.extract_dist_one_replica()
	print(e,len(e))
#	print(len(d))

#	print(a.n_atoms, b.n_atoms, c.n_atoms)



	
