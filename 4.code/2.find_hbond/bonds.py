'''O código abaixo analisa um sistema molecular (por exemplo, uma proteína ou complexo molecular) e identifica:

Ligações entre átomos pesados (não hidrogênio).

Ligações de hidrogênio (H-bonds) no sistema.

Ligações de hidrogênio envolvendo átomos de hidrogênio reativos (especificados pelo usuário).

Ele também calcula distâncias entre átomos e gera rótulos (labels) para identificar essas ligações.
'''
from select_and_extract import Select_Extract
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

class Bonds(Select_Extract):

	def _is_hbond(self, donor, hydrogen, acceptor):
		"""Check if Hydrogen bonded using the default mdanalysis.HBA metric."""
		dh_dist   = self.get_dist(donor, hydrogen)
		da_dist   = self.get_dist(donor, acceptor)
		dha_angle = self.get_angle(donor, hydrogen, acceptor)
		if dh_dist <= 1.2 and da_dist <= 3.0 and dha_angle >= 110:
			if dha_angle >= 150:
				print(f"####Found dha_angle >150 match: {donor} - {acceptor}")
				return True
		else:
			return False

	def hvychembonds_label(self):
		"""
		Find all heavy atoms pairin bond should be the same for all systems and all replica_id
		"""
		heavy_atoms = self.qmhvy_sel()
		all_bonds = heavy_atoms.bonds
		add_to_the_list_value = 0  # for y72 structure, we add the TYR.71.CZ:TYR.71.OH to the last part of the label list

		bonds_label = []
		for i in range(0, heavy_atoms.n_atoms-1):
			for j in range(i+1, heavy_atoms.n_atoms):

				atom_i = heavy_atoms[i]
				atom_j = heavy_atoms[j]

				for bond in all_bonds:
					if atom_i in bond and atom_j in bond:
						# l = f"{atom_i.residue.resname}.{atom_i.residue.resid}.{atom_i.name}:" \
						# 	f"{atom_j.residue.resname}.{atom_j.residue.resid}.{atom_j.name}"
						l = f"{atom_i.residue.resid}.{atom_i.name}:" \
							f"{atom_j.residue.resid}.{atom_j.name}"

						# if l in ["TYR.71.CZ:TYR.71.OH"]:
						if l in ["71.CZ:71.OH"]:
							add_to_the_list_value = add_to_the_list_value + 1
							continue
						# remove the pseudo-chembond
					    # if l not in ["IMI.291.N4:IMI.291.C7", "IMI.291.N4:IMI.291.C7"]:
						if l not in ["291.N4:291.C7", "291.N4:291.C7"]:
							bonds_label.append(l)

		# Manually add some heavy-atom-bonds not in the topology
		# bonds_label.append("SER.69.OG:IMI.291.C7")
		# bonds_label.append("IMI.291.C7:TIP3.292.OH2")
		bonds_label.append("69.OG:291.C7")
		bonds_label.append("291.C7:292.OH2")

		if add_to_the_list_value > 0:
			# bonds_label.append("TYR.71.CZ:TYR.71.OH")
			bonds_label.append("71.CZ:71.OH")

		return bonds_label

	def extract_dist_one_replica(self):
		"""
		based on the bonds label we get, extract the bond distance from the bonds label
		"""
		bonds_label = self.hvychembonds_label()
		dist_row = []

		for bond in bonds_label:
			atomlabel_i = bond.split(':')[0]
			atomlabel_j = bond.split(':')[1]

			atom_i = (self.u).select_atoms(f"segid Q{self.replica_id} and resname {atomlabel_i.split('.')[0]}" \
										   f"                         and resid   {atomlabel_i.split('.')[1]}" \
										   f"                         and name    {atomlabel_i.split('.')[2]}" )
			atom_j = (self.u).select_atoms(f"segid Q{self.replica_id} and resname {atomlabel_j.split('.')[0]}" \
										   f"                         and resid   {atomlabel_j.split('.')[1]}" \
										   f"                         and name    {atomlabel_j.split('.')[2]}" )

			if atom_i.n_atoms != 1 or atom_j.n_atoms != 1:raise ValueError(f"Bad selection on reactive bonds")




			d = self.get_dist(atom_i[0], atom_j[0])
			dist_row.append(d)
		return dist_row

	def hbond_labels(self):
		"""
		extracting all Hbonds in the qm regions for each replica
		hbonds = HBA(
		  universe,
		  donors_sel,
		  hydrogen_sel,
		  acceptor_sel,
		  d_h_cutoff,  # D-H
		  d_a_cutoof,  # D-A
		  d_h_a_angle_cutoff, #D-H-A
		  update_selections
		 )
		Output:
		results = [
		  [
		     <frame>,
		     <donor index (0-based)>,
		     <hydrogen index (0-based)>,
		     <acceptor index (0-based)>,
		     <distance>
		     <angle>
		  ],
		 ...
		]

		"""
		hvy_atoms = self.qmhvy_sel(return_str = True)
		hdy_atoms = self.qmhdy_sel(return_str = True)
		# The Hbonds only exists in the QM region.
		hba = HBA(
			universe = self.u,
			donors_sel = hvy_atoms,
			hydrogens_sel = hdy_atoms,
			acceptors_sel = hvy_atoms,
                        d_a_cutoff = 2.8,
                        d_h_a_angle_cutoff = 150,
			)
		hba.run()
		hbonds_tmp = hba.hbonds

		hlabels = []

		for hb in hbonds_tmp:
			donor_hvyatom = (self.u).atoms[int(hb[1])]
			hydrogen_atom = (self.u).atoms[int(hb[2])]
			accep_hvyatom = (self.u).atoms[int(hb[3])]

			# exclude intra-residue hydrogen bond
			if (donor_hvyatom.residue.resname == accep_hvyatom.residue.resname and \
				donor_hvyatom.residue.resid == accep_hvyatom.residue.resid):
				continue
			dha_l = f"{donor_hvyatom.residue.resname}.{donor_hvyatom.residue.resid}.{donor_hvyatom.name}:" \
				   f"{hydrogen_atom.residue.resname}.{hydrogen_atom.residue.resid}.{hydrogen_atom.name}:" \
			       f"{accep_hvyatom.residue.resname}.{accep_hvyatom.residue.resid}.{accep_hvyatom.name}"
			# dimiss the Hbonds between TYR.71.OH:TYR.71.HH:GLU:165.OE2, add it manualy in the last step(merge all )
			if dha_l == "TYR.71.OH:TYR.71.HH:GLU.165.OE2":
				continue
			hlabels.append(dha_l)

		# fix possiable Hbonds concerning IMI-N4-H --acceptor
		imi_n4 = (self.u).select_atoms(f"segid Q{self.replica_id} and resname IMI and resid 291 and name  N4")
		s70_hg = (self.u).select_atoms(f"segid Q{self.replica_id} and resname SER and resid  69 and name HG1")

		if imi_n4.n_atoms != 1 or s70_hg.n_atoms != 1 : raise ValueError(f"Bad selection on IMI N-H atoms")

		donor = imi_n4[0]
		hydro = s70_hg[0]
		candid_acceptors = (self.u).select_atoms(
			f"(around 3.0 (segid Q{self.replica_id} and resname IMI and resid 291 and name N4)) and (group qmhvyatom )",
			qmhvyatom = self.qmhvy_sel()
		)
#		print(candid_acceptors)
		if candid_acceptors.n_atoms == 0 : raise ValueError(f"Bad selection on IMI N-H atoms")

		for accep in candid_acceptors:
			if self._is_hbond(donor, hydro, accep) == True:
				da_l = f"{donor.residue.resname}.{donor.residue.resid}.{donor.name}:" \
					   f"{hydro.residue.resname}.{hydro.residue.resid}.{hydro.name}:" \
     				   f"{accep.residue.resname}.{accep.residue.resid}.{accep.name}"
				ad_l = f"{accep.residue.resname}.{accep.residue.resid}.{accep.name}:" \
					   f"{hydro.residue.resname}.{hydro.residue.resid}.{hydro.name}:" \
					   f"{donor.residue.resname}.{donor.residue.resid}.{donor.name}"
				if not (da_l in hlabels or ad_l in hlabels):
					hlabels.append(da_l)
		# print(len(hlabels))
		return hlabels


	def rxhbond_labels(self):
		"""
		Find all possibale hydrogen bond related with rxhydrogen
		"""
		hvy_atoms = self.qmhvy_sel(return_str = True)
		hdy_atoms = self.qmrxhdy_sel(return_str = True)
#		print(hvy_atoms)
		# The Hbonds only exists in the QM region.
		hba = HBA(
			universe = self.u,
			donors_sel = hvy_atoms,
			hydrogens_sel = hdy_atoms,
			acceptors_sel = hvy_atoms,
			)
		hba.run()
		rxhbonds_tmp = hba.hbonds

		rxhlabels = []
		for hb in rxhbonds_tmp:
			donor_hvyatom = (self.u).atoms[int(hb[1])]
			hydrogen_atom = (self.u).atoms[int(hb[2])]
			accep_hvyatom = (self.u).atoms[int(hb[3])]

			dh_l = f"{donor_hvyatom.residue.resname}.{donor_hvyatom.residue.resid}.{donor_hvyatom.name}:" \
				   f"{hydrogen_atom.residue.resname}.{hydrogen_atom.residue.resid}.{hydrogen_atom.name}"
			rxhlabels.append(dh_l)

			ah_l = f"{accep_hvyatom.residue.resname}.{accep_hvyatom.residue.resid}.{accep_hvyatom.name}:" \
			f"{hydrogen_atom.residue.resname}.{hydrogen_atom.residue.resid}.{hydrogen_atom.name}"
			rxhlabels.append(ah_l)
		return rxhlabels



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
	X = Bonds(_basedir, path_id, replica_id, NPath, NReplicas, sysname, whichstate)
	#	ene_pathway, barrier_energy, reverse_barrier  = X.load_ene()
	#	U = X.load_unverise()
	c = X.qmneighbor_sel(qmatomsel)
	b = X.qmhdy_sel(True)
	a = X.qmhvy_sel(True)
	d = X.hvychembonds_label()
	print(d,len(d))
#	e = X.extract_dist_one_replica()
#	print(e,len(e))
#	f = X.rxhbond_labels()#
#	print(f)

#	print(a.n_atoms, b.n_atoms, c.n_atoms)
