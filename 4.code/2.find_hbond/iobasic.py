# Load the basic info of system
# Chao Yin, Jan 10, 2022

import os
import warnings
warnings.filterwarnings("ignore")

# Third party package
import numpy as np
from MDAnalysis import Universe

class Base_Info():
	def __init__(self, _base_dir, _path_id, _replica_id, NPath, NReplicas, sysname, whichstate,  Single_path = False, Single_replica = False):
		self.base_dir = _base_dir           # basic directory
		self.path_id = _path_id
		self.replica_id = _replica_id
		self.NPath = NPath
		self.NReplicas = NReplicas
		self.sysname = sysname             # protein name, kpc_wt or kpc_y72
		self.whichstate = whichstate       # structure state, d1 or d2
		self.dir_psf, self.dir_cor, self.dir_charge, self.dir_ene  = Base_Info.get_path(self)
		self.Single_path = Single_path     # Analysis one Path containing NReplicas if Single Path is true, if false, analyze 1 to NReplicas
		self.Single_replica = Single_replica # Analyze one replicas if Single replicas is true, if false, analyze NReplicas

	def get_path(self):
		if self.whichstate == 'd2' and self.sysname == 'kpc_wt':
			_proj_path = f"{self.base_dir}/0.kpc_wt_d2"
		elif self.whichstate == 'd1' and self.sysname == 'kpc_wt':
			_proj_path = f"{self.base_dir}/1.kpc_wt_d1"
		elif self.whichstate == 'd2' and self.sysname == 'kpc_y72':
			_proj_path = f"{self.base_dir}/2.kpc_y72_d2"
		elif self.whichstate == 'd1' and self.sysname == 'kpc_y72':
			_proj_path = f"{self.base_dir}/3.kpc_y72_d1"

		struc_path = f"{_proj_path}/6.rpath/path_opt"
		data_path  = f"{_proj_path}/7.sp/0.b3lyp_d3/charge_ene"
		_psf_file  = f"{struc_path}/{self.sysname}.path_f{self.path_id}.psf"
		_cor_file  = f"{struc_path}/{self.sysname}.path_f{self.path_id}.cor"
		_chrg_file = f"{data_path}/charges_all_paths.npz"
		_ene_file  = f"{data_path}/enes_all_paths.npz"

		return _psf_file, _cor_file, _chrg_file, _ene_file

	def load_unverise(self):
		"""load the system, which includes NReplicas coordinates, and return as a mda.universe object"""
		u = Universe(self.dir_psf, self.dir_cor, topology_format = 'PSF', format = 'CRD')
		# print(u)
		return u

	def load_ene(self):
		"""
        load the pathway energy, if self.Single_path is true, return the energy for one pathway.
		If flase, return the energy profile for NPath pathways
	    !!!! Attention: enegy profile is a npz file containing 3 arraies, which are 'enes_pathway', 'ene_barrier' and 'reverse_enes',
		return 'enes_pathway', 'ene_barrier' and 'reverse_enes'
		"""
		energy_profile = np.load(self.dir_ene)
		if self.Single_path == True:
			return (energy_profile['enes_pathway'][(self.NPath - 1)],
		            energy_profile['ene_barrier'][(self.NPath - 1)] ,
				    energy_profile['reverse_enes'][(self.NPath - 1)] ) # return  pathway energy profile
		else :
			return (energy_profile['enes_pathway'] ,
		            energy_profile['ene_barrier'] ,
				    energy_profile['reverse_enes'] )

	def load_chrg(self):
		"""
		load the charge profile for the pathway
		self.Single_path = false, self.Single_replica = false, dimension = (NPath * NReplicas * NQMAtoms)
		self.Single_path = true, self.Single_replica = false, dimension for the charge profile = (NReplicas * NQMatoms)
		self.Single_path = true, self.Single_replace = ture, dimension for the charge profile = (NQMAtoms)
		return two np.array, atom_label, charge_all_path
		!!! ATTENTION, here we read charge npz file, which contains two arraies, 'atom_label' and 'charge_all_path'
		"""
		chrg_profile = np.load(self.dir_charge)
		if (self.Single_path == False and self.Single_replica == False):
			return (chrg_profile['atom_label'] ,
	                chrg_profile['charge_all_path'])
		elif (self.Single_path == True and self.Single_replica == False):
			return (chrg_profile['atom_label'],
		            chrg_profile['charge_all_path'][NPath - 1] )
		elif (self.Single_path == True and self.Single_replica == False):
			return (chrg_profile['atom_label'],
		            chrg_profile['charge_all_path'][NPath - 1][NReplicas - 1] )

	def get_dist(self, atom1, atom2):
		"""
		get the distance bewteen two atoms, note here, atom1 and atom2 must be atom, not AtomGroup
		"""
		position_atom1 = atom1.position
		position_atom2 = atom2.position

		dist = np.linalg.norm(position_atom1 - position_atom2)
		return dist

	def get_angle(self, atom1, atom2, atom3):
		""" compute the D-H-A angle"""
		d = atom1.position
		h = atom2.position
		a = atom3.position

		hd = d - h
		ha = a - h

		cosine_angle = np.dot(hd, ha) / (np.linalg.norm(hd) * np.linalg.norm(ha))
		angle = np.arccos(cosine_angle)
		return np.abs(np.degrees(angle))



if __name__ == "__main__":
	path_id = 1
	replica_id = 1
	sysname = 'kpc_wt'
	whichstate = 'd2'
	_basedir = "/users/chaoy/scratch/0.proj_kpc/1.sample"
	NPath = 200
	NReplicas = 36
	X = Base_Info(_basedir, path_id, replica_id, NPath, NReplicas, sysname, whichstate)
	ene_pathway, barrier_energy, reverse_barrier  = X.load_ene()
	print(ene_pathway.shape)
	_, chrg = X.load_chrg()
	print(chrg.shape)
