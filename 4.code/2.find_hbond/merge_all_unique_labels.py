from multiprocessing import Pool
from extract_on_batch import batch_extract

import numpy as np

class merge_unique_labels():
	def __init__(self, _basedir, NPath, NReplicas):
		self.base_dir = _basedir
		self.NPath = NPath
		self.NReplicas = NReplicas
		self.wt_d1 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_wt", "d1")
		self.wt_d2 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_wt", "d2")
		self.y72_d1 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_y72", "d1")
		self.y72_d2 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_y72", "d2")


	def merge_hbonds(self):
		"""
		Merge all hydrogen bonds in all structure(all paths and all replicas),  here we should add the Hydrogen bonds TYR.71.OH:TYR.71.HH:GLU.165.OE2
		"""
		unique_hbonds = []
		hbonds_mat = []

		with Pool(processes = 4) as pool:
			# wt_d1 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_wt", "d1")
			# wt_d2 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_wt", "d2")
			# y72_d1 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_y72", "d1")
			# y72_d2 = batch_extract(self.base_dir,self.NPath, self.NReplicas, "kpc_y72", "d2")

			wt_d1_hbonds = pool.apply_async((self.wt_d1).extract_hbonds_on_one_structure, args =())
			wt_d2_hbonds = pool.apply_async((self.wt_d2).extract_hbonds_on_one_structure, args =())
			y72_d1_hbonds = pool.apply_async((self.y72_d1).extract_hbonds_on_one_structure, args =())
			y72_d2_hbonds = pool.apply_async((self.y72_d2).extract_hbonds_on_one_structure, args =())

			wt_d1_hbonds_labels = wt_d1_hbonds.get()
			wt_d2_hbonds_labels = wt_d2_hbonds.get()
			y72_d1_hbonds_labels = y72_d1_hbonds.get()
			y72_d2_hbonds_labels = y72_d2_hbonds.get()

		# merge the unique hydrogen bonds labels
		hbonds_mat.append(wt_d1_hbonds_labels)
		hbonds_mat.append(wt_d2_hbonds_labels)
		hbonds_mat.append(y72_d1_hbonds_labels)
		hbonds_mat.append(y72_d2_hbonds_labels)
		for i in range(len(hbonds_mat)):
			for j in range(len(hbonds_mat[i])):
				if hbonds_mat[i][j] not in unique_hbonds:
					unique_hbonds.append(hbonds_mat[i][j])
		# add the hydrogen bond (TYR.71.OH:TYR.71.HH:GLU.165.OE2)
		unique_hbonds.append("TYR.71.OH:TYR.71.HH:GLU.165.OE2")
		return unique_hbonds


	def merge_rxhbonds(self):
		"""
		Merge all reaction hydrogen bonds in all structure(all paths and all replicas)
		"""
		rxhbonds_mat = []
		unique_rxhbonds = []
		with Pool(processes = 4) as pool:
			wt_d1_rxhbonds = pool.apply_async((self.wt_d1).extract_rxhbonds_on_one_structure, args =())
			wt_d2_rxhbonds = pool.apply_async((self.wt_d2).extract_rxhbonds_on_one_structure, args =())
			y72_d1_rxhbonds = pool.apply_async((self.y72_d1).extract_rxhbonds_on_one_structure, args =())
			y72_d2_rxhbonds = pool.apply_async((self.y72_d2).extract_rxhbonds_on_one_structure, args =())

			wt_d1_rxhbonds_labels = wt_d1_rxhbonds.get()
			wt_d2_rxhbonds_labels = wt_d2_rxhbonds.get()
			y72_d1_rxhbonds_labels = y72_d1_rxhbonds.get()
			y72_d2_rxhbonds_labels = y72_d2_rxhbonds.get()

		rxhbonds_mat.append(wt_d1_rxhbonds_labels)
		rxhbonds_mat.append(wt_d2_rxhbonds_labels)
		rxhbonds_mat.append(y72_d1_rxhbonds_labels)
		rxhbonds_mat.append(y72_d2_rxhbonds_labels)
		for i in range(len(rxhbonds_mat)):
			for j in range(len(rxhbonds_mat[i])):
				if rxhbonds_mat[i][j] not in unique_rxhbonds:
					unique_rxhbonds.append(rxhbonds_mat[i][j])
		# add the hydrogen bond (TYR.71.OH:TYR.71.HH:GLU.165.OE2)
		return unique_rxhbonds



	def merge_chembonds(self):
		"""
		Merge all heavy atoms chemical bonds
		"""
		chembonds_mat = []
		unique_chembonds = []
		wt_d1_chembonds = (self.wt_d1).extract_chembonds_on_one_structure()
		wt_d2_chembonds = (self.wt_d2).extract_chembonds_on_one_structure()
		y72_d1_chembonds = (self.y72_d1).extract_chembonds_on_one_structure()
		y72_d2_chembonds = (self.y72_d2).extract_chembonds_on_one_structure()

		chembonds_mat.append(wt_d1_chembonds)
		chembonds_mat.append(wt_d2_chembonds)
		chembonds_mat.append(y72_d1_chembonds)
		chembonds_mat.append(y72_d2_chembonds)

		# merge all unique chem bonds label
		for i in range(len(chembonds_mat)):
			for j in range(len(chembonds_mat[i])):
				if chembonds_mat[i][j] not in unique_chembonds:
					unique_chembonds.append(chembonds_mat[i][j])
		return unique_chembonds



if __name__ == "__main__":
	_basedir = "/users/chaoy/scratch/0.proj_kpc/1.sample"
	NPath = 10
	NReplicas = 1
	Z = merge_unique_labels(_basedir, NPath, NReplicas)

	# Hbonds labels
	unique_hbonds = Z.merge_hbonds()

	# rxhbonds labels
	#unique_rxhbonds = Z.merge_rxhbonds()

	# chembonds labels
	#unique_chembonds = Z.merge_chembonds()
	print(unique_hbonds, len(unique_hbonds))
	#print(unique_rxhbonds, len(unique_rxhbonds))
	#print( unique_chembonds, len(unique_chembonds))
	#np.savez("./labels/labels", unique_hbonds=unique_hbonds, unique_rxhbonds = unique_rxhbonds, unique_chembonds = unique_chembonds)
