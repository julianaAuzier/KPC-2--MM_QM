from merge_all_unique_labels import merge_unique_labels
import numpy as np

if __name__ == "__main__":
	dir = "/users/chaoy/scratch/0.proj_kpc/1.sample"
	NPath = 200
	NReplicas = 1# 36
	XZ = merge_unique_labels(dir, NPath, NReplicas)

	# Hbonds labels
	unique_hbonds = XZ.merge_hbonds()

	print(unique_hbonds, len(unique_hbonds))
	np.save("./reactant_label/hbonds_labels_reactant_da_ct.npy", unique_hbonds)
