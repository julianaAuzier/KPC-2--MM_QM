from bonds import Bonds

class batch_extract():
    def __init__(self, base_dir, NPath, NReplicas, sysname, whichstate):
        self.base_dir = base_dir
        self.NPath = NPath
        self.NReplicas = NReplicas
        self.sysname =  sysname
        self.whichstate = whichstate

    def extract_hbonds_on_one_structure(self):
        hbonds_cub = []
        for i in  range(1, self.NPath+1):
            hbonds_mat = []
            for j in range(1, self.NReplicas+1):
                X = Bonds(self.base_dir, i, j, self.NPath, self.NReplicas, self.sysname, self.whichstate)
                hbondslabel = X.hbond_labels()
                hbonds_mat.append(hbondslabel)
            hbonds_cub.append(hbonds_mat)
        uniqe_hbonds = []
        for i in range(len(hbonds_cub)):
            for j in range(len(hbonds_cub[i])):
                for k in range(len(hbonds_cub[i][j])):
                    if hbonds_cub[i][j][k] not in uniqe_hbonds:
                        uniqe_hbonds.append(hbonds_cub[i][j][k])



        return uniqe_hbonds

    def extract_rxhbonds_on_one_structure(self):
        rxhbonds_cub = []
        for i in range(1, self.NPath+1):
            rxhbonds_mat = []
            for j in range(1, self.NReplicas+1):
                X = Bonds(self.base_dir, i, j, self.NPath, self.NReplicas, self.sysname, self.whichstate)
                rxhbondlabels = X.rxhbond_labels()
                rxhbonds_mat.append(rxhbondlabels)
            rxhbonds_cub.append(rxhbonds_mat)

        unique_rxhbonds = []
        for i in range(len(rxhbonds_cub)):
            for j in range(len(rxhbonds_cub[i])):
                for k in range(len(rxhbonds_cub[i][j])):
                    if rxhbonds_cub[i][j][k] not in unique_rxhbonds:
                        unique_rxhbonds.append(rxhbonds_cub[i][j][k])


        return unique_rxhbonds


    def extract_chembonds_on_one_structure(self):
        chembonds_cub = []
        for i in range(1, self.NPath+1):
            chembonds_mat = []
            for j in range(1, self.NReplicas+1):
                X = Bonds(self.base_dir, i, j, self.NPath, self.NReplicas, self.sysname, self.whichstate)
                chembondslabels = X.hvychembonds_label()
                chembonds_mat.append(chembondslabels)
            chembonds_cub.append(chembonds_mat)
        unique_chembonds = []
        for i in range(len(chembonds_cub)):
            for j in range(len(chembonds_cub[i])):
                for k in range(len(chembonds_cub[i][j])):
                    if chembonds_cub[i][j][k] not in uniqe_chembonds:
                        uniqe_chembonds.append(chembonds_cub[i][j][k])
        return unique_chembonds




if __name__ == "__main__":
    _basedir = "/users/chaoy/scratch/0.proj_kpc/1.sample"
    NPath = 10
    NReplicas = 1
    sysname = 'kpc_wt'
    whichstate = 'd1'
    X = batch_extract(_basedir, NPath, NReplicas, sysname, whichstate)
    hbondslables = X.extract_hbonds_on_one_structure()
    print(hbondslables)
    #rxhbondslabel = X.extract_rxhbonds_on_one_structure()
    #print(rxhbondslabel)
    #chemlabels = X.extract_chembonds_on_one_structure()
    #print(chemlabels)
