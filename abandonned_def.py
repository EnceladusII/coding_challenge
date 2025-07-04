# TESTING ZONE

def list_sdf(path=str):
    """
        Function for append in a list all sdf files names in a folder
    """
    sdf_list = []
    for root, folder, files in os.walk(path):
        for fname in files:
            if fname.endswith('.sdf'):
                sdf_list.append(fname)
    return sdf_list

def list_pdb(path=str):
    """
        Function for append in a  list all pdb files names in a folder
    """
    pdb_list = []
    for root, folder, files in os.walk(path):
        for fname in files:
            if fname.endswith('.pdb'):
                pdb_list.append(fname)
    return pdb_list

def calculate_rmsd(mol1_path=str,mol2_path=str):
    """
        Return calculate rmsd between two different ligands
    """
    mol1=io.loadmol(mol1_path)
    mol2=io.loadmol(mol2_path)

    mol1.strip()
    mol2.strip()

    rmsd=rmsdwrapper(mol1,mol2)
    return rmsd

def compare_ligands(refpath=str,lig1path=str,lig2path=str):
    """
        Return calculate rmsd between two different ligands compared to a ref 
    """
    rmsd1=calculate_rmsd(refpath,lig1path)
    rmsd2=calculate_rmsd(refpath,lig2path)
    return rmsd1[0],rmsd2[0]

def prediction_methods_df(dico=dict, refpath=str, met1name=str, met1path=str, met2name=str, met2path=str):
    """
    
    """
    results=[]
    for comp,sdfname in dico.items():
        rmsd1,rmsd2=compare_ligands(refpath+sdfname, met1path+sdfname, met2path+sdfname)
        results.append({
            'Complex ID' : comp,
            'Ligand ID' : sdfname.split('.sdf')[0],
            f'Ref/{met1name} RMSD' : rmsd1,
            f'Ref/{met2name} RMSD' : rmsd2
        })
    return pd.DataFrame(results)