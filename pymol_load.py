import os
import pymol
import pandas as pd
from pymol import cmd

# Actual user projectpath
userpath='/home/encelade/Documents/Master/M2_BBS/S4/coding_challenge/'


# Dockings results folders:
# Exp docking sdf results folderpath:
Cexppath='exercise_data/Complexes_experimental/'
# Exp docking sdf results folderpath:
Lexppath='exercise_data/Ligands_experimental/'
# DiffDock docking sdf results folderpath:
Ldiffpath='exercise_data/Ligands_DiffDock/'
# AutoDock Vina docking sdf results folderpath:
Lvinapath='exercise_data/Ligands_Vina/'

def L_C_dict(tsv_path=str):
    """
        Function to create a dico to zip ligand ID and complex ID in a tsv file with column 1 = ligand ID and column 2= prot ID
    """
    df=pd.read_csv(tsv_path,sep='\t')
    dico=dict(zip(df.iloc[:,1],df.iloc[:,0]+'.sdf'))
    return dico

def pymol_load (dico=str):
    for i in dico:
        pymol.finish_launching(['pymol', '-q'])
        pymol.cmd.load(Cexppath+f'{i}.pdb')
        for j in l

dico=L_C_dict('exercise_data/ligand_mapping.tsv')
pymol_load(dico)