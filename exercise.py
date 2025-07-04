import os
import spyrmsd
import pandas as pd
import matplotlib.pyplot as plt
from spyrmsd import io, rmsd
from spyrmsd.rmsd import rmsdwrapper
from spyrmsd.molecule import Molecule
from rdkit import Chem
from rdkit.Chem import AllChem

# Actual user projectpath
userpath='/home/encelade/Documents/Master/M2_BBS/S4/coding_challenge/'


# Dockings results folders:
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

def calculate_symcorrected_rmsd(ref_path=str,met_path=str):
    """
        Return calculate rmsd between two different ligands
    """
    ref=io.loadmol(ref_path)
    met=io.loadmol(met_path)
    ref.strip()
    met.strip()

    RMSD=rmsd.symmrmsd(
        ref.coordinates,
        met.coordinates,
        ref.atomicnums,
        met.atomicnums,
        ref.adjacency_matrix,
        met.adjacency_matrix,
        minimize=False     
    )
    return RMSD

def compare_symcorr(refpath=str,lig1path=str,lig2path=str):
    """
        Return calculate rmsd between two different ligands compared to a ref 
    """
    rmsd1=calculate_symcorrected_rmsd(refpath,lig1path)
    rmsd2=calculate_symcorrected_rmsd(refpath,lig2path)
    return rmsd1,rmsd2

def prediction_symcorr_df(dico=dict, refpath=str, met1name=str, met1path=str, met2name=str, met2path=str):
    """
        Return DataFrame of the compare_symcorr() for each comparison between docking method and ref
    """
    results=[]
    for comp,sdfname in dico.items():
        rmsd1,rmsd2=compare_symcorr(refpath+sdfname, met1path+sdfname, met2path+sdfname)
        results.append({
            'Complex ID' : comp,
            'Ligand ID' : sdfname.split('.sdf')[0],
            f'Ref/{met1name} SymCorrRMSD' : rmsd1,
            f'Ref/{met2name} SymmCorrRMSD' : rmsd2
        })
    return pd.DataFrame(results)

def plot_df1(df=pd.DataFrame):
    """
        Plot the DataFrame from prediction_symcorr_df()
    """
    plt.figure(figsize=(12,6))
    plt.scatter(df[df.columns[0]], df[df.columns[2]], 
            color='blue', label=df.columns[2], alpha=0.7, marker='o')
    plt.scatter(df[df.columns[0]], df[str(df.columns[3])], 
            color='red', label=df.columns[3], alpha=0.7, marker='x')
    # Docking is considered successful when SymCorrRMSD <Å :
    plt.axhline(y=2, color='orange', linestyle='--', label='Considered docking success limit')

    plt.title('RMSD Values for Complexes from Docking Tools', fontsize=16)
    plt.xlabel(df.columns[0], fontsize=14)
    plt.ylabel('RMSD (Å)', fontsize=14)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

def plot_df(df=pd.DataFrame):
    """
    Plot the DataFrame with enhanced visibility and emphasis on key elements.
    """
    plt.figure(figsize=(10, 5))  # Taille un peu plus large pour une meilleure lisibilité

    # Scatter plots avec marqueurs plus grands et plus distincts
    plt.scatter(df[df.columns[0]], df[df.columns[2]], 
                color='blue', label=f"**{df.columns[2]}**", alpha=0.8, s=50, marker='o', edgecolor='black', linewidth=0.8)
    plt.scatter(df[df.columns[0]], df[df.columns[3]], 
                color='red', label=f"**{df.columns[3]}**", alpha=0.8, s=50, marker='x', linewidth=0.8)

    # Ligne limite avec style renforcé
    plt.axhline(y=2, color='orange', linestyle='--', label='**Docking success limit**', linewidth=1.5)

    # Titre en gras
    plt.title('**RMSD Values for Complexes from Docking Tools**', fontsize=16, weight='bold', pad=15)

    # Axes avec textes en gras et police plus grande
    plt.xlabel(f"**{df.columns[0]}**", fontsize=14, weight='bold')
    plt.ylabel('**RMSD (Å)**', fontsize=14, weight='bold')

    # Axes et ticks renforcés
    plt.xticks(fontsize=12, weight='bold')
    plt.yticks(fontsize=12, weight='bold')

    # Grille renforcée
    plt.grid(alpha=0.5, linewidth=0.8, linestyle='--')

    # Légende mise en évidence et déplacée
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12, frameon=True, shadow=True, edgecolor='black')

    # Bordures des axes plus épaisses
    ax = plt.gca()
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    # Ajustement final
    plt.tight_layout()
    plt.show()





dico=L_C_dict('exercise_data/ligand_mapping.tsv')
symcorr=prediction_symcorr_df(dico,Lexppath, 'Vina', Lvinapath,'DiffDock', Ldiffpath)
print(symcorr)
plot_df(symcorr)
