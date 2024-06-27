import biotite.database.rcsb as rcsb
import os
from biotite.structure.io.pdbx import get_structure
import biotite.structure.io.pdbx as pdbx
import biotite
import numpy as np

CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'
mapping = {'Aba': 'A', 'Ace': 'X', 'Acr': 'X', 'Ala': 'A', 'Aly': 'K', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cas': 'C',
           'Ccs': 'C', 'Cme': 'C', 'Csd': 'C', 'Cso': 'C', 'Csx': 'C', 'Cys': 'C', 'Dal': 'A', 'Dbb': 'T', 'Dbu': 'T',
           'Dha': 'S', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'Glz': 'G', 'His': 'H', 'Hse': 'S', 'Ile': 'I', 'Leu': 'L',
           'Llp': 'K', 'Lys': 'K', 'Men': 'N', 'Met': 'M', 'Mly': 'K', 'Mse': 'M', 'Nh2': 'X', 'Nle': 'L', 'Ocs': 'C',
           'Pca': 'E', 'Phe': 'F', 'Pro': 'P', 'Ptr': 'Y', 'Sep': 'S', 'Ser': 'S', 'Thr': 'T', 'Tih': 'A', 'Tpo': 'T',
           'Trp': 'W', 'Tyr': 'Y', 'Unk': 'X', 'Val': 'V', 'Ycm': 'C', 'Sec': 'U', 'Pyl': 'O', 'Mhs': 'H', 'Snm': 'S',
           'Mis': 'S', 'Seb': 'S', 'Hic': 'H', 'Fme': 'M', 'Asb': 'D', 'Sah': 'C', 'Smc': 'C', 'Tpq': 'Y', 'Onl': 'X',
           'Tox': 'W', '5x8': 'X', 'Ddz': 'A'}


def three_to_one(three_letter_code):
    if three_letter_code[0].upper() + three_letter_code[1:].lower() not in mapping:
        return 'X'
    return mapping[three_letter_code[0].upper() + three_letter_code[1:].lower()]

def get_c_alphas(input_path):

    counter = 0
    ids = [filename[:-4] for filename in os.listdir(input_path)]
    results = {}

    for id in ids:
        print(f'{counter} / {len(ids)}: {id}.txt')
        counter += 1

        pdb_id = id[:4]
        chain_id = id[4:]
        cif_file_path = rcsb.fetch(pdb_id, "cif", CIF_FILES_PATH)

        cif_file = pdbx.CIFFile.read(cif_file_path)

        structure = get_structure(cif_file, model=1)
        all_atoms = structure[(structure.chain_id == chain_id) & (biotite.structure.filter_peptide_backbone(structure))]

        # the following code taken from the biotite source code:
        # https://github.com/biotite-dev/biotite/blob/v0.41.0/src/biotite/structure/residues.py#L22
        chain_id_changes = (all_atoms.chain_id[1:] != all_atoms.chain_id[:-1])
        res_id_changes   = (all_atoms.res_id[1:]   != all_atoms.res_id[:-1]  )
        ins_code_changes = (all_atoms.ins_code[1:] != all_atoms.ins_code[:-1])
        res_name_changes = (all_atoms.res_name[1:] != all_atoms.res_name[:-1])
        residue_change_mask = (
            chain_id_changes |
            res_id_changes |
            ins_code_changes |
            res_name_changes
        )

        residue_starts = np.where(residue_change_mask)[0] + 1
        residue_starts = np.concatenate(([0], residue_starts))

        # take the C_alphas if possible (now we have the N atoms, however, regularly the C_alphas are located right behind the N atoms, so increase the indices by 1)
        if len(all_atoms) > residue_starts[-1] + 1:
            residue_starts = residue_starts + 1

        c_alphas = all_atoms[residue_starts]

        results[id] = c_alphas

    return results

