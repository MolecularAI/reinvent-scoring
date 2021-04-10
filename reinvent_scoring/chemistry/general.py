from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol
from rdkit.DataStructs.cDataStructs import UIntSparseIntVect


class GeneralChemistry:

    @staticmethod
    def smiles_to_mols_and_indices(query_smiles: List[str]) -> Tuple[List[Mol], List[int]]:
        mols = [Chem.MolFromSmiles(smile) for smile in query_smiles]
        valid_mask = [mol is not None for mol in mols]
        valid_idxs = [idx for idx, is_valid in enumerate(valid_mask) if is_valid]
        valid_mols = [mols[idx] for idx in valid_idxs]
        return valid_mols, valid_idxs

    @staticmethod
    def mols_to_fingerprints(molecules: List[Mol], radius=3, use_counts=True, use_features=True) \
            -> List[UIntSparseIntVect]:
        fingerprints = [AllChem.GetMorganFingerprint(mol, radius, useCounts=use_counts, useFeatures=use_features) for
                        mol in molecules]
        return fingerprints

    @staticmethod
    def smiles_to_mols(query_smiles: List[str]) -> List[Mol]:
        mols = [Chem.MolFromSmiles(smile) for smile in query_smiles]
        valid_mask = [mol is not None for mol in mols]
        valid_idxs = [idx for idx, is_valid in enumerate(valid_mask) if is_valid]
        valid_mols = [mols[idx] for idx in valid_idxs]
        return valid_mols

    def smiles_to_fingerprints(self, query_smiles: List[str], radius=3, use_counts=True, use_features=True) -> List[
        UIntSparseIntVect]:
        mols = self.smiles_to_mols(query_smiles)
        fingerprints = self.mols_to_fingerprints(mols, radius=radius, use_counts=use_counts, use_features=use_features)
        return fingerprints

    def smile_to_mol(self, smile: str) -> Mol:
        """
        Creates a Mol object from a SMILES string.
        :param smile: SMILES string.
        :return: A Mol object or None if it's not valid.
        """
        if smile:
            return MolFromSmiles(smile)

    def mols_to_smiles(self, molecules: List[Mol]) -> List[str]:
        """This method assumes that all molecules are valid."""
        valid_smiles = [MolToSmiles(mol, isomericSmiles=False) for mol in molecules]
        return valid_smiles
