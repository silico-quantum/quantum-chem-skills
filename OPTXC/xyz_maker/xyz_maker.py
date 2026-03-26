import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os

class XYZGenerator:
    def __init__(self, output_dir='xyz_molecules'):
        """
        XYZ File Generator
        
        Parameters:
        output_dir: Directory path for saving XYZ files
        """
        self.output_dir = output_dir
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
    
    def generate_xyz(self, smiles, index=0):
        """
        Generate XYZ file for a single molecule
        
        Parameters:
        smiles: SMILES string
        index: Molecule index (for filename generation)
        
        Returns:
        Absolute path to the generated XYZ file
        """
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        state = AllChem.EmbedMolecule(mol, useRandomCoords=True)
        if state == -1:
            raise RuntimeError(f"Failed to generate 3D coordinates for SMILES: {smiles}")
        
        # Optimize molecular geometry using MMFF force field
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Generate filename based on index or SMILES hash
        if index is None:
            # Use SMILES hash as filename if index not provided
            filename = f"molecule_{abs(hash(smiles))}.xyz"
        else:
            filename = f"job_{index+1}.xyz"
        
        # Save XYZ file
        xyz_path = os.path.join(self.output_dir, filename)
        Chem.MolToXYZFile(mol, xyz_path)
        
        return os.path.abspath(xyz_path)
    
    def generate_xyz_batch(self, smiles_list):
        """
        Batch generate XYZ files for multiple molecules
        
        Parameters:
        smiles_list: List of SMILES strings
        
        Returns:
        List of absolute paths to generated XYZ files
        """
        results = []
        
        for i, smiles in enumerate(smiles_list):
            try:
                # Generate XYZ file for current SMILES
                xyz_path = self.generate_xyz(smiles, i)
                results.append(xyz_path)
                print(f"Saved XYZ file: {xyz_path}")
            except Exception as e:
                # Handle errors while keeping process running
                print(f"Error processing SMILES {smiles}: {str(e)}")
                results.append(None)
        
        return results

    def generate_xyz_from_csv(self, csv_path, smiles_column='SMILES'):
        """
        Generate XYZ files from SMILES in a CSV file
        
        Parameters:
        csv_path: Path to input CSV file
        smiles_column: Column name containing SMILES strings
        
        Returns:
        List of absolute paths to generated XYZ files
        """
        try:
            # Read input CSV file
            df = pd.read_csv(csv_path)
            if smiles_column not in df.columns:
                raise ValueError(f"Column not found in CSV file: {smiles_column}")
            
            # Extract SMILES list and generate XYZ files
            smiles_list = df[smiles_column].tolist()
            return self.generate_xyz_batch(smiles_list)
        except Exception as e:
            # Handle file reading errors
            print(f"Error reading CSV file: {str(e)}")
            return []