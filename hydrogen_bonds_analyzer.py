import os
import json
import logging
import math
from pathlib import Path
from typing import Dict, Set, Tuple, List, Optional
from dataclasses import dataclass
from collections import defaultdict

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import cKDTree

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class HydrogenBondConfig:
    """Configuration for hydrogen bond detection."""
    distance_threshold: float = 3.5  # Angstroms
    angle_threshold: float = 150.0   # Degrees


class TopologyLoader:
    """Handles loading and processing protein topology data."""
    
    @staticmethod
    def load(topology_path: str) -> Dict:
        """
        Load topology data from JSON file.
        
        Args:
            topology_path: Path to topology JSON file
            
        Returns:
            Dictionary with residue topology information
            
        Raises:
            FileNotFoundError: If topology file not found
            json.JSONDecodeError: If JSON is malformed
        """
        if not os.path.exists(topology_path):
            raise FileNotFoundError(f"Topology file not found: {topology_path}")
        
        logger.info(f"Loading topology from {topology_path}")
        try:
            with open(topology_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse JSON: {e}")
            raise
        
        # Build adjacency dictionaries for efficient bond lookups
        for residue_name, residue_data in data.items():
            adjacency = defaultdict(set)
            for u, v in residue_data.get('bonds', []):
                adjacency[u].add(v)
                adjacency[v].add(u)
            residue_data['adjancey_bonds'] = dict(adjacency)
        
        logger.info(f"Loaded topology for {len(data)} residue types")
        return data


class PDBLoader:
    """Handles loading PDB structure files."""
    
    @staticmethod
    def load(pdb_path: str) -> pd.DataFrame:
        """
        Load atomic coordinates from PDB file.
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            DataFrame with atom information
            
        Raises:
            FileNotFoundError: If PDB file not found
        """
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        
        logger.info(f"Loading PDB from {pdb_path}")
        column_names = [
            'type', 'atom_id', 'name', 'residue', 'chain',
            'residue_id', 'x', 'y', 'z', 'par1', 'par2', 'element'
        ]
        atoms_df = pd.read_csv(
            pdb_path,
            sep=r'\s+',
            header=None,
            names=column_names
        )
        
        # Drop unnecessary columns
        atoms_df = atoms_df.drop(columns=['type', 'chain', 'par1', 'par2', 'element'])
        logger.info(f"Loaded {len(atoms_df)} atoms")
        return atoms_df


class GraphBuilder:
    """Builds protein structure graphs from topology and coordinates."""
    
    def __init__(self, atoms_df: pd.DataFrame, topology: Dict):
        """
        Initialize graph builder.
        
        Args:
            atoms_df: DataFrame with atom information
            topology: Dictionary with topology data
        """
        self.atoms_df = atoms_df
        self.topology = topology
        self.G_atoms = nx.Graph()
        self.G_residues = nx.Graph()
    
    def build(self) -> Tuple[nx.Graph, nx.Graph]:
        """
        Build atom and residue graphs with covalent bonds.
        
        Returns:
            Tuple of (atom_graph, residue_graph)
        """
        logger.info("Building molecular graphs...")
        
        for res_id in self.atoms_df['residue_id'].unique():
            self._process_residue(res_id)
        
        logger.info(f"Graphs built: {self.G_atoms.number_of_nodes()} atoms, "
                   f"{self.G_residues.number_of_nodes()} residues")
        return self.G_atoms, self.G_residues
    
    def _process_residue(self, res_id: int) -> None:
        """Process a single residue and add it to graphs."""
        res_df = self.atoms_df[self.atoms_df['residue_id'] == res_id]
        res_name = res_df['residue'].iat[0]
        
        if res_name not in self.topology:
            logger.warning(f"Unknown residue type: {res_name}")
            return
        
        res_topology = self.topology[res_name]
        existing_bonds = res_topology.get('adjancey_bonds', {})
        
        # Build donor/acceptor role mappings
        donor_role = {pair[0]: 'donor' for pair in res_topology.get('donors', [])}
        acceptor_role = {pair[0]: 'acceptor' for pair in res_topology.get('acceptors', [])}
        
        # Add nodes and intra-residue covalent bonds
        self._add_atoms_and_bonds(res_df, existing_bonds, donor_role, acceptor_role)
        
        # Add inter-residue bonds (C-N peptide bonds)
        self._add_peptide_bonds(res_id, res_df, res_topology)
        
        # Add residue node
        self.G_residues.add_node(res_id, kind='covalent')
    
    def _add_atoms_and_bonds(self, res_df: pd.DataFrame, existing_bonds: Dict,
                             donor_role: Dict, acceptor_role: Dict) -> None:
        """Add atoms and intra-residue covalent bonds."""
        res_atoms = res_df['atom_id'].values
        
        # Add nodes
        for atom_id in res_atoms:
            atom_name = self.atoms_df.at[atom_id - 1, 'name']
            role = donor_role.get(atom_name) or acceptor_role.get(atom_name) or 'none'
            self.G_atoms.add_node(int(atom_id), role=role)
        
        # Add covalent bonds
        atom_names = {aid: self.atoms_df.at[aid - 1, 'name'] for aid in res_atoms}
        for idx, atom_id_1 in enumerate(res_atoms):
            for atom_id_2 in res_atoms[idx + 1:]:
                name_1 = atom_names[atom_id_1]
                name_2 = atom_names[atom_id_2]
                if name_2 in existing_bonds.get(name_1, []):
                    self.G_atoms.add_edge(int(atom_id_1), int(atom_id_2), kind='covalent')
    
    def _add_peptide_bonds(self, res_id: int, res_df: pd.DataFrame,
                          res_topology: Dict) -> None:
        """Add C-N peptide bonds connecting consecutive residues."""
        if '+N' not in res_topology.get('adjancey_bonds', {}).get('C', []):
            return
        
        c_atoms = res_df[res_df['name'] == 'C']['atom_id']
        if c_atoms.empty:
            return
        
        next_res = self.atoms_df[self.atoms_df['residue_id'] == res_id + 1]
        n_atoms = next_res[next_res['name'] == 'N']['atom_id']
        
        for c_atom in c_atoms:
            for n_atom in n_atoms:
                c_atom_int = int(c_atom)
                n_atom_int = int(n_atom)
                self.G_atoms.add_node(c_atom_int, role='none')
                self.G_atoms.add_node(n_atom_int, role='none')
                self.G_atoms.add_edge(c_atom_int, n_atom_int, kind='covalent')
                self.G_residues.add_edge(res_id, res_id + 1, kind='covalent')


class HydrogenBondDetector:
    """Detects hydrogen bonds based on geometric criteria."""
    
    def __init__(self, atoms_df: pd.DataFrame, G_atoms: nx.Graph,
                 G_residues: nx.Graph, config: Optional[HydrogenBondConfig] = None):
        """
        Initialize hydrogen bond detector.
        
        Args:
            atoms_df: DataFrame with atom coordinates
            G_atoms: Atom-level graph
            G_residues: Residue-level graph
            config: Detection configuration
        """
        self.atoms_df = atoms_df
        self.G_atoms = G_atoms
        self.G_residues = G_residues
        self.config = config or HydrogenBondConfig()
        self.coords = atoms_df[['x', 'y', 'z']].to_numpy()
        self.id_to_idx = {int(aid): idx for idx, aid in enumerate(atoms_df['atom_id'])}
        self.kd_tree = cKDTree(self.coords)
    
    def detect(self) -> Tuple[nx.Graph, nx.Graph]:
        """
        Detect hydrogen bonds in the structure.
        
        Returns:
            Tuple of (updated_atom_graph, updated_residue_graph)
        """
        logger.info("Detecting hydrogen bonds...")
        
        # Get donor and acceptor atoms
        donors = [
            int(aid) for aid, props in self.G_atoms.nodes(data=True)
            if props.get('role') == 'donor'
        ]
        acceptors = {
            int(aid) for aid, props in self.G_atoms.nodes(data=True)
            if props.get('role') == 'acceptor'
        }
        
        logger.info(f"Found {len(donors)} donors and {len(acceptors)} acceptors")
        
        hbond_count = 0
        for donor_id in donors:
            di = self.id_to_idx[donor_id]
            d_coord = self.coords[di]
            
            # Find all atoms within distance threshold
            neighbor_indices = self.kd_tree.query_ball_point(
                d_coord, r=self.config.distance_threshold
            )
            
            for nj in neighbor_indices:
                acceptor_id = int(self.atoms_df.iloc[nj]['atom_id'])
                if acceptor_id not in acceptors:
                    continue
                
                # Check angle criterion
                if self._check_angle(donor_id, acceptor_id, di, nj):
                    self._add_hydrogen_bond(donor_id, acceptor_id, di, nj)
                    hbond_count += 1
        
        logger.info(f"Detected {hbond_count} hydrogen bonds")
        return self.G_atoms, self.G_residues
    
    def _check_angle(self, donor_id: int, acceptor_id: int, di: int, ai: int) -> bool:
        """Check if donor-acceptor angle meets threshold."""
        d_coord = self.coords[di]
        a_coord = self.coords[ai]
        
        # Find reference atom bonded to donor
        neighbors = list(self.G_atoms.neighbors(donor_id))
        if not neighbors:
            return False
        
        ref_id = neighbors[0]
        ri = self.id_to_idx[ref_id]
        r_coord = self.coords[ri]
        
        # Calculate angle
        v1 = a_coord - d_coord  # donor -> acceptor
        v2 = r_coord - d_coord  # donor -> reference
        
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        theta = math.degrees(math.acos(np.clip(cos_theta, -1.0, 1.0)))
        
        return theta >= self.config.angle_threshold
    
    def _add_hydrogen_bond(self, donor_id: int, acceptor_id: int,
                           di: int, ai: int) -> None:
        """Add hydrogen bond to both atom and residue graphs."""
        self.G_atoms.add_edge(donor_id, acceptor_id, kind='hydrogen')
        
        donor_res_id = int(self.atoms_df.iloc[di]['residue_id'])
        acceptor_res_id = int(self.atoms_df.iloc[ai]['residue_id'])
        self.G_residues.add_edge(donor_res_id, acceptor_res_id, kind='hydrogen')


class HydrogenBondComparator:
    """Compares hydrogen bonding patterns between structures."""
    
    def __init__(self, atoms_df_1: pd.DataFrame, atoms_df_2: pd.DataFrame,
                 G1_atoms: nx.Graph, G2_atoms: nx.Graph):
        """
        Initialize comparator.
        
        Args:
            atoms_df_1: Reference structure atoms
            atoms_df_2: Query structure atoms
            G1_atoms: Reference structure atom graph
            G2_atoms: Query structure atom graph
        """
        self.atoms_df_1 = atoms_df_1
        self.atoms_df_2 = atoms_df_2
        self.G1_atoms = G1_atoms
        self.G2_atoms = G2_atoms
        self.comparison_df = pd.DataFrame()
    
    def compare(self, exclude_water: bool = True) -> pd.DataFrame:
        """
        Compare hydrogen bonds between two structures.
        
        Args:
            exclude_water: Whether to exclude water (HOH) residues
            
        Returns:
            DataFrame with comparison results
        """
        logger.info("Comparing hydrogen bonds between structures...")
        
        # Get hydrogen bonds from both structures
        h_bonds_1 = self._extract_hydrogen_bonds(self.G1_atoms)
        h_bonds_2 = self._extract_hydrogen_bonds(self.G2_atoms)
        
        # Merge structures for atom mapping
        merged = pd.merge(
            self.atoms_df_1, self.atoms_df_2,
            on=['residue_id', 'residue', 'name'],
            how='outer',
            indicator=True,
            suffixes=('_old', '_new')
        )
        
        # Process bonds from structure 1
        self._process_reference_bonds(h_bonds_1, h_bonds_2, merged)
        
        # Process bonds only in structure 2
        self._process_query_bonds(h_bonds_2, merged)
        
        # Post-process: remove water if requested
        if exclude_water:
            self._remove_water_bonds()
        
        logger.info(f"Comparison complete: {len(self.comparison_df)} bond interactions")
        return self.comparison_df
    
    @staticmethod
    def _extract_hydrogen_bonds(G_atoms: nx.Graph) -> Set[Tuple[int, int]]:
        """Extract all hydrogen bonds from atom graph."""
        return {
            tuple(sorted((int(u), int(v))))
            for u, v, attrs in G_atoms.edges(data=True)
            if attrs.get('kind') == 'hydrogen'
        }
    
    def _process_reference_bonds(self, h_bonds_1: Set, h_bonds_2: Set,
                                 merged: pd.DataFrame) -> None:
        """Process bonds from reference structure."""
        h_bonds_2_copy = h_bonds_2.copy()
        
        for donor_id, acceptor_id in h_bonds_1:
            row = {
                'atom_id_old_donor': donor_id,
                'atom_id_old_acceptor': acceptor_id,
            }
            
            # Get reference structure bond info
            donor_info = self._get_atom_info(self.atoms_df_1, donor_id, '_old')
            acceptor_info = self._get_atom_info(self.atoms_df_1, acceptor_id, '_old')
            row.update(donor_info)
            row.update(acceptor_info)
            
            # Get mapping to query structure
            try:
                mut_donor_id = int(merged.loc[merged['atom_id_old'] == donor_id, 'atom_id_new'].iloc[0])
                mut_acceptor_id = int(merged.loc[merged['atom_id_old'] == acceptor_id, 'atom_id_new'].iloc[0])
            except (ValueError, IndexError):
                row['status'] = 'lost'
                self.comparison_df = pd.concat([self.comparison_df, pd.DataFrame([row])],
                                              ignore_index=True)
                continue
            
            row['atom_id_new_donor'] = mut_donor_id
            row['atom_id_new_acceptor'] = mut_acceptor_id
            
            # Check if bond exists in query structure
            if (mut_donor_id, mut_acceptor_id) in h_bonds_2_copy:
                row['status'] = 'conserved'
                h_bonds_2_copy.remove((mut_donor_id, mut_acceptor_id))
            else:
                row['status'] = 'lost'
            
            # Add query structure bond info if conserved
            if row['status'] == 'conserved':
                mut_donor_info = self._get_atom_info(self.atoms_df_2, mut_donor_id, '_new')
                mut_acceptor_info = self._get_atom_info(self.atoms_df_2, mut_acceptor_id, '_new')
                row.update(mut_donor_info)
                row.update(mut_acceptor_info)
            
            self.comparison_df = pd.concat([self.comparison_df, pd.DataFrame([row])],
                                          ignore_index=True)
        
        # Process bonds only in query structure
        for donor_id, acceptor_id in h_bonds_2_copy:
            row = {
                'atom_id_new_donor': donor_id,
                'atom_id_new_acceptor': acceptor_id,
                'status': 'gained'
            }
            
            donor_info = self._get_atom_info(self.atoms_df_2, donor_id, '_new')
            acceptor_info = self._get_atom_info(self.atoms_df_2, acceptor_id, '_new')
            row.update(donor_info)
            row.update(acceptor_info)
            
            self.comparison_df = pd.concat([self.comparison_df, pd.DataFrame([row])],
                                          ignore_index=True)
    
    def _process_query_bonds(self, h_bonds_2: Set, merged: pd.DataFrame) -> None:
        """Process bonds only in query structure (already handled in _process_reference_bonds)."""
        pass
    
    @staticmethod
    def _get_atom_info(atoms_df: pd.DataFrame, atom_id: int, suffix: str) -> Dict:
        """Extract atom information from DataFrame."""
        try:
            atom_row = atoms_df[atoms_df['atom_id'] == atom_id].iloc[0]
            return {
                f'name{suffix}': atom_row['name'],
                f'residue{suffix}': atom_row['residue'],
                f'residue_id{suffix}': atom_row['residue_id']
            }
        except IndexError:
            return {}
    
    def _remove_water_bonds(self) -> None:
        """Remove hydrogen bonds involving water molecules."""
        water_residues = ['HOH', 'WAT']
        mask = (
            self.comparison_df['residue_old'].isin(water_residues) |
            self.comparison_df['residue_new'].isin(water_residues)
        )
        self.comparison_df = self.comparison_df[~mask].reset_index(drop=True)


class PipelineOrchestrator:
    """Orchestrates the complete analysis pipeline."""
    
    def __init__(self, config: HydrogenBondConfig):
        """Initialize pipeline."""
        self.config = config
    
    def analyze_structure(self, pdb_path: str, topology_path: str,
                         output_dir: str) -> Tuple[nx.Graph, nx.Graph]:
        """
        Analyze a single protein structure.
        
        Args:
            pdb_path: Path to PDB file
            topology_path: Path to topology file
            output_dir: Directory for output files
            
        Returns:
            Tuple of (atom_graph, residue_graph)
        """
        # Load data
        topology = TopologyLoader.load(topology_path)
        atoms_df = PDBLoader.load(pdb_path)
        
        # Build graphs
        builder = GraphBuilder(atoms_df, topology)
        G_atoms, G_residues = builder.build()
        
        # Detect hydrogen bonds
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues, self.config)
        G_atoms, G_residues = detector.detect()
        
        # Save graphs
        self._save_graphs(G_atoms, G_residues, output_dir)
        
        return G_atoms, G_residues
    
    def compare_structures(self, pdb_1: str, pdb_2: str, topology_path: str,
                          output_dir: str) -> pd.DataFrame:
        """
        Compare hydrogen bonding between two structures.
        
        Args:
            pdb_1: Path to first PDB file
            pdb_2: Path to second PDB file
            topology_path: Path to topology file
            output_dir: Directory for output files
            
        Returns:
            DataFrame with comparison results
        """
        logger.info("Starting structure comparison analysis...")
        
        # Analyze both structures
        G1_atoms, _ = self.analyze_structure(pdb_1, topology_path, output_dir)
        G2_atoms, _ = self.analyze_structure(pdb_2, topology_path, output_dir)
        
        # Load atom data
        atoms_df_1 = PDBLoader.load(pdb_1)
        atoms_df_2 = PDBLoader.load(pdb_2)
        
        # Compare
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        comparison_df = comparator.compare(exclude_water=True)
        
        # Save results
        output_path = os.path.join(output_dir, 'hydrogen_bonds_comparison.csv')
        comparison_df.to_csv(output_path, index=False)
        logger.info(f"Comparison results saved to {output_path}")
        
        return comparison_df
    
    @staticmethod
    def _save_graphs(G_atoms: nx.Graph, G_residues: nx.Graph,
                    output_dir: str) -> None:
        """Save graphs to GraphML files."""
        os.makedirs(output_dir, exist_ok=True)
        
        atoms_path = os.path.join(output_dir, 'atoms.graphml')
        residues_path = os.path.join(output_dir, 'residues.graphml')
        
        nx.write_graphml(G_atoms, atoms_path)
        nx.write_graphml(G_residues, residues_path)
        
        logger.info(f"Graphs saved to {output_dir}")


def main():
    """Main entry point."""
    import sys
    from dotenv import load_dotenv
    
    # Load environment variables
    load_dotenv()
    
    # Configuration
    topology_path = os.getenv('TOPOLOGY_FILE', './data/topology_complete.json')
    output_dir = os.getenv('OUTPUT_DIR', './output')
    
    config = HydrogenBondConfig(
        distance_threshold=float(os.getenv('HB_DISTANCE_THRESHOLD', 3.5)),
        angle_threshold=float(os.getenv('HB_ANGLE_THRESHOLD', 150.0))
    )
    
    orchestrator = PipelineOrchestrator(config)
    
    if len(sys.argv) == 2:
        # Single structure analysis
        pdb_path = sys.argv[1]
        orchestrator.analyze_structure(pdb_path, topology_path, output_dir)
    elif len(sys.argv) == 3:
        # Comparison analysis
        pdb_1 = sys.argv[1]
        pdb_2 = sys.argv[2]
        orchestrator.compare_structures(pdb_1, pdb_2, topology_path, output_dir)
    else:
        print("Usage:")
        print("  Single structure: python hydrogen_bonds_analyzer.py <pdb_file>")
        print("  Compare structures: python hydrogen_bonds_analyzer.py <pdb_1> <pdb_2>")


if __name__ == '__main__':
    main()
