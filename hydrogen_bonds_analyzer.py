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


@dataclass
class ChargeInteractionConfig:
    """Configuration for charge-charge interaction detection."""
    dist_threshold: float = 4.5  # Angstroms
    charge_threshold: float = 0.3  # Minimum absolute charge to consider


class TopologyLoader:
    """Handles loading and processing protein topology data."""
    
    @staticmethod
    def load(topology_path: str) -> Dict:
        """Load topology data from JSON file."""
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

    @staticmethod
    def get_charged_atoms_for_residue(topology: Dict, residue_name: str,
                                      charge_threshold: float = 0.3) -> Dict[str, float]:
        """
        Extract charged atoms from a residue definition in topology.
        
        Args:
            topology: Topology dictionary
            residue_name: Name of residue (e.g., 'ARG', 'ASP')
            charge_threshold: Minimum absolute charge to consider
            
        Returns:
            Dictionary mapping atom names to charges
        """
        if residue_name not in topology:
            return {}
        
        res_data = topology[residue_name]
        charged_atoms = {}
        
        for atom_name, atom_data in res_data.get('atoms', {}).items():
            charge = atom_data.get('charge', 0.0)
            if abs(charge) >= charge_threshold:
                charged_atoms[atom_name] = charge
        
        return charged_atoms

    @staticmethod
    def get_sidechain_atoms_for_residue(topology: Dict, residue_name: str) -> Set[str]:
        """
        Extract side-chain atoms from a residue definition.
        
        Backbone atoms (N, CA, C, O, H) are excluded.
        """
        if residue_name not in topology:
            return set()
        
        backbone = {'N', 'CA', 'C', 'O', 'H', 'HN', 'OXT', 'HO'}
        res_data = topology[residue_name]
        side_chain = {
            name for name in res_data.get('atoms', {}).keys()
            if name not in backbone
        }
        
        return side_chain

    @staticmethod
    def get_charged_sidechain_atoms(topology: Dict, residue_name: str,
                                   charge_threshold: float = 0.3) -> Dict[str, float]:
        """
        Extract charged side-chain atoms from a residue definition.
        Combines charged atoms and side-chain atom filtering.
        """
        if residue_name not in topology:
            return {}
        
        charged = TopologyLoader.get_charged_atoms_for_residue(topology, residue_name, charge_threshold)
        side_chain = TopologyLoader.get_sidechain_atoms_for_residue(topology, residue_name)
        
        # Return only charged atoms that are in side chain
        return {name: charge for name, charge in charged.items() if name in side_chain}


class PDBLoader:
    """Handles loading PDB structure files."""
    
    @staticmethod
    def load(pdb_path: str) -> pd.DataFrame:
        """Load atomic coordinates from PDB file."""
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
        """Initialize graph builder."""
        self.atoms_df = atoms_df
        self.topology = topology
        self.G_atoms = nx.Graph()
        self.G_residues = nx.Graph()
    
    def build(self) -> Tuple[nx.Graph, nx.Graph]:
        """Build atom and residue graphs with covalent bonds."""
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
        """Initialize hydrogen bond detector."""
        self.atoms_df = atoms_df
        self.G_atoms = G_atoms
        self.G_residues = G_residues
        self.config = config or HydrogenBondConfig()
        self.coords = atoms_df[['x', 'y', 'z']].to_numpy()
        self.id_to_idx = {int(aid): idx for idx, aid in enumerate(atoms_df['atom_id'])}
        self.kd_tree = cKDTree(self.coords)
    
    def detect(self) -> Tuple[nx.Graph, nx.Graph]:
        """Detect hydrogen bonds in the structure."""
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
        """Initialize comparator."""
        self.atoms_df_1 = atoms_df_1
        self.atoms_df_2 = atoms_df_2
        self.G1_atoms = G1_atoms
        self.G2_atoms = G2_atoms
        self.comparison_df = pd.DataFrame()
    
    def compare(self, exclude_water: bool = True) -> pd.DataFrame:
        """Compare hydrogen bonds between two structures."""
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
        
        # Process bonds only in query structure (gained bonds)
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


class ChargeInteractionDetector:
    """Detects charge-charge interactions based on distance and topology charges."""
    
    def __init__(self, atomsdf: pd.DataFrame, Gatoms: nx.Graph, topology: Dict,
                 config: Optional[ChargeInteractionConfig] = None):
        """
        Initialize charge interaction detector.
        
        Args:
            atomsdf: DataFrame with atom coordinates
            Gatoms: Atom-level graph
            topology: Topology dictionary with charge information
            config: Detection configuration
        """
        self.atomsdf = atomsdf
        self.Gatoms = Gatoms
        self.topology = topology
        self.config = config or ChargeInteractionConfig()
        self.coords = atomsdf[['x', 'y', 'z']].to_numpy()
        self.charged_atoms = self._get_charged_atoms_from_topology()
        if self.charged_atoms:
            self.kdtree = cKDTree(self.coords)
        else:
            self.kdtree = None

    def _get_charged_atoms_from_topology(self) -> Dict[int, float]:
        """Extract charged atoms from topology for this structure."""
        charged = {}
        
        for _, row in self.atomsdf.iterrows():
            resname = row['residue']
            atomname = row['name']
            atomid = int(row['atom_id'])
            
            # Get charged side-chain atoms from topology
            charged_sidechain = TopologyLoader.get_charged_sidechain_atoms(
                self.topology, resname, self.config.charge_threshold
            )
            
            if atomname in charged_sidechain:
                charged[atomid] = charged_sidechain[atomname]
        
        logger.info(f"Found {len(charged)} charged atoms in structure")
        return charged

    def detect(self) -> nx.Graph:
        """Annotate Gatoms graph with charge-charge interactions."""
        if not self.charged_atoms or self.kdtree is None:
            logger.info("No charged atoms found or KDTree not available")
            return self.Gatoms
        
        atom_ids = list(self.charged_atoms.keys())
        charges = np.array([self.charged_atoms[aid] for aid in atom_ids])
        points = np.array([
            self.atomsdf[self.atomsdf['atom_id'] == aid][['x', 'y', 'z']].values[0]
            for aid in atom_ids
        ])
        
        tree = cKDTree(points)
        charge_count = 0
        
        for i, aid1 in enumerate(atom_ids):
            charge1 = charges[i]
            neighbors = tree.query_ball_point(points[i], self.config.dist_threshold)
            
            for j in neighbors:
                aid2 = atom_ids[j]
                charge2 = charges[j]
                
                if aid1 >= aid2:  # Avoid duplicates and self-interactions
                    continue
                
                # Opposite charges (electrostatic attraction)
                if charge1 * charge2 < 0:
                    distance = np.linalg.norm(points[i] - points[j])
                    self.Gatoms.add_edge(
                        aid1, aid2, 
                        kind='charge', 
                        charge1=charge1,
                        charge2=charge2,
                        distance=distance
                    )
                    charge_count += 1
        
        logger.info(f"Detected {charge_count} charge-charge interactions")
        return self.Gatoms


class ChargeInteractionComparator:
    """Compares charge-charge interactions between structures."""
    
    def __init__(self, atomsdf1: pd.DataFrame, atomsdf2: pd.DataFrame,
                 G1atoms: nx.Graph, G2atoms: nx.Graph):
        """Initialize comparator."""
        self.atomsdf1 = atomsdf1
        self.atomsdf2 = atomsdf2
        self.G1atoms = G1atoms
        self.G2atoms = G2atoms
        self.comparison_df = pd.DataFrame()

    @staticmethod
    def extract_charge_edges(G: nx.Graph) -> Set[Tuple[int, int]]:
        """Extract all charge interaction edges from graph."""
        return {
            tuple(sorted((u, v))) 
            for u, v, d in G.edges(data=True) 
            if d.get('kind') == 'charge'
        }

    def compare(self, exclude_water: bool = True) -> pd.DataFrame:
        """
        Compare charge interactions between two proteins.
        
        Output columns mirror hydrogen bond comparison for consistency:
        - atom_id_old_1/2: Old structure atom IDs (first and second charged atom)
        - atom_id_new_1/2: New structure atom IDs
        - name_old_1/2: Atom names in old structure
        - residue_old_1/2: Residue names in old structure
        - residue_id_old_1/2: Residue indices in old structure
        - charge1/2_old: Charges from old structure
        - name_new_1/2: Atom names in new structure (if conserved/gained)
        - residue_new_1/2: Residue names in new structure (if conserved/gained)
        - residue_id_new_1/2: Residue indices in new structure (if conserved/gained)
        - charge1/2_new: Charges from new structure (if conserved/gained)
        - status: conserved | lost | gained
        """
        logger.info("Comparing charge-charge interactions...")
        
        charges1 = self.extract_charge_edges(self.G1atoms)
        charges2 = self.extract_charge_edges(self.G2atoms)

        # Atom-level mapping
        merged = pd.merge(
            self.atomsdf1, self.atomsdf2, 
            on=['residue_id', 'residue', 'name'], 
            how='outer', 
            indicator=True, 
            suffixes=('_old', '_new')
        )

        comparison = []
        
        # Process charge interactions from structure 1
        for aid1, aid2 in charges1:
            row = {}
            
            # Old structure: atom IDs
            row['atom_id_old_1'] = aid1
            row['atom_id_old_2'] = aid2
            
            # Old structure: atom 1 info
            info1_old = self._get_atom_info(self.atomsdf1, aid1)
            row['name_old_1'] = info1_old.get('name')
            row['residue_old_1'] = info1_old.get('residue')
            row['residue_id_old_1'] = info1_old.get('residue_id')
            
            # Old structure: atom 2 info
            info2_old = self._get_atom_info(self.atomsdf1, aid2)
            row['name_old_2'] = info2_old.get('name')
            row['residue_old_2'] = info2_old.get('residue')
            row['residue_id_old_2'] = info2_old.get('residue_id')
            
            # Get charges from graph
            edge_data = self.G1atoms.get_edge_data(aid1, aid2)
            row['charge1_old'] = edge_data.get('charge1') if edge_data else None
            row['charge2_old'] = edge_data.get('charge2') if edge_data else None
            row['distance_old'] = edge_data.get('distance') if edge_data else None
            
            # Try to map to new structure
            try:
                mut_aid1 = int(merged.loc[merged['atom_id_old'] == aid1, 'atom_id_new'].values[0])
                mut_aid2 = int(merged.loc[merged['atom_id_old'] == aid2, 'atom_id_new'].values[0])
            except (ValueError, IndexError):
                status = "lost"
                mut_aid1, mut_aid2 = None, None
            else:
                if (mut_aid1, mut_aid2) in charges2 or (mut_aid2, mut_aid1) in charges2:
                    status = "conserved"
                    charges2.discard((mut_aid1, mut_aid2))
                    charges2.discard((mut_aid2, mut_aid1))
                else:
                    status = "lost"
            
            row['atom_id_new_1'] = mut_aid1
            row['atom_id_new_2'] = mut_aid2
            row['status'] = status
            
            # If conserved or lost, still try to get new structure info
            if mut_aid1 and mut_aid2:
                info1_new = self._get_atom_info(self.atomsdf2, mut_aid1)
                row['name_new_1'] = info1_new.get('name')
                row['residue_new_1'] = info1_new.get('residue')
                row['residue_id_new_1'] = info1_new.get('residue_id')
                
                info2_new = self._get_atom_info(self.atomsdf2, mut_aid2)
                row['name_new_2'] = info2_new.get('name')
                row['residue_new_2'] = info2_new.get('residue')
                row['residue_id_new_2'] = info2_new.get('residue_id')
                
                edge_data_new = self.G2atoms.get_edge_data(mut_aid1, mut_aid2)
                row['charge1_new'] = edge_data_new.get('charge1') if edge_data_new else None
                row['charge2_new'] = edge_data_new.get('charge2') if edge_data_new else None
                row['distance_new'] = edge_data_new.get('distance') if edge_data_new else None
            else:
                row['name_new_1'] = None
                row['residue_new_1'] = None
                row['residue_id_new_1'] = None
                row['name_new_2'] = None
                row['residue_new_2'] = None
                row['residue_id_new_2'] = None
                row['charge1_new'] = None
                row['charge2_new'] = None
                row['distance_new'] = None
            
            comparison.append(row)
        
        # Process charge interactions only in structure 2 (gained)
        for aid1, aid2 in charges2:
            row = {}
            
            # No old structure info for gained interactions
            row['atom_id_old_1'] = None
            row['atom_id_old_2'] = None
            row['name_old_1'] = None
            row['residue_old_1'] = None
            row['residue_id_old_1'] = None
            row['name_old_2'] = None
            row['residue_old_2'] = None
            row['residue_id_old_2'] = None
            row['charge1_old'] = None
            row['charge2_old'] = None
            row['distance_old'] = None
            
            # New structure info
            row['atom_id_new_1'] = aid1
            row['atom_id_new_2'] = aid2
            
            info1_new = self._get_atom_info(self.atomsdf2, aid1)
            row['name_new_1'] = info1_new.get('name')
            row['residue_new_1'] = info1_new.get('residue')
            row['residue_id_new_1'] = info1_new.get('residue_id')
            
            info2_new = self._get_atom_info(self.atomsdf2, aid2)
            row['name_new_2'] = info2_new.get('name')
            row['residue_new_2'] = info2_new.get('residue')
            row['residue_id_new_2'] = info2_new.get('residue_id')
            
            edge_data_new = self.G2atoms.get_edge_data(aid1, aid2)
            row['charge1_new'] = edge_data_new.get('charge1') if edge_data_new else None
            row['charge2_new'] = edge_data_new.get('charge2') if edge_data_new else None
            row['distance_new'] = edge_data_new.get('distance') if edge_data_new else None
            
            row['status'] = 'gained'
            
            comparison.append(row)
        
        # Create DataFrame with organized columns
        if comparison:
            self.comparison_df = pd.DataFrame(comparison)
            
            # Reorder columns for easy analysis
            col_order = [
                'atom_id_old_1', 'atom_id_old_2',
                'name_old_1', 'name_old_2',
                'residue_old_1', 'residue_old_2',
                'residue_id_old_1', 'residue_id_old_2',
                'charge1_old', 'charge2_old',
                'distance_old',
                'atom_id_new_1', 'atom_id_new_2',
                'name_new_1', 'name_new_2',
                'residue_new_1', 'residue_new_2',
                'residue_id_new_1', 'residue_id_new_2',
                'charge1_new', 'charge2_new',
                'distance_new',
                'status'
            ]
            self.comparison_df = self.comparison_df[col_order]
        
        # Post-process: remove water if requested
        if exclude_water:
            self._remove_water_interactions()
        
        logger.info(f"Charge comparison complete: {len(self.comparison_df)} interactions")
        return self.comparison_df

    @staticmethod
    def _get_atom_info(df: pd.DataFrame, atomid: int) -> Dict:
        """Extract atom information from DataFrame."""
        try:
            row = df[df['atom_id'] == atomid].iloc[0]
            return {
                'name': row['name'],
                'residue': row['residue'],
                'residue_id': row['residue_id']
            }
        except IndexError:
            return {'name': None, 'residue': None, 'residue_id': None}

    def _remove_water_interactions(self) -> None:
        """Remove charge interactions involving water molecules."""
        water_residues = ['HOH', 'WAT']
        mask = (
            self.comparison_df['residue_old_1'].isin(water_residues) |
            self.comparison_df['residue_old_2'].isin(water_residues) |
            self.comparison_df['residue_new_1'].isin(water_residues) |
            self.comparison_df['residue_new_2'].isin(water_residues)
        )
        self.comparison_df = self.comparison_df[~mask].reset_index(drop=True)


class PipelineOrchestrator:
    """Orchestrates the complete analysis pipeline."""
    
    def __init__(self, hb_config: Optional[HydrogenBondConfig] = None,
                 charge_config: Optional[ChargeInteractionConfig] = None):
        """Initialize pipeline with configurations."""
        self.hb_config = hb_config or HydrogenBondConfig()
        self.charge_config = charge_config or ChargeInteractionConfig()
    
    def analyze_structure(self, pdb_path: str, topology_path: str,
                         output_dir: str) -> Tuple[nx.Graph, nx.Graph, pd.DataFrame]:
        """
        Analyze a single protein structure.
        
        Returns:
            Tuple of (atom_graph, residue_graph, atoms_df)
        """
        # Load data
        topology = TopologyLoader.load(topology_path)
        atoms_df = PDBLoader.load(pdb_path)
        
        # Build graphs
        builder = GraphBuilder(atoms_df, topology)
        G_atoms, G_residues = builder.build()
        
        # Detect hydrogen bonds
        hb_detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues, self.hb_config)
        G_atoms, G_residues = hb_detector.detect()

        # Detect charge-charge interactions (now uses topology)
        charge_detector = ChargeInteractionDetector(atoms_df, G_atoms, topology, self.charge_config)
        G_atoms = charge_detector.detect()
        
        # Save graphs
        self._save_graphs(G_atoms, G_residues, output_dir)
        
        return G_atoms, G_residues, atoms_df
    
    def compare_structures(self, pdb_1: str, pdb_2: str, topology_path: str,
                          output_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Compare hydrogen bonding and charge interactions between two structures.
        
        Returns:
            Tuple of (hb_comparison_df, charge_comparison_df)
        """
        logger.info("Starting structure comparison analysis...")
        
        # Analyze both structures
        G1_atoms, _, atoms_df_1 = self.analyze_structure(pdb_1, topology_path, output_dir)
        G2_atoms, _, atoms_df_2 = self.analyze_structure(pdb_2, topology_path, output_dir)

        # Hydrogen bond comparison
        hb_comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        hb_comparison_df = hb_comparator.compare()

        # Charge interaction comparison
        charge_comparator = ChargeInteractionComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        charge_comparison_df = charge_comparator.compare()
        
        # Save results
        hb_comparison_df.to_csv(f"{output_dir}/hydrogen_bonds_comparison.csv", index=False)
        charge_comparison_df.to_csv(f"{output_dir}/charge_interactions_comparison.csv", index=False)
        logger.info(f"Comparison results saved to {output_dir}")
        
        return hb_comparison_df, charge_comparison_df
    
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
    
    hb_config = HydrogenBondConfig(
        distance_threshold=float(os.getenv('HB_DISTANCE_THRESHOLD', 3.5)),
        angle_threshold=float(os.getenv('HB_ANGLE_THRESHOLD', 150.0))
    )
    
    charge_config = ChargeInteractionConfig(
        dist_threshold=float(os.getenv('CHARGE_DIST_THRESHOLD', 4.5)),
        charge_threshold=float(os.getenv('CHARGE_THRESHOLD', 0.3))
    )
    
    orchestrator = PipelineOrchestrator(hb_config, charge_config)
    
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
        print("  Single structure: python script.py <pdb_file>")
        print("  Compare structures: python script.py <pdb_1> <pdb_2>")


if __name__ == '__main__':
    main()
