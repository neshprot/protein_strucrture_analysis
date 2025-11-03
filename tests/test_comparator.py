"""
Исправленные юнит-тесты для модуля HydrogenBondComparator.

Проверяет корректность сравнения водородных связей между структурами.
"""

import pytest
import sys
from pathlib import Path

import pandas as pd
import networkx as nx

# Добавляем корневую директорию проекта в path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hydrogen_bonds_analyzer import HydrogenBondComparator


class TestHydrogenBondComparator:
    """Набор тестов для класса HydrogenBondComparator."""
    
    @pytest.fixture
    def simple_atoms_dataframes(self):
        """Создаёт две простые структуры атомов для сравнения."""
        atoms_df_1 = pd.DataFrame({
            'atom_id': [1, 2, 3, 4],
            'name': ['N', 'CA', 'O', 'OG'],
            'residue': ['ALA', 'ALA', 'ALA', 'SER'],
            'residue_id': [1, 1, 1, 2],
            'x': [0.0, 1.0, 2.0, 3.0],
            'y': [0.0, 0.0, 0.0, 1.0],
            'z': [0.0, 0.0, 0.0, 0.0]
        })
        
        atoms_df_2 = pd.DataFrame({
            'atom_id': [1, 2, 3, 4],
            'name': ['N', 'CA', 'O', 'OG'],
            'residue': ['ALA', 'ALA', 'ALA', 'SER'],
            'residue_id': [1, 1, 1, 2],
            'x': [0.1, 1.1, 2.1, 3.5],
            'y': [0.1, 0.1, 0.1, 1.2],
            'z': [0.0, 0.0, 0.0, 0.0]
        })
        
        return atoms_df_1, atoms_df_2
    
    @pytest.fixture
    def simple_graphs(self):
        """Создаёт простые графы для тестирования."""
        # Граф 1: водородная связь 1-3
        G1_atoms = nx.Graph()
        G1_atoms.add_node(1, role='donor')
        G1_atoms.add_node(2, role='none')
        G1_atoms.add_node(3, role='acceptor')
        G1_atoms.add_node(4, role='acceptor')
        G1_atoms.add_edge(1, 3, kind='hydrogen')
        
        # Граф 2: водородная связь 1-3 сохранена, добавлена 1-4
        G2_atoms = nx.Graph()
        G2_atoms.add_node(1, role='donor')
        G2_atoms.add_node(2, role='none')
        G2_atoms.add_node(3, role='acceptor')
        G2_atoms.add_node(4, role='acceptor')
        G2_atoms.add_edge(1, 3, kind='hydrogen')
        G2_atoms.add_edge(1, 4, kind='hydrogen')
        
        return G1_atoms, G2_atoms
    
    def test_comparator_initialization(self, simple_atoms_dataframes, simple_graphs):
        """Проверяет инициализацию компаратора."""
        atoms_df_1, atoms_df_2 = simple_atoms_dataframes
        G1_atoms, G2_atoms = simple_graphs
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        
        assert comparator.atoms_df_1 is not None
        assert comparator.atoms_df_2 is not None
        assert len(comparator.comparison_df) == 0
    
    def test_conserved_bond_detected(self, simple_atoms_dataframes, simple_graphs):
        """Проверяет обнаружение сохранённой водородной связи."""
        atoms_df_1, atoms_df_2 = simple_atoms_dataframes
        G1_atoms, G2_atoms = simple_graphs
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=False)
        
        # Проверяем, что есть строка со статусом 'conserved'
        conserved = result_df[result_df['status'] == 'conserved']
        assert len(conserved) > 0
    
    def test_gained_bond_detected(self, simple_atoms_dataframes, simple_graphs):
        """Проверяет обнаружение приобретённой водородной связи."""
        atoms_df_1, atoms_df_2 = simple_atoms_dataframes
        G1_atoms, G2_atoms = simple_graphs
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=False)
        
        # Проверяем, что есть строка со статусом 'gained'
        gained = result_df[result_df['status'] == 'gained']
        assert len(gained) > 0
    
    def test_lost_bond_detected(self):
        """Проверяет обнаружение потерянной водородной связи."""
        atoms_df_1 = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 2.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        atoms_df_2 = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 10.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G1_atoms = nx.Graph()
        G1_atoms.add_node(1, role='donor')
        G1_atoms.add_node(2, role='none')
        G1_atoms.add_node(3, role='acceptor')
        G1_atoms.add_edge(1, 3, kind='hydrogen')
        
        G2_atoms = nx.Graph()
        G2_atoms.add_node(1, role='donor')
        G2_atoms.add_node(2, role='none')
        G2_atoms.add_node(3, role='acceptor')
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=False)
        
        # Проверяем, что есть строка со статусом 'lost'
        lost = result_df[result_df['status'] == 'lost']
        assert len(lost) > 0
    
    def test_water_exclusion(self):
        """Проверяет исключение водородных связей с водой."""
        atoms_df_1 = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'HOH'],
            'residue_id': [1, 1, 3],
            'x': [0.0, 1.0, 2.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        atoms_df_2 = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'HOH'],
            'residue_id': [1, 1, 3],
            'x': [0.0, 1.0, 2.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G1_atoms = nx.Graph()
        G1_atoms.add_node(1, role='donor')
        G1_atoms.add_node(2, role='none')
        G1_atoms.add_node(3, role='acceptor')
        G1_atoms.add_edge(1, 3, kind='hydrogen')
        
        G2_atoms = nx.Graph()
        G2_atoms.add_node(1, role='donor')
        G2_atoms.add_node(2, role='none')
        G2_atoms.add_node(3, role='acceptor')
        G2_atoms.add_edge(1, 3, kind='hydrogen')
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=True)
        
        # После исключения воды, должно быть меньше или вообще 0 связей
        assert len(result_df) == 0
    
    def test_output_columns_present(self, simple_atoms_dataframes, simple_graphs):
        """Проверяет наличие всех необходимых колонок в результате."""
        atoms_df_1, atoms_df_2 = simple_atoms_dataframes
        G1_atoms, G2_atoms = simple_graphs
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=False)
        
        # Проверяем наличие ключевых колонок
        expected_columns = ['status']
        for col in expected_columns:
            assert col in result_df.columns
    
    def test_status_values(self, simple_atoms_dataframes, simple_graphs):
        """Проверяет валидные значения статуса."""
        atoms_df_1, atoms_df_2 = simple_atoms_dataframes
        G1_atoms, G2_atoms = simple_graphs
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=False)
        
        # Проверяем, что все значения статуса из допустимого набора
        valid_statuses = {'conserved', 'lost', 'gained', 'only_wt', 'only_mut'}
        for status in result_df['status'].unique():
            assert status in valid_statuses
    
    def test_empty_graphs_comparison(self):
        """Проверяет сравнение пустых графов."""
        atoms_df_1 = pd.DataFrame({
            'atom_id': [1, 2],
            'name': ['N', 'O'],
            'residue': ['ALA', 'ALA'],
            'residue_id': [1, 1],
            'x': [0.0, 2.0],
            'y': [0.0, 0.0],
            'z': [0.0, 0.0]
        })
        
        atoms_df_2 = pd.DataFrame({
            'atom_id': [1, 2],
            'name': ['N', 'O'],
            'residue': ['ALA', 'ALA'],
            'residue_id': [1, 1],
            'x': [0.0, 2.0],
            'y': [0.0, 0.0],
            'z': [0.0, 0.0]
        })
        
        # Пустые графы (без H-bonds)
        G1_atoms = nx.Graph()
        G1_atoms.add_node(1, role='none')
        G1_atoms.add_node(2, role='none')
        
        G2_atoms = nx.Graph()
        G2_atoms.add_node(1, role='none')
        G2_atoms.add_node(2, role='none')
        
        comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
        result_df = comparator.compare(exclude_water=False)
        
        # Результат должен быть пустым
        assert len(result_df) == 0
    
    def test_extract_hydrogen_bonds(self):
        """Проверяет извлечение водородных связей из графа."""
        G = nx.Graph()
        G.add_edge(1, 2, kind='hydrogen')
        G.add_edge(2, 3, kind='covalent')
        G.add_edge(3, 4, kind='hydrogen')
        
        h_bonds = HydrogenBondComparator._extract_hydrogen_bonds(G)
        
        assert len(h_bonds) == 2
        assert (1, 2) in h_bonds
        assert (3, 4) in h_bonds
        assert (2, 3) not in h_bonds
