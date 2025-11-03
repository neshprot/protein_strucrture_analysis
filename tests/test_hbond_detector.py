"""
Юнит-тесты для модуля HydrogenBondDetector.

Проверяет корректность обнаружения водородных связей.
"""

import pytest
import numpy as np
import pandas as pd
import networkx as nx
from hydrogen_bonds_analyzer import (
    HydrogenBondDetector,
    HydrogenBondConfig,
    GraphBuilder
)


class TestHydrogenBondConfig:
    """Тесты для конфигурации детекции."""
    
    def test_default_config(self):
        """Проверяет значения конфигурации по умолчанию."""
        config = HydrogenBondConfig()
        
        assert config.distance_threshold == 3.5
        assert config.angle_threshold == 150.0
    
    def test_custom_config(self):
        """Проверяет создание пользовательской конфигурации."""
        config = HydrogenBondConfig(
            distance_threshold=3.2,
            angle_threshold=160.0
        )
        
        assert config.distance_threshold == 3.2
        assert config.angle_threshold == 160.0


class TestHydrogenBondDetector:
    """Набор тестов для класса HydrogenBondDetector."""
    
    @pytest.fixture
    def simple_structure(self):
        """Создаёт простую структуру для тестирования.
        
        Симулирует две молекулы с потенциальной водородной связью:
        - Атом 1 (N): донор при (0, 0, 0)
        - Атом 2 (reference): при (1, 0, 0)
        - Атом 3 (O): акцептор при (3, 0, 0) - должна быть H-bond с углом 180°
        """
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        # Создаём простой граф
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        return atoms_df, G_atoms, G_residues
    
    def test_detector_initialization(self, simple_structure):
        """Проверяет инициализацию детектора."""
        atoms_df, G_atoms, G_residues = simple_structure
        config = HydrogenBondConfig()
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues, config)
        
        assert detector.config.distance_threshold == 3.5
        assert detector.config.angle_threshold == 150.0
        assert len(detector.atoms_df) == 3
    
    def test_kdtree_built(self, simple_structure):
        """Проверяет построение KD-дерева."""
        atoms_df, G_atoms, G_residues = simple_structure
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        
        assert detector.kd_tree is not None
        assert len(detector.id_to_idx) == 3
    
    def test_detect_ideal_hydrogen_bond(self, simple_structure):
        """Проверяет обнаружение идеальной водородной связи."""
        atoms_df, G_atoms, G_residues = simple_structure
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, G_residues_out = detector.detect()
        
        # Проверяем, что водородная связь найдена
        assert G_atoms_out.has_edge(1, 3)
        edge_data = G_atoms_out[1][3]
        assert edge_data['kind'] == 'hydrogen'
    
    def test_distance_threshold_respected(self):
        """Проверяет соблюдение порога расстояния."""
        # Создаём две точки на расстоянии > 3.5 Å
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 5.0],  # Расстояние N-O = 5 Å
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(
            atoms_df, G_atoms, G_residues,
            HydrogenBondConfig(distance_threshold=3.5)
        )
        G_atoms_out, _ = detector.detect()
        
        # Связь не должна быть обнаружена
        assert not G_atoms_out.has_edge(1, 3) or G_atoms_out[1][3]['kind'] != 'hydrogen'
    
    def test_angle_threshold_respected(self):
        """Проверяет соблюдение углового порога."""
        # Создаём структуру с плохим углом (90° вместо 180°)
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 0.0],  # Угол 90°, а не 180°
            'y': [0.0, 0.0, 3.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(
            atoms_df, G_atoms, G_residues,
            HydrogenBondConfig(angle_threshold=150.0)
        )
        G_atoms_out, _ = detector.detect()
        
        # Связь не должна быть обнаружена из-за плохого угла
        assert not G_atoms_out.has_edge(1, 3) or G_atoms_out[1][3]['kind'] != 'hydrogen'
    
    def test_only_donors_and_acceptors_considered(self):
        """Проверяет, что учитываются только доноры и акцепторы."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['C', 'CA', 'C'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='none')  # Ни донор, ни акцептор
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='none')
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Никаких водородных связей не должно быть
        h_bonds = [
            (u, v) for u, v, d in G_atoms_out.edges(data=True)
            if d['kind'] == 'hydrogen'
        ]
        assert len(h_bonds) == 0
    
    def test_hydrogen_bonds_in_residue_graph(self, simple_structure):
        """Проверяет добавление водородных связей в граф остатков."""
        atoms_df, G_atoms, G_residues = simple_structure
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        _, G_residues_out = detector.detect()
        
        # Проверяем связь между остатками
        if G_residues_out.has_edge(1, 2):
            edge_data = G_residues_out[1][2]
            assert edge_data['kind'] == 'hydrogen'
    
    def test_no_hydrogen_bonds_without_reference_atom(self):
        """Проверяет, что H-связь не создаётся без atom соседа донора."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 3],
            'name': ['N', 'O'],
            'residue': ['ALA', 'ALA'],
            'residue_id': [1, 2],
            'x': [0.0, 3.0],
            'y': [0.0, 0.0],
            'z': [0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')  # Нет соседей
        G_atoms.add_node(3, role='acceptor')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Связь не должна быть обнаружена
        assert not G_atoms_out.has_edge(1, 3) or G_atoms_out[1][3]['kind'] != 'hydrogen'
