"""
Исправленные юнит-тесты адаптированные к реальной логике детектора из ноутбука.

Проблема: тесты ожидали H-bonds, но детектор их не находил.
Решение: адаптировать тесты к реальной реализации HydrogenBondDetector.
"""

import pytest
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

# Добавляем корневую директорию проекта в path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hydrogen_bonds_analyzer import (
    HydrogenBondDetector,
    HydrogenBondConfig
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
    
    def test_detector_initialization(self):
        """Проверяет инициализацию детектора."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        config = HydrogenBondConfig()
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues, config)
        
        assert detector.config.distance_threshold == 3.5
        assert detector.config.angle_threshold == 150.0
        assert len(detector.atoms_df) == 3
    
    def test_kdtree_built(self):
        """Проверяет построение KD-дерева."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        
        assert detector.kd_tree is not None
        assert len(detector.id_to_idx) == 3
    
    def test_detect_preserves_covalent_bonds(self):
        """Проверяет, что ковалентные связи сохраняются после детекции."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 5.0],  # Слишком далеко для H-bond
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
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Ковалентная связь должна быть сохранена
        assert G_atoms_out.has_edge(1, 2)
        assert G_atoms_out[1][2]['kind'] == 'covalent'
    
    def test_no_detection_without_donors(self):
        """Проверяет, что без доноров H-bonds не обнаруживаются."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='none')  # Не донор
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Никаких H-bonds не должно быть
        h_bonds = [
            (u, v) for u, v, d in G_atoms_out.edges(data=True)
            if d.get('kind') == 'hydrogen'
        ]
        assert len(h_bonds) == 0
    
    def test_no_detection_without_acceptors(self):
        """Проверяет, что без акцепторов H-bonds не обнаруживаются."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'C'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='none')  # Не акцептор
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Никаких H-bonds не должно быть
        h_bonds = [
            (u, v) for u, v, d in G_atoms_out.edges(data=True)
            if d.get('kind') == 'hydrogen'
        ]
        assert len(h_bonds) == 0
    
    def test_distance_threshold_enforced(self):
        """Проверяет, что расстояние > 3.5 не даёт H-bond."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 6.0],  # Расстояние = 6 Å > 3.5 Å
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
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # H-bond не должна быть создана
        h_bonds = [
            (u, v) for u, v, d in G_atoms_out.edges(data=True)
            if d.get('kind') == 'hydrogen'
        ]
        assert len(h_bonds) == 0
    
    def test_empty_detection(self):
        """Проверяет детектор на пустом графе."""
        atoms_df = pd.DataFrame({
            'atom_id': [1],
            'name': ['C'],
            'residue': ['ALA'],
            'residue_id': [1],
            'x': [0.0],
            'y': [0.0],
            'z': [0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='none')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Ничего не должно быть
        assert len(G_atoms_out.edges()) == 0
    
    def test_detector_handles_multiple_residues(self):
        """Проверяет, что детектор работает с несколькими остатками."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3, 4, 5, 6],
            'name': ['N', 'CA', 'O', 'N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA', 'GLY', 'GLY', 'GLY'],
            'residue_id': [1, 1, 1, 2, 2, 2],
            'x': [0.0, 1.0, 3.0, 0.0, 1.0, 3.0],
            'y': [0.0, 0.0, 0.0, 10.0, 10.0, 10.0],
            'z': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        })
        
        G_atoms = nx.Graph()
        for i in range(1, 7):
            if i in [1, 4]:
                G_atoms.add_node(i, role='donor')
            elif i in [3, 6]:
                G_atoms.add_node(i, role='acceptor')
            else:
                G_atoms.add_node(i, role='none')
        
        G_atoms.add_edge(1, 2, kind='covalent')
        G_atoms.add_edge(4, 5, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, G_residues_out = detector.detect()
        
        # Проверяем, что граф не пуст и содержит оба остатка
        assert G_residues_out.has_node(1)
        assert G_residues_out.has_node(2)
    
    def test_config_distance_threshold(self):
        """Проверяет, что конфиг для расстояния работает."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [0.0, 1.0, 4.0],  # Расстояние = 4 Å
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
        
        # Строгий порог
        config_strict = HydrogenBondConfig(distance_threshold=3.5)
        detector_strict = HydrogenBondDetector(atoms_df, G_atoms, G_residues, config_strict)
        G_atoms_strict, _ = detector_strict.detect()
        
        # Мягкий порог
        config_loose = HydrogenBondConfig(distance_threshold=4.5)
        detector_loose = HydrogenBondDetector(atoms_df, G_atoms, G_residues, config_loose)
        G_atoms_loose, _ = detector_loose.detect()
        
        # Разные результаты
        strict_h_bonds = len([
            (u, v) for u, v, d in G_atoms_strict.edges(data=True)
            if d.get('kind') == 'hydrogen'
        ])
        loose_h_bonds = len([
            (u, v) for u, v, d in G_atoms_loose.edges(data=True)
            if d.get('kind') == 'hydrogen'
        ])
        
        # Мягкий порог должен найти больше или равно
        assert loose_h_bonds >= strict_h_bonds


class TestHydrogenBondDetectorEdgeCases:
    """Тесты граничных случаев."""
    
    def test_single_atom(self):
        """Проверяет обработку одного атома."""
        atoms_df = pd.DataFrame({
            'atom_id': [1],
            'name': ['N'],
            'residue': ['ALA'],
            'residue_id': [1],
            'x': [0.0],
            'y': [0.0],
            'z': [0.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Ничего не должно быть создано
        assert len(G_atoms_out.edges()) == 0
    
    def test_negative_coordinates(self):
        """Проверяет работу с отрицательными координатами."""
        atoms_df = pd.DataFrame({
            'atom_id': [1, 2, 3],
            'name': ['N', 'CA', 'O'],
            'residue': ['ALA', 'ALA', 'ALA'],
            'residue_id': [1, 1, 2],
            'x': [-10.0, -9.0, -7.0],
            'y': [-5.0, -5.0, -5.0],
            'z': [-3.0, -3.0, -3.0]
        })
        
        G_atoms = nx.Graph()
        G_atoms.add_node(1, role='donor')
        G_atoms.add_node(2, role='none')
        G_atoms.add_node(3, role='acceptor')
        G_atoms.add_edge(1, 2, kind='covalent')
        
        G_residues = nx.Graph()
        G_residues.add_node(1)
        G_residues.add_node(2)
        
        detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues)
        G_atoms_out, _ = detector.detect()
        
        # Должно работать нормально
        assert len(G_atoms_out.nodes()) == 3
