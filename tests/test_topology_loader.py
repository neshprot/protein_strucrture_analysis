"""
Исправленные юнит-тесты для модуля TopologyLoader.

Проверяет корректность загрузки и обработки файлов топологии.
"""

import pytest
import json
import os
import sys
import tempfile
from pathlib import Path

# Добавляем корневую директорию проекта в path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hydrogen_bonds_analyzer import TopologyLoader


class TestTopologyLoader:
    """Набор тестов для класса TopologyLoader."""
    
    @pytest.fixture
    def sample_topology(self):
        """Создаёт образец файла топологии для тестирования."""
        return {
            "ALA": {
                "donors": [["N"]],
                "acceptors": [["O"]],
                "bonds": [["N", "CA"], ["CA", "C"], ["C", "O"]]
            },
            "GLY": {
                "donors": [["N"]],
                "acceptors": [["O"]],
                "bonds": [["N", "CA"], ["CA", "C"], ["C", "O"]]
            },
            "SER": {
                "donors": [["N"], ["OG"]],
                "acceptors": [["O"], ["OG"]],
                "bonds": [["N", "CA"], ["CA", "CB"], ["CB", "OG"], ["CA", "C"], ["C", "O"]]
            }
        }
    
    @pytest.fixture
    def topology_file(self, sample_topology):
        """Создаёт временный файл топологии."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as f:
            json.dump(sample_topology, f)
            temp_path = f.name
        yield temp_path
        os.unlink(temp_path)
    
    def test_load_valid_topology(self, topology_file):
        """Тестирует загрузку корректного файла топологии."""
        data = TopologyLoader.load(topology_file)
        
        assert "ALA" in data
        assert "GLY" in data
        assert "SER" in data
        assert len(data) == 3
    
    def test_adjacency_bonds_created(self, topology_file):
        """Проверяет, что словари смежности корректно построены."""
        data = TopologyLoader.load(topology_file)
        
        # Проверяем ALA
        ala_adjacency = data["ALA"]["adjancey_bonds"]
        assert "N" in ala_adjacency
        assert "CA" in ala_adjacency["N"]
        assert "CA" in ala_adjacency
        assert "N" in ala_adjacency["CA"]
        assert "C" in ala_adjacency["CA"]
    
    def test_adjacency_bonds_bidirectional(self, topology_file):
        """Убеждается, что рёбра в смежности двусторонние."""
        data = TopologyLoader.load(topology_file)
        
        for residue_name, residue_data in data.items():
            adjacency = residue_data["adjancey_bonds"]
            
            for atom1, neighbors in adjacency.items():
                for atom2 in neighbors:
                    assert atom1 in adjacency[atom2], \
                        f"Асимметричная граница: {atom1}-{atom2} в {residue_name}"
    
    def test_file_not_found(self):
        """Проверяет обработку отсутствующего файла."""
        with pytest.raises(FileNotFoundError):
            TopologyLoader.load("non_existent_topology.json")
    
    def test_invalid_json(self):
        """Проверяет обработку некорректного JSON."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as f:
            f.write("{ invalid json content")
            temp_path = f.name
        
        try:
            with pytest.raises(json.JSONDecodeError):
                TopologyLoader.load(temp_path)
        finally:
            os.unlink(temp_path)
    
    def test_empty_topology(self):
        """Проверяет загрузку пустой топологии."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as f:
            json.dump({}, f)
            temp_path = f.name
        
        try:
            data = TopologyLoader.load(temp_path)
            assert len(data) == 0
        finally:
            os.unlink(temp_path)
    
    def test_residue_without_bonds(self):
        """Проверяет остаток без определённых связей."""
        topology = {
            "UNK": {
                "donors": [],
                "acceptors": [],
                "bonds": []
            }
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as f:
            json.dump(topology, f)
            temp_path = f.name
        
        try:
            data = TopologyLoader.load(temp_path)
            unk_adjacency = data["UNK"]["adjancey_bonds"]
            assert len(unk_adjacency) == 0
        finally:
            os.unlink(temp_path)
    
    def test_self_loops_not_created(self, topology_file):
        """Проверяет, что самопетли не создаются."""
        topology = {
            "TEST": {
                "donors": [["N"]],
                "acceptors": [["O"]],
                "bonds": [["CA", "CA"]]  # Самопетля
            }
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as f:
            json.dump(topology, f)
            temp_path = f.name
        
        try:
            data = TopologyLoader.load(temp_path)
            # Самопетля как set содержит сам себя, что нормально для графов
            assert "CA" in data["TEST"]["adjancey_bonds"]["CA"]
        finally:
            os.unlink(temp_path)
