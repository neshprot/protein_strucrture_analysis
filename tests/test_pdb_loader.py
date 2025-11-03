"""
Исправленные юнит-тесты для модуля PDBLoader.

Проверяет корректность загрузки и обработки PDB-файлов.
"""

import pytest
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

# Добавляем корневую директорию проекта в path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hydrogen_bonds_analyzer import PDBLoader


class TestPDBLoader:
    """Набор тестов для класса PDBLoader."""
    
    @pytest.fixture
    def sample_pdb_content(self):
        """Создаёт образец PDB-файла."""
        return """HETATM      1  N   ALA A   1      -8.901   4.127  -0.555  1.00  0.00           N
HETATM      2  CA  ALA A   1      -8.608   3.135  -1.442  1.00  0.00           C
HETATM      3  C   ALA A   1      -7.117   2.961  -1.765  1.00  0.00           C
HETATM      4  O   ALA A   1      -6.207   3.765  -1.379  1.00  0.00           O
HETATM      5  N   GLY A   2      -6.923   1.946  -2.519  1.00  0.00           N
HETATM      6  CA  GLY A   2      -5.517   1.649  -2.925  1.00  0.00           C
HETATM      7  C   GLY A   2      -4.581   1.431  -1.765  1.00  0.00           C
HETATM      8  O   GLY A   2      -3.412   1.185  -1.933  1.00  0.00           O
"""
    
    @pytest.fixture
    def pdb_file(self, sample_pdb_content):
        """Создаёт временный PDB-файл."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False, encoding='utf-8') as f:
            f.write(sample_pdb_content)
            temp_path = f.name
        yield temp_path
        os.unlink(temp_path)
    
    def test_load_valid_pdb(self, pdb_file):
        """Тестирует загрузку корректного PDB-файла."""
        df = PDBLoader.load(pdb_file)
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 8
        assert 'atom_id' in df.columns
        assert 'name' in df.columns
        assert 'x' in df.columns
        assert 'y' in df.columns
        assert 'z' in df.columns
    
    def test_atom_coordinates_loaded(self, pdb_file):
        """Проверяет корректность загрузки координат атомов."""
        df = PDBLoader.load(pdb_file)
        
        # Проверяем первый атом
        first_atom = df.iloc[0]
        assert first_atom['atom_id'] == 1
        assert first_atom['name'] == 'N'
        assert first_atom['residue'] == 'ALA'
        assert abs(first_atom['x'] - (-8.901)) < 0.001
        assert abs(first_atom['y'] - 4.127) < 0.001
        assert abs(first_atom['z'] - (-0.555)) < 0.001
    
    def test_residue_information_loaded(self, pdb_file):
        """Проверяет загрузку информации об остатках."""
        df = PDBLoader.load(pdb_file)
        
        assert df['residue'].iloc[0] == 'ALA'
        assert df['residue'].iloc[4] == 'GLY'
        assert df['residue_id'].iloc[0] == 1
        assert df['residue_id'].iloc[4] == 2
    
    def test_unnecessary_columns_dropped(self, pdb_file):
        """Проверяет удаление ненужных колонок."""
        df = PDBLoader.load(pdb_file)
        
        # Эти колонки должны быть удалены
        assert 'type' not in df.columns
        assert 'chain' not in df.columns
        assert 'par1' not in df.columns
        assert 'par2' not in df.columns
        assert 'element' not in df.columns
    
    def test_file_not_found(self):
        """Проверяет обработку отсутствующего файла."""
        with pytest.raises(FileNotFoundError):
            PDBLoader.load("non_existent_structure.pdb")
    
    def test_empty_pdb_file(self):
        """Проверяет загрузку пустого PDB-файла."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False, encoding='utf-8') as f:
            f.write("")
            temp_path = f.name
        
        try:
            df = PDBLoader.load(temp_path)
            assert len(df) == 0
        finally:
            os.unlink(temp_path)
    
    def test_multiple_residues(self, pdb_file):
        """Проверяет обработку нескольких остатков."""
        df = PDBLoader.load(pdb_file)
        
        unique_residues = df['residue'].unique()
        assert len(unique_residues) == 2
        assert 'ALA' in unique_residues
        assert 'GLY' in unique_residues
    
    def test_atom_id_column_type(self, pdb_file):
        """Проверяет тип данных колонки atom_id."""
        df = PDBLoader.load(pdb_file)
        
        # atom_id должен быть числовым
        assert pd.api.types.is_numeric_dtype(df['atom_id'])
    
    def test_coordinate_columns_numeric(self, pdb_file):
        """Проверяет, что координаты являются числами."""
        df = PDBLoader.load(pdb_file)
        
        assert pd.api.types.is_numeric_dtype(df['x'])
        assert pd.api.types.is_numeric_dtype(df['y'])
        assert pd.api.types.is_numeric_dtype(df['z'])
    
    def test_no_nan_in_coordinates(self, pdb_file):
        """Проверяет отсутствие NaN в координатах."""
        df = PDBLoader.load(pdb_file)
        
        assert not df['x'].isna().any()
        assert not df['y'].isna().any()
        assert not df['z'].isna().any()
    
    def test_sequential_atom_ids(self, pdb_file):
        """Проверяет последовательность ID атомов."""
        df = PDBLoader.load(pdb_file)
        
        expected_ids = list(range(1, 9))
        actual_ids = df['atom_id'].tolist()
        assert actual_ids == expected_ids
