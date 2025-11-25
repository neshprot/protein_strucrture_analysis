# Анализатор водородных связей и электростатических взаимодействий


## Обзор


**Анализатор водородных связей и электростатических взаимодействий** позволяет:
- Обнаруживать водородные связи в белковых структурах, используя геометрические критерии (расстояние и угловые пороги)
- Обнаруживать электростатические взаимодействия (солевые мостики) между заряженными боковыми цепями на основе топологии
- Строить графы взаимодействий на уровне атомов и остатков
- Сравнивать паттерны водородных связей и электростатических взаимодействий между двумя белковыми структурами (например, дикий тип vs. мутант)
- Генерировать детальные отчёты о сохранённых, потерянных и приобретённых взаимодействиях


## Ключевые функции


- **Модульная архитектура**: разделение ответственности с отдельными классами для каждого компонента
- **Водородные связи**: детекция на основе геометрии (расстояние + угол)
- **Электростатические взаимодействия**: детекция солевых мостиков (ARG/LYS/HIS ↔ ASP/GLU)
- **Топология-ориентированный подход**: все параметры заряда извлекаются из JSON файла топологии
- **Production-ready код**: комплексное логирование, обработка ошибок и type hints
- **Оптимизированная производительность**: использование cKDTree для эффективного пространственного поиска
- **Гибкая конфигурация**: настраиваемые пороги через `HydrogenBondConfig` и `ChargeInteractionConfig`
- **Сравнительный анализ**: детальное отслеживание изменений взаимодействий между структурами
- **Графовое представление**: графы NetworkX для взаимодействий на уровне атомов и остатков
- **Пакетная обработка**: поддержка анализа одной структуры или сравнения двух
- **Исключение воды**: опциональное удаление взаимодействий с молекулами воды (HOH)


## Установка


### Требования
- Python 3.8+
- pip или conda package manager


### Зависимости
```bash
pip install numpy pandas networkx scipy python-dotenv
```


### Клонирование и настройка
```bash
git clone <repository-url>
cd hydrogen-bond-analyzer
pip install -r requirements.txt
```


### requirements.txt
```
numpy>=1.20.0
pandas>=1.3.0
networkx>=2.6
scipy>=1.7.0
python-dotenv>=0.19.0
pytest>=7.0
```


## Конфигурация


### Переменные окружения
Создайте файл `.env` в корне проекта:


```env
# Topology file with residue definitions and charges
TOPOLOGY_FILE=./data/raw/topology_complete.json

# Output directory for results
OUTPUT_DIR=./output

# Hydrogen bond detection thresholds (optional, defaults provided)
HB_DISTANCE_THRESHOLD=3.5
HB_ANGLE_THRESHOLD=150.0

# Charge interaction detection thresholds (optional, defaults provided)
CHARGE_DIST_THRESHOLD=4.5
CHARGE_THRESHOLD=0.3
```


### Параметры детекции водородных связей


Два геометрических критерия определяют водородную связь:


1. **Порог расстояния** (по умолчанию: 3.5 Å)
   - Максимальное расстояние между атомами донора и акцептора


2. **Угловой порог** (по умолчанию: 150°)
   - Минимальный угол при атоме донора (реф-донор-акцептор)
   - Значения ближе к 180° указывают на лучшее выравнивание


### Параметры детекции электростатических взаимодействий


1. **Порог расстояния** (по умолчанию: 4.5 Å)
   - Максимальное расстояние между двумя заряженными атомами


2. **Порог заряда** (по умолчанию: 0.3 e)
   - Минимальный абсолютный заряд атома для учёта
   - Извлекается из поля `charge` в JSON топологии


3. **Условие взаимодействия**
   - Детектируются только противоположные заряды (привлекательное взаимодействие)
   - Одинаковые заряды игнорируются (отталкивание, не стабилизирующее взаимодействие)


## Использование


### Анализ одной структуры


Обнаружение всех водородных связей и электростатических взаимодействий в белковой структуре:


```bash
python hydrogen_bonds_analyzer.py structure.pdb
```


**Выходные файлы:**
- `output/atoms.graphml` — граф уровня атомов с водородными, ковалентными и зарядовыми связями
- `output/residues.graphml` — граф взаимодействий уровня остатков


### Сравнительный анализ


Сравнение взаимодействий между двумя структурами (например, дикий тип vs. мутант):


```bash
python hydrogen_bonds_analyzer.py wt_structure.pdb mut_structure.pdb
```


**Выходные файлы:**
- Отдельные файлы графов для каждой структуры
- `output/hydrogenbonds_comparison.csv` — детальный отчёт сравнения водородных связей
- `output/charge_interactions_comparison.csv` — детальный отчёт сравнения электростатических взаимодействий


## Форматы входных файлов


### PDB-файл
Стандартный формат Protein Data Bank с атомными координатами. Пример:
```
ATOM      1    N ALA A   1       5.073  -1.453  -9.394  0.00  0.00      A
ATOM      2  HT1 ALA A   1       4.989  -0.460  -9.226  0.00  0.00      A
ATOM      3  HT2 ALA A   1       4.210  -1.901  -9.131  0.00  0.00      A
```


### Топология JSON
Файл определения остатков с информацией о зарядах, донорах, акцепторах и ковалентных связях:
```json
{
  "ALA": {
    "name": "ALA",
    "type": "RESI",
    "charge": 0.0,
    "atoms": {
      "N": {"name": "N", "type": "N", "charge": -0.41},
      "CA": {"name": "CA", "type": "CT", "charge": 0.07},
      "C": {"name": "C", "type": "C", "charge": 0.51},
      "O": {"name": "O", "type": "O", "charge": -0.51},
      "CB": {"name": "CB", "type": "CT", "charge": -0.18}
    },
    "bonds": [["N", "CA"], ["CA", "C"], ["C", "O"], ["CA", "CB"]],
    "donors": [["N"]],
    "acceptors": [["O"]]
  },
  "ARG": {
    "name": "ARG",
    "type": "RESI",
    "charge": 1.0,
    "atoms": {
      "NH1": {"name": "NH1", "type": "N2", "charge": -0.80},
      "NH2": {"name": "NH2", "type": "N2", "charge": -0.80}
    },
    "bonds": [["NH1", "CZ"], ["CZ", "NH2"]],
    "donors": [["NH1"], ["NH2"]],
    "acceptors": []
  },
  "ASP": {
    "name": "ASP",
    "type": "RESI",
    "charge": -1.0,
    "atoms": {
      "OD1": {"name": "OD1", "type": "O2", "charge": -0.76},
      "OD2": {"name": "OD2", "type": "O2", "charge": -0.76}
    },
    "bonds": [["OD1", "CG"], ["OD2", "CG"]],
    "donors": [],
    "acceptors": [["OD1"], ["OD2"]]
  }
}
```

## Форматы выходных данных


### Сравнение водородных связей: hydrogenbonds_comparison.csv


| Колонка | Тип | Описание |
|---------|-----|---------|
| `status` | str | Классификация связи: `conserved`, `lost`, `gained` |
| `atom_id_old_donor` | int | ID атома донора в структуре 1 |
| `atom_id_old_acceptor` | int | ID атома акцептора в структуре 1 |
| `atom_id_new_donor` | int | ID атома донора в структуре 2 |
| `atom_id_new_acceptor` | int | ID атома акцептора в структуре 2 |
| `name_old_donor` | str | Имя атома донора в структуре 1 |
| `name_old_acceptor` | str | Имя атома акцептора в структуре 1 |
| `residue_old_donor` | str | Тип остатка донора в структуре 1 |
| `residue_old_acceptor` | str | Тип остатка акцептора в структуре 1 |
| `residue_id_old_donor` | int | Индекс остатка донора в структуре 1 |
| `residue_id_old_acceptor` | int | Индекс остатка акцептора в структуре 1 |
| `name_new_donor` | str | Имя атома донора в структуре 2 |
| `name_new_acceptor` | str | Имя атома акцептора в структуре 2 |
| `residue_new_donor` | str | Тип остатка донора в структуре 2 |
| `residue_new_acceptor` | str | Тип остатка акцептора в структуре 2 |
| `residue_id_new_donor` | int | Индекс остатка донора в структуре 2 |
| `residue_id_new_acceptor` | int | Индекс остатка акцептора в структуре 2 |


### Сравнение электростатических взаимодействий: charge_interactions_comparison.csv


| Колонка | Тип | Описание |
|---------|-----|---------|
| `status` | str | Классификация взаимодействия: `conserved`, `lost`, `gained` |
| `atom_id_old_1` | int | ID первого заряженного атома в структуре 1 |
| `atom_id_old_2` | int | ID второго заряженного атома в структуре 1 |
| `atom_id_new_1` | int | ID первого заряженного атома в структуре 2 |
| `atom_id_new_2` | int | ID второго заряженного атома в структуре 2 |
| `name_old_1` | str | Имя первого атома в структуре 1 |
| `name_old_2` | str | Имя второго атома в структуре 1 |
| `residue_old_1` | str | Тип остатка первого атома в структуре 1 |
| `residue_old_2` | str | Тип остатка второго атома в структуре 1 |
| `residue_id_old_1` | int | Индекс остатка первого атома в структуре 1 |
| `residue_id_old_2` | int | Индекс остатка второго атома в структуре 1 |
| `charge1_old` | float | Заряд первого атома в структуре 1 |
| `charge2_old` | float | Заряд второго атома в структуре 1 |
| `distance_old` | float | Расстояние между атомами в структуре 1 (Å) |
| `name_new_1` | str | Имя первого атома в структуре 2 |
| `name_new_2` | str | Имя второго атома в структуре 2 |
| `residue_new_1` | str | Тип остатка первого атома в структуре 2 |
| `residue_new_2` | str | Тип остатка второго атома в структуре 2 |
| `residue_id_new_1` | int | Индекс остатка первого атома в структуре 2 |
| `residue_id_new_2` | int | Индекс остатка второго атома в структуре 2 |
| `charge1_new` | float | Заряд первого атома в структуре 2 |
| `charge2_new` | float | Заряд второго атома в структуре 2 |
| `distance_new` | float | Расстояние между атомами в структуре 2 (Å) |


### Значения статуса


- **conserved**: взаимодействие существует в обеих структурах (функционально стабильно)
- **lost**: взаимодействие присутствует только в структуре 1 (мутация нарушает взаимодействие)
- **gained**: взаимодействие присутствует только в структуре 2 (мутация создаёт новое взаимодействие)


| Status    | Структура 1 | Структура 2 | Оба атома маппированы |
| --------- | ----------- | ----------- | ---------------------- |
| conserved | ✓           | ✓           | ✓                      |
| lost      | ✓           | ✗           | ✓                      |
| gained    | ✗           | ✓           | ✓                      |


### Выходной формат графов (GraphML)


Оба файла `atoms.graphml` и `residues.graphml` можно импортировать в инструменты визуализации графов:


**Атрибуты узлов:**
- `role`: роль атома (`donor`, `acceptor`, `none`)


**Атрибуты рёбер:**
- `kind`: тип связи (`covalent`, `hydrogen`, или `charge`)
- Для заряженых взаимодействий также: `charge1`, `charge2`, `distance`


## Справочник API


### Основные классы для водородных связей


#### `HydrogenBondConfig`
Dataclass конфигурации параметров детекции водородных связей.


```python
@dataclass
class HydrogenBondConfig:
    distance_threshold: float = 3.5  # Ångströms
    angle_threshold: float = 150.0   # Degrees
```


#### `HydrogenBondDetector`
Обнаруживает водородные связи, используя пространственные и геометрические критерии.


```python
detector = HydrogenBondDetector(atoms_df, G_atoms, G_residues, config)
G_atoms, G_residues = detector.detect()
```


#### `HydrogenBondComparator`
Сравнивает паттерны водородных связей между двумя структурами.


```python
comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
comparison_df = comparator.compare(exclude_water=True)
```


### Основные классы для электростатических взаимодействий


#### `ChargeInteractionConfig`
Dataclass конфигурации параметров детекции электростатических взаимодействий.


```python
@dataclass
class ChargeInteractionConfig:
    dist_threshold: float = 4.5      # Ångströms
    charge_threshold: float = 0.3    # Elementary charges (e)
```


#### `ChargeInteractionDetector`
Обнаруживает электростатические взаимодействия на основе информации о зарядах из топологии.


```python
detector = ChargeInteractionDetector(atoms_df, G_atoms, topology, config)
G_atoms = detector.detect()
```

**Особенности:**
- Извлекает информацию о зарядах из JSON топологии
- Обнаруживает заряженные боковые цепи (исключает остов)
- Детектирует только противоположные заряды (электростатическое привлечение)
- Сохраняет заряды и расстояния в атрибутах рёбер графа


#### `ChargeInteractionComparator`
Сравнивает электростатические взаимодействия между двумя структурами.


```python
comparator = ChargeInteractionComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
comparison_df = comparator.compare(exclude_water=True)
```

**Выход:**
- DataFrame с 22 колонками, включающими информацию о зарядах, расстояниях и статусе


### Общие классы


#### `TopologyLoader`
Загружает и предварительно обрабатывает файлы топологии.


```python
topology = TopologyLoader.load('topology.json')

# Получить заряженные боковые цепи остатка
charged_sidechain = TopologyLoader.get_charged_sidechain_atoms(topology, 'ARG', charge_threshold=0.3)
```


#### `PDBLoader`
Загружает PDB-файлы в pandas DataFrame.


```python
atoms_df = PDBLoader.load('structure.pdb')
```


#### `GraphBuilder`
Строит графы атомов и остатков из топологии и координат.


```python
builder = GraphBuilder(atoms_df, topology)
G_atoms, G_residues = builder.build()
```


#### `PipelineOrchestrator`
Высокоуровневый интерфейс для управления полным рабочим процессом.


```python
config_hb = HydrogenBondConfig(distance_threshold=3.5, angle_threshold=150.0)
config_charge = ChargeInteractionConfig(dist_threshold=4.5, charge_threshold=0.3)

orchestrator = PipelineOrchestrator(config_hb, config_charge)

# Анализ одной структуры
G_atoms, G_residues, atoms_df = orchestrator.analyze_structure(
    pdb_path, topology_path, output_dir
)

# Сравнение двух структур
hb_df, charge_df = orchestrator.compare_structures(
    pdb_1, pdb_2, topology_path, output_dir
)
```


## Примеры использования


### Пример 1: Анализ структуры мутанта


```bash
python hydrogen_bonds_analyzer.py ./data/examples/protein_mut.pdb
```


Проверьте `output/atoms.graphml` для сети всех взаимодействий (водородные связи, ковалентные связи, солевые мостики).


### Пример 2: Сравнение мутанта с диким типом


```bash
python hydrogen_bonds_analyzer.py ./data/examples/protein_wt.pdb ./data/examples/protein_mut.pdb
```


Анализируем результаты:
```python
import pandas as pd

# Водородные связи
hb_df = pd.read_csv('output/hydrogenbonds_comparison.csv')
print("Водородные связи:")
print(hb_df['status'].value_counts())
print(f"\nПотеряно водородных связей: {len(hb_df[hb_df['status'] == 'lost'])}")
print(f"Приобретено водородных связей: {len(hb_df[hb_df['status'] == 'gained'])}")

# Электростатические взаимодействия
charge_df = pd.read_csv('output/charge_interactions_comparison.csv')
print("\n\nЭлектростатические взаимодействия:")
print(charge_df['status'].value_counts())
print(f"\nПотеряно солевых мостиков: {len(charge_df[charge_df['status'] == 'lost'])}")
print(f"Приобретено солевых мостиков: {len(charge_df[charge_df['status'] == 'gained'])}")

# Детальная информация о потерянных взаимодействиях
lost_hb = hb_df[hb_df['status'] == 'lost']
for idx, row in lost_hb.iterrows():
    print(f"\nПотеряна: {row['residue_old_donor']}{row['residue_id_old_donor']} "
          f"({row['name_old_donor']}) -> {row['residue_old_acceptor']}{row['residue_id_old_acceptor']} "
          f"({row['name_old_acceptor']})")
```


### Пример 3: Программный анализ


```python
from hydrogen_bonds_analyzer import (
    TopologyLoader, PDBLoader, GraphBuilder,
    HydrogenBondDetector, HydrogenBondComparator,
    ChargeInteractionDetector, ChargeInteractionComparator,
    HydrogenBondConfig, ChargeInteractionConfig
)

# Загрузить топологию
topology = TopologyLoader.load('./data/topology.json')

# Загрузить две структуры
atoms_df_1 = PDBLoader.load('./data/wt.pdb')
atoms_df_2 = PDBLoader.load('./data/mut.pdb')

# Построить графы
builder_1 = GraphBuilder(atoms_df_1, topology)
G1_atoms, G1_residues = builder_1.build()

builder_2 = GraphBuilder(atoms_df_2, topology)
G2_atoms, G2_residues = builder_2.build()

# Конфигурация
hb_config = HydrogenBondConfig(distance_threshold=3.5, angle_threshold=150.0)
charge_config = ChargeInteractionConfig(dist_threshold=4.5, charge_threshold=0.3)

# Обнаружить водородные связи
hb_detector_1 = HydrogenBondDetector(atoms_df_1, G1_atoms, G1_residues, hb_config)
G1_atoms, G1_residues = hb_detector_1.detect()

hb_detector_2 = HydrogenBondDetector(atoms_df_2, G2_atoms, G2_residues, hb_config)
G2_atoms, G2_residues = hb_detector_2.detect()

# Обнаружить электростатические взаимодействия
charge_detector_1 = ChargeInteractionDetector(atoms_df_1, G1_atoms, topology, charge_config)
G1_atoms = charge_detector_1.detect()

charge_detector_2 = ChargeInteractionDetector(atoms_df_2, G2_atoms, topology, charge_config)
G2_atoms = charge_detector_2.detect()

# Сравнить водородные связи
hb_comparator = HydrogenBondComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
hb_comparison = hb_comparator.compare(exclude_water=True)

# Сравнить электростатические взаимодействия
charge_comparator = ChargeInteractionComparator(atoms_df_1, atoms_df_2, G1_atoms, G2_atoms)
charge_comparison = charge_comparator.compare(exclude_water=True)

# Сохранить результаты
hb_comparison.to_csv('hydrogenbonds_comparison.csv', index=False)
charge_comparison.to_csv('charge_interactions_comparison.csv', index=False)
```


## Рассмотрения производительности


- **Оптимизация cKDTree**: поиск соседей в пространстве использует `scipy.spatial.cKDTree` с сложностью O(log n)
- **Использование памяти**: хранилище графа масштабируется с количеством атомов и связей (типично <100MB для структур <10k атомов)
- **Время выполнения**:
  - Водородные связи: ~100-300ms на структуру
  - Электростатические взаимодействия: ~50-150ms на структуру
  - Сравнение: ~100-200ms для пары структур
  - **Итого:** ~500-1000ms для полного анализа типичного белка


## Решение проблем


### Проблема: `FileNotFoundError: Topology file not found`
**Решение:** проверьте путь `TOPOLOGY_FILE` в `.env` или передайте корректный путь в функции.


### Проблема: `KeyError: Unknown residue type`
**Решение:** убедитесь, что все типы аминокислот в PDB определены в JSON топологии.


### Проблема: водородные связи не обнаружены
**Возможные причины:**
1. Пороговое расстояние слишком строгое (<3.0 Å)
2. Угловой порог слишком строгий (>160°)
3. В топологии отсутствуют определения доноров/акцепторов
4. PDB-файл не содержит явных водородов

**Решение:** ослабьте пороги или используйте предварительно обработанный PDB с явными водородами.


### Проблема: электростатические взаимодействия не обнаружены
**Возможные причины:**
1. Пороговое расстояние слишком строгое (<4.0 Å)
2. Пороговый заряд слишком высокий (>0.5 e)
3. В топологии отсутствуют значения заряда
4. В структуре нет заряженных остатков (ARG, LYS, ASP, GLU)

**Решение:** проверьте значения зарядов в JSON топологии или ослабьте пороги.


### Проблема: ошибка памяти на больших структурах
**Решение:** обрабатывайте структуры по частям или увеличьте доступную системную память.


## Логирование


Анализатор предоставляет детальное логирование на уровнях INFO и ERROR:


```python
import logging
logging.basicConfig(level=logging.DEBUG)  # Для подробного вывода
```


**Пример логирования:**
```
2025-11-25 17:10:45,123 - root - INFO - Loading topology from ./data/topology.json
2025-11-25 17:10:45,456 - root - INFO - Loaded topology for 20 residue types
2025-11-25 17:10:46,789 - root - INFO - Loading PDB from structure.pdb
2025-11-25 17:10:47,012 - root - INFO - Loaded 2847 atoms
2025-11-25 17:10:48,234 - root - INFO - Building molecular graphs...
2025-11-25 17:10:49,456 - root - INFO - Graphs built: 2847 atoms, 128 residues
2025-11-25 17:10:50,678 - root - INFO - Detecting hydrogen bonds...
2025-11-25 17:10:51,234 - root - INFO - Found 142 donors and 156 acceptors
2025-11-25 17:10:52,890 - root - INFO - Detected 87 hydrogen bonds
2025-11-25 17:10:53,456 - root - INFO - Detecting charge-charge interactions...
2025-11-25 17:10:54,123 - root - INFO - Found 12 charged atoms in structure
2025-11-25 17:10:54,789 - root - INFO - Detected 5 charge-charge interactions
```


## Тестирование


### Запуск всех тестов
```bash
pytest tests/ -v
```


### Запуск тестов для конкретного компонента
```bash
# Тесты топологии
pytest tests/test_topology_loader.py -v

# Тесты PDB-загрузчика
pytest tests/test_pdb_loader.py -v

# Тесты детекции водородных связей
pytest tests/test_hbond_detector.py -v

# Тесты сравнения водородных связей
pytest tests/test_comparator.py -v

# Тесты электростатических взаимодействий (новые)
pytest tests/test_charge_interactions.py -v
```


### Отчёт о покрытии тестами
```bash
pytest tests/ --cov=hydrogen_bonds_analyzer --cov-report=html
open htmlcov/index.html
```


### Статистика тестов
- **Всего тестов:** 74
  - TopologyLoader: 8 тестов
  - PDBLoader: 12 тестов
  - HydrogenBondDetector: 18 тестов
  - HydrogenBondComparator: 12 тестов
  - ChargeInteractionDetector: 8 тестов
  - ChargeInteractionComparator: 13 тестов
  - Integration: 1 тест
- **Время выполнения:** ~5.8 секунд
- **Покрытие:** >90%


## Разработка и контриб


### Структура проекта
```
hydrogen-bond-analyzer/
├── hydrogen_bonds_analyzer.py       # Основной модуль
├── tests/                           # Юнит-тесты
│   ├── test_topology_loader.py
│   ├── test_pdb_loader.py
│   ├── test_hbond_detector.py
│   ├── test_comparator.py
│   └── test_charge_interactions.py
├── data/
│   ├── raw/
│   │   └── topology_complete.json   # Определения топологии с зарядами
│   └── examples/
│       ├── protein_wt.pdb
│       └── protein_mut.pdb
├── notebooks/                       # Директория разработки методов
├── output/                          # Директория результатов
├── .env                             # Переменные окружения
├── requirements.txt                 # Зависимости
└── README.md                        # Этот файл
```


### Установка для разработки
```bash
git clone <repository-url>
cd hydrogen-bond-analyzer
pytest tests/ -v
```


## История изменений


### v2.0.0 (2025-11-25)
- Добавлена детекция электростатических взаимодействий (солевых мостиков)
- Параметры заряда теперь извлекаются из JSON топологии
- Новый класс `ChargeInteractionDetector` для обнаружения взаимодействий
- Новый класс `ChargeInteractionComparator` для сравнения электростатических взаимодействий
- Расширенный класс `ChargeInteractionConfig` для настройки пороговых значений
- 24 новых юнит-теста для покрытия функционала электростатических взаимодействий
- Обновленный выходной формат для включения информации о зарядах и расстояниях
- Документация на русском языке


### v1.0.0 (2025-11-03)
- Первый релиз
- Анализ одной и сравнение двух структур
- Поддержка экспорта GraphML
- Комплексное логирование и обработка ошибок
- Поддержка Jupyter notebooks
