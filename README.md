# Анализатор водородных связей

## Обзор

**Анализатор водородных связей** позволяет:
- Обнаруживать водородные связи в белковых структурах, используя геометрические критерии (расстояние и угловые пороги)
- Строить графы взаимодействий на уровне атомов и остатков
- Сравнивать паттерны водородных связей между двумя белковыми структурами (например, дикий тип vs. мутант)
- Генерировать детальные отчёты о сохранённых, потерянных и приобретённых взаимодействиях

## Ключевые функции

- **Модульная архитектура**: разделение ответственности с отдельными классами для каждого компонента
- **Production-ready код**: комплексное логирование, обработка ошибок и type hints
- **Оптимизированная производительность**: использование cKDTree для эффективного пространственного поиска
- **Гибкая конфигурация**: настраиваемые пороги расстояния и угла через `HydrogenBondConfig`
- **Сравнительный анализ**: детальное отслеживание изменений водородных связей между структурами
- **Графовое представление**: графы NetworkX для взаимодействий на уровне атомов и остатков
- **Пакетная обработка**: поддержка анализа одной структуры или сравнения двух

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
```

## Конфигурация

### Переменные окружения
Создайте файл `.env` в корне проекта:

```env
# Topology file with residue definitions
TOPOLOGY_FILE=./../data/raw/topology_complete.json

# Output directory for results
OUTPUT_DIR=./output

# Detection thresholds (optional, defaults provided)
HB_DISTANCE_THRESHOLD=3.5
HB_ANGLE_THRESHOLD=150.0
```

### Параметры детекции водородных связей

Два геометрических критерия определяют водородную связь:

1. **Порог расстояния** (по умолчанию: 3.5 Å)
- Максимальное расстояние между атомами донора и акцептора

2. **Угловой порог** (по умолчанию: 150°)
- Минимальный угол при атоме донора (реф-донор-акцептор)
- Значения ближе к 180° указывают на лучшее выравнивание

## Использование

### Анализ одной структуры

Обнаружение всех водородных связей в белковой структуре:

```bash
python hydrogen_bonds_analyzer.py structure.pdb
```

**Выходные файлы:**
- `output/atoms.graphml` — граф уровня атомов с водородными и ковалентными связями
- `output/residues.graphml` — граф взаимодействий уровня остатков

### Сравнительный анализ

Сравнение водородных связей между двумя структурами (например, дикий тип vs. мутант):

```bash
python hydrogen_bonds_analyzer.py wt_structure.pdb mut_structure.pdb
```

**Выходные файлы:**
- Отдельные файлы графов для каждой структуры
- `output/hydrogen_bonds_comparison.csv` — детальный отчёт сравнения

## Форматы входных файлов

### PDB-файл
Стандартный формат Protein Data Bank с атомными координатами. Пример:
```
ATOM      1    N ALA A   1       5.073  -1.453  -9.394  0.00  0.00      A
ATOM      2  HT1 ALA A   1       4.989  -0.460  -9.226  0.00  0.00      A
ATOM      3  HT2 ALA A   1       4.210  -1.901  -9.131  0.00  0.00      A
```

### Топология JSON
Файл определения остатков, указывающий доноры, акцепторы и ковалентные связи:
```json
{
"ALA": {
"donors": [["N"]],
"acceptors": [["O"]],
"bonds": [["N", "CA"], ["CA", "C"], ["C", "O"]]
},
"GLY": {
"donors": [["N"]],
"acceptors": [["O"]],
"bonds": [["N", "CA"], ["CA", "C"], ["C", "O"]]
}
}
```

## Формат выходных данных

### Колонки CSV при сравнении

| Колонка | Тип | Описание |
|---------|-----|---------|
| `status` | str | Классификация связи: `conserved`, `lost`, `gained`, `only_wt`, `only_mut` |
| `atom_id_old_donor` | int | ID атома донора в структуре 1 |
| `atom_id_old_acceptor` | int | ID атома акцептора в структуре 1 |
| `atom_id_new_donor` | int | ID атома донора в структуре 2 |
| `atom_id_new_acceptor` | int | ID атома акцептора в структуре 2 |
| `name_old` | str | Имя атома в структуре 1 |
| `name_new` | str | Имя атома в структуре 2 |
| `residue_old` | str | Тип остатка в структуре 1 |
| `residue_new` | str | Тип остатка в структуре 2 |
| `residue_id_old` | int | Индекс остатка в структуре 1 |
| `residue_id_new` | int | Индекс остатка в структуре 2 |

### Значения статуса

- **conserved**: водородная связь существует в обеих структурах (функционально стабильна)
- **lost**: связь присутствует только в структуре 1 (мутация нарушает взаимодействие)
- **gained**: связь присутствует только в структуре 2 (мутация создаёт новое взаимодействие)
- **only_wt**: связь в структуре 1, но атом отсутствует в структуре 2 (удаление/замена атома)
- **only_mut**: связь в структуре 2, но атом отсутствует в структуре 1 (вставка/замена атома)

| Status    | Структура 1 | Структура 2 | Оба атома маппированы  |
| --------- | ----------- | ----------- | ---------------------- |
| conserved | ✓           | ✓           | ✓                      |
| lost      | ✓           | ✗           | ✓                      |
| gained    | ✗           | ✓           | ✓                      |
| only_wt   | ✓           | -           | ✗ (атом не существует) |
| only_mut  | -           | ✓           | ✗ (атом не существует) |

### Выходной формат графов (GraphML)

Оба файла `atoms.graphml` и `residues.graphml` можно импортировать в инструменты визуализации графов:

**Атрибуты узлов:**
- `role`: роль атома (`donor`, `acceptor`, `none`)

**Атрибуты рёбер:**
- `kind`: тип связи (`covalent` или `hydrogen`)

## Справочник API

### Основные классы

#### `HydrogenBondConfig`
Dataclass конфигурации параметров детекции.

```python
@dataclass
class HydrogenBondConfig:
distance_threshold: float = 3.5
angle_threshold: float = 150.0
```

#### `TopologyLoader`
Загружает и предварительно обрабатывает файлы топологии.

```python
topology = TopologyLoader.load('topology.json')
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

#### `PipelineOrchestrator`
Высокоуровневый интерфейс для управления полным рабочим процессом.

```python
orchestrator = PipelineOrchestrator(config)
G_atoms, G_residues = orchestrator.analyze_structure(pdb_path, topology_path, output_dir)
comparison_df = orchestrator.compare_structures(pdb_1, pdb_2, topology_path, output_dir)
```

## Примеры использования

### Пример 1: Анализ структуры дикого типа

```bash
python hydrogen_bonds_analyzer.py ./data/examples/protein_mut.pdb ./data/examples/protein_wt.pdb
```

Проверьте `output/atoms.graphml` для сети водородных связей.

### Пример 2: Сравнение мутанта с диким типом

```bash
python hydrogen_bonds_analyzer.py ./data/examples/protein_mut.pdb ./data/examples/protein_wt.pdb
```

Анализируем `output/hydrogen_bonds_comparison.csv`:
```python
import pandas as pd

df = pd.read_csv('output/hydrogen_bonds_comparison.csv')

# Считаем распределение статусов
print(df['status'].value_counts())

# Просмотрим потерянные взаимодействия
lost = df[df['status'] == 'lost']
print(f"Потеряно взаимодействий: {len(lost)}")
```

## Рассмотрения производительности

- **Оптимизация cKDTree**: поиск соседей в пространстве использует `scipy.spatial.cKDTree` с сложностью O(log n)
- **Использование памяти**: хранилище графа масштабируется с количеством атомов и связей (типично <100MB для структур <10k атомов)
- **Время выполнения**: ~100-500ms на структуру для типичных белков, в зависимости от размера и настроек пороговых значений

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
2025-11-03 17:10:45,123 - __main__ - INFO - Loading topology from ./data/topology.json
2025-11-03 17:10:45,456 - __main__ - INFO - Loaded topology for 20 residue types
2025-11-03 17:10:46,789 - __main__ - INFO - Loading PDB from structure.pdb
2025-11-03 17:10:47,012 - __main__ - INFO - Loaded 2847 atoms
```

## Разработка и контриб

### Структура проекта
```
hydrogen-bond-analyzer/
├── hydrogen_bonds_analyzer.py    # Основной модуль
├── tests/                        # Юнит-тесты
│   ├── test_topology_loader.py
│   ├── test_hbond_detector.py
│   └── test_comparator.py
├── data/
│   ├── raw/
│   │   ├── topology_complete.json    # Определения топологии
│   └── examples/
│       ├── protein_wt.pdb
│       └── protein_mut.pdb
├── notebooks/                    # Директория разработки методов
├── output/                       # Директория результатов
├── .env                          # Переменные окружения
├── requirements.txt              # Зависимости
└── README.md                     # Этот файл
```

### Запуск тестов
```bash
pytest tests/ -v
```

## История изменений

### v1.0.0 (2025-11-03)
- Первый релиз
- Анализ одной и сравнение двух структур
- Поддержка экспорта GraphML
- Комплексное логирование и обработка ошибок
