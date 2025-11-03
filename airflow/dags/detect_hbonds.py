
from datetime import datetime
from airflow import DAG
from airflow.operators.python import PythonOperator
import os

import sys
sys.path.insert(0, os.path.join(os.getenv("AIRFLOW_HOME"), ".."))

import pandas as pd
import networkx as nx
from src.config import get_source_dir
from src.io_utils import write_graphml
from src.hydrogen_bonds import detect_hydrogen_bonds

default_args = {
    "owner": "airflow",
    "start_date": datetime(2025, 9, 24),
    "depends_on_past": False,
}

def task_detect_hbonds(**kwargs):
    src = get_source_dir()
    atoms_path = os.path.join(src, "data", "processed", "atoms.parquet")
    atoms_g_path = os.path.join(src, "data", "graphs", "atoms.graphml")
    res_g_path = os.path.join(src, "data", "graphs", "residues.graphml")

    atoms_df = pd.read_parquet(atoms_path)
    G_atoms = nx.read_graphml(atoms_g_path)
    G_residues = nx.read_graphml(res_g_path)

    G_atoms_updated, G_res_updated = detect_hydrogen_bonds(atoms_df, G_atoms, G_residues)

    write_graphml(G_atoms_updated, os.path.join(src, "data", "graphs", "atoms_hbonds.graphml"))
    write_graphml(G_res_updated, os.path.join(src, "data", "graphs", "residues_hbonds.graphml"))

with DAG(
    dag_id="detect_hbonds",
    default_args=default_args,
    schedule_interval=None,
    catchup=False,
    tags=["protein", "dvc"],
) as dag:
    detect_hbonds = PythonOperator(
        task_id="detect_hydrogen_bonds",
        python_callable=task_detect_hbonds,
        provide_context=True,
    )
