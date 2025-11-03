
from datetime import datetime
from airflow import DAG
from airflow.operators.python import PythonOperator
import os

import sys
sys.path.insert(0, os.path.join(os.getenv("AIRFLOW_HOME"), ".."))

import pandas as pd
import json
from src.config import get_source_dir
from src.graph_builder import build_graphs
from src.io_utils import write_graphml

default_args = {
    "owner": "airflow",
    "start_date": datetime(2025, 9, 24),
    "depends_on_past": False,
}

def task_build_graphs(**kwargs):
    src = get_source_dir()
    atoms_path = os.path.join(src, "data", "processed", "atoms.parquet")
    topo_path = os.path.join(src, "data", "processed", "topology.json")

    atoms_df = pd.read_parquet(atoms_path)
    with open(topo_path, "r") as f:
        topology = json.load(f)

    G_atoms, G_residues = build_graphs(atoms_df, topology)

    out_dir = os.path.join(src, "data", "graphs")
    os.makedirs(out_dir, exist_ok=True)
    write_graphml(G_atoms, os.path.join(out_dir, "atoms.graphml"))
    write_graphml(G_residues, os.path.join(out_dir, "residues.graphml"))

with DAG(
    dag_id="build_graphs",
    default_args=default_args,
    schedule_interval=None,
    catchup=False,
    tags=["protein", "dvc"],
) as dag:
    build_graphs_task = PythonOperator(
        task_id="build_graphs",
        python_callable=task_build_graphs,
        provide_context=True,
    )
