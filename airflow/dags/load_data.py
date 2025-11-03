
from datetime import datetime
from airflow import DAG
from airflow.operators.python import PythonOperator
import os

# Ensure src on PYTHONPATH
import sys
sys.path.insert(0, os.path.join(os.getenv("AIRFLOW_HOME"), ".."))

from src.config import get_source_dir
from src.data_loader import load_atoms, load_topology

default_args = {
    "owner": "airflow",
    "start_date": datetime(2025, 9, 24),
    "depends_on_past": False,
}

def task_load_data(**kwargs):
    src = get_source_dir()
    pdb_path = kwargs["dag_run"].conf.get("pdb_path")
    topo_path = os.path.join(src, "data", "raw", "topology_complete.json")
    atoms_df = load_atoms(pdb_path)
    topo = load_topology(topo_path)
    # Save intermediate artifacts
    atoms_df.to_parquet(os.path.join(src, "data", "processed", "atoms.parquet"), index=False)
    with open(os.path.join(src, "data", "processed", "topology.json"), "w") as f:
        import json
        json.dump(topo, f)

with DAG(
    dag_id="load_data",
    default_args=default_args,
    schedule_interval=None,
    catchup=False,
    tags=["protein", "dvc"],
) as dag:
    load_data = PythonOperator(
        task_id="load_and_process_data",
        python_callable=task_load_data,
        provide_context=True,
    )
