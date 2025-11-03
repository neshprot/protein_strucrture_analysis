# airflow/dags/protein_pipeline_master.py

from datetime import datetime
import os
import glob

from airflow import DAG
from airflow.operators.trigger_dagrun import TriggerDagRunOperator

# Ensure src package is on path if your load_data, etc. DAGs import from src
import sys
sys.path.insert(0, os.path.join(os.getenv("AIRFLOW_HOME"), ".."))

from src.config import get_source_dir

default_args = {
    "owner": "airflow",
    "start_date": datetime(2025, 9, 24),
    "depends_on_past": False,
}

with DAG(
    dag_id="protein_pipeline_master",
    default_args=default_args,
    schedule_interval=None,
    catchup=False,
    tags=["protein", "dvc"],
) as dag:

    src = get_source_dir()
    pdb_dir = os.path.join(src, "data", "pdbs")
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

    for pdb in pdb_files:
        base = os.path.splitext(os.path.basename(pdb))[0]

        # 1. Trigger load_data, passing pdb_path
        t1 = TriggerDagRunOperator(
            task_id=f"load_data_{base}",
            trigger_dag_id="load_data",
            conf={"pdb_path": pdb},
        )

        # 2. Trigger build_graphs, downstream of load_data
        t2 = TriggerDagRunOperator(
            task_id=f"build_graphs_{base}",
            trigger_dag_id="build_graphs",
            conf={"pdb_path": pdb},
        )

        # 3. Trigger detect_hbonds, downstream of build_graphs
        t3 = TriggerDagRunOperator(
            task_id=f"detect_hbonds_{base}",
            trigger_dag_id="detect_hbonds",
            conf={"pdb_path": pdb},
        )

        t1 >> t2 >> t3
