
# Укажите версию Airflow, которую хотите установить
AIRFLOW_VERSION=2.7.3

# Определите вашу версию Python
PYTHON_VERSION="$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')"

# Сформируйте URL для constraint-файла
CONSTRAINT_URL="https://raw.githubusercontent.com/apache/airflow/constraints-${AIRFLOW_VERSION}/constraints-${PYTHON_VERSION}.txt"

# Установите Airflow с помощью uv, используя constraints
uv pip install "apache-airflow==${AIRFLOW_VERSION}" --constraint "${CONSTRAINT_URL}"
