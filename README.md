# 3d-scere

## Setup your environment

Clone the repository:
```
git clone https://github.com/data-fun/3d-scere.git
```

Move to the new directory:
```
cd 3d-scere
```

Create a [conda](https://docs.conda.io/en/latest/miniconda.html) environment:
```
conda env create -f dashboard_conda_env.yml
```

Load the `dashboard` conda environment:
```
conda activate dashboard
```

## Required static data

Download the SQL database (~340 Mb):
```
wget -O static/SCERE.db https://zenodo.org/record/5526011/files/SCERE.db

Download the distance matrix (~550 Mb):
```
wget -O static/3D_distances.parquet.gzip https://zenodo.org/record/5526011/files/3D_distances.parquet.gzip
```

## Run the dashboard

```
make run
```

then open your web browser on <http://127.0.0.1:8050/>


## Run the dashboard on a server

```
make run-gunicorn
```

then open your web browser on <http://127.0.0.1:8000/>
