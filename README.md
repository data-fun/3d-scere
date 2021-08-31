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

Download the SQL database:
```
wget -P ./static https://zenodo.org/blablablab/SCERE.db
```

Download the distances matrix:
```
wget -P ./static https://zenodo.org/blablablab/adjacency_matrix_V4.parquet.gzip
```

## Run the dashboard

```
python app.py
```

then open your web browser on <http://127.0.0.1:8050/>
