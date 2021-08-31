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

Optionally, install DB Browser for SQLite with `apt`:
```
grep -vE '^#' binder/apt.txt | xargs sudo apt install -y
```
