# 3d-scere

## What is 3d-scere ?

3D-Scere is an open-source tool for interactive visualization and exploration. Source code is available here. The tool is also freely usable online at https://3d-scere.ijm.fr/. It allows the visualization of any list of genes in the context of the 3D model of S. cerevisiae genome. Further information can easily be added, like functional annotations (GO terms) or gene expression measurements. Qualitative or quantitative functional properties are highlighted in the large-scale 3D context of the genome with only a few mouse clicks.

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
```

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

## Test the dashboard with example data

Use the files in example data folder.

Example file for the "GO term projection" and "3D distances histogram and network" tabs:
gene_list_example_UPC2_38_targets.csv : Targets list of the UPC2 trancription factor, extracted from the [supplementary data](https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-020-74043-7/MediaObjects/41598_2020_74043_MOESM1_ESM.zip) of Monteiro et al., Assessing regulatory features of the current transcriptional network of Saccharomyces cerevisiae. Sci Rep. 2020 Dec;10(1):17744.


Example file for the "Quantitative variable projection" tab:
quantitative_variables_example.csv: microarray data (heatshock 25°C to 37°C) from 2010.Gasch00_HS25-37_formated.flt.knn.avg.pcl downloaded from [SGD archives](http://sgd-archive.yeastgenome.org/expression/microarray/Gasch_2000_PMID_11102521/).

