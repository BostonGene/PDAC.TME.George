[![N|Solid](https://bostongene.com/wp-content/uploads/2022/01/bg-logo.svg)](https://bostongene.com/)

# PDAC Tumor Microenvironment Classification
This repository contains the code and resources for classifying the tumor microenvironment (TME) of pancreatic ductal adenocarcinoma (PDAC) based on functional gene signatures (FGEs) calculated on bulk RNA-Seq or microarray gene expressions data.

## Overview
PDAC presents a complex tumor environment, which has historically posed challenges in the development of reliable predictive biomarkers for targeted therapies and immunomodulation. To address this, we've implemented a classification approach based on transcriptomic profiling of the TME.

Four TME subtypes were developed: 
* **Immune Enriched (IE)**
* **Immune Enriched, Fibrotic (IE/F)**
* **Fibrotic (F)**
* **Immune Depleted (D)**

### Publication
Developed classification is described in details in [DOI link - TODO]  

![plot](/img/abstract.png)

## Features
This repository contains:
+ Data
    * PDAC Meta-Cohort: annotation & calculated FGEs for all samples 
    * ICGC PACA-CA (test cohort): annotation & expressions
+ Genesets 
    * **genesets/PAAD_genesets.gmt** - table of all used genesets
+ Code
    * **utils/** - scripts for data preprocessing, ssGSEA score calculation, median scaling and plots
    * Notebook **MFP_TME_classificataion.ipynb** with a classification of  ICGC PACA-CA samples into four TME subtypes 
    using supervised clustering

## Getting Started
### Prerequisites
+ **Python 3.x** - Python 3.9 is advised
+ Required Python libraries: please see requirements.txt
### Installation
Clone the repository:
```
git clone https://github.com/BostonGene/PDAC.TME.George.git
```
Navigate to the repository directory:
```
cd PDAC.TME.George
```
Install required dependencies:  
```
python3.9 -m venv venv  
source venv/bin/activate  
venv/bin/python3.9 -m pip install --upgrade pip  
pip install --no-deps -r requirements.txt  
jupyter nbextension enable --py widgetsnbextension  
python -m ipykernel install --user --name=pdac_tme_venv 
``` 

## License
Licensed by BostonGene Licence -- for more info please see **LICENSE** file  
Please direct any inquiries concerning usage to askusepermission@bostongene.com.  
© 2023 BostonGene Corporation.
