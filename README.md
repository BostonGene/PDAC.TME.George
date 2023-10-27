[![N|Solid](https://bostongene.com/wp-content/uploads/2022/01/bg-logo.svg)](https://bostongene.com/)

# PDAC Tumor Microenvironment Classification

This repository contains the code and resources for classifying the tumor microenvironment (TME) of pancreatic ductal adenocarcinoma (PDAC) based on functional gene signatures (FGEs).

## Overview
PDAC presents a complex tumor environment, which has historically posed challenges in the development of reliable predictive biomarkers for targeted therapies and immunomodulation. To address this, we've implemented a classification approach based on transcriptomic profiling of the TME.

Developed classification is described in [Paper's DOI link]

## Features
This repository contains:
- Data
-- Meta-Cohort: annotation & calculated FGEs
-- ICGC PACA-CA (test cohort): annotation & expressions
Code presented in "MFP_TME_classificataion.ipynb" include:
- Data Processing: Scripts to preprocess and clean the PDAC datasets.
-- Calculation of FGEs
- Classification: Algorithms to classify TME into four distinct subtypes:
-- Immune enriched (IE)
-- Immune enriched, fibrotic (IE/F)
-- Fibrotic (F)
-- Immune depleted (D)

## Getting Started
### Prerequisites
Python 3.x
Required Python libraries: [List down libraries used, e.g., pandas, numpy, scikit-learn, etc.]
### Installation
Clone the repository:
git clone [repository-link]

Navigate to the repository directory:
bash
Copy code
cd [repository-name]
Install required dependencies:
bash
Copy code
pip install -r requirements.txt

## License
[]

## Acknowledgments
[]

© 2023 BostonGene Corporation.