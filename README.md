# Brain Structural-Functional Coupling in Trigeminal Neuralgia

This repository contains the analysis code and sample data for the study investigating the macroscopic structural-functional (S-F) coupling alterations and their microscopic transcriptomic underpinnings in patients with Trigeminal Neuralgia (TN) before and following surgical treatment.

## 🧠 Overview
In this study, we utilized a graph harmonic model to quantify S-F coupling from multimodal MRI (dMRI and rs-fMRI). Furthermore, we bridged the macroscopic neuroimaging phenotypes (treatment-induced recovery) with microscopic gene expression profiles using Partial Least Squares (PLS) regression and rigorous spatial permutation testing (spin-tests).

## 📂 Repository Structure
* `/01_SF_Coupling`: Scripts for constructing structural/functional networks and calculating the S-F coupling metric based on the graph harmonic model.
* `/02_Transcriptomic_PLS`: Code for integrating AHBA gene expression data, performing PLS regression, and conducting spatial permutation tests (spin-tests).
* `/Data_Sample`: Anonymized, region-level derived data (Schaefer 400 parcellation) required to run the demonstration scripts.
* `[Update Pending] /Statistical_Analysis`: *The scripts for basic statistical analyses (e.g., longitudinal comparisons, FDR correction) are currently being organized and cleaned for better readability, and will be updated to this repository soon.*

## 🛠️ Dependencies
To run the scripts in this repository, the following software and toolboxes are recommended:
* MATLAB (R2021a or later)
* [PANDA Toolbox](https://www.nitrc.org/projects/panda/) (for dMRI preprocessing)
* ENIGMA Toolbox (for spatial permutation/spin-tests)

## ⚠️ Data Availability
Due to patient privacy and ethical restrictions, the raw multimodal MRI datasets and identifiable clinical information are not publicly available here. We provide anonymized, parcellated group-level summary maps (e.g., unthresholded t-maps) in the `/Data_Sample` folder to ensure code reproducibility. 

## ✉️ Contact
For any questions or code-related issues, please open an issue in this repository or contact the corresponding author.