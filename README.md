QSP-mets: Quantitative systems pharmacology model for immunotherapy in metastatic triple negative breast cancer

This repository contains the QSP platform for metastatic TNBC developed by the Popel Systems Biology Laboratory. 
https://popellab.johnshopkins.edu/. 

This QSP model is an extension of the QSPIO-TNBC model (https://github.com/HanwenWang95/QSPIO-TNBC) with metastatic tumors. 

Contents of the repository
- Folder "model" contains different modules to generate the simbiology model
- Folder "@struct" contains some functions required to generate the simbiology model
- Folder "parameters" contains parameter files
- Folder "postprocessing" contains scripts for analyzing simulation results
- Folder "scripts" contains scripts to run simulations
- Folder "scripts-mets" contains new scripts adapted to analyze/visualize simulation results with multiple tumors
- Folder "visualization" contains scripts for visualization

To simulate the QSP model once:
1) use the script immune_oncology_model_TNBC.m in the folder "scripts". 

To run the in-silico clinical trial simulations: 
1) use the script PSA_script_TNBC.m to run simulations in a single core. Set the number of virtual patients in this script.

For setting up multiple job arrays:
1) use the script arch_generatemodel.m to generate/setup the simbiology model. 
Set the total number of virtual patients and number of virtual patients to be simulated per job in this script.
2) use the script arch_simulate.m to simulate the model. This function uses the argument job number.
3) Collect the results from all jobs in a single folder and use the script arch_analyze.m to analyze the results.

Reference: Arulraj, T., Wang, H., Emens, L.A., Santa-Maria, C.A. and Popel, A.S., 2023. 
A transcriptome-informed QSP model of metastatic triple-negative breast cancer identifies predictive biomarkers for PD-1 inhibition. Science Advances, 9(26), p.eadg0289.
