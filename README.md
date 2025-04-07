# This repository contains the raw data, MATLAB scripts, and supporting files associated with the study:
**Lotfi, M. et al. Measuring Alignment of Structural Proteins in Engineered Tissue Constructs Using Polarized Raman Spectroscopy. PLOS ONE (2025)**

To reproduce the analysis:
1. Open `PRS_alignment_analysis.m` in MATLAB.
2. Load raw spectra for a specific region (e.g., `corner_edge_spectra`).
3. Run the script to preprocess spectra (background subtraction, smoothing, SNV normalization).
4. Project spectra using `master_loading_function.mat` to obtain PC1 scores.
5. Fit PC1 scores to a sine curve using `fit_sine_to_PC1.m` to extract alignment parameters.

For questions, please contact:
**Malisa Sarntinoranont**  
Email: msarnt@ufl.edu  
Department of Mechanical and Aerospace Engineering, University of Florida
