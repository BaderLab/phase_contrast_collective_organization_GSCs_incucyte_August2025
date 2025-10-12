------------------------------------------------------------------------------------------------------------------------
# Glioblastoma stem cells show transcriptionally correlated spatial organization
------------------------------------------------------------------------------------------------------------------------

##  Contributors :

This project was led by researchers at the **Bader Lab** (University of Toronto) focused on spatial and single-cell systems biology.

AUTHORSHIP :
*** Shamini Ayyadhury ***, Patty Sachamitr, Michelle M. Kushida, Nicole I Park, Fiona J. Coutinho, Owen Whitley, Panagiotis Prinos, Cheryl H. Arrowsmith, Peter B. Dirks, Trevor J. Pugh, Gary D. Bader

----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------

### Overview

1. Glioblastoma stem cells (GSCs) are an important GBM model system. 
2. In culture, these cells form spatial structures that share morphological aspects with their source tumors. 
3. We collected 17,000 phase-contrast images of 15 patient-derived GSC lines growing to confluence. 
4. We find that GSCs grow in characteristic multicellular patterns depending on their transcriptional state. 
5. Interpretable computer vision algorithms identified specific image features that predict transcriptional state across multiple cell confluency levels. 
5. This relationship will be useful in developing GSC screens where image features can be used to identify how GSC biology changes in response to perturbations simply by imaging cultured cells on plates.

----------------------------------------------------------------------------------------------------------------


## Repository Structure

```
├── CellProfilerPipelines/                  # Image segmentation & feature extraction pipelines
├── code/                                   # Python scripts for analysis and preprocessing
├── ilastik_training_scripts/               # Scripts for training Ilastik models
├── ilastik_training_scripts_validation/    # Scripts for validating trained Ilastik models
├── renv.lock                               # R environment snapshot
└── README.md                               # This file
```

## Getting Started

1. Install dependencies:
    - CellProfiler (https://cellprofiler.org/)
    - Ilastik (https://www.ilastik.org/)
    - RStudio(2024.12.0+467) and R(R version 4.4.2 (2024-10-31)) was used in for the analysis of the code in this paper
        - R with packages and dependecies are defined in `renv.lock`

2. Data Input:
    - GSC cells seeded in 384 well plates were imaged using the incucyte imager(https://www.sartorius.com/en/products/live-cell-analysis/incucyte) across different samples at multiple time-points
    - Images were saved as jpeg formats but .tiff and.png is accepted

3. Run Pipelines:
    - BINARY MASK TRAINING : Train Ilastik models using scripts in `ilastik_training_scripts/`. This will allow you to differentiate between foreground (cell-regions) and background (plate/accellular regions)
        * Collect all binary masks together with the original phase-contrast images and save them into folders. 
        * Save the ilastik training file (.ilp) so that you can generate further images headless without going through the training process
    - FEATURE EXTRACTION : Apply CellProfiler pipelines from `CellProfilerPipelines/`. Images were loaded onto Cellprofiler software, together with the ilastik generated masks. The background masks are used to ensure taht feature extraction is only retained on the region of the images with cells. 
    - CSV FILE : Cellprofiler feature extractions can be saved as a csv file. Values from the output are used for downstream analysis of features using R scripts in code/
    
4. Validation pipelines :
    - For the validation cohort, the above steps 1 to 3 are repeated. 


## Goals

1. To determine how GSC populations `self-organize` over time.
2. Correlate image features across samples and confluencey groups with bulk gene expression generated from matched samples. 



## License
    - MIT LICENSE
