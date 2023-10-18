# cellstruct-reproducibility
This repository contains the codes for reproducing the single-cell analysis in the paper "cellstruct: Metrics scores to quantify the biological preservation between two embeddings".

There are three parts of analysis:
1. Human liver cell atlas
   * Seurat object can be downloaded from [https://zenodo.org/records/7770308/files/ref.Rds?download=1]
   * The script used to generate most analysis/figures is humanLiver.R
   * downsampling_humanLiver.R, compiling_downsampling.R, and plot_downsamplingFigure.R are used to test the stability of tuneUMAP function by varying dataset size

2. COVID-19 peripheral immune atlas
   * Seurat object can be downloaded from [https://s3.us-west-2.amazonaws.com/atlas.fredhutch.org/data/hutch/covid19/downloads/wilk_2020_processed.HDF5]
   * The scripts used are wilk2020.R (main analysis) and enrichR.R (enrichment analysis of differentially expressed genes)

3. Simulated datasets from scDEED
   * Zenodo: [https://zenodo.org/records/7216361/files/Simulated_Data.zip?download=1] (20 simulated datasets are used for evaluating KNN and KNC metrics; only simulated dataset 1 is used for the detailed comparison between scDEED and cellstruct)
   * The scripts used are run_scDEED.R and edited_scDEED.R (manually edit scDEED to disable hyperparameter optimization, since we only used the default hyperparameters embeddings)
