
We used BARseq3 pipeline to analyze Forg, Mouse and Songbird virus screening data and BARseq2 pipeline to analyze Marmoset and Macaque virus screening data. The pipeline was included in the folder of each species.

For BARseq3 pipeline, install the environment using environment.yml file. We used cellpose 4 (Pachitariu et al., 2025, bioRxiv) to perform cell segmentation and Batch_process_win.py to get gene expression matrix. For barcode, we only used the pipeline to do round alignment and cell segmentation in barcode module. We used barcode_last_step.ipynb to call barcode sequence and combine gene and barcode matrice. The virus_and_cell_type.py script was used to define cell types.

For BARseq2 general pipeline, please refer to https://github.com/chenxy877/barseq_utils.git. Species-specific downstream analysis scripts are located in the corresponding species folders. 
