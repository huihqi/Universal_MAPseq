import os
import pandas as pd
from datetime import datetime
#from path import path
from starmap.config import Config
from starmap.pipeline import Pipeline

# Base/default config parameters
base_config = {
    "filepath_dapi_all_z": "",
    "filepath_nissl": "",
    "filepath_dapi": "",
    "round_index": [1,2,3],
    "radius": 5,
    "fdr_thresh": 0.10,
    "expand_pixel": 0,
    "ifov": 0,
    "igenes": ["Slc17a6", "Gad1", "Atoh1", "Sox14"],
    "SpGene": 'Finished',# 'ChAlign', 'RoAlign', 'CbCorrection', 'Normalization', 'BardensrReg', 'GeneCalling', 'Finished'
    "EpGene": "Finished",#"Finished", # 'PosReg', 'CellSeg', 'ExpMask', 'Gene2Mask', 'RmvOverlap', 'Finished'
    "SpMask": "Gene2Mask",
    "EpMask": "Finished",
    "grid_stitch": False,
    "bardensr_reg": False,
    "find_param": False,
    "align_method": 'OpenCV',
    "round_num": 3,
    "thresh_refined": 0.78, #0.60 gave good results
    "bardensr_patch_size": (1000, 1000),
    "bardensr_overlap": (100, 100),
    "bardensr_niter": 10,
    "round_align": 1,
    "dapi_round": 1
}
'''

'''
# Define your list of sample-specific configurations
sample_subsets = [
    {
        "name": "Brain1",
        "filepath_homedir": r"T:\Huihui\Bird\HVC12152025\Analysis\B2\Brain1\\",
        "filepath_rawdata": r"T:\Huihui\Bird\HVC12152025\Analysis\B2\Brain1\\",
        "filepath_codebook": r"T:\Huihui\Bird\HVC12152025\Gene_List.csv",
        "fov_num": 20,
        "fov_align": 0,
        "fov_minmax": [0]
    }
    # Add more samples as needed...
]
if __name__ == '__main__': 
# Process all samples
    for subset in sample_subsets:
        print(f"\n🧠 Processing sample: {subset['name']}")

        # Merge base config with sample-specific config
        config_params = {**subset, **base_config}
        print(config_params)
        config = Config(**config_params)

        # Run pipeline
        pipeline = Pipeline(config)
        pipeline.run()

        # Run remove_overlap if required
        if config.CpMask.get('RmvOverlap') and getattr(config, 'remove_overlap', False):
            dim_arg = str(config.dim) if hasattr(config, 'dim') else ""
            if config.grid_stitch:
                os.system(f'python remove_overlap_mp.py {config.filepath_homedir} position_reg {dim_arg}')
            else:
                os.system(f'python remove_overlap_mp.py {config.filepath_homedir} position_corr {dim_arg}')
