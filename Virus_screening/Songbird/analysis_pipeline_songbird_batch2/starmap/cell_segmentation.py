"""
find cell masks, expand cell masks
map gene calls to cell masks
remove repeats in overlap
convert gene result dataframe into anndata
"""

import numpy as np
import pandas as pd
import anndata as ad
from cellpose import models
from cellpose import denoise, io
from skimage.segmentation import expand_labels
import copy
from datetime import datetime
# print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"")

# import cv2
# import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
import multiprocessing as mp
from multiprocessing import Process, Pool
import warnings
warnings.simplefilter('ignore', pd.errors.SettingWithCopyWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
import skimage as sk
# from skimage import io
# from skimage.color import rgb2gray
# from skimage.data import stereo_motorcycle, vortex
# from skimage.transform import warp
# from skimage.registration import optical_flow_tvl1, optical_flow_ilk
# from skimage.segmentation import expand_labels, watershed
# from skimage.color import label2rgb
# import IPython.display
# import SimpleITK as sitk
# import time, os, sys
# from urllib.parse import urlparse
# from cellpose import utils, io, plot
# from scipy.spatial import KDTree
import numpy as np
from scipy import ndimage as ndi
from skimage import exposure
from cellpose import models


def normalize_percentile(img, p_low=1, p_high=99.8):
    img = img.astype(np.float32)
    low = np.percentile(img, p_low)
    high = np.percentile(img, p_high)

    if high <= low:
        return np.zeros_like(img, dtype=np.float32)

    img = (img - low) / (high - low)
    img = np.clip(img, 0, 1)
    return img


def preprocess_nissl(img, background_sigma=40, clahe_clip=0.01, gamma=0.8):
    # robust normalization
    img = normalize_percentile(img, p_low=1, p_high=99.8)

    # remove slow background
    background = ndi.gaussian_filter(img, sigma=background_sigma)
    img = img - background
    img = np.clip(img, 0, None)

    if np.max(img) > 0:
        img = img / np.max(img)

    # enhance local contrast
    img = exposure.equalize_adapthist(img, clip_limit=clahe_clip)

    # enhance weak cells
    img = np.power(img, gamma)

    img = np.clip(img, 0, 1).astype(np.float32)
    return img


def find_masks(
    dapi_images,
    nissl_images,
    cellpose_model,
    cellpose_channel,
    cellpose_diameter,
    use_cpsam_default=False,
    cpsam_flow_threshold=0.7,
    cpsam_cellprob_threshold=-1.1,
    cpsam_min_size=40,
    cpsam_p_low=0.2,
    cpsam_p_high=99.3
):
    # 0-255 intensity values

    import numpy as np
    from cellpose import models

   
    # increase contrast
  
    def preprocess_for_cpsam(img, p_low=1, p_high=99.8):
        img = img.astype(np.float32)

        vmin = np.percentile(img, p_low)
        vmax = np.percentile(img, p_high)

        if vmax <= vmin:
            return img.astype(np.uint8)

        img = np.clip(img, vmin, vmax)
        img = (img - vmin) / (vmax - vmin + 1e-8)
        img = (img * 255).astype(np.uint8)

        return img

    all_masks = []

  
    # Cellpose 4.1.1 / Cellpose-SAM
  
    if cellpose_model in ["cpsam", "cellpose_sam", "cellpose-SAM", "Cellpose-SAM"]:

        # For cpsam, just use whichever image input is available.
        # In pipeline, dapi_images are Nissl images.
        if dapi_images is not None:
            input_images = dapi_images
        elif nissl_images is not None:
            input_images = nissl_images
        else:
            raise ValueError("No image input found. Both dapi_images and nissl_images are None.")

        model = models.CellposeModel(
            gpu=True,
            pretrained_model="cpsam"
        )

        for i in range(len(input_images)):

            img = input_images[i]

            # Make raw ND2 image brighter
            img = preprocess_for_cpsam(
                img,
                p_low=cpsam_p_low,
                p_high=cpsam_p_high
            )

            masks, flows, styles = model.eval(
                img,
                flow_threshold=cpsam_flow_threshold,
                cellprob_threshold=cpsam_cellprob_threshold,
                min_size=cpsam_min_size
            )

            all_masks.append(masks)

        return all_masks

    
    # Other Cellpose models, for example cyto3 / nuclei
  

    if (dapi_images is not None) and (nissl_images is not None):

        merged_images = []

        for i in range(len(dapi_images)):

            nissl_max = np.max(nissl_images[i])
            dapi_max = np.max(dapi_images[i])

            if nissl_max == 0:
                nissl_norm = nissl_images[i]
            else:
                nissl_norm = nissl_images[i] / nissl_max

            if dapi_max == 0:
                dapi_norm = dapi_images[i]
            else:
                dapi_norm = dapi_images[i] / dapi_max

            merged_images.append(
                np.stack(
                    [
                        np.zeros(dapi_images[i].shape),
                        nissl_norm,
                        dapi_norm
                    ],
                    axis=2
                )
            )

    elif (dapi_images is not None) and (nissl_images is None):

        merged_images = dapi_images

    elif (dapi_images is None) and (nissl_images is not None):

        merged_images = nissl_images

    else:
        raise ValueError("No image input found. Both dapi_images and nissl_images are None.")

    model = models.CellposeModel(
        gpu=True,
        pretrained_model=cellpose_model
    )

    for i in range(len(merged_images)):

        masks, flows, styles = model.eval(
            merged_images[i],
            diameter=cellpose_diameter,
            channels=cellpose_channel,
            flow_threshold=2.5,
            cellprob_threshold=-5,
            min_size=30
        )

        all_masks.append(masks)

    return all_masks

# Chatgpt version find_masks_better
from cellpose import models
import numpy as np

def _p01(x):
    
    lo, hi = np.percentile(x, 1.0), np.percentile(x, 99.0)
    if hi <= lo: hi = lo + 1.0
    y = (x - lo) / (hi - lo)
    return np.clip(y, 0, 1)

def find_masks_better(
    dapi_images,                
    nissl_images=None,          
    cellpose_model='nuclei',     
    cellpose_channel=(0,0),      
    cellpose_diameter=24,       
    flow_threshold=0.2,          
    cellprob_threshold=0.0,     
    normalize='p99',             
    use_gpu=True                 
):
    
    if isinstance(dapi_images, np.ndarray) and dapi_images.ndim == 3:
        imgs = [dapi_images[i] for i in range(dapi_images.shape[0])]
    else:
        imgs = list(dapi_images)

    
    model = models.Cellpose(gpu=use_gpu, model_type=str(cellpose_model))

    masks_all = []
    for idx, im in enumerate(imgs):
        im = im.astype(np.float32)
        im_norm = _p01(im) if normalize == 'p99' else (im / (im.max() + 1e-8))
        im_u8 = (im_norm * 255).astype(np.uint8)

       
        masks, flows, styles, diams = model.eval(
            im_u8,
            diameter=float(cellpose_diameter),
            channels=list(cellpose_channel) if isinstance(cellpose_channel, (list,tuple)) else [0,0],
            flow_threshold=float(flow_threshold),
            cellprob_threshold=float(cellprob_threshold)
        )
        masks_all.append(masks)

       

    return masks_all

#denoise images
def denoise_segment_images(dapi_all_z_images, diameter=6/0.16, channels=[0,0], model_type="cyto3", restore_type="deblur_cyto3", nframes=1, int_thresh = 0.9999): #0-255 intensity values
    dapi_all_z = dapi_all_z_images

    dapi_merged = []
    for i in range(len(dapi_all_z)):
        print("Merging top 8 z images for FOV " +str(i) + "...")
        intensity = []
        for j in range(len(dapi_all_z[i])):
            intensity.append(np.mean(dapi_all_z[i][j].ravel()))
        intensity = intensity - np.min(intensity)
        z_pos = []
        for k in range(len(intensity)):
            if intensity[k]>int_thresh*np.max(intensity):
                if len(z_pos) <nframes:
                    z_pos.append(k)
        dapi_merged.append(np.sum(dapi_all_z[i][z_pos], axis = 0))
    print("Denoising images and segmenting cells...")


    model = denoise.CellposeDenoiseModel(gpu=True, model_type=model_type, restore_type=restore_type)
    masks, flows, styles, dapi_merged_dn = model.eval(dapi_merged, diameter=diameter, channels=channels)
    
    return masks # , flows, styles, dapi_merged_dn

def expand_masks(masks, n):
    expanded_masks= []
    for i in range(len(masks)):
        expanded = expand_labels(masks[i], distance=n)
        expanded_masks.append(expanded)
    return expanded_masks

def find_cell_center(image, cell_value):
    
    indices = np.where(image == cell_value)
    center_x = np.mean(indices[0]) # original 
    center_y = np.mean(indices[1]) #original  
    return center_x, center_y
def map_gene2cell(genes_result_global, expanded_masks, position_reg, dim, xdim, ydim):
    genes_result_global = genes_result_global.copy()
    genes_result_global['FOV'] = genes_result_global['FOV'].astype(int)
    fovs = np.unique(genes_result_global['FOV'])
    genes_result_cell = []

    for i, fov in enumerate(fovs):
        print("Mapping FOV " + str(fov))

        mask = expanded_masks[i]

        genes_result_fov = genes_result_global[genes_result_global['FOV'] == fov].copy()

        # m1=row, m2=col
        rr = genes_result_fov['m1'].astype(int).to_numpy()
        cc = genes_result_fov['m2'].astype(int).to_numpy()

        valid = (
            (rr >= 0) & (rr < mask.shape[0]) &
            (cc >= 0) & (cc < mask.shape[1])
        )
        genes_result_fov = genes_result_fov.loc[valid].copy()
        rr = rr[valid]
        cc = cc[valid]

        genes_result_fov["cell_ID"] = mask[rr, cc]
        genes_result_fov = genes_result_fov[genes_result_fov["cell_ID"] != 0].copy()

        props = sk.measure.regionprops(mask)

        label_to_centroid = {
            prop.label: prop.centroid for prop in props
        }

        # centroid: (row, col) = (y, x)
        genes_result_fov["cell_center_X"] = genes_result_fov["cell_ID"].map(
            lambda cid: label_to_centroid[cid][1] if cid in label_to_centroid else np.nan
        )
        genes_result_fov["cell_center_Y"] = genes_result_fov["cell_ID"].map(
            lambda cid: label_to_centroid[cid][0] if cid in label_to_centroid else np.nan
        )

        genes_result_fov["cell_center_globX"] = (
            genes_result_fov["cell_center_X"] + position_reg.iloc[i]['X'] - 2304/2
        )
        genes_result_fov["cell_center_globY"] = (
            genes_result_fov["cell_center_Y"] - position_reg.iloc[i]['Y'] + np.max(position_reg['Y'])
        )

        genes_result_cell.append(genes_result_fov)

    genes_result_cell = pd.concat(genes_result_cell, ignore_index=True)
    print("Creating Cell IDs...")

    genes_result_cell['cell_number'] = (
        genes_result_cell['cell_ID'].astype(str) + '_' + genes_result_cell['FOV'].astype(str)
    )
    genes_result_cell['cell_globID'] = genes_result_cell.groupby(['FOV', 'cell_ID']).ngroup()

    return genes_result_cell

# find neighbor fovs (only with larger index)
def find_neighb_fovs(m, position_reg, dimensions):
    neighbor = []
    for n in range(len(position_reg)):
        if (abs(position_reg.loc[n, 'X']-position_reg.loc[m, 'X']) < dimensions) & (abs(position_reg.loc[n, 'Y']-position_reg.loc[m, 'Y']) < dimensions) & (n > m): 
            y_m, x_m = position_reg.loc[m, 'X'], position_reg.loc[m, 'Y']
            y_n, x_n = position_reg.loc[n, 'X'], position_reg.loc[n, 'Y']
            if (x_m-x_n > dimensions*0.1) & (y_m-y_n > dimensions*0.1):
                neighbor.append((n, "top_right", "bottom_left"))
            elif (x_m-x_n > dimensions*0.1) & (abs(y_n-y_m) < dimensions*0.1):
                neighbor.append((n, "right", "left"))
            elif (x_m-x_n > dimensions*0.1) & (y_n-y_m > dimensions*0.1):
                neighbor.append((n, "bottom_right", "top_left"))
            elif (abs(x_n-x_m) < dimensions*0.1) & (y_n-y_m > dimensions*0.1):
                neighbor.append((n, "bottom", "top"))
            elif (x_n-x_m > dimensions*0.1) & (y_n-y_m > dimensions*0.1):
                neighbor.append((n, "bottom_left", "top_right"))
            elif (x_n-x_m > dimensions*0.1) & (abs(y_m-y_n) < dimensions*0.1):
                neighbor.append((n, "left", "right"))
            elif (x_n-x_m > dimensions*0.1) & (y_m-y_n > dimensions*0.1):
                neighbor.append((n, "top_left", "bottom_right"))
            elif (abs(x_m-x_n) < dimensions*0.1) & (y_m-y_n > dimensions*0.1):
                neighbor.append((n, "top", "bottom"))
    return neighbor

# find cell in the boundaries with genecalls
def find_cell_boundary(fov, mask, cells_mapped, boundaryType, dimensions, overlap=0.15):
    boundary = int(dimensions * overlap)
    xTemp, yTemp = np.where(mask != 0)
    indexTemp = list(zip(xTemp, yTemp))
    cells_boundary = []
    for (x, y) in indexTemp:
        if boundaryType == 'top_right':
            if (x<=boundary) & (y>=(dimensions-boundary)):
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'right':
            if y>=(dimensions-boundary):
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'bottom_right':
            if (x>=(dimensions-boundary)) & (y>=(dimensions-boundary)):
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'bottom':
            if x>=(dimensions-boundary):
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'bottom_left':
            if (x>=(dimensions-boundary)) & (y<=boundary):
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'left':
            if y<=boundary:
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'top_left':
            if (x<=boundary) & (y<=boundary):
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
        elif boundaryType == 'top':
            if x<=boundary:
                cells_boundary.append(str(mask[x, y])+'_'+str(fov))
    cells_boundary = np.unique(cells_boundary)
    cells_ret = list(set(cells_boundary).intersection(set(cells_mapped)))
    cells_result = [int(i.split('_')[0]) for i in cells_ret]
    return sorted(cells_result)

# check overlapped cells between two fovs
def overlap_condition_check(m, n, maskM, maskN, cellM, cellN, position_reg, dimensions):
    xpos_corr, ypos_corr = [int(i) for i in position_reg['X']], [int(i) for i in position_reg['Y']]
    width = np.max(xpos_corr)-np.min(xpos_corr)+int(dimensions)
    height = np.max(ypos_corr)-np.min(ypos_corr)+int(dimensions)
    half = int(dimensions/2)
    
    xindexM, yindexM = np.where(maskM == cellM) # np array
    xposM = xpos_corr[m] - half + xindexM
    yposM = height - ypos_corr[m] - half + yindexM
    posM = list(zip(xposM, yposM))
    xindexN, yindexN = np.where(maskN == cellN)
    xposN = xpos_corr[n] - half + xindexN
    yposN = height - ypos_corr[n] - half + yindexN
    posN = list(zip(xposN, yposN))
    ret = list(set(posM).intersection(set(posN)))

    # start conditions
    if len(posM) == 0 or len(posN) == 0:
        return False
    else:
        ratio = len(posM)/len(posN)
        if ratio > 0.5 and ratio < 2: 
            if len(ret)/min(len(posM), len(posN)) >= 0.5:# (len(ret)/len(posM) > 0.5) & (len(ret)/len(posN) > 0.5)
                return True
            else:
                return False
        else:
            if len(ret)/min(len(posM), len(posN)) >= 0.5:
                return True
            else:
                return False

# delete the cell with less genes
def overlap_condition_remove(m, n, cellM, cellN, gene_trimmed):
    gene_cellM = np.where(gene_trimmed['cell_number']==str(cellM)+'_'+str(m))[0]
    gene_cellN = np.where(gene_trimmed['cell_number']==str(cellN)+'_'+str(n))[0]
    if len(gene_cellM) > len(gene_cellN): # delete N
        # maskN[np.where(maskN==cellN)] = 0
        cell_num = str(cellN)+'_'+str(n)
    else:                                 # delete M
        # maskM[np.where(maskM==cellM)] = 0
        cell_num = str(cellM)+'_'+str(m)
    return cell_num

def overlap_rmv_cells_and_calls(gene_mapped, expanded_masks, position_reg, dimensions):
    fov_num = len(expanded_masks)
    cell_deleted = []
    for m in range(fov_num):
        neighb_fovs = find_neighb_fovs(m, position_reg, dimensions)
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"main fov "+str(m)+":", neighb_fovs)
        for (n, bTypeM, bTypeN) in neighb_fovs:
            print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"    check neighbor fov "+str(n))
            maskM, maskN = expanded_masks[m], expanded_masks[n]
            cells_mapped = np.unique(gene_mapped['cell_number'])
            cells_resultM = find_cell_boundary(m, maskM, cells_mapped, bTypeM, dimensions)
            cells_resultN = find_cell_boundary(n, maskN, cells_mapped, bTypeN, dimensions)
            for cellN in cells_resultN:
                for cellM in cells_resultM:
                    if overlap_condition_check(m, n, maskM, maskN, cellM, cellN, position_reg, dimensions):
                        # print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"        repeated cells: ("+str(cellM)+"_"+str(m)+", "+str(cellN)+"_"+str(n)+") delete: "+cell_num)
                        cell_deleted.append(overlap_condition_remove(m, n, cellM, cellN, gene_mapped))
    trimmed_masks = copy.deepcopy(expanded_masks)
    for i in cell_deleted:
        cell_number, fov = int(i.split('_')[0]), int(i.split('_')[1])
        trimmed_masks[fov][np.where(trimmed_masks[fov]==cell_number)] = 0
    gene_trimmed = gene_mapped[~gene_mapped['cell_number'].isin(cell_deleted)]
    return gene_trimmed, trimmed_masks, cell_deleted

def overlap_rmv_subset(q_cell_deleted, fov, gene_mapped, mask_mapped, position_reg, dimensions):
    cell_deleted = []
    m = fov
    neighb_fovs = find_neighb_fovs(m, position_reg, dimensions)
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"main fov "+str(m)+":", neighb_fovs)
    for (n, bTypeM, bTypeN) in neighb_fovs:
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"    check neighbor fov "+str(n))
        maskM, maskN = mask_mapped[m], mask_mapped[n]
        cells_mapped = np.unique(gene_mapped['cell_number'])
        cells_resultM = find_cell_boundary(m, maskM, cells_mapped, bTypeM, dimensions)
        cells_resultN = find_cell_boundary(n, maskN, cells_mapped, bTypeN, dimensions)
        for cellN in cells_resultN:
            for cellM in cells_resultM:
                if overlap_condition_check(m, n, maskM, maskN, cellM, cellN, position_reg, dimensions):
                    # print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"        repeated cells: ("+str(cellM)+"_"+str(m)+", "+str(cellN)+"_"+str(n)+") delete: "+cell_num)
                    cell_deleted.append(overlap_condition_remove(m, n, cellM, cellN, gene_mapped))
    q_cell_deleted.put(cell_deleted)

    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"    finish FOV: "+str(fov)+",", str(q_cell_deleted.qsize())+"/"+str(len(position_reg))+' fovs done')

def overlap_rmv_mp(gene_mapped, mask_mapped, position_reg, dim):
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"start multiprocessing; "+"number of cpu: "+str(mp.cpu_count()))
    p = Pool(processes=mp.cpu_count()) # cannot write Queue()
    q_cell_deleted = mp.Manager().Queue()
    for fov in range(len(mask_mapped)):
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"start multi process FOV: "+str(fov))
        # temp_gene = gene_mapped[gene_mapped['FOV']==fov]
        # temp_mask = mask_mapped[fov]
        p.apply_async(overlap_rmv_subset, args=(q_cell_deleted, fov, gene_mapped, mask_mapped, position_reg, dim))
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"Wait for all subprocesses done ...")
    p.close()
    p.join()
    print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "        "+"All subprocesses done")

    cell_deleted = []
    for i in range(q_cell_deleted.qsize()):
        cell_deleted.extend(q_cell_deleted.get())
    mask_trimmed = copy.deepcopy(mask_mapped)
    for i in cell_deleted:
        cell_number, fov = int(i.split('_')[0]), int(i.split('_')[1])
        mask_trimmed[fov][np.where(mask_trimmed[fov]==cell_number)] = 0
    gene_trimmed = gene_mapped[~gene_mapped['cell_number'].isin(cell_deleted)]
    return gene_trimmed, mask_trimmed, cell_deleted

def gene_expression_matrix(masks, spots, dim):
    spots_data = pd.DataFrame(columns = ['Cell number', 'Spot number', 'Spot X Pos', 'Spot Y Pos', 'Gene'])
    cells_data = pd.DataFrame(columns = ['Cell number', 'Cell Max X', 'Cell Max Y', 'Cell Min X', 'Cell Min Y'])
    
    num_cells = np.unique(masks)
    
    spots_cells = []
    spots_temp = []
    all_spots = []
    for f in range(len(np.unique(spots['FOV']))):
        num_cells = np.unique(masks[f])
        num_cells = [num_cells != 0]

        for i in range(len(num_cells)):
            masks_sc = np.where(masks[f] == i)
            cell_position = list(zip(masks_sc[0],masks_sc[1]))
            spots_temp = []
            
            points = list(zip((spots[spots['FOV']==f]['m2']),(dim-spots[spots['FOV']==f]['m1'])))
            spots_cells = list(set(points) & set(cell_position))
            spots_cells_ind = np.where(np.array(spots_cells) == np.array(points))[0]
            if len(spots_cells_ind) != 0:
                spots_temp = pd.DataFrame()
                cells_temp = pd.DataFrame()
                spots_temp['Cell number'] = [f+i]*len(spots_cells_ind)
                spots_temp['FOV'] = [f]*len(spots_cells_ind)
                cells_temp['Cell number'] =  [f+i]*len(spots_cells_ind)

                for m in spots_cells_ind:
                    spots_temp['Spot X Pos'] = spots[spots['FOV']==f].iloc[m]['m2']
                    spots_temp['Spot Y Pos'] = dim-spots[spots['FOV']==f].iloc[m]['m1']
                    spots_temp['Gene'] = spots[spots['FOV']==f].iloc[m]['Names']
                cells_temp['Cell Max X'] =  np.max(masks_sc[0])
                cells_temp['Cell Max Y'] =  np.max(masks_sc[1])
                cells_temp['Cell Min Y'] =  np.min(masks_sc[1])
                cells_temp['Cell Min X'] =  np.min(masks_sc[0])
                spots_data = pd.concat([spots_data, spots_temp])
                cells_data = pd.concat([cells_data, cells_temp])
    return spots_data, cells_data

def convert_to_annData(gene_mat_df):
    genes = np.unique(gene_mat_df['Names'])
    cells = np.unique(gene_mat_df['cell_number'])

    gene_counts = []
    coordinates = []
    for i in range(len(cells)):
        temp = gene_mat_df[gene_mat_df['cell_number']==cells[i]]
        gene_counts.append([len(np.array(temp[temp['Names']==gene])) for gene in genes])
        coordinates.append([np.array(temp['cell_center_globX'])[0], np.array(temp['cell_center_globY'])[0]])
    gene_expression_matrix = pd.DataFrame(gene_counts, columns = genes, index = cells)
    adata = ad.AnnData(gene_expression_matrix)
    adata.obs_names = cells
    adata.var_names = genes
    adata.obsm['spatial'] = np.array(coordinates)
    #adata.uns['spatial'] = np.array(coordinates)
    return adata

