import skimage.io as io 
import numpy as np 
from cellpose import models 
import czifile
import matplotlib.pyplot as plt
import glob
import os 
from tqdm import tqdm 
from skimage.segmentation import expand_labels
import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation
from skimage.filters import threshold_otsu as otsu
from skimage.filters import threshold_triangle as triangle
from skimage.filters import try_all_threshold
import math
from skimage.measure import regionprops
diameters = {"25x": 200}

def circ(props, masks):

    circs = []
    for i in range(len(props)):
        bbox = props[i].bbox
        centroid =  list(props[i].centroid)
        area = props[i].area
        test = (masks == i+1)[bbox[0]:bbox[2], bbox[1]:bbox[3]]
        
        centroid[0] = centroid[0] - bbox[0]
        centroid[1] = centroid[1] - bbox[1]
        
        
        dist = []
        for x in range(test.shape[0]):
            for y in range(test.shape[1]): 
                if test[x, y] == 1: 
                    dist.append(np.linalg.norm([np.array(centroid) - np.array([x, y])]))
        
        max_dist = int(np.max(dist))
        
        area_circle = np.pi * max_dist **2
        
        circ = area / area_circle
        circs.append(circ)
    
    return circs


df_inj = pd.DataFrame(columns = ["mouse_id", "slice_id", "antibody", "construct", "on_pixels", "off_pixels", "s2n_on", "s2n_off", "circularity" ,"colocs_mean", "colocs_max", "cy5_pixels", "background"])
df_dapi = pd.DataFrame(columns = ["mouse_id", "slice_id", "antibody", "construct", "on_pixels", "off_pixels", "s2n_on", "s2n_off", "cy5_likelihood", "inj_likelihood", "cy5_pixels", "background"])


# ---- User settings: set these to your own folders ----
folder = "path/to/KO_images"        # renamed .tif stacks of the KO experiment
folder_fig = "path/to/figures"
folder_masks = os.path.join(folder, "masks")
os.makedirs(folder_masks, exist_ok = True)
files = glob.glob(os.path.join(folder, "*.tif"))


if os.path.exists(os.path.join(folder, "results_inj.csv")):
    df_inj = pd.read_csv(os.path.join(folder, "results_inj.csv"), index_col = 0)
if os.path.exists(os.path.join(folder, "results_dapi.csv")):
    df_dapi = pd.read_csv(os.path.join(folder, "results_dapi.csv"), index_col = 0)


df_dapi = df_dapi.drop_duplicates()
df_inj = df_inj.drop_duplicates()

for file in tqdm(files): 
    filename = file.split("/")[-1]
    mouse_id = filename.split("_")[1][-3:]
    slice_id = filename.split("_")[3].strip("Slice")
    antibody = filename.split("_")[2]
    construct = filename.split("_")[1][:-3]
    
    
    
    
    
    if antibody == "2CNTL": 
        continue
    
    if list(df_inj.loc[(df_inj.antibody == antibody) & (df_inj.construct == construct) & (df_inj.mouse_id == int(mouse_id)) & (df_inj.slice_id == int(slice_id))].values):
        print("Already done")
        continue
    else: 
        img = io.imread(file)
        cy5 = img[:, :, :, 0]
        t4 = img[:, :, :, 1]
        inj = img[:, :, :, 2]
        dapi = img[:, :, :, 3]
        
        
        cy5_mip = np.max(cy5, axis = 0)
        t4_mip = np.max(t4, axis = 0)
        inj_mip = np.max(inj, axis = 0)
        dapi_mip = np.max(dapi, axis = 0)
        
        
        
        model = models.Cellpose(gpu=False, model_type='cyto2')
        masks_inj, _, _, _ = model.eval(inj_mip, diameter=200, channels=[0,0],flow_threshold=0.4, do_3D=False)
        
        
        model_2 = models.Cellpose(gpu=False, model_type = "cyto2")
        masks_dapi, _, _, _ = model.eval(dapi_mip, diameter=200, channels=[0,0],flow_threshold=0.4, do_3D=False)
        
        props_inj = regionprops(masks_inj)
        props_dapi = regionprops(masks_dapi)
        
        
        background = np.min(t4_mip[masks_inj == 0])
        background_cy5_inj = np.min(cy5_mip[masks_inj == 0])
        background_cy5_dapi = np.min(cy5_mip[masks_dapi == 0])
        
        circularity = circ(props_inj, masks_inj)
    
        
        t4_thr = triangle(t4_mip)
        cy5_thr = triangle(cy5_mip)
        
       
        
        
        on_pixels = []
        off_pixels = []
        s2n_on = []
        s2n_off = []
        colocs_mean = []
        colocs_max = []
        cy5_pixels = []
        for cell in tqdm(range(1, np.max(masks_inj) + 1)): 
            
            on_cell = []
            off_cell = []
            coloc = []
            cy5_data = []
            for slice_ in range(t4.shape[0]): 
                t4_slice = t4[slice_]
                cy5_slice = cy5[slice_]
    
                on = t4_slice[(masks_inj == cell) & (t4_slice > t4_thr)]
                off = t4_slice[(masks_inj == cell) & (t4_slice < t4_thr)]
                c = np.corrcoef(t4_slice[masks_inj == cell], cy5_slice[masks_inj == cell])[0,1]
                coloc.append(c)
    
                on_cell.extend(on)
                off_cell.extend(off)
                cy5_data.extend(cy5_slice[masks_inj == cell])
    
            on_pixels.append(np.sum(on_cell) / props_inj[cell-1].area)
            off_pixels.append(np.sum(off_cell) / props_inj[cell-1].area)
            s2n_on.append(np.mean(on_cell))
            s2n_off.append(np.mean(off_cell))
            colocs_mean.append(np.mean(coloc))
            colocs_max.append(np.max(coloc))
            cy5_pixels.append(np.median(cy5_data /background_cy5_inj))
    
        df_inner_inj = pd.DataFrame([on_pixels, off_pixels, s2n_on, s2n_off], index = ["on_pixels", "off_pixels", "s2n_on", "s2n_off"]).T
        df_inner_inj["antibody"] = antibody
        df_inner_inj["construct"] = construct
        df_inner_inj["slice_id"] = int(slice_id)
        df_inner_inj["mouse_id"] = int(mouse_id)
        df_inner_inj["background"] = background
        df_inner_inj["circularity"] = circularity
        df_inner_inj["colocs_mean"] = colocs_mean
        df_inner_inj["colocs_max"] = colocs_max
        df_inner_inj["cy5_pixels"] = cy5_pixels
    
    
        df_inj = df_inj.merge(df_inner_inj, on = list(df_inj.columns), how = "outer")
    
        df_inj.to_csv(os.path.join(folder, "results_inj.csv"))
    
        cy5_likelihood = []
        inj_likelihood = []
        cy5_mean_intensity = []
        on_pixels_dapi = []
        off_pixels_dapi = []
        s2n_on_dapi = []
        s2n_off_dapi = []
        for cell in tqdm(range(1, np.max(masks_dapi) + 1)): 
            cy5_on = cy5_mip[(masks_dapi == cell) & (cy5_mip > cy5_thr)]
            cy5_off = cy5_mip[(masks_dapi == cell) & (cy5_mip < cy5_thr)]
            cy5_ratio = cy5_on.shape[0] / (cy5_off.shape[0] + cy5_on.shape[0])
            #GRAP threshold for considering gfap plus set to 0.2
            cy5_likelihood.append(cy5_ratio)
            
            inj_on = inj_mip[(masks_dapi == cell) & (masks_inj > 0)]
            inj_off = inj_mip[(masks_dapi == cell) & (masks_inj == 0)]
            inj_ratio = inj_on.shape[0] / (inj_off.shape[0] + inj_on.shape[0])
            #GRAP threshold for considering gfap plus set to 0.2
            inj_likelihood.append(inj_ratio)
            
            
            
            #calculating the gfap intesity of each dapi cell 
            cy5_int = []
            on_cell = []
            off_cell = []
            for slice_ in range(cy5.shape[0]): 
                
                cy5_slice = cy5[slice_]
                t4_slice = t4[slice_]
    
                cy5_int.append(np.mean(cy5_slice[masks_dapi == cell]))
                on = t4_slice[(masks_dapi == cell) & (t4_slice > t4_thr)]
                off = t4_slice[(masks_dapi == cell) & (t4_slice < t4_thr)]
                on_cell.extend(on)
                off_cell.extend(off)
            cy5_mean_intensity.append(np.median(cy5_int) / background_cy5_dapi)
            on_pixels_dapi.append(np.sum(on_cell) / props_dapi[cell-1].area)
            off_pixels_dapi.append(np.sum(off_cell) / props_dapi[cell-1].area)
            s2n_on_dapi.append(np.mean(on_cell))
            s2n_off_dapi.append(np.mean(off_cell))
            
        df_inner_dapi = pd.DataFrame([cy5_likelihood,inj_likelihood, cy5_mean_intensity, on_pixels_dapi, off_pixels_dapi, s2n_on_dapi, s2n_off_dapi], index = ["cy5_likelihood","inj_likelihood","cy5_pixels","on_pixels", "off_pixels", "s2n_on", "s2n_off"]).T
        df_inner_dapi["antibody"] = antibody
        df_inner_dapi["construct"] = construct
        df_inner_dapi["slice_id"] = int(slice_id)
        df_inner_dapi["mouse_id"] = int(mouse_id)
        df_inner_dapi["background"] = background
    
        df_dapi = df_dapi.merge(df_inner_dapi, on = list(df_dapi.columns), how = "outer")
        
        df_dapi.to_csv(os.path.join(folder, "results_dapi.csv"))
    
        np.save(os.path.join(folder_masks, filename.rstrip(".tif") + "_injmask.npy"), masks_inj)
        np.save(os.path.join(folder_masks, filename.rstrip(".tif") + "_dapimask.npy"), masks_dapi)



t4 = df_inj[df_inj.antibody != "2CNTL"]
t4.construct = t4.construct.str.replace("CNTL", "WT").str.replace("KO", "ACSL4 Cat KO")

t4.sort_values("construct",ascending = False,  inplace = True)

N = t4.groupby("construct").mouse_id.agg("unique")
n = t4.groupby("construct").mouse_id.agg("count")
mean = t4.groupby("construct").s2n_on.agg("mean").astype(int)
median = t4.groupby("construct").s2n_on.agg("median").astype(int)

sns.set_palette(["#835396", "#ba8e47"])

fig,ax = plt.subplots()
sns.violinplot(data = t4, x = "construct", y = "s2n_on", ax = ax)
test_results = add_stat_annotation(ax, data = t4.groupby(["mouse_id", "slice_id", "construct"]).agg("mean").reset_index(), x = "construct", y= "s2n_on", test = "Mann-Whitney", box_pairs = [("WT", "ACSL4 Cat KO")], verbose = 1, loc = "inside")
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
ax.text(0.1, 215, f"N={len(N.WT)}")
ax.text(0.1, 200, f"n={n.WT}")
ax.text(0.1, 185, f"mean={mean.WT}")
ax.text(0.1, 170, f"median={median.WT}")


ax.text(1.1, 215, f"N={len(N['ACSL4 Cat KO'])}")
ax.text(1.1, 200, f"n={n['ACSL4 Cat KO']}")
ax.text(1.1, 185, f"mean={mean['ACSL4 Cat KO']}")
ax.text(1.1, 170, f"median={median['ACSL4 Cat KO']}")


sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, 't4_KO_on.eps'), format = "eps")

N = t4.groupby("construct").mouse_id.agg("unique")
n = t4.groupby("construct").mouse_id.agg("count")
mean = t4.groupby("construct").s2n_off.agg("mean").astype(int)
median = t4.groupby("construct").s2n_off.agg("median").astype(int)


fig,ax = plt.subplots()
sns.violinplot(data = t4, x = "construct", y = "s2n_off", ax = ax)
test_results = add_stat_annotation(ax, data = t4.groupby(["mouse_id", "slice_id", "construct"]).agg("mean").reset_index(), x = "construct", y= "s2n_off", test = "Mann-Whitney", box_pairs = [("WT", "ACSL4 Cat KO")], verbose = 1, loc = "inside")
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
ax.text(0.1, 50, f"N={len(N.WT)}")
ax.text(0.1, 47, f"n={n.WT}")
ax.text(0.1, 44, f"mean={mean.WT}")
ax.text(0.1, 41, f"median={median.WT}")



ax.text(1.1, 50, f"N={len(N['ACSL4 Cat KO'])}")
ax.text(1.1, 47, f"n={n['ACSL4 Cat KO']}")
ax.text(1.1, 44, f"mean={mean['ACSL4 Cat KO']}")
ax.text(1.1, 41, f"median={median['ACSL4 Cat KO']}")


sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, 't4_KO_off.eps'), format = "eps")
