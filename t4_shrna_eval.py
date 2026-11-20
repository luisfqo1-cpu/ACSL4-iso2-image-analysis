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
def return_pairs(labels): 
    pairs = []
    
    for i in range(len(labels)): 
        rest = list(range(len(labels)))
        rest.remove(i)
        for j in rest: 
            
            if sum([True if labels[i] in k and labels[j] in k else False for k in pairs]) < 1:
                pairs.append((labels[i], labels[j]))

    return pairs

df_inj = pd.DataFrame(columns = ["mouse_id", "slice_id", "antibody", "construct", "on_pixels", "off_pixels", "s2n_on", "s2n_off", "circularity" ,"colocs_mean", "colocs_max", "cy5_pixels", "background"])
df_dapi = pd.DataFrame(columns = ["mouse_id", "slice_id", "antibody", "construct", "on_pixels", "off_pixels", "s2n_on", "s2n_off", "cy5_likelihood", "inj_likelihood", "cy5_pixels", "background"])

# ---- User settings: set these to your own folders ----
folder = "path/to/shRNA_everything"          # renamed .tif stacks (output of renames_shRNA.py)
folder_fig = "path/to/figures"
folder_cb1_magl = "path/to/CB1_MAGL_images"   # CB1/MAGL co-staining stacks
folder_neun_cb1 = "path/to/NeuN_CB1_images"   # same folder as in t4_neun_cb1.py

folder_masks = os.path.join(folder, "masks")
os.makedirs(folder_masks, exist_ok = True)
files = glob.glob(os.path.join(folder, "*.tif"))



if os.path.exists(os.path.join(folder, "results_inj_feb.csv")):
    df_inj = pd.read_csv(os.path.join(folder, "results_inj_feb.csv"), index_col = 0)
if os.path.exists(os.path.join(folder, "results_dapi.csv")):
    df_dapi = pd.read_csv(os.path.join(folder, "results_dapi.csv"), index_col = 0)


df_dapi = df_dapi.drop_duplicates()
df_inj = df_inj.drop_duplicates()


for file in tqdm(files): 
    filename = file.split("/")[-1]
    mouse_id = filename.split("_")[0]
    slice_id = filename.split("_")[1]
    construct = filename.split("_")[2]
    antibody = filename.split("_")[3].rstrip(".tif")
    
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
    
        df_inj.to_csv(os.path.join(folder, "results_inj_feb.csv"))
    
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




#Secondary controls 


pairs = [('2CNTL', 'MAGL'),
 ('2CNTL', 'T4'),
 ('2CNTL', 'GFAP'),
 ('2CNTL', 'PSD95'),
 ('2CNTL', 'SYT'),]


sec = df_inj[df_inj.construct != "shRNA"]
sec.antibody = sec.antibody.replace("T4-MAGL", "MAGL").replace("T4-PSD95", "PSD95").replace("T4-SYT", "SYT")

fig, ax = plt.subplots()
sns.boxplot(data = sec, y = "s2n_on", x = "antibody", showfliers = False)
test_results = add_stat_annotation(ax, data = sec, x = "antibody", y= "s2n_on", test = "Mann-Whitney", box_pairs = pairs, verbose = 1, loc = "inside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15, rotation = 30)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "onpixels_seccontrol.eps"))


fig, ax = plt.subplots()

sns.boxplot(data = sec, y = "s2n_off", x = "antibody", showfliers = False)
test_results = add_stat_annotation(ax, data = sec, x = "antibody", y= "s2n_off", test = "Mann-Whitney", box_pairs = pairs, verbose = 1, loc = "inside")
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15, rotation = 30)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "offpixels_seccontrol.eps"))



t4 = df_inj[df_inj.antibody != "2CNTL"]
t4.construct = t4.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")

t4 = t4.groupby(["mouse_id", "slice_id", "construct"]).agg("mean").reset_index()

t4.sort_values("construct", inplace = True, ascending = False)

N = t4.groupby("construct").mouse_id.agg("unique")
n = t4.groupby("construct").mouse_id.agg("count")
mean = t4.groupby("construct").s2n_on.agg("mean").astype(int)
median = t4.groupby("construct").s2n_on.agg("median").astype(int)

sns.set_palette(["#835396", "#ba8e47"])

fig,ax = plt.subplots()
sns.violinplot(data = t4, x = "construct", y = "s2n_on", ax = ax)
sns.swarmplot(data = t4, x = "construct", y = "s2n_on", ax = ax, color = "k", alpha = 0.8)

test_results = add_stat_annotation(ax, data = t4, x = "construct", y= "s2n_on", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 0, loc = "outside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
ax.text(0.2, 100, f"N={len(N.Scrambled)}", fontsize = 15)
ax.text(0.2, 90, f"n={n.Scrambled}", fontsize = 15)
ax.text(0.2, 80, f"mean={mean.Scrambled}", fontsize = 15)
ax.text(0.2, 70, f"median={median.Scrambled}", fontsize = 15)


ax.text(1.2, 100, f"N={len(N['ACSL4 Iso 2'])}", fontsize = 15)
ax.text(1.2, 90, f"n={n['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.2, 80, f"mean={mean['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.2, 70, f"median={median['ACSL4 Iso 2']}", fontsize = 15)


sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, 't4_on.pdf'))

N = t4.groupby("construct").mouse_id.agg("unique")
n = t4.groupby("construct").mouse_id.agg("count")
mean = t4.groupby("construct").s2n_off.agg("mean").astype(int)
median = t4.groupby("construct").s2n_off.agg("median").astype(int)


fig,ax = plt.subplots()
sns.violinplot(data = t4, x = "construct", y = "s2n_off", ax = ax)
sns.swarmplot(data = t4, x = "construct", y = "s2n_off", ax = ax, color = "k", alpha = 0.8)

test_results = add_stat_annotation(ax, data = t4, x = "construct", y= "s2n_off", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 0, loc = "outside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
ax.text(0.2, 30, f"N={len(N.Scrambled)}", fontsize = 15)
ax.text(0.2, 27, f"n={n.Scrambled}", fontsize = 15)
ax.text(0.2, 24, f"mean={mean.Scrambled}", fontsize = 15)
ax.text(0.2, 21, f"median={median.Scrambled}", fontsize = 15)


ax.text(1.2, 30, f"N={len(N['ACSL4 Iso 2'])}", fontsize = 15)
ax.text(1.2, 27, f"n={n['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.2, 24, f"mean={mean['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.2, 21, f"median={median['ACSL4 Iso 2']}", fontsize = 15)


sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, 't4_off.pdf'))




gfp_cells = df_dapi[df_dapi.antibody != "2CNTL"]
gfp_cells.inj_likelihood = gfp_cells.inj_likelihood > 0.5
gfp_cells.cy5_likelihood = gfp_cells.cy5_likelihood > 0.2

gfp_rate = pd.DataFrame(gfp_cells.groupby(["construct", "mouse_id"]).inj_likelihood.agg("sum") / gfp_cells.groupby(["construct", "mouse_id"]).inj_likelihood.agg("count") * 100)
gfp_rate.reset_index(inplace =True)
gfp_rate.construct = gfp_rate.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "T4 Iso 2")


fig,ax = plt.subplots()
sns.swarmplot(data = gfp_rate, x = "construct", y = "inj_likelihood", ax = ax, s = 10)
test_results = add_stat_annotation(ax, data = gfp_rate, x = "construct", y= "inj_likelihood", test = "Mann-Whitney", box_pairs = [("Scrambled", "T4 Iso 2")], verbose = 0, loc = "inside")
ax.set_ylabel("GFP+ %", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()

plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "GFP_percentage.pdf"))



gfap = df_dapi[df_dapi.antibody == "GFAP"]
gfap.construct = gfap.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")

gfap.cy5_likelihood = gfap.cy5_likelihood > 0.2
gfap.inj_likelihood = gfap.inj_likelihood > 0.5

gfap = gfap.sort_values(by ="construct", ascending = False)




cell_type = []

for i in range(gfap.shape[0]): 
    
    if gfap.cy5_likelihood.iloc[i] == True: 
        cell_type.append("GFAP +")
    elif gfap.inj_likelihood.iloc[i] == True: 
        cell_type.append("GFP +")
    else: 
        cell_type.append("GFAP-/GFP-")

gfap["Cell Type"] = cell_type



gfap_rate = pd.DataFrame(gfap.groupby(["construct", "slice_id"]).cy5_likelihood.agg("sum") / gfap.groupby(["construct", "slice_id"]).cy5_likelihood.agg("count") * 100)
gfap_rate.cy5_likelihood = gfap_rate.cy5_likelihood.astype(int)
gfap_rate.reset_index(inplace =True)

gfap_rate.construct = gfap_rate.construct.str.replace("T4", "ACSL4")

fig,ax = plt.subplots()
sns.swarmplot(data = gfap_rate, x = "construct", y = "cy5_likelihood", ax = ax, s = 10)
test_results = add_stat_annotation(ax, data = gfap_rate, x = "construct", y= "cy5_likelihood", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 0, loc = "inside", fontsize = 15)
ax.set_ylabel("GFAP+ %", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
plt.tight_layout()





box_pairs = [
    (("Scrambled", "GFP +"),("Scrambled", "GFAP +")),
    (("Scrambled", "GFP +"),("Scrambled", "GFAP-/GFP-")),
    (("Scrambled", "GFAP +"),("Scrambled", "GFAP-/GFP-")),
    (("Scrambled", "GFAP +"),("ACSL4 Iso 2", "GFAP +")),

    (("ACSL4 Iso 2", "GFP +"),("ACSL4 Iso 2", "GFAP +")),
    (("ACSL4 Iso 2", "GFP +"),("ACSL4 Iso 2", "GFAP-/GFP-")),
    (("ACSL4 Iso 2", "GFAP +"),("ACSL4 Iso 2", "GFAP-/GFP-")),]


box_pairs = [
    (("Scrambled", "GFP +"),("ACSL4 Iso 2", "GFP +")),
    (("Scrambled", "GFAP +"),("ACSL4 Iso 2", "GFAP +")),
    (("Scrambled", "GFAP-/GFP-"),("ACSL4 Iso 2", "GFAP-/GFP-"))]
 
fig,ax = plt.subplots()
sns.boxplot(data = gfap, x = "construct",hue = "Cell Type", y = "s2n_on", ax = ax, showfliers = False, palette = ["green", "magenta", "cyan"])
test_results = add_stat_annotation(ax, data = gfap,hue  = "Cell Type", x = "construct", y= "s2n_on", test = "Mann-Whitney", box_pairs = box_pairs, verbose = 0, loc = "inside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "gfap_t4_on.pdf"))



fig,ax = plt.subplots()
sns.boxplot(data = gfap, x = "construct",hue = "Cell Type", y = "s2n_off", ax = ax, showfliers = False, palette = ["green", "magenta", "cyan"])
test_results = add_stat_annotation(ax, data = gfap,hue  = "Cell Type", x = "construct", y= "s2n_off", test = "Mann-Whitney", box_pairs = box_pairs, verbose = 0, loc = "inside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "gfap_t4_off.pdf"))




img = io.imread(os.path.join(folder, "20_4_CNTL_T4.tif"))
mask = np.load(os.path.join(folder, "masks/20_4_CNTL_T4_injmask.npy"))

cy5 = img[:, :, :, 0]
t4 = img[:, :, :, 1]
inj = img[:, :, :, 2]
dapi = img[:, :, :, 3]


cy5_mip = np.max(cy5, axis = 0)
t4_mip = np.max(t4, axis = 0)
inj_mip = np.max(inj, axis = 0)
dapi_mip = np.max(dapi, axis = 0)

t4_thr = triangle(t4_mip)


np.save(os.path.join(folder, "example_on_pixels.npy"), (t4_mip > t4_thr))
np.save(os.path.join(folder, "example_off_pixels.npy"), ((mask != 0) &( t4_mip < t4_thr)))



img = io.imread(os.path.join(folder, "13_1_shRNA_T4.tif"))
mask = np.load(os.path.join(folder, "masks/13_1_shRNA_T4_injmask.npy"))

cy5 = img[:, :, :, 0]
t4 = img[:, :, :, 1]
inj = img[:, :, :, 2]
dapi = img[:, :, :, 3]


cy5_mip = np.max(cy5, axis = 0)
t4_mip = np.max(t4, axis = 0)
inj_mip = np.max(inj, axis = 0)
dapi_mip = np.max(dapi, axis = 0)

t4_thr = triangle(t4_mip)


np.save(os.path.join(folder, "example_on_pixels_2.npy"), (t4_mip > t4_thr))
np.save(os.path.join(folder, "example_off_pixels_2.npy"), ((mask != 0) &( t4_mip < t4_thr)))









sns.set_palette(["#835396","#ba8e47"])

#MAGL
magl = df_inj.loc[df_inj.antibody == "T4-MAGL"]
magl.construct = magl.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")
magl = magl.groupby(["mouse_id", "slice_id", "construct"]).agg("mean").reset_index()

magl.sort_values("construct", inplace = True, ascending = False)


fig,ax = plt.subplots()
sns.violinplot(data = magl, x = "construct", y = "cy5_pixels", ax = ax, s = 10)
sns.swarmplot(data = magl, x = "construct", y = "cy5_pixels", ax = ax, alpha = 0.8, color = "k")

test_results = add_stat_annotation(ax, data = magl, x = "construct", y= "cy5_pixels", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()

N = magl.groupby("construct").mouse_id.agg("unique")
n = magl.groupby("construct").mouse_id.agg("count")
mean = magl.groupby("construct").cy5_pixels.agg("mean").astype(int)
median = magl.groupby("construct").cy5_pixels.agg("median").astype(int)

ax.text(0.1, 30, f"N={len(N.Scrambled)}", fontsize = 15)
ax.text(0.1, 26, f"n={n.Scrambled}", fontsize = 15)
ax.text(0.1, 22, f"mean={mean.Scrambled}", fontsize = 15)
ax.text(0.1, 18, f"median={median.Scrambled}", fontsize = 15)


ax.text(1.1, 30, f"N={len(N['ACSL4 Iso 2'])}", fontsize = 15)
ax.text(1.1, 26, f"n={n['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.1, 22, f"mean={mean['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.1, 18, f"median={median['ACSL4 Iso 2']}", fontsize = 15)

plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "magl.pdf"))



fig,ax = plt.subplots()
sns.violinplot(data = magl, x = "construct", y = "colocs_mean", ax = ax, s = 10)
sns.swarmplot(data = magl, x = "construct", y = "colocs_mean", ax = ax, alpha = 0.8, color ="k")

test_results = add_stat_annotation(ax, data = magl, x = "construct", y= "colocs_mean", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside", fontsize = 15)
ax.set_ylabel("Mean Colocalization", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "magl_coloc.pdf"))


sns.set_palette(["#835396","#ba8e47"])

fig, ax = plt.subplots()
magl = df_inj.loc[df_inj.antibody == "T4-MAGL"]
magl.construct = magl.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")

magl.sort_values("construct", inplace = True, ascending = False)

h =sns.jointplot(data = magl, y = "s2n_on", x = "cy5_pixels",hue = "construct", alpha = .8, ax = ax, s= 40, linewidth = 0)
h.ax_joint.set_ylabel("ACSL4 On Pixels (AU)", fontsize = 15)
h.ax_joint.set_xlabel("MAGL Intensity (AU)", fontsize = 15)
h.refline(x=20, y=100, marginal = False)
h.plot_joint(sns.kdeplot, zorder=0, levels=2)
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "magl_scatter.pdf"))


#SYT staining
sns.set_palette(["#835396","#ba8e47"])

syt = df_inj.loc[df_inj.antibody == "T4-SYT"]
syt.construct = syt.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")
syt = syt.groupby(["mouse_id", "slice_id", "construct"]).agg("mean").reset_index()

syt.sort_values("construct", inplace = True, ascending = False)


fig,ax = plt.subplots()
sns.violinplot(data = syt, x = "construct", y = "cy5_pixels", ax = ax, s = 10)
sns.swarmplot(data = syt, x = "construct", y = "cy5_pixels", ax = ax, alpha = 0.8, color ="k")

test_results = add_stat_annotation(ax, data = syt, x = "construct", y= "cy5_pixels", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()

N = syt.groupby("construct").mouse_id.agg("unique")
n = syt.groupby("construct").mouse_id.agg("count")
mean = syt.groupby("construct").cy5_pixels.agg("mean").astype(int)
median = syt.groupby("construct").cy5_pixels.agg("median").astype(int)

ax.text(0.1, 40, f"N={len(N.Scrambled)}", fontsize = 15)
ax.text(0.1, 36, f"n={n.Scrambled}", fontsize = 15)
ax.text(0.1, 32, f"mean={mean.Scrambled}", fontsize = 15)
ax.text(0.1, 28, f"median={median.Scrambled}", fontsize = 15)


ax.text(1.1, 40, f"N={len(N['ACSL4 Iso 2'])}", fontsize = 15)
ax.text(1.1, 36, f"n={n['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.1, 32, f"mean={mean['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.1, 28, f"median={median['ACSL4 Iso 2']}", fontsize = 15)


plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "syt.pdf"))


fig,ax = plt.subplots()
sns.violinplot(data = syt, x = "construct", y = "colocs_mean", ax = ax, s = 10)
sns.swarmplot(data = syt, x = "construct", y = "colocs_mean", ax = ax, alpha = 0.8, color = "k")

test_results = add_stat_annotation(ax, data = syt, x = "construct", y= "colocs_mean", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside", fontsize = 15)
ax.set_ylabel("Mean Colocalization", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()

plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "syt_coloc.pdf"))
#Between cells colocalization for SYT

syt = df_inj.loc[df_inj.antibody == "T4-SYT"]
syt.construct = syt.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")
syt.sort_values("construct", inplace = True, ascending = False)

fig, ax = plt.subplots()

h =sns.jointplot(data = syt, y = "s2n_on", x = "cy5_pixels",hue = "construct", alpha = .8, ax = ax, s= 40, linewidth = 0)
h.ax_joint.set_ylabel("ACSL4 On Pixels (AU)", fontsize = 15)
h.ax_joint.set_xlabel("SYT Intensity (AU)", fontsize = 15)
#h.refline(x=20, y=100, marginal = False)
h.plot_joint(sns.kdeplot, zorder=0, levels=2)
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "syt_scatter.pdf"))



#PSD95 staining
sns.set_palette(["#835396", "#ba8e47"])

psd95 = df_inj.loc[df_inj.antibody == "T4-PSD95"]
psd95 = psd95.loc[~((psd95.mouse_id == 13) & (psd95.slice_id == 4)), :]
psd95.construct = psd95.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")
psd95 = psd95.groupby(["mouse_id", "slice_id", "construct"]).agg("mean").reset_index()

psd95.sort_values("construct", inplace = True, ascending = False)


fig,ax = plt.subplots()
sns.violinplot(data = psd95, x = "construct", y = "cy5_pixels", ax = ax, s = 10)
sns.swarmplot(data = psd95, x = "construct", y = "cy5_pixels", ax = ax, alpha = 0.8, color = "k")

test_results = add_stat_annotation(ax, data = psd95, x = "construct", y= "cy5_pixels", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside", fontsize = 15)
ax.set_ylabel("Pixel Intensity (AU)", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()

N = psd95.groupby("construct").mouse_id.agg("unique")
n = psd95.groupby("construct").mouse_id.agg("count")
mean = psd95.groupby("construct").cy5_pixels.agg("mean").astype(int)
median = psd95.groupby("construct").cy5_pixels.agg("median").astype(int)

ax.text(0.1, 100, f"N={len(N.Scrambled)}", fontsize = 15)
ax.text(0.1, 90, f"n={n.Scrambled}", fontsize = 15)
ax.text(0.1, 80, f"mean={mean.Scrambled}", fontsize = 15)
ax.text(0.1, 70, f"median={median.Scrambled}", fontsize = 15)


ax.text(1.1, 100, f"N={len(N['ACSL4 Iso 2'])}", fontsize = 15)
ax.text(1.1, 90, f"n={n['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.1, 80, f"mean={mean['ACSL4 Iso 2']}", fontsize = 15)
ax.text(1.1, 70, f"median={median['ACSL4 Iso 2']}", fontsize = 15)


plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "psd95.pdf"))

fig,ax = plt.subplots()
sns.violinplot(data = psd95, x = "construct", y = "colocs_mean", ax = ax, s = 10)
sns.swarmplot(data = psd95, x = "construct", y = "colocs_mean", ax = ax, alpha = 0.8, color = "k")

test_results = add_stat_annotation(ax, data = psd95, x = "construct", y= "colocs_mean", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside", fontsize = 15)
ax.set_ylabel("Mean Colocalization", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "psd95_coloc.pdf"))

#Between cells colocalization for PSD95
fig, ax = plt.subplots()
psd95 = df_inj.loc[df_inj.antibody == "T4-PSD95"]
psd95 = psd95.loc[~((psd95.mouse_id == 13) & (psd95.slice_id == 4)), :]
psd95.construct = psd95.construct.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")

psd95.sort_values("construct", inplace = True, ascending = False)


h =sns.jointplot(data = psd95, y = "s2n_on", x = "cy5_pixels",hue = "construct", alpha = .8, ax = ax, s= 40, linewidth = 0)
h.ax_joint.set_ylabel("ACSL4 On Pixels (AU)", fontsize = 15)
h.ax_joint.set_xlabel("PSD95 Intensity (AU)", fontsize = 15)
#h.refline(x=20, y=100, marginal = False)
h.plot_joint(sns.kdeplot, zorder=0, levels=2)
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "psd95_scatter.pdf"))

#colocalizations of MAGL and CB1

df_cb1_magl = pd.DataFrame(columns = ["mouse_id", "magl_on", "magl_off", "colocs_mean", "colocs_max", "cb1_pixels"])


files = glob.glob(os.path.join(folder_cb1_magl, "*.tif"))

for file in tqdm(files): 
    img = io.imread(file)
    magl = img[:, :, :, 0]
    cb1 = img[:,:,:, 2]
    dapi = img[:, :, :, 3]
    
    mouse_id = file.split("/")[-1].split("_")[2]
    
    
    magl_mip = np.max(magl, axis = 0)
    cb1_mip = np.max(cb1, axis = 0)
    dapi_mip = np.max(dapi, axis = 0)
    
    model = models.Cellpose(gpu=False, model_type='cyto2')
    masks, _, _, _ = model.eval(dapi_mip, diameter=200, channels=[0,0],flow_threshold=0.4, do_3D=False)
    
    props = regionprops(masks)
    background_cb1 = np.min(cb1_mip[masks == 0])

    magl_thr = triangle(magl_mip)
    cb1_thr = triangle(cb1_mip)

    magl_on = []
    magl_off = []
    colocs_mean = []
    colocs_max = []
    cb1_pixels = []
    for cell in tqdm(range(1, np.max(masks) + 1)): 
        
        on_cell = []
        off_cell = []
        coloc = []
        cb1_data = []
        for slice_ in range(magl.shape[0]): 
            magl_slice = magl[slice_]
            cb1_slice = cb1[slice_]
    
            on = magl_slice[(masks == cell) & (magl_slice > magl_thr)]
            off = magl_slice[(masks == cell) & (magl_slice < magl_thr)]
            c = np.corrcoef(magl_slice[masks == cell], cb1_slice[masks == cell])[0,1]
            coloc.append(c)
    
            on_cell.extend(on)
            off_cell.extend(off)
            cb1_data.extend(cb1_slice[masks == cell])
    
        magl_on.append(np.mean(on_cell))
        magl_off.append(np.mean(off_cell))
        colocs_mean.append(np.mean(coloc))
        colocs_max.append(np.max(coloc))
        cb1_pixels.append(np.median(cb1_data /background_cb1))

    df_inner = pd.DataFrame([magl_on, magl_off, colocs_mean, colocs_max, cb1_pixels], index = ["magl_on", "magl_off", "colocs_mean", "colocs_max", "cb1_pixels"]).T
    df_inner["mouse_id"] = mouse_id
    
    df_cb1_magl = df_cb1_magl.merge(df_inner, on = list(df_cb1_magl.columns), how = "outer")
    
    df_cb1_magl.to_csv(os.path.join(folder_cb1_magl, "cb1_magl.csv"))



plt.hist(df_cb1_magl.colocs_mean, bins = 40, edgecolor = "k")


#Between cells colocalization for CB1 and MAGL
fig, ax = plt.subplots()

h =sns.jointplot(data = df_cb1_magl, y = "magl_on", x = "cb1_pixels", alpha = .2, ax = ax, s= 40)
h.ax_joint.set_ylabel("MAGL Intensity (AU)", fontsize = 15)
h.ax_joint.set_xlabel("CB1 Intensity (AU)", fontsize = 15)
#h.refline(x=20, y=100, marginal = False)
h.plot_joint(sns.kdeplot, zorder=0, levels=2)
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "cb1_MAGL_scatter.pdf"))



#Colocalization of all control cells (non-SHRNA)

colocs_t4_cb1 = pd.read_csv(os.path.join(folder_neun_cb1, "results_inj.csv"), index_col = 0)
colocs_t4_cb1 = colocs_t4_cb1.loc[colocs_t4_cb1.region == "CA3"]
colocs_t4_cb1 = pd.DataFrame(colocs_t4_cb1["colocs_mean"])

colocs_t4_cb1["antibody"] = "ACSL4-CB1"




colocs_all = df_inj.loc[(df_inj.construct == "CNTL") & (df_inj.antibody != "2CNTL") & (df_inj.antibody != "T4")][["antibody", "colocs_mean"]]

colocs_all.antibody = colocs_all.antibody.str.replace("GFAP", "T4-GFAP")
colocs_all.antibody = colocs_all.antibody.str.replace("T4", "ACSL4")

colocs_magl_cb1 = pd.DataFrame(df_cb1_magl["colocs_mean"])
colocs_magl_cb1["antibody"] = "CB1-MAGL"



colocs_all = colocs_all.merge(colocs_magl_cb1, how = "outer", on = ["antibody", "colocs_mean"])
colocs_all = colocs_all.merge(colocs_t4_cb1, how = "outer", on = ["antibody", "colocs_mean"])




fig, ax = plt.subplots()

order = colocs_all.groupby("antibody").colocs_mean.agg("mean").sort_values().index.unique()

pairs = return_pairs(colocs_all.antibody.unique())


sns.violinplot(data = colocs_all, x = "antibody", y = "colocs_mean", ax= ax, order = order)
ax.set_ylabel("Mean Colocalization", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15, rotation = 15)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "all_colocalizations.pdf"))








