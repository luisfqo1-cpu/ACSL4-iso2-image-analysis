import pandas as pd 
import os 
import glob 
from tqdm import tqdm 
import matplotlib.pyplot as plt 
import matplotlib.colors as clr
import numpy as np 
from skimage.draw import polygon2mask
from skimage.measure import label, regionprops_table, regionprops
from imagetools.util import calculate_lag, calculate_bounding_rectangle, calculate_props, calculate_rate, fit_dim, generate_polygon_mask, crop_with_polygon_mask, df_remove_suffixes, return_mask, return_positive_negative_cells
from skimage import io
import seaborn as sns 
from statannot import add_stat_annotation

# ---- User settings: set these to your own folders ----
folder_data = "path/to/circularity_csvs"   # contains circ.csv and circ_ko.csv
folder_fig = "path/to/figures"


def measure_circularity(props, mask):

    circs = []
    for i in range(len(props)):
        bbox = props[i].bbox
        centroid =  list(props[i].centroid)
        area = props[i].area
        test = (mask == i+1)[bbox[0]:bbox[2], bbox[1]:bbox[3]]
        
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



circ = pd.read_csv(os.path.join(folder_data, "circ.csv"), index_col = 0)
circs = []
for i in tqdm(circ.polygon_coords): 
    polygon = [[int(float(j)) for j in k.split(",")] for k in i.split(";")]
    polygon = [j[2:] for j in polygon]
    rectangle = calculate_bounding_rectangle(polygon)
    mask = polygon2mask((4000,4000),polygon)
    mask2 = mask[rectangle[0][0]:rectangle[2][0], rectangle[0][1]: rectangle[2][1]]
    cell_label = mask.astype(int)
    props = regionprops(cell_label)
    circularity = measure_circularity(props, cell_label)
    circs.extend(circularity)
    
    
circ["circularity"] = circs


circ.condition = circ.condition.str.replace("CNTL", "Scrambled").str.replace("shRNA", "ACSL4 Iso 2")
sns.set_palette(["#835396", "#ba8e47"])

fig, ax = plt.subplots()
sns.swarmplot(data = circ, y = "circularity", x="condition",ax = ax)
test_results = add_stat_annotation(ax, data = circ, x = "condition", y= "circularity", test = "Mann-Whitney", box_pairs = [("Scrambled", "ACSL4 Iso 2")], verbose = 1, loc = "outside")
ax.set_ylabel("Area of Cell/Area of Circle", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)

sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, 'circularity.eps'), format = "eps")



circ_control = circ[circ.condition == "Scrambled"]
circ_sh = circ[circ.condition == "ACSL4 Iso 2"]


fig, ax = plt.subplots(5,5, figsize = (10,10))

counter = 0

for i in range(5): 
    for j in range(5):
        polygon = [[int(float(j)) for j in k.split(",")] for k in circ_control.polygon_coords.iloc[counter].split(";")]
        polygon = [j[2:] for j in polygon]
        rectangle = calculate_bounding_rectangle(polygon)
        mask = polygon2mask((4000,4000),polygon)
        mask2 = mask[rectangle[0][0]:rectangle[2][0], rectangle[0][1]: rectangle[2][1]]
        
        ax[i,j].imshow(mask2, cmap = clr.ListedColormap(['white', '#835396']))
        ax[i,j].axis("off")
        counter += 1
        
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "control_cell_examples.eps"))

fig, ax = plt.subplots(5,5, figsize = (10,10))

counter = 0

for i in range(5): 
    for j in range(5):
        polygon = [[int(float(j)) for j in k.split(",")] for k in circ_sh.polygon_coords.iloc[counter].split(";")]
        polygon = [j[2:] for j in polygon]
        rectangle = calculate_bounding_rectangle(polygon)
        mask = polygon2mask((4000,4000),polygon)
        mask2 = mask[rectangle[0][0]:rectangle[2][0], rectangle[0][1]: rectangle[2][1]]
        
        ax[i,j].imshow(mask2, cmap = clr.ListedColormap(['white', '#ba8e47']))
        ax[i,j].axis("off")
        counter += 1
        
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "sh_cell_examples.eps"))



circ = pd.read_csv(os.path.join(folder_data, "circ_ko.csv"), index_col = 0)
circs = []
for i in tqdm(circ.polygon_coords): 
    polygon = [[int(float(j)) for j in k.split(",")] for k in i.split(";")]
    polygon = [j[2:] for j in polygon]
    rectangle = calculate_bounding_rectangle(polygon)
    mask = polygon2mask((4000,4000),polygon)
    mask2 = mask[rectangle[0][0]:rectangle[2][0], rectangle[0][1]: rectangle[2][1]]
    cell_label = mask.astype(int)
    props = regionprops(cell_label)
    circularity = measure_circularity(props, cell_label)
    circs.extend(circularity)
    
    
circ["circularity"] = circs


circ.condition = circ.condition.str.replace("CNTL", "WT").str.replace("KO", "ACSL4 Cat KO")
sns.set_palette(["#835396", "#ba8e47"])

fig, ax = plt.subplots()
sns.swarmplot(data = circ, y = "circularity", x="condition",ax = ax)
test_results = add_stat_annotation(ax, data = circ, x = "condition", y= "circularity", test = "Mann-Whitney", box_pairs = [("WT", "ACSL4 Cat KO")], verbose = 1, loc = "outside")
ax.set_ylabel("Area of Cell/Area of Circle", fontsize = 15)
ax.set_xlabel("")
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='x', labelsize=15)

sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, 'circularity_ko.eps'), format = "eps")



circ_control = circ[circ.condition == "WT"]
circ_sh = circ[circ.condition == "ACSL4 Cat KO"]


fig, ax = plt.subplots(5,5, figsize = (10,10))

counter = 0

for i in range(5): 
    for j in range(5):
        polygon = [[int(float(j)) for j in k.split(",")] for k in circ_control.polygon_coords.iloc[counter].split(";")]
        polygon = [j[2:] for j in polygon]
        rectangle = calculate_bounding_rectangle(polygon)
        mask = polygon2mask((4000,4000),polygon)
        mask2 = mask[rectangle[0][0]:rectangle[2][0], rectangle[0][1]: rectangle[2][1]]
        
        ax[i,j].imshow(mask2, cmap = clr.ListedColormap(['white', '#835396']))
        ax[i,j].axis("off")
        counter += 1
        
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "control_cell_examples_ko.eps"))

fig, ax = plt.subplots(5,5, figsize = (10,10))

counter = 0

for i in range(5): 
    for j in range(5):
        polygon = [[int(float(j)) for j in k.split(",")] for k in circ_sh.polygon_coords.iloc[counter].split(";")]
        polygon = [j[2:] for j in polygon]
        rectangle = calculate_bounding_rectangle(polygon)
        mask = polygon2mask((4000,4000),polygon)
        mask2 = mask[rectangle[0][0]:rectangle[2][0], rectangle[0][1]: rectangle[2][1]]
        
        ax[i,j].imshow(mask2, cmap = clr.ListedColormap(['white', '#ba8e47']))
        ax[i,j].axis("off")
        counter += 1
        
plt.tight_layout()
plt.savefig(os.path.join(folder_fig, "sh_cell_examples_ko.eps"))
