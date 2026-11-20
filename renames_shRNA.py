import os 
import glob
from tqdm import tqdm 
# ---- User settings: set these to your own folders ----
# Folder that receives all renamed files
folder_export = "path/to/shRNA_everything"

# Source folders with the deconvolved .tif stacks of each staining
folder_t4 = "path/to/raw/T4"
folder_magl = "path/to/raw/MAGL"
folder_gfap = "path/to/raw/GFAP"
folder_syt = "path/to/raw/SYT"

folder_psd_feb = "path/to/raw/PSD95"



for file in tqdm(glob.glob(os.path.join(folder_t4, "*.tif"))):
    construct = file.split("/")[-1].split("_")[3][:-2]
    antibody = file.split("/")[-1].split("_")[1]
    mouse_id = file.split("/")[-1].split("_")[3][-2:]
    slice_id = file.split("/")[-1].split("_")[2].strip("Slice")
    
    
    filename_new = mouse_id + "_" + slice_id + "_"  +construct + "_" + antibody + ".tif"

    os.rename(file, os.path.join(folder_export, filename_new))
    
    
    

for file in tqdm(glob.glob(os.path.join(folder_magl, "*.tif"))):
    construct = file.split("/")[-1].split("_")[3][:-2]
    antibody = file.split("/")[-1].split("_")[1]
    mouse_id = file.split("/")[-1].split("_")[3][-2:]
    slice_id = file.split("/")[-1].split("_")[2].strip("Slice")
    
    
    filename_new = mouse_id + "_" + slice_id + "_"  +construct + "_" + antibody + ".tif"

    os.rename(file, os.path.join(folder_export, filename_new))
    
    

for file in tqdm(glob.glob(os.path.join(folder_gfap, "*.tif"))):
    construct = file.split("/")[-1].split("_")[3][:-2]
    antibody = file.split("/")[-1].split("_")[4].split("-")[0]
    mouse_id = file.split("/")[-1].split("_")[3][-2:]
    slice_id = file.split("/")[-1].split("_")[4].split("Slice")[1].strip("-OME")
    
    
    filename_new = mouse_id + "_" + slice_id + "_"  +construct + "_" + antibody + ".tif"
    print(filename_new)
    os.rename(file, os.path.join(folder_export, filename_new))
    
    
    
for file in tqdm(glob.glob(os.path.join(folder_syt, "*.tif"))):
    construct = file.split("/")[-1].split("_")[3][:-2]
    antibody = file.split("/")[-1].split("_")[1]
    mouse_id = file.split("/")[-1].split("_")[3][-2:]
    slice_id = file.split("/")[-1].split("_")[2].strip("Slice")
    
    
    filename_new = mouse_id + "_" + slice_id + "_"  +construct + "_" + antibody + ".tif"
    print(filename_new)
    os.rename(file, os.path.join(folder_export, filename_new))
    
    
    

for file in tqdm(glob.glob(os.path.join(folder_psd_feb, "*.tif"))):
    construct = file.split("/")[-1].split("_")[3][:-2]
    antibody = file.split("/")[-1].split("_")[1]
    mouse_id = file.split("/")[-1].split("_")[3][-2:]
    slice_id = file.split("/")[-1].split("_")[2].strip("Slice")
    
    
    filename_new = mouse_id + "_" + slice_id + "_"  +construct + "_" + antibody + ".tif"
    print(filename_new)
    os.rename(file, os.path.join(folder_export, filename_new))
    
    
    
    
    
    
