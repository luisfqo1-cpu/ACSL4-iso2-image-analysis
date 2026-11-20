# ACSL4-iso2-image-analysis

Image analysis pipeline for: *A druggable moonlighting function of ACSL4 couples membrane arachidonic acid metabolism to endocannabinoid transport*

## ACSL4 (T4) immunofluorescence quantification in mouse brain slices

Python scripts we used to quantify antibody signal in confocal z-stacks of mouse brain slices. They measure **ACSL4** expression (called **T4** in the code and file names) and compare it with other markers: GFAP, MAGL, synaptotagmin (SYT), PSD95, CB1 and NeuN. Two genetic manipulations are analysed:

| Experiment | Control label in files | Treatment label in files | Label in figures |
|---|---|---|---|
| shRNA knock-down | `CNTL` | `shRNA` | Scrambled vs. **ACSL4 Iso 2** |
| Conditional knock-out (Sox2-Cre) | `CNTL` | `KO` | WT vs. **ACSL4 Cat KO** |

`2CNTL` marks secondary-antibody-only controls, which are used to estimate background staining.

---

## Repository contents

| Script | Purpose |
|---|---|
| `renames_shRNA.py` | Collects the deconvolved `.tif` stacks from several acquisition folders and renames them to a single scheme, `<mouse>_<slice>_<construct>_<antibody>.tif`, in one export folder. |
| `t4_shrna_eval.py` | Main analysis for the shRNA experiment: segmentation, per-cell quantification, secondary-control comparison, GFP/GFAP percentages, MAGL/SYT/PSD95 intensity and colocalization, CB1–MAGL colocalization and a summary plot of all colocalizations. |
| `t4_KO_eval.py` | The same per-cell pipeline applied to the ACSL4 conditional knock-out (WT vs. KO) with NeuN staining. |
| `t4_neun_cb1.py` | ACSL4 vs. CB1 in NeuN⁺ cells, split by hippocampal region (e.g. CA3), using precomputed masks. |
| `circularity.py` | Cell-shape (circularity) analysis from polygon outlines exported to CSV (`circ.csv`, `circ_ko.csv`), plus 5×5 galleries of example cell outlines for each condition. |
| `imagetools/` | Our helper package (`imagetools.util`) and napari plugin for interactive ROI-based quantification. See [`imagetools/README.md`](imagetools/README.md). A newer, extended version of the plugin is maintained in [luisfqo1-cpu/plugin2025](https://github.com/luisfqo1-cpu/plugin2025). |

---

## Input data

### Image format
Each `.tif` is a 4-D confocal stack with shape **(z, y, x, channel)** and this channel order:

| Index | Variable | Content |
|---|---|---|
| 0 | `cy5` | Co-stained marker (GFAP, MAGL, SYT, PSD95, CB1, …) |
| 1 | `t4` | ACSL4 |
| 2 | `inj` | Injection / reporter channel (GFP from the viral construct, or NeuN) |
| 3 | `dapi` | Nuclei |

### File naming
- After `renames_shRNA.py`: `<mouse_id>_<slice_id>_<construct>_<antibody>.tif`, e.g. `20_4_CNTL_T4.tif`, `13_1_shRNA_T4-MAGL.tif`
- KO dataset: `<prefix>_<construct><mouse_id(3 digits)>_<antibody>_Slice<slice_id>_….tif`
- NeuN/CB1 dataset: `<region>_…_…_<mouse_id>_….tif`, with masks saved as `<mouse_id>_<region>_masks_neun.npy` / `_masks_dapi.npy`

### Circularity CSVs
`circ.csv` / `circ_ko.csv` are exported from napari with the **polygons** widget of [plugin2025](https://github.com/luisfqo1-cpu/plugin2025), after drawing cell outlines on images named `<mouse>_<slice>_<condition>_<antibody>`. They need a `condition` column and a `polygon_coords` column containing `;`-separated vertices, each a comma-separated list whose last two values are the (row, col) coordinates.

---

## Analysis pipeline

1. **Maximum-intensity projection (MIP)** of each channel along z.
2. **Segmentation** with [Cellpose](https://github.com/MouseLand/cellpose) (`cyto2`, diameter = 200 px, flow threshold = 0.4):
   - `masks_inj`: transduced (GFP⁺) or NeuN⁺ cells
   - `masks_dapi`: all nuclei
   Masks are saved as `.npy` in a `masks/` subfolder so they can be reused.
3. **Thresholding** of the ACSL4 and Cy5 MIPs with the triangle method (`skimage.filters.threshold_triangle`).
4. **Per-cell measurements**, computed over every z-slice inside each mask:

   | Column | Meaning |
   |---|---|
   | `s2n_on` | Mean ACSL4 intensity of pixels **above** threshold ("on pixels") |
   | `s2n_off` | Mean ACSL4 intensity of pixels **below** threshold ("off pixels") |
   | `on_pixels` / `off_pixels` | Summed on/off intensity divided by cell area |
   | `cy5_pixels` | Median Cy5 intensity normalised to the background (minimum Cy5 value outside the cells) |
   | `colocs_mean` / `colocs_max` | Pearson correlation between ACSL4 and Cy5 in the cell, mean / max over z-slices |
   | `circularity` | Cell area / area of the circle whose radius is the largest centroid-to-pixel distance |
   | `background` | Minimum ACSL4 MIP value outside the cells |

5. **DAPI-based cell classification** (`results_dapi.csv`):
   - `inj_likelihood`: fraction of the nucleus overlapping a GFP⁺ mask. A cell is **GFP⁺** if > 0.5.
   - `cy5_likelihood`: fraction of nuclear pixels above the Cy5 threshold. A cell is **GFAP⁺** if > 0.2.
6. **Statistics and plots**: two-sided Mann–Whitney U tests through `statannot`. For the main comparisons, per-cell values are first averaged per *(mouse, slice)*, so each data point is one slice. N = number of mice, n = number of slices. Figures are written as `.pdf` / `.eps`.

Results are appended to `results_inj*.csv` / `results_dapi.csv` after every image. Images that already appear in those files are skipped, so an interrupted run can be resumed.

---

## Output figures (selection)

| File | Content |
|---|---|
| `onpixels_seccontrol.eps`, `offpixels_seccontrol.eps` | Each antibody vs. secondary-only control |
| `t4_on.pdf`, `t4_off.pdf` | ACSL4 on/off-pixel intensity, Scrambled vs. ACSL4 Iso 2 |
| `t4_KO_on.eps`, `t4_KO_off.eps` | The same comparison for WT vs. ACSL4 Cat KO |
| `GFP_percentage.pdf` | Percentage of GFP⁺ cells per mouse |
| `gfap_t4_on.pdf`, `gfap_t4_off.pdf` | ACSL4 intensity by cell type (GFAP⁺ / GFP⁺ / double negative) |
| `magl.pdf`, `syt.pdf`, `psd95.pdf` (+ `_coloc`, `_scatter`) | Marker intensity, ACSL4 colocalization and per-cell scatter plots |
| `cb1_MAGL_scatter.pdf`, `cb1_scatter.eps` | CB1 vs. MAGL and CB1 vs. ACSL4 |
| `all_colocalizations.pdf` | Mean colocalization for every marker pair |
| `circularity*.eps`, `*_cell_examples*.eps` | Cell-shape comparison and example outlines |

---

## Requirements

- Python 3.8 – 3.11
- All dependencies, with the version ranges the scripts need, are listed in [`requirements.txt`](requirements.txt). This includes the `imagetools` helper package from this repository, so run the command from the repository root:

```bash
pip install -r requirements.txt
```

- To also use the napari plugin, install `pip install "./imagetools[napari]"` (see [`imagetools/README.md`](imagetools/README.md)).

The scripts run Cellpose on the CPU (`gpu=False`). Set `gpu=True` if a CUDA GPU is available, which is much faster.

---

## Usage

The scripts are written as analysis notebooks: run them top to bottom, or cell by cell in Spyder/VS Code.

1. Open each script and set the folders in the **User settings** block at the top (`folder`, `folder_fig`, `folder_export`, …) to your own locations. The defaults are placeholders such as `"path/to/figures"`.
2. `renames_shRNA.py` **moves** files (`os.rename`) instead of copying them. Back up the raw data before running it.
3. Run `t4_shrna_eval.py` / `t4_KO_eval.py`. If `results_inj*.csv` / `results_dapi.csv` already exist in the data folder they are loaded and images already in them are skipped. Otherwise the analysis starts from scratch. The `masks/` subfolder is created automatically.
4. Run `t4_neun_cb1.py` before the last section of `t4_shrna_eval.py`, which reads its `results_inj.csv` for the ACSL4–CB1 colocalization.
5. Run `circularity.py` once `circ.csv` / `circ_ko.csv` have been placed in `folder_data`.

---

## Citation

If you use this code, please cite: *A druggable moonlighting function of ACSL4 couples membrane arachidonic acid metabolism to endocannabinoid transport* (full reference to be added on publication).

## License

This project is released under the MIT License (see also [`LICENSE`](LICENSE)). This includes the `imagetools/` package (see [`imagetools/LICENSE`](imagetools/LICENSE)).

```
MIT License

Copyright (c) 2026 University of Bern (Gertsch lab)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
