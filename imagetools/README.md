# imagetools

Helper functions and a [napari](https://napari.org) plugin we used to quantify antibody staining in brain slices. The analysis scripts in the repository root use `imagetools.util` (for example `circularity.py`). The napari plugin offers the same tools interactively.

> **Newer version:** the plugin is developed further in [luisfqo1-cpu/plugin2025](https://github.com/luisfqo1-cpu/plugin2025), which adds widgets for two-channel and time-series pixel intensity, colocalization, polygon ROIs, zebrafish images and denoising, plus a `.czi` reader. The copy here is the version used for the manuscript and is kept unchanged so the analysis can be reproduced. Its `util.py` matches plugin2025 as of commit [`191162f`](https://github.com/luisfqo1-cpu/plugin2025/tree/191162f8e6c3686c7f55bd0cb728d34196a0cd9c) (27 Oct 2025).

## Contents

| Path | What it is |
|---|---|
| `src/imagetools/util.py` | Helper functions: ROI cropping (`calculate_bounding_rectangle`, `generate_polygon_mask`, `crop_with_polygon_mask`, `fit_dim`), segmentation (`return_mask` with StarDist or Cellpose, chunked segmentation of large images), per-cell measurements (`calculate_props`, `calculate_rate`, `return_positive_negative_cells`), image registration (`calculate_lag`, `register_image`) and atlas slice matching (`match_atlas`, BrainGlobe) |
| `src/imagetools_napari/` | napari plugin (**ImageTools**): a reader for `.npy` / `.nd2` files and the widgets listed below |
| `scripts/nd2_npy_converter.py` | Converts every `.nd2` file in a folder (recursively) to `.npy` |
| `examples/brainreg.json` | Example brainreg registration settings (example mouse atlas, 100 µm) |
| `examples/target_rate_example.csv` | Example output of *Calculate Target Rate* for hippocampal ROIs |

### napari widgets

| Widget | What it does |
|---|---|
| Register Images | Aligns a target image to a reference by cross-correlating binarised crops inside the drawn ROIs, then adds the shifted image |
| Adjust Contrast | Sets an image's contrast limits to `[0, max]` |
| Slider | Interactive threshold: shows the selected image above the slider value as a new layer |
| Stardist Segmentation | Segments cells inside rectangular ROIs with StarDist or Cellpose (`method` = `stardist` / `cellpose`) |
| Shapes2Labels | Converts a Shapes layer into a Labels layer |
| Calculate Target Rate | For each ROI, the fraction of reference cells (labels) overlapping the target signal by more than `coloc_thr` |
| Calculate Target Intensity | For each ROI, the mean (or another `regionprops` property) of the target image over the reference cells |
| Calculate Antibody Intensity | Splits the reference cells in each ROI into target-positive and target-negative, then measures antibody intensity in each group |
| Draw Contours | Draws ROI outlines filled with the rates from a saved CSV |

The measurement widgets save their results as `<save>.csv` in the folder napari was started from. Running them again with the same name merges the new rows into the existing file.

## Installation

From the repository root:

```bash
pip install ./imagetools                 # helper functions only (enough for the analysis scripts)
pip install "./imagetools[napari]"       # + napari plugin and .nd2 reader
pip install "./imagetools[all]"          # + StarDist/Cellpose segmentation and BrainGlobe atlas matching
```

Segmentation and atlas matching are optional extras, so `imagetools.util` imports without them. `return_mask` and `match_atlas` need `[segmentation]` and `[atlas]` respectively.

After installing with `[napari]`, the widgets appear under **Plugins → ImageTools** in napari.

## Tests

```bash
pip install "./imagetools[testing]"
pytest imagetools/src/imagetools_napari/_tests
```

## License

`imagetools` is distributed under the MIT license (see [`LICENSE`](LICENSE)), like the rest of this repository.
