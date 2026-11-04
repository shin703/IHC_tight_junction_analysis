# Mica → DAPI / ZO-1 Junction Analysis GUI

A single desktop GUI that connects two steps of a Caco-2 / ZO-1 image-analysis workflow:

1. **Mica multipage TIFF → channel-wise Maximum Intensity Projection**
2. **DAPI StarDist nuclei review + ZO-1 junction analysis**

The application uses one window with two tabs, so preprocessing and analysis can be performed from the same GUI. The DAPI / ZO-1 filename postfixes are editable in the analysis tab, so datasets do not need to use fixed `_DAPI` / `_ZO1` names.

---

## Japanese / English GUI

The integrated GUI can be switched between **日本語** and **English** without restarting the application.

Use the language buttons at the top of the main window:

```text
Language / 言語:  [日本語] [English]
```

The selected language is applied to:

- Mica MAX Projection tab
- DAPI / ZO-1 Analysis tab
- StarDist nuclei-review window
- ZO-1 preprocessing QC window
- dialog messages and file-dialog titles
- analysis status labels and table headings

The selected language is remembered for the next launch. Scientific names, file postfixes, sample names and output filenames are not translated.

---

## Workflow

```text
Mica multipage Z-stack TIFF
        |
        v
[1. Mica MAX Projection]
        |
        |  Ch1 -> <sample>_DAPI.tif
        |  Ch2 -> <sample>_ZO1.tif
        v
Output folder
        |
        |  "Open output folder in analysis tab"
        v
[2. DAPI / ZO-1 Analysis]
        |
        +--> StarDist nuclei detection
        +--> manual nuclei review / correction
        +--> ZO-1 preprocessing
        +--> nucleus-pair junction analysis
        +--> sample / group CSV summaries
```

---

## Tab 1 — Mica MAX Projection

The Mica TIFF is assumed to contain interleaved channel pages:

```text
2 channels:
Z1-Ch1, Z1-Ch2, Z2-Ch1, Z2-Ch2, ...

3 channels:
Z1-Ch1, Z1-Ch2, Z1-Ch3, Z2-Ch1, ...

4 channels:
Z1-Ch1, Z1-Ch2, Z1-Ch3, Z1-Ch4, Z2-Ch1, ...
```

For every channel, the program creates a **Maximum Intensity Projection along Z**.

### Default postfixes

For a 2-channel DAPI / ZO-1 dataset:

```text
Ch1: _DAPI
Ch2: _ZO1
```

For input:

```text
con1.tif
```

the output becomes:

```text
con1_DAPI.tif
con1_ZO1.tif
```

These filenames are directly compatible with the analysis tab.

After projection, click:

```text
出力フォルダを解析タブで開く
```

to move directly to Tab 2 and scan the generated images.

---

## Tab 2 — DAPI / ZO-1 Analysis

The analysis tab pairs DAPI and ZO-1 TIFF files using **user-editable postfixes**.

Default values are:

```text
DAPI postfix: _DAPI
ZO-1 postfix: _ZO1
```

With these defaults, the expected filenames are:

```text
<sample>_DAPI.tif
<sample>_ZO1.tif
```

Example:

```text
con1_DAPI.tif
con1_ZO1.tif
con2_DAPI.tif
con2_ZO1.tif
```

You can change the postfixes directly in the analysis tab. For example:

```text
DAPI postfix: _DAPI_MAX
ZO-1 postfix: _ZO1_MAX
```

will pair:

```text
con1_DAPI_MAX.tif
con1_ZO1_MAX.tif
```

Use **再スキャン / Rescan** after changing the postfixes. The values are saved in the global settings and restored at the next launch.

When you move from the MAX Projection tab using **出力フォルダを解析タブで開く**, the Ch1 and Ch2 postfix values are copied automatically into the DAPI and ZO-1 postfix fields. You can also copy them manually using **MAX設定から取得**.

### DAPI nuclei

Initial nuclei detection uses:

```text
StarDist 2D_versatile_fluo
```

Detected nuclei can then be manually reviewed and edited.

Main review controls include:

```text
Left drag          Create / add to a nucleus
Ctrl + left drag   Eraser
Double click       Select / deselect a nucleus
Right drag         Pan
Mouse wheel        Zoom
Delete             Delete selected nucleus
Ctrl + Z           Undo
[ / ]              Brush size
S                   Save
N                   Save and move to next sample
```

---

## ZO-1 junction analysis

The workflow includes:

```text
Raw ZO-1
  -> rolling-ball-style background subtraction
  -> white top-hat
  -> fixed raw-scale threshold
  -> small-object removal
  -> skeletonization
```

The analysis then creates candidate neighboring nucleus pairs and evaluates ZO-1 signal within the pair corridor.

Main outputs include junction coverage and lane-based ZO-1 intensity metrics.

---

## Output structure

For flat-file input:

```text
ROOT/
├── con1_DAPI.tif
├── con1_ZO1.tif
├── con2_DAPI.tif
├── con2_ZO1.tif
└── ZO1_junction_analysis/
    ├── con1/
    └── con2/
```

Per-sample outputs can include:

```text
dapi_nuclei_annotation.npz
dapi_nuclei_labels.tif
dapi_nuclei_binary.tif
dapi_nuclei_review_overlay.png

junction_pair_analysis.csv
junction_analysis_summary.csv
junction_analysis_overlay.png

zo1_background_subtracted.tif
zo1_white_tophat.tif
zo1_binary_mask.tif
zo1_skeleton.tif

zo1_raw_plus_skeleton_overlay.png
zo1_tophat_plus_skeleton_overlay.png
zo1_skeleton_plus_nuclei_overlay.png
```

Root-level summary files include:

```text
ZO1_junction_summary.csv
ZO1_junction_summary_by_group.csv
zo1_sample_groups.json
```

---

## Installation

Python environment:

```bash
pip install -r requirements.txt
```

Run:

```bash
python mica_zo1_junction_gui.py
```

---

## Main dependencies

```text
numpy
opencv-python
Pillow
tifffile
scipy
scikit-image
tensorflow
stardist
csbdeep
imagecodecs
```

The code uses `numpy==1.26.4` to maintain compatibility with the current analysis environment.

---

## Windows EXE

A PyInstaller configuration is included:

```text
Mica_ZO1_Junction_GUI.spec
build_exe.bat
```

On Windows, double-click:

```text
build_exe.bat
```

The executable is created under:

```text
dist/Mica_ZO1_Junction_GUI.exe
```

### Note about EXE size

Because the application includes TensorFlow, StarDist and scientific Python libraries, the executable can be large. Building from a clean virtual environment is recommended.

For GitHub distribution, keep `dist/` and `.exe` files out of the repository and attach the final EXE to **GitHub Releases**.

---

## Repository structure

```text
Mica-ZO1-Junction-GUI/
├── mica_zo1_junction_gui.py
├── requirements.txt
├── Mica_ZO1_Junction_GUI.spec
├── build_exe.bat
├── .gitignore
└── README.md
```

---

## Important analysis notes

- The MAX-projection tab does not alter channel pixel values except for the Z-axis maximum operation.
- DAPI / ZO-1 image dimensions must match for analysis.
- Batch StarDist results are saved separately from manually reviewed annotations so an automatic draft is not confused with a human-reviewed result.
- Group statistics are summarized at the **sample/image level**, rather than treating individual nucleus pairs as independent biological replicates.
- Analysis settings are shared across parent folders using the existing global-settings behavior.

---

## Recommended first test

Before processing a full experiment:

1. Select one Mica TIFF.
2. Use 2 channels with `_DAPI` and `_ZO1`.
3. Create MAX projections.
4. Click **出力フォルダを解析タブで開く**.
5. Confirm that the sample appears as a DAPI / ZO-1 pair.
6. Run StarDist on the sample.
7. Review and save nuclei.
8. Run ZO-1 junction analysis.
9. Check the overlay images and CSV output.

This verifies the complete preprocessing-to-analysis pipeline on one sample before batch processing.
