# 3D_GFAP_Puromycin_Analysis_Workflow
This repository contains a Fiji workflow for processing  `.lif` files and performing 3D analysis of nuclei, virtual cyto, and cell.  The workflow is designed for multichannel 3D images containing DAPI, GFAP, and puromycin. Nuclear instance masks are generated with Cellpose and are then matched to the corresponding images during the 3D analysis.

## Workflow overview

```text
Leica .lif files
      |
      v
1. splitFiles.groovy
   Export every LIF series as an individual TIFF hyperstack
      |
      v
Exported 3D TIFF images
      |
      +------------------------------+
      |                              |
      v                              |
2. Nuclear segmentation             |
   Cellpose or equivalent            |
      |                              |
      v                              v
Label masks: *_cp_masks.tif     Raw multichannel TIFFs
      |                              |
      +---------------+--------------+
                      |
                      v
3. 3DAnalysis_4markers_cell_nuclei_cyto_accelerated.groovy
   - Nuclear measurements
   - Calibrated 3D virtual-cytoplasm expansion
   - Whole-cell measurements
   - GFAP and puromycin classification
   - Per-nucleus and per-image CSV export
```

> **Important:** The repository contains the Groovy scripts for TIFF extraction and 3D measurement. Nuclear segmentation is an intermediate step and must be completed before running the 3D analysis.

## Repository contents

```text
.
├── splitFiles.groovy
├── 3DAnalysis_4markers_cell_nuclei_cyto_accelerated.groovy
└── README.md
```

### `splitFiles.groovy`

Opens all image files found in the configured input directory using Bio-Formats. For each Leica LIF file, all series are opened as composite `XYCZT` hyperstacks and saved individually as TIFF images.

### `3DAnalysis_4markers_cell_nuclei_cyto_accelerated.groovy`

Matches each raw TIFF image to a labelled nuclear mask, separates the fluorescence channels, creates calibrated 3D virtual-cytoplasm regions around each nucleus, measures GFAP and puromycin intensity, classifies objects, and writes results at nucleus and image level.

## Requirements

Run the scripts in **Fiji/ImageJ** with the following components available:

- Fiji/ImageJ
- Groovy script support
- Bio-Formats
- 3D ImageJ Suite / MCIB3D
- A compatible Java installation supplied with Fiji
- Cellpose, or another method that produces labelled 3D nuclear masks

The nuclear masks must be label images in which:

- background pixels/voxels have value `0`;
- each nucleus has its own positive integer label;
- mask dimensions match the corresponding raw image in X, Y, and Z.

## Recommended directory structure

```text
project/
├── lif/
│   ├── experiment_01.lif
│   └── experiment_02.lif
├── images/
│   ├── image_01.tif
│   └── image_02.tif
├── nuclei/
│   ├── image_01_cp_masks.tif
│   └── image_02_cp_masks.tif
└── csv/
```

## Step 1: Export LIF series as TIFF images

Open `splitFiles.groovy` in the Fiji Script Editor and edit these variables:

```groovy
def inputFiles = new File("/path/to/project/lif")
def outputDir = new File("/path/to/project/images")
```

The script:

1. lists the files in `inputFiles`;
2. opens each file with Bio-Formats;
3. loads all series as `XYCZT` hyperstacks;
4. preserves the multichannel composite structure;
5. saves every series as an individual TIFF in `outputDir`.

Create the output directory before running the script. Review the exported filenames because they are used to match raw images and nuclear masks later.

## Step 2: Generate nuclear masks

Segment nuclei from the DAPI channel using Cellpose or an equivalent instance-segmentation method.

For every raw image, save its mask in the `nuclei` directory using the following naming convention:

```text
Raw image:     image_01.tif
Nuclear mask:  image_01_cp_masks.tif
```

For `.tiff` inputs, the expected mask name is generated in the same way:

```text
Raw image:     image_01.tiff
Nuclear mask:  image_01_cp_masks.tif
```

Do not save binary masks if several nuclei are present. The analysis requires an instance-labelled mask with a different integer value for each nucleus.

## Step 3: Run the 3D analysis

Open `3DAnalysis_4markers_cell_nuclei_cyto_accelerated.groovy` in the Fiji Script Editor.

### 3.1 Configure directories

```groovy
final File inputFilesNuclei = new File(
        "/path/to/project/nuclei"
)

final File inputFilesRaw = new File(
        "/path/to/project/images"
)

final File outputDir = new File(
        "/path/to/project/csv"
)
```

The script checks that the raw-image and nuclear-mask directories exist. It creates the CSV output directory when possible.

### 3.2 Configure channel indices

```groovy
final int DAPI_CHANNEL = 3
final int PUROMYCIN_CHANNEL = 1
final int GFAP_CHANNEL = 0
```

Channel indices are **zero-based**:

- Fiji channel 1 = index `0`
- Fiji channel 2 = index `1`
- Fiji channel 3 = index `2`
- Fiji channel 4 = index `3`

Adjust these values to match the acquisition channel order. The current configuration assumes GFAP in Fiji channel 1, puromycin in Fiji channel 2, and DAPI in Fiji channel 4.

Although the script filename includes `4markers`, the current implementation explicitly measures GFAP and puromycin and uses DAPI-associated nuclear masks. Any additional channel is not measured unless the script is extended.

### 3.3 Configure the virtual cytoplasm

```groovy
final float CYTO_DILATION_XY_PX = 5.0f
final boolean EXCLUDE_ALL_NUCLEI_FROM_CYTO = true
```

`CYTO_DILATION_XY_PX` defines the XY expansion radius around each nucleus. The script converts this value into a physical distance using the image calibration and calculates the corresponding Y and Z radii. This produces a calibrated anisotropic 3D dilation.

When `EXCLUDE_ALL_NUCLEI_FROM_CYTO` is `true`, voxels belonging to any segmented nucleus are removed from every virtual-cytoplasm region.

The virtual cytoplasm is defined as:

```text
virtual cytoplasm = dilated nuclear object - original nucleus
```

The complete virtual cell is calculated without creating an additional mask:

```text
cell = nucleus + virtual cytoplasm
```

### 3.4 Configure positivity thresholds

```groovy
final double GFAP_THRESHOLD_K = 0.5d
final double PURO_THRESHOLD_K = 0.5d
```

For each image and marker, the threshold is calculated from the distribution of whole-cell mean intensities:

```text
threshold = mean(cell means) + K x standard deviation(cell means)
```

The default value is `K = 0.5` for both GFAP and puromycin.

### 3.5 Choose the compartment used for classification

```groovy
final boolean CLASSIFY_USING_CELL_MEAN = false
```

- `false`: thresholds are estimated from whole-cell means, but positivity is assigned using nuclear mean intensity.
- `true`: thresholds are estimated from whole-cell means and positivity is assigned using whole-cell mean intensity.

This distinction is important when interpreting the positivity columns.

### 3.6 Configure progressive saving

```groovy
final boolean SAVE_AFTER_EACH_IMAGE = true
```

When enabled, both result tables are saved after each successfully processed image. This reduces the risk of losing completed measurements if a later image fails.

## Measurements

For each nucleus, the script calculates:

- nucleus ID and centroid coordinates;
- nuclear Z range;
- nuclear volume in voxels;
- GFAP mean and integrated intensity in the nucleus;
- puromycin mean and integrated intensity in the nucleus;
- virtual-cytoplasm volume;
- GFAP mean and integrated intensity in the virtual cytoplasm;
- puromycin mean and integrated intensity in the virtual cytoplasm;
- virtual-cell volume;
- GFAP mean and integrated intensity in the virtual cell;
- puromycin mean and integrated intensity in the virtual cell;
- marker-specific thresholds;
- GFAP and puromycin positivity;
- combined DAPI/GFAP/puromycin status;
- virtual-cytoplasm dilation parameters.

DAPI positivity is currently assigned as `true` for every segmented nuclear object. Therefore, `DAPI+` means that the object originates from the supplied nuclear label mask, not that a separate DAPI intensity threshold was applied.

## Outputs

The analysis produces two CSV files.

### Per-nucleus output

```text
3DAnalysis_perNucleus_GFAP_Puromycin_virtualCyto3D.csv
```

Contains one row per segmented nucleus, including spatial coordinates, compartment-specific measurements, thresholds, and classification status.

### Per-image output

```text
3DAnalysis_perImage_GFAP_Puromycin_virtualCyto3D.csv
```

Contains one row per successfully analysed image, including:

- total number of nuclei;
- GFAP and puromycin thresholds;
- number of GFAP-positive objects;
- number of GFAP/puromycin double-positive objects;
- measurements for DAPI+/GFAP+ objects;
- measurements for DAPI+/GFAP+/puromycin+ objects;
- whole-cell summaries;
- virtual-cytoplasm summaries;
- dilation parameters.

Mean columns ending in `_Avg` are averages across the relevant objects. Integrated-intensity columns ending in `_Sum` are sums across the relevant objects.

## Input matching and validation

The analysis accepts `.tif` and `.tiff` raw images. It ignores `.DS_Store` files and files containing `_seg.npy`.

For every raw image, the script automatically searches for a mask with the `_cp_masks.tif` suffix. Images are skipped if:

- the matching mask is missing;
- the raw image or mask cannot be opened;
- requested channel indices are unavailable;
- raw channels and masks have different XYZ dimensions;
- the mask contains no labelled nuclei.

Errors and progress information are written to the Fiji Log window. Processing continues with the next image whenever possible.

## Running the scripts in Fiji

1. Start Fiji.
2. Open **File > New > Script**.
3. Select **Language > Groovy**.
4. Open or paste the required script.
5. Edit the directories and configuration variables at the top.
6. Click **Run**.
7. Monitor the Fiji Log window for progress, warnings, and output paths.

Run the scripts in this order:

1. `splitFiles.groovy`
2. nuclear segmentation with Cellpose
3. `3DAnalysis_4markers_cell_nuclei_cyto_accelerated.groovy`

## Quality-control recommendations

Before analysing the full dataset:

1. test the workflow on one or two representative images;
2. verify the raw-to-mask filename matching;
3. confirm channel assignments in Fiji;
4. check that image calibration is correct in X, Y, and Z;
5. overlay the nuclear masks on the raw images;
6. inspect the selected virtual-cytoplasm radius;
7. validate GFAP and puromycin classifications visually;
8. review the per-nucleus CSV for implausible volumes or intensities.

## Known assumptions and limitations

- The workflow assumes one labelled nuclear mask per raw TIFF image.
- The mask and raw image must have identical XYZ dimensions.
- Thresholds are calculated independently for each image, so positivity is image-relative.
- DAPI is not independently thresholded.
- Virtual cytoplasm is an approximation based on nuclear dilation and is not a membrane-based cell segmentation.
- Virtual cytoplasm excludes nuclear voxels when configured, but neighbouring cytoplasmic expansions may still overlap because unique ownership between nearby cells is not explicitly resolved.
- Measurements are exported in voxel and raw-intensity units. Their biological interpretation depends on acquisition consistency and image calibration.
- The current script measures GFAP and puromycin only.



