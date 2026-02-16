<div align="center">

# StatCraft

**Second-Level Neuroimaging Analysis Tool**

[Features](#features) | [Installation](#installation) | [Quick Start](#quick-start) | [Configuration](#configuration) | [Outputs](#outputs) | [License](#license)  | [Citation](#citation) | [Acknowledgments](#acknowledgments)

</div>

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL%203.0-blue.svg)](https://www.gnu.org/licenses/agpl-3.0.html)

CLI tool for second-level neuroimaging analysis, supporting group-level comparisons, method comparisons, and statistical inference on brain images (e.g., fMRI, PET) or connectivity matrices. Uses flexible pattern matching for file discovery, leverages Nilearn, and produces reproducible, interpretable results with minimal user input.

## Features

- **Multiple Analysis Types**
  - One-sample t-tests
  - Two-sample t-tests (group comparisons)
  - Paired t-tests (within-subject comparisons)
  - General Linear Model (GLM) with custom design matrices

- **Flexible File Discovery**
  - Pattern-based file matching with glob patterns
  - Participant filtering by subject ID
  - Works with derivatives from any first-level analysis
  - Supports BIDS-like and custom file naming

- **Statistical Inference**
  - Uncorrected thresholding
  - FDR (False Discovery Rate) correction
  - FWER (Family-Wise Error Rate) via Bonferroni
  - Permutation-based FWER correction

- **Anatomical Annotation of Clusters**
  - Harvard-Oxford atlas (default)
  - AAL, Destrieux, Schaefer atlases
  - Custom atlas support

- **Comprehensive Reporting**
  - HTML reports with visualizations
  - Design matrix plots
  - Activation maps (thresholded and unthresholded)
  - Cluster tables with anatomical labels

## Installation

```bash
git clone https://github.com/arovai/statCraft.git
cd statCraft
pip install -e .
```

## Quick Start

### Command Line Interface

#### Usage Patterns

StatCraft uses flexible pattern matching for file discovery. The idea is that you point to one or more directories and define filenaming patterns using wildcards, and this define the data on which the second-level analysis will be performed.

This allows the user to use this tool in a variety of first-level tools without any particular contraints on the output structure.

Here are the typical usage, depending on the analysis type:

**One-sample t-tests**

```bash
statcraft <INPUD_DIR> <OUTPUT_DIR> --analysis-type one-sample --pattern "PATTERN"
```
The pipeline will search for files matching `PATTERN` in the `<INPUT_DIR>` (exploring subfolders), perform the one-sample t-test, and save the results in `<OUTPUT_DIR>`.

**Two-sample t-tests**

The logic is similar as above, except that now one must specify two patterns to define the two groups to compare. We can also assign these groups custom names to ease the output filenaming:
```bash
statcraft <INPUD_DIR> <OUTPUT_DIR> --analysis-type two-samples --patterns "Group1=PATTERN1 Group2=PATTERN2"
```
The strings `Group1` and `Group2` are artitraty and are used in the output filenaming and in the analysis report.

**Paired t-tests**

This case is similar to the two-samples case except that one must provide a key to pair the data across the two groups. We do this by using the `--pair-by` argument:
```bash
statcraft <INPUD_DIR> <OUTPUT_DIR> --analysis-type paired --patterns "Group1=PATTERN1 Group2=PATTERN2" --pair-by "sub"
```
The pairing is done by using the value given in `--pair-by` using a BIDS-like `key-value` structure: for instance, for `--pair-by "sub"`, files with `sub-001_*` will be paired together.
*Note:* `--pair-by` is optional, it's default value being `sub`.

**General Linear Model (GLM)**

For GLM analysis one must provide a design matrix. This is done by passing the `--participants-file`, which follows the structure of the `participants.tsv` file in a BIDS directory:

```
participant_id  sex age (other columns)
sub-001 F 42  ...
sub-CTL1 F 93  ...
sub-abcd M 17  ...
```
By default, a design matrix is build using all the columns of this file. If only a subset of columns should be included, this can be achieved using the `--regressors` option. Moreover, columns that should be treated as categorical can be specified with `--categorical-regressors` (dummy-coded in the analysis). Finally, the constrast to compute is defined using the `--contrasts` argument.

Here is a complete example featuring this functionality:
```bash
statcraft <INPUT_DIR> <OUTPUT_DIR> --analysis-type glm --participant-files <PATH_TO_PARTICIPANTS.TSV> --regressors sex age IQ --categorical-regressors sex --contrasts age M-F
```
In this example: colums `sex`, `age` and `IQ` are used to build the design matrix, with `sex` being treated as a categorical variable, and two contrasts are computed: the effet of `age` as well as the difference between `M` and `F` (assuming those are the labels in original `participants.tsv` file).

*Notes*: 
- By default, non-categorical variables are z-scored. To control this, one can use `--no-standardize-regressors`, e.g. `--no-standardize-regressors IQ`
- An intercept is also added in the design matrix. Contrasts using the intercept can be build using the word `mean`, e.g. `--contrats mean`. To remove the interecept, use `--no-intercept` option.
- Non-trivial constrasts can be build using simple operands, e.g. `--contrasts 0.5*M+0.5*F-mean`

**File Discovery**

All file discovery uses glob patterns specified with `--pattern`. Include task, session, space, and other filters directly in your pattern:
```bash
# Pattern with embedded filters
--pattern '*task-rest*ses-01*space-MNI152*stat-effect*.nii.gz'
```

#### Examples

```bash
# One-sample t-test with pattern matching
statcraft /data/dataset /data/output group \
    --derivatives /data/derivatives/fmriprep \
    --analysis-type one-sample \
    --pattern '*task-nback*space-MNI152*stat-effect*.nii.gz'

# One-sample with specific participants
statcraft /data/dataset /data/output group \
    --derivatives /data/derivatives/fmriprep \
    --analysis-type one-sample \
    --participant-label 01 02 03 \
    --pattern '*task-rest*.nii.gz'

# Derivatives-only with explicit participants file
statcraft /data/derivatives/fmriprep /data/output group \
    --participants-file /data/rawdata/participants.tsv \
    --analysis-type one-sample \
    --pattern '*task-nback*stat-effect*.nii.gz'

# Two-sample t-test (group comparison)
statcraft /data/dataset /data/output group \
    --derivatives /data/derivatives/fmriprep \
    --analysis-type two-sample \
    --group-column group \
    --pattern '*stat-effect*.nii.gz'

# Paired t-test with sample patterns
statcraft /data/dataset /data/output group \
    --derivatives /data/derivatives/fmriprep \
    --analysis-type paired \
    --pair-by sub \
    --patterns 'pre=*ses-pre*.nii.gz post=*ses-post*.nii.gz' \
    --contrast 'post-pre'

# Method comparison (two-sample with patterns)
statcraft /data/dataset /data/output group \
    --derivatives /data/derivatives/cvrmap \
    --analysis-type two-sample \
    --patterns 'GS=**/*GS*cvr*.nii.gz SSS=**/*SSS*cvr*.nii.gz' \
    --contrast 'GS-SSS'

# Using a configuration file
statcraft /data/dataset /data/output group \
    --derivatives /data/derivatives/fmriprep \
    --config config.yaml
```

## Configuration

StatCraft can be configured via command-line arguments, configuration files (YAML/JSON), or Python dictionaries.

### Generate Default Configuration

```bash
statcraft --init-config config.yaml
```

### Example Configuration (YAML)

```yaml
# Analysis type
analysis_type: glm

# Participant filtering (optional)
participant_label: null  # Or: ["01", "02", "03"]

# File pattern with embedded filters
# Include task, session, space directly in the pattern
file_pattern: "**/*task-nback*ses-baseline*space-MNI152*stat-effect*.nii.gz"

# Design matrix columns (from participants.tsv)
design_matrix:
  columns:
    - age
    - group
  add_intercept: true
  categorical_columns:
    - group

# Contrasts to compute
contrasts:
  - age                    # Effect of age
  - patients - controls    # Group difference

# Paired test settings (for paired analysis)
paired_test:
  pair_by: sub
  sample_patterns:
    pre: "**/*ses-pre*.nii.gz"
    post: "**/*ses-post*.nii.gz"

# Statistical inference
inference:
  alpha_corrected: 0.05       # Significance for corrected thresholds
  alpha_uncorrected: 0.001    # Cluster-forming threshold
  cluster_threshold: 10
  corrections:
    - uncorrected
    - fdr
    - bonferroni

# Atlas for cluster annotation
atlas: harvard_oxford

# Output settings
output:
  generate_report: true
  report_filename: report.html
```

## Outputs

StatCraft generates:

1. **Statistical Maps** (`*.nii.gz`)
   - T-statistic maps
   - P-value maps
   - Effect size maps
   - Thresholded maps (uncorrected, FDR, Bonferroni)

2. **Cluster Tables** (`*.tsv`)
   - Peak coordinates (MNI)
   - Cluster sizes
   - Peak statistics
   - Anatomical labels

3. **HTML Report** (`report.html`)
   - Methodology description
   - Design matrix visualization
   - Activation maps
   - Glass brain views
   - Cluster tables

4. **Configuration** (`config.yaml`)
   - Complete configuration for reproducibility

## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Citation

If you use StatCraft in your research, please cite:

```bibtex
@software{statcraft,
  title = {StatCraft: Second-Level Neuroimaging Analysis Tool},
  author = {StatCraft Contributors},
  year = {2024},
  url = {https://github.com/arovai/StatCraft}
}
```

## Acknowledgments

- [Nilearn](https://nilearn.github.io/) - Core neuroimaging functionality
- [PyBIDS](https://bids-standard.github.io/pybids/) - BIDS data handling
- [NiBabel](https://nipy.org/nibabel/) - NIfTI file handling
