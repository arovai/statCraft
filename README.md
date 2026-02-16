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
  - Works with derivatives from fMRIPrep, SPM, FSL, etc.
  - Supports BIDS-like and custom file naming

- **Statistical Inference**
  - Uncorrected thresholding (default p < 0.001)
  - FDR (False Discovery Rate) correction
  - FWER (Family-Wise Error Rate) via Bonferroni
  - Permutation-based FWER correction

- **Anatomical Annotation**
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

StatCraft uses flexible pattern matching for file discovery. You can organize your data in two ways:

**Pattern 1: Dataset with Separate Derivatives (Recommended)**
```bash
statcraft <INPUT_DIR> <OUTPUT_DIR> group --derivatives <DERIVATIVES_PATH> [options]
```
Use this when you have a dataset root with participants.tsv and processed derivatives (e.g., fMRIPrep output).

**Pattern 2: Derivatives-Only Directory**
```bash
statcraft <DERIVATIVES_PATH> <OUTPUT_DIR> group [--participants-file <PATH>] [options]
```
Use this when analyzing data directly from a derivatives folder. You may need to specify the `--participants-file` path if it's not in the derivatives folder.

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
