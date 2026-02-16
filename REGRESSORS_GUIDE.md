# Regressors & Contrasts Guide

This guide covers how to specify regressors (covariates) and contrasts in StatCraft for both NIfTI and connectivity analyses.

## Overview

When running a `glm` analysis type, you can include columns from `participants.tsv` as regressors in the design matrix. StatCraft handles:

- **Continuous regressors** (e.g., age, IQ): z-scored by default (can be disabled)
- **Categorical regressors** (e.g., sex, group): dummy-coded into indicator columns

---

## CLI Options

### `--regressors`

Specify one or more column names from `participants.tsv` to include in the design matrix. By default, all regressors are treated as **continuous** and z-scored.

```bash
--regressors age sex IQ
```

For columns with **non-numerical values** (e.g., M/F for sex, treatment codes), use `--categorical-regressors`.

### `--categorical-regressors`

Explicitly declare which regressors are categorical (dummy-coded). Supports both numerical and non-numerical values.

```bash
--categorical-regressors sex group treatment
```

**Examples:**
- `sex` with values: M, F (non-numerical) ✓
- `group` with values: control, patient (non-numerical) ✓
- `treatment` with values: 0, 1 (numerical but treat as categorical) ✓

If a column contains non-numerical values and is NOT marked as categorical, StatCraft will raise an error suggesting to use `--categorical-regressors`.

### `--no-standardize-regressors`

Disable z-scoring for specific continuous regressors. By default, continuous regressors are standardized to mean=0, std=1. Specify column names to keep in original units.

```bash
--no-standardize-regressors age IQ
```

This is useful when coefficients need to be interpretable in original units (e.g., effect per year of age).

> **Note:** This is independent of `--zscore`, which z-scores brain maps at the subject level.

### `--contrast`

Specify the contrast expression to test. Contrast expressions reference design matrix column names.

```bash
--contrast 'age'           # Effect of age
--contrast 'sex_M-sex_F'   # Males vs Females
```

---

## Design Matrix Column Naming

Understanding how columns are named in the design matrix is key to writing contrasts:

| Regressor type | Input column | Design matrix columns | Notes |
|---|---|---|---|
| Continuous (numerical) | `age` | `age` | z-scored by default unless in `--no-standardize-regressors` |
| Continuous (numerical) | `IQ` | `IQ` | z-scored by default |
| Categorical (non-numerical) | `sex` (M, F) | `sex_M`, `sex_F` | Requires `--categorical-regressors sex` |
| Categorical (non-numerical) | `group` (HC, PAT) | `group_HC`, `group_PAT` | Requires `--categorical-regressors group` |
| Categorical (numerical) | `treatment` (0, 1) | `treatment_0`, `treatment_1` | Requires `--categorical-regressors treatment` |
| Intercept | (automatic) | `intercept` | Added by default |

So with `--regressors age sex --categorical-regressors sex`, the design matrix columns would be:

```
intercept | age | sex_M | sex_F
```

### Handling Non-Numerical Categorical Values

When a column contains non-numerical values (strings), it **must** be marked as categorical:

**Example: sex column with M/F values**

```bash
statcraft ... --regressors age sex --categorical-regressors sex
```

This creates:
- `sex_M`: 1 if sex=M, 0 otherwise
- `sex_F`: 1 if sex=F, 0 otherwise

**Without `--categorical-regressors sex`:** Error!
```
ValueError: Column 'sex' contains non-numerical values and cannot be treated 
as continuous. Found values: ['M', 'F']. Please use --categorical-regressors sex 
to treat it as categorical.
```

---

## Contrast Specification

Contrasts are mathematical expressions that test hypotheses about effects. For categorical variables, StatCraft supports **both original values and dummy-coded column names**.

### Contrast Notation

When you use `--categorical-regressors sex` with values M and F, StatCraft creates dummy columns: `sex_M` and `sex_F`.

You can write contrasts using **either notation**:

| Original Values | Dummy-Coded Columns | Both work? |
|---|---|---|
| `M-F` | `sex_M-sex_F` | ✓ Yes (equivalent) |
| `0.5*M+0.5*F` | `0.5*sex_M+0.5*sex_F` | ✓ Yes (equivalent) |

The first notation is closer to what you see in `participants.tsv`, making it more intuitive for users.

### Examples

| Goal | Contrast expression | Explanation |
|---|---|---|
| Effect of age | `age` | Tests whether age has a significant effect |
| Males vs Females | `M-F` or `sex_M-sex_F` | Difference between male and female groups (both notations work) |
| Mean across groups | `0.5*M+0.5*F` or `0.5*sex_M+0.5*sex_F` | Average effect across both sexes (both notations work) |
| Intercept (mean) | `intercept` | Tests whether the mean effect differs from 0 |

### Default contrast

If no `--contrast` is specified, StatCraft defaults to testing the intercept (`mean_effect` = `[1, 0, 0, ...]`).

---

## Full CLI Examples

### 1. Effect of age on connectivity (NIfTI)

```bash
statcraft /data/bids /data/output \
  -d /data/derivatives/fmriprep \
  -t glm \
  --regressors age \
  --contrast age \
  -p '*task-rest*bold*.nii.gz'
```

### 2. Sex differences in connectivity matrices

```bash
statcraft /data/bids /data/output \
  -d /data/derivatives \
  --data-type connectivity \
  -t glm \
  --regressors sex \
  --categorical-regressors sex \
  --contrast 'M-F' \
  -p '**/*_connmat.npy'
```

Note: You can also use `--contrast 'sex_M-sex_F'` (dummy-coded names) instead of `'M-F'` (original values). Both work the same way.

### 3. Age effect controlling for sex

```bash
statcraft /data/bids /data/output \
  -d /data/derivatives \
  --data-type connectivity \
  -t glm \
  --regressors age --regressors sex \
  --categorical-regressors sex \
  --contrast age \
  -p '**/*_connmat.npy'
```

### 4. Multiple regressors without standardization

```bash
statcraft /data/bids /data/output \
  -d /data/derivatives/fmriprep \
  -t glm \
  --regressors age --regressors IQ \
  --no-standardize-regressors \
  --contrast age \
  -p '*task-rest*bold*.nii.gz'
```

### 5. Simple one-sample test (no regressors)

```bash
statcraft /data/bids /data/output \
  -d /data/derivatives \
  --data-type connectivity \
  -t one-sample \
  -p '**/*_connmat.npy'
```

---

## YAML Configuration

All CLI options can also be set in a config file:

```yaml
analysis_type: glm

design_matrix:
  columns:
    - age
    - sex
  categorical_columns:
    - sex
  standardize_continuous: true   # set to false to disable z-scoring
  add_intercept: true

# Contrast expression or list of contrasts
contrast: "age"

# Or multiple contrasts:
# contrasts:
#   - "age"
#   - expression: "sex_M-sex_F"
#     name: "sexDifference"
```

---

## Output Naming

Output files include the contrast name using BIDS-compatible naming:

```
output_dir/
├── data/
│   ├── {bids_prefix}_contrast-effectOfAge_stat-tstat.npy
│   ├── {bids_prefix}_contrast-effectOfAge_stat-pval.npy
│   ├── {bids_prefix}_contrast-effectOfAge_stat-tstat.json
│   ├── {bids_prefix}_contrast-effectOfAge_stat-fdr_threshold.npy
│   └── {bids_prefix}_contrast-effectOfAge_stat-bonferroni_threshold.npy
├── tables/
│   ├── {bids_prefix}_contrast-effectOfAge_stat-fdr_edges.tsv
│   └── {bids_prefix}_contrast-effectOfAge_stat-bonferroni_edges.tsv
└── {bids_prefix}_contrast-effectOfAge_report.html
```

The contrast name is auto-generated from the expression:
- `age` → `effectOfAge`
- `sex_M-sex_F` → `sex_MVersusSex_f`
- Or provide a custom name via the YAML config's `name` field
