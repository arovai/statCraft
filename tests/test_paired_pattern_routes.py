"""Regression tests exercise discovery -> pairing -> arithmetic -> GLM -> inference."""

import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm, ttest_1samp

from statcraft.pipeline import StatCraftPipeline


def save_image(path, values, translation=0):
    affine = np.eye(4)
    affine[0, 3] = translation
    nib.save(nib.Nifti1Image(np.asarray(values, dtype=float), affine), path)
    return path


@pytest.fixture
def cohort(tmp_path):
    data = tmp_path / "inputs"
    data.mkdir()
    (data / "dataset_description.json").write_text(json.dumps({
        "Name": "Synthetic paired comparisons", "BIDSVersion": "1.8.0"}))
    for subject in range(1, 5):
        for offset, name in enumerate(("A", "B", "C")):
            values = np.random.default_rng(subject * 10 + offset).normal(size=(3, 3, 3))
            values += subject * 0.2 + offset
            # A-B and A+B-2*C are exactly zero here, and must remain in the mask.
            values[0, 0, 0] = 5
            save_image(data / f"sub-{subject:02}_ses-01_run-1_desc-{name}_cvr.nii.gz", values)
    mask = np.ones((3, 3, 3))
    mask[2, 2, 2] = 0
    mask_path = save_image(tmp_path / "domain.nii.gz", mask)
    return data, mask_path


def pipeline_for(cohort, tmp_path, combine=False, **overrides):
    data, mask_path = cohort
    patterns = {name: f"*desc-{name}_*" for name in (("A", "B", "C") if combine else ("A", "B"))}
    config = {
        "analysis_type": "paired",
        "paired_test": {"sample_patterns": patterns, "pair_by": "sub,ses,run",
                        "combine": "A+B-2*C" if combine else None},
        "inference_mask": str(mask_path),
        "output": {"generate_report": False},
        "inference": {"corrections": [], "cluster_threshold": 0},
        "verbose": 0,
    }
    config.update(overrides)
    pipeline = StatCraftPipeline(data, tmp_path / "out", [], config=config)
    # Annotation may download atlases; it does not alter any tested numeric path.
    pipeline.annotate_clusters = lambda: {}
    return pipeline, patterns


def run_route(pipeline, patterns, combine=False, **kwargs):
    return pipeline.run(sample_patterns=patterns, **kwargs)


@pytest.mark.parametrize("combine", [False, True])
def test_actual_routes_reject_duplicate_complete_keys(cohort, tmp_path, combine):
    data, _ = cohort
    original = next(data.glob("sub-01*desc-A*"))
    save_image(data / original.name.replace("_cvr", "_variant-extra_cvr"), np.ones((3, 3, 3)))
    pipeline, patterns = pipeline_for(cohort, tmp_path, combine)
    with pytest.raises(ValueError, match="Duplicate pairing key"):
        run_route(pipeline, patterns, combine)


@pytest.mark.parametrize("combine", [False, True])
def test_actual_routes_exclude_whole_geometry_invalid_pair(cohort, tmp_path, combine):
    data, _ = cohort
    bad_path = next(data.glob("sub-04*desc-B*"))
    save_image(bad_path, np.full((3, 3, 3), 1e6), translation=20)
    pipeline, patterns = pipeline_for(cohort, tmp_path, combine)
    result = run_route(pipeline, patterns, combine)
    assert len(pipeline._diff_images) == 3
    assert all(info["subject"] != "04" for info in result["images"])
    assert len(pipeline._invalid_images) == (3 if combine else 2)
    saved_exclusions = result["saved_files"]["excluded_files"].read_text()
    assert "Affine mismatch" in saved_exclusions
    assert "sub-04" in saved_exclusions


@pytest.mark.parametrize("combine", [False, True])
def test_actual_routes_finite_mask_keeps_zero_differences_and_matches_reference(cohort, tmp_path, combine):
    data, _ = cohort
    corrupted = next(data.glob("sub-03*desc-B*"))
    values = nib.load(corrupted).get_fdata()
    values[1, 1, 1] = np.nan
    save_image(corrupted, values)
    pipeline, patterns = pipeline_for(cohort, tmp_path, combine)
    result = run_route(pipeline, patterns, combine)
    mask = pipeline.glm.mask.get_fdata().astype(bool)
    assert mask[0, 0, 0]  # zero differences are valid
    assert not mask[1, 1, 1]  # nonfinite in any original map
    assert not mask[2, 2, 2]  # predetermined anatomical exclusion
    assert int(mask.sum()) == 25
    assert np.array_equal(pipeline.glm.model.masker_.mask_img_.get_fdata(), mask)
    assert all(np.all(np.isfinite(img.get_fdata())) for img in pipeline._diff_images)
    expected = []
    for subject in range(1, 5):
        a = nib.load(next(data.glob(f"sub-{subject:02}*desc-A*"))).get_fdata()
        b = nib.load(next(data.glob(f"sub-{subject:02}*desc-B*"))).get_fdata()
        if combine:
            c = nib.load(next(data.glob(f"sub-{subject:02}*desc-C*"))).get_fdata()
            expected.append((a + b - 2*c)[mask])
        else:
            expected.append((a - b)[mask])
    expected = np.stack(expected)
    np.testing.assert_allclose(np.stack([img.get_fdata()[mask] for img in pipeline._diff_images]), expected)
    reference = ttest_1samp(expected, 0, axis=0).statistic
    contrast_name = next(iter(pipeline.glm.results))
    actual = pipeline.glm.results[contrast_name]["maps"]["stat"].get_fdata()[mask]
    np.testing.assert_allclose(actual[1:], reference[1:], rtol=1e-5)
    saved = result["saved_files"]
    np.testing.assert_array_equal(nib.load(saved["inference_mask"]).get_fdata(), mask)
    metadata = json.loads(saved["inference_mask_metadata"].read_text())
    assert metadata["VoxelCount"] == 25
    assert metadata["ZeroValuesAreValid"]
    assert not metadata["NormalizationMasksDefineInference"]


@pytest.mark.parametrize("combine", [False, True])
def test_actual_routes_reject_foreign_normalization_grid(cohort, tmp_path, combine):
    foreign = save_image(tmp_path / "foreign.nii.gz", np.ones((3, 3, 3)), translation=20)
    pipeline, patterns = pipeline_for(cohort, tmp_path, combine)
    with pytest.raises(RuntimeError, match="shape/affine mismatch"):
        run_route(pipeline, patterns, combine, scaling=str(foreign))


@pytest.mark.parametrize("combine", [False, True])
def test_normalization_roi_does_not_define_test_domain(cohort, tmp_path, combine):
    roi = np.zeros((3, 3, 3))
    roi[0, 0, 0] = 1
    roi_path = save_image(tmp_path / "normalization_roi.nii.gz", roi)
    pipeline, patterns = pipeline_for(cohort, tmp_path, combine)
    run_route(pipeline, patterns, combine, scaling=str(roi_path))
    assert int(pipeline.glm.mask.get_fdata().sum()) == 26


def test_permutation_and_all_parametric_families_receive_identical_mask(cohort, tmp_path, monkeypatch):
    import statcraft.core.inference as inference_module

    threshold_masks, permutation_masks = [], []
    real_threshold = inference_module.threshold_stats_img
    real_permutation = inference_module.non_parametric_inference

    def threshold(*args, **kwargs):
        threshold_masks.append(kwargs.get("mask_img"))
        return real_threshold(*args, **kwargs)

    def permutation(*args, **kwargs):
        permutation_masks.append(kwargs.get("mask"))
        return real_permutation(*args, **kwargs)

    monkeypatch.setattr(inference_module, "threshold_stats_img", threshold)
    monkeypatch.setattr(inference_module, "non_parametric_inference", permutation)
    pipeline, patterns = pipeline_for(cohort, tmp_path, inference={
        "corrections": ["uncorrected", "fdr", "bonferroni"],
        "run_permutation": True, "n_permutations": 16, "cluster_threshold": 0,
    }, random_state=7)
    run_route(pipeline, patterns)
    assert len(threshold_masks) == 3
    assert len(permutation_masks) == 1
    assert all(mask is pipeline.glm.mask for mask in threshold_masks + permutation_masks)
    assert "A_minus_B_perm" in pipeline.inference._inference_results
    np.testing.assert_allclose(pipeline.inference.threshold_values["A_minus_B"]["uncorrected"],
                               norm.isf(0.001 / 2))


def test_missing_complete_entity_is_error(cohort, tmp_path):
    pipeline, patterns = pipeline_for(cohort, tmp_path,
                                     paired_test={"sample_patterns": {"A": "*desc-A*", "B": "*desc-B*"},
                                                  "pair_by": "sub,ses,echo"})
    with pytest.raises(ValueError, match="missing pairing entities.*echo"):
        run_route(pipeline, patterns)


def test_simple_ttest_rejects_repeated_subject_runs(cohort, tmp_path):
    data, _ = cohort
    for name in ("A", "B"):
        image = nib.load(next(data.glob(f"sub-01*desc-{name}*")))
        save_image(data / f"sub-01_ses-01_run-2_desc-{name}_cvr.nii.gz", image.get_fdata())
    pipeline, patterns = pipeline_for(cohort, tmp_path)
    with pytest.raises(ValueError, match="Repeated observations"):
        run_route(pipeline, patterns)


@pytest.mark.parametrize("combine", [False, True])
def test_actual_routes_reject_inference_mask_geometry(cohort, tmp_path, combine):
    foreign = save_image(tmp_path / "foreign_domain.nii.gz", np.ones((3, 3, 3)), translation=0.01)
    pipeline, patterns = pipeline_for(cohort, tmp_path, combine, inference_mask=str(foreign))
    with pytest.raises(ValueError, match="Inference mask.*shape/affine mismatch"):
        run_route(pipeline, patterns, combine)


def test_actual_route_rejects_zscore_mask_geometry(cohort, tmp_path):
    foreign = save_image(tmp_path / "foreign_zscore.nii.gz", np.ones((3, 3, 3)), translation=10)
    pipeline, patterns = pipeline_for(cohort, tmp_path)
    with pytest.raises(RuntimeError, match="shape/affine mismatch"):
        run_route(pipeline, patterns, zscore=True, mask=str(foreign))


def test_participant_selection_is_applied_before_pairing(cohort, tmp_path):
    pipeline, patterns = pipeline_for(cohort, tmp_path)
    result = pipeline.run(sample_patterns=patterns, participant_label=["01", "02"])
    assert {info['subject'] for info in result['images']} == {"01", "02"}
    assert len(pipeline._diff_images) == 2


def test_legacy_condition_route_uses_shared_mask_and_excludes_geometry_partner(cohort, tmp_path):
    data, mask = cohort
    for path in data.glob("*desc-C*"):
        path.unlink()
    bad_path = next(data.glob("sub-04*desc-B*"))
    save_image(bad_path, np.ones((3, 3, 3)), translation=20)
    pipeline = StatCraftPipeline(data, tmp_path / 'legacy_output', [], config={
        'analysis_type': 'paired', 'paired_test': {
            'pair_by': 'description', 'condition1': 'A', 'condition2': 'B'},
        'inference_mask': str(mask), 'output': {'generate_report': False},
        'inference': {'corrections': []}, 'verbose': 0,
    })
    pipeline.annotate_clusters = lambda: {}
    result = pipeline.run(pattern='*cvr.nii.gz')
    assert len(pipeline._diff_images) == 3
    assert all(info['subject'] != '04' for info in result['images'])
    assert int(pipeline.glm.mask.get_fdata().sum()) == 26
    assert len(pipeline._invalid_images) == 2


def test_cli_example_runs_full_pattern_path(cohort, tmp_path, monkeypatch):
    import sys
    from statcraft.cli import main

    data, mask = cohort
    monkeypatch.setattr(StatCraftPipeline, 'annotate_clusters', lambda self: {})
    monkeypatch.setattr(StatCraftPipeline, 'generate_report', lambda self: tmp_path / 'report.html')
    monkeypatch.setattr(sys, 'argv', [
        'statcraft', str(data), str(tmp_path / 'cli_output'), 'group',
        '--analysis-type', 'paired', '--patterns', 'A=*desc-A_* B=*desc-B_*',
        '--pair-by', 'sub,ses,run', '--inference-mask', str(mask),
    ])
    main()
    saved_mask = list((tmp_path / 'cli_output' / 'maps').glob('*desc-inference_mask.nii.gz'))
    assert len(saved_mask) == 1
    assert int(nib.load(saved_mask[0]).get_fdata().sum()) == 26
