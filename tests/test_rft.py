"""Software integration checks; these do not establish scientific FWE calibration."""

import json

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from scipy import ndimage, special, stats

from statcraft.config import Config
from statcraft.core.glm import SecondLevelGLM
from statcraft.core.inference import StatisticalInference
from statcraft.core.rft import prepare_context
from statcraft.pipeline import StatCraftPipeline
import test_inference_brain_masks

inputs = test_inference_brain_masks.inputs


@pytest.mark.parametrize("design_kind", ["one-sample", "two-sample", "glm", "paired"])
@pytest.mark.parametrize("two_sided", [True, False])
def test_map_rft_matches_direct_backend_and_fitted_residuals(
    design_kind, two_sided, tmp_path
):
    pyrft = pytest.importorskip("pyrft")
    rng = np.random.default_rng(41)
    shape = (9, 10, 11)
    mask = nib.Nifti1Image(np.ones(shape, dtype=np.uint8), np.eye(4))
    mask.header.set_xyzt_units("mm")
    images = []
    for _ in range(12):
        data = ndimage.gaussian_filter(rng.normal(size=shape), 1.3) + 0.3
        img = nib.Nifti1Image(data, np.eye(4))
        img.header.set_xyzt_units("mm")
        images.append(img)
    glm = SecondLevelGLM(mask=mask, smoothing_fwhm=1.5, retain_residuals=True)
    if design_kind == "one-sample":
        glm.one_sample_ttest(images, contrast_name="test")
    elif design_kind == "paired":
        glm.paired_ttest(images[:6], images[6:], contrast_name="test")
    elif design_kind == "two-sample":
        glm.two_sample_ttest(images[:6], images[6:], contrast_name="test")
    else:
        glm.fit(images, pd.DataFrame({"intercept": np.ones(12), "age": np.arange(12)}))
        glm.compute_contrast("intercept", contrast_name="test")
    context = prepare_context(glm.model)
    n = 6 if design_kind == "paired" else 12
    rank = 2 if design_kind in ("glm", "two-sample") else 1
    assert context["df"] == n - rank
    # Independent projection of exactly the fitted (including paired difference) inputs.
    y = glm.model.masker_.transform(glm._images)
    x = glm._design_matrix.to_numpy()
    residuals = y - x @ np.linalg.pinv(x) @ y
    expected_smoothness = pyrft.estimate_smoothness_residuals(
        glm.model.masker_.inverse_transform(residuals).get_fdata(),
        mask.get_fdata(),
        residual_df=n - rank,
    )
    np.testing.assert_allclose(
        context["smoothness"]["fwhm_vox"], expected_smoothness["fwhm_vox"]
    )
    inference = StatisticalInference(
        cluster_threshold=0, two_sided=two_sided, mask=mask
    )
    output, _ = inference.threshold_rft(glm.results["test"], context, "test")
    direct = pyrft.threshold_map(
        glm.results["test"]["maps"]["stat"],
        mask,
        stat="T",
        df=n - rank,
        alpha=0.05,
        fwhm=expected_smoothness["fwhm_vox"],
        tail="two-sided" if two_sided else "positive",
    )
    expected = np.where(
        direct.thresholded_map.get_fdata() != 0,
        glm.results["test"]["maps"]["z_score"].get_fdata(),
        0,
    )
    np.testing.assert_allclose(output.get_fdata(), expected)
    zcut = -special.ndtri_exp(stats.t.logsf(direct.threshold, n - rank))
    assert inference.threshold_values["test"]["rft"] == pytest.approx(zcut)
    json.dumps(inference.rft_metadata["test"], allow_nan=False)


def test_connectivity_rejects_rft_before_loading_data(tmp_path):
    (tmp_path / "participants.tsv").write_text("participant_id\nsub-01\n")
    pipeline = StatCraftPipeline(
        tmp_path, tmp_path / "out", [], config={"inference": {"corrections": ["rft"]}}
    )
    with pytest.raises(ValueError, match="only for volumetric maps"):
        pipeline.run_connectivity_analysis()


def test_rft_alpha_validation():
    with pytest.raises(ValueError, match="alpha_corrected <= 0.1"):
        Config(inference={"corrections": ["rft"], "alpha_corrected": 0.2})


@pytest.mark.parametrize(
    "route", ["one-sample", "two-sample", "paired", "combine", "legacy"]
)
def test_pipeline_rft_routes_and_sidecars(inputs, tmp_path, route):
    pytest.importorskip("pyrft")
    # Use enough residual degrees of freedom for the supported T-field inversion.
    with (inputs / "participants.tsv").open("a") as stream:
        for subject in range(5, 13):
            stream.write(f"sub-{subject:02}\n")
            for condition in ("A", "B"):
                values = np.random.default_rng(subject * 3 + (condition == "B")).normal(
                    size=(4, 4, 4)
                )
                nib.save(
                    nib.Nifti1Image(values, np.eye(4)),
                    inputs / f"sub-{subject:02}_desc-{condition}_stat.nii.gz",
                )
            nib.save(
                nib.Nifti1Image(np.ones((4, 4, 4)), np.eye(4)),
                inputs / f"sub-{subject:02}_desc-brain_mask.nii.gz",
            )
    patterns = {"A": "*desc-A_stat.nii.gz", "B": "*desc-B_stat.nii.gz"}
    paired = {"sample_patterns": patterns, "pair_by": "sub"}
    if route == "combine":
        paired["combine"] = "A-B"
    if route == "legacy":
        for path in list(inputs.glob("*_stat.nii.gz")):
            path.rename(path.with_name(path.name.replace("_desc-", "_ses-")))
        paired = {"pair_by": "session", "condition1": "A", "condition2": "B"}
    config = {
        "analysis_type": "paired"
        if route in ("paired", "combine", "legacy")
        else route,
        "paired_test": paired,
        "inference": {"corrections": ["rft"], "cluster_threshold": 0},
        "output": {"generate_report": False},
        "verbose": 0,
    }
    if route == "two-sample":
        config["sample_patterns"] = patterns
    pipeline = StatCraftPipeline(inputs, tmp_path / "out", [], config=config)
    pipeline.annotate_clusters = lambda: {}
    if route in ("two-sample", "paired", "combine"):
        result = pipeline.run(sample_patterns=patterns)
    else:
        result = pipeline.run(
            pattern="*_stat.nii.gz" if route == "legacy" else patterns["A"]
        )
    sidecars = [
        path
        for key, path in result["saved_files"].items()
        if key.endswith("_rft_map_json")
    ]
    assert sidecars
    for path in sidecars:
        metadata = json.loads(path.read_text())
        assert metadata["StatisticType"] == "Z"
        assert metadata["RFT"]["smoothness"]["source"] == "residuals"
        assert metadata["ActualZThreshold"] == metadata["RFT"]["z_threshold"]


def test_rejects_f_contrasts_and_missing_maps():
    from statcraft.core.rft import threshold_contrast

    with pytest.raises(ValueError, match="T contrasts only"):
        threshold_contrast({"stat_type": "F"}, {}, 0.05, True)
    with pytest.raises(ValueError, match="both native T and Z"):
        threshold_contrast({"stat_type": "t", "maps": {}}, {}, 0.05, True)


def test_missing_optional_backend_has_installation_guidance(monkeypatch):
    import sys
    from statcraft.core.rft import require_pyrft

    monkeypatch.setitem(sys.modules, "pyrft", None)
    with pytest.raises(ImportError, match=r"\[rft\]"):
        require_pyrft()
