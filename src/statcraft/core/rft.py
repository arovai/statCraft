"""Residual-based, voxel-level RFT inference for volumetric group models.

The backend estimates stationary equivalent smoothness from full OLS residuals.
No statistic-map smoothness approximation or smoothing-kernel substitution is used.
"""

import numpy as np
from scipy import special, stats


def require_pyrft():
    """Load the optional backend without affecting ordinary or table analyses."""
    try:
        import pyrft
    except ImportError as exc:
        raise ImportError(
            "RFT requires the optional backend: python -m pip install -e '.[rft]'"
        ) from exc
    return pyrft


def prepare_context(model):
    """Estimate smoothness on the fitted mask from retained Nilearn OLS residuals."""
    backend = require_pyrft()
    if model is None or model.minimize_memory or not model.results_:
        raise ValueError("RFT requires a fitted map model with retained residuals")
    labels = np.asarray(model.labels_)
    n = len(model.design_matrix_)
    residuals = np.empty((n, labels.size), dtype=float)
    assigned = np.zeros(labels.size, dtype=bool)
    dfs = []
    for key, result in model.results_.items():
        selected = labels == key
        values = np.asarray(result.whitened_residuals)
        if values.shape != (n, int(selected.sum())):
            raise ValueError("Unexpected fitted residual dimensions")
        residuals[:, selected] = values
        assigned[selected] = True
        dfs.append(float(result.df_residuals))
    if (
        not assigned.all()
        or not np.allclose(dfs, dfs[0])
        or not np.isfinite(dfs[0])
        or dfs[0] <= 1
    ):
        raise ValueError("RFT requires complete residuals and a common residual df > 1")
    mask = model.masker_.mask_img_
    residual_img = model.masker_.inverse_transform(residuals)
    smoothness = backend.estimate_smoothness_residuals(
        residual_img.get_fdata(), mask.get_fdata(), residual_df=dfs[0], n_observations=n
    )
    return {
        "mask": mask,
        "df": dfs[0],
        "smoothness": smoothness,
        "n_observations": n,
        "design_rank": int(np.linalg.matrix_rank(model.design_matrix_)),
    }


def threshold_contrast(result, context, alpha, two_sided):
    """Threshold native T statistics and preserve surviving Z scores for output."""
    if result["stat_type"].lower() != "t":
        raise ValueError("RFT currently supports T contrasts only, not F contrasts")
    maps = result["maps"]
    if not {"stat", "z_score"}.issubset(maps):
        raise ValueError("RFT requires both native T and Z maps")
    t_img, z_img = maps["stat"], maps["z_score"]
    if t_img.shape != z_img.shape or not np.allclose(t_img.affine, z_img.affine):
        raise ValueError("Native T and Z maps must share the same grid")
    outcome = require_pyrft().threshold_map(
        t_img,
        context["mask"],
        stat="T",
        df=context["df"],
        alpha=alpha,
        tail="two-sided" if two_sided else "positive",
        method="rft",
        fwhm=context["smoothness"]["fwhm_vox"],
    )
    z_threshold = float(
        -special.ndtri_exp(stats.t.logsf(outcome.threshold, context["df"]))
    )
    metadata = outcome.metadata
    metadata.update(
        smoothness=context["smoothness"],
        z_threshold=z_threshold,
        output_statistic="Z",
        inference_scope="voxel-level FWE per contrast",
        residual_preparation={
            "source": "Nilearn fitted OLS whitened_residuals",
            "standardization": "residual / sqrt(RSS / residual_df)",
            "n_observations": context["n_observations"],
            "design_rank": context["design_rank"],
        },
    )
    from nilearn.image import new_img_like

    selected = outcome.thresholded_map.get_fdata() != 0
    data = z_img.get_fdata()
    if not np.isfinite(data[context["mask"].get_fdata() != 0]).all():
        raise ValueError("RFT requires finite in-mask Z scores")
    return new_img_like(z_img, np.where(selected, data, 0)), z_threshold, metadata
