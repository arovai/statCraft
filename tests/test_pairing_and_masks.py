import nibabel as nib
import numpy as np
import pytest

from statcraft.core.data_loader import DataLoader


def image_info(condition, path):
    return {
        "subject": "01",
        "path": str(path),
        "entities": {"session": condition},
    }


def test_ordinary_pairing_rejects_duplicate_condition(tmp_path):
    loader = object.__new__(DataLoader)
    images = [
        image_info("pre", tmp_path / "a.nii.gz"),
        image_info("pre", tmp_path / "b.nii.gz"),
        image_info("post", tmp_path / "c.nii.gz"),
    ]
    with pytest.raises(ValueError, match="Duplicate pairing key"):
        loader.organize_paired_data(images, "session", "pre", "post")


def test_shared_mask_is_true_intersection(tmp_path):
    from statcraft.pipeline import StatCraftPipeline

    first = tmp_path / "first.nii.gz"
    second = tmp_path / "second.nii.gz"
    nib.save(nib.Nifti1Image(np.array([[[1, 1], [0, 0]]], dtype=np.uint8), np.eye(4)), first)
    nib.save(nib.Nifti1Image(np.array([[[1, 0], [1, 0]]], dtype=np.uint8), np.eye(4)), second)
    pipeline = object.__new__(StatCraftPipeline)
    pipeline._zscore_masks_used = {"a": str(first), "b": str(second)}
    result = pipeline._compute_intersection_mask().get_fdata()
    assert int(result.sum()) == 1
