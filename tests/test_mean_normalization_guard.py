"""Cancellation safeguards must be independent of CVR units and preserve signs."""

import nibabel as nib
import numpy as np
import pytest

from statcraft.pipeline import StatCraftPipeline


def images(values):
    data = np.asarray(values, dtype=float).reshape(1, 1, -1)
    return (nib.Nifti1Image(data, np.eye(4)),
            nib.Nifti1Image(np.ones(data.shape), np.eye(4)))


@pytest.mark.parametrize('scale', [1e-9, 1, 1e9])
def test_stable_signed_maps_keep_negative_voxels_and_record_guard(scale):
    image, mask = images(scale * np.array([2, 2, -1]))
    qc = {}
    result = StatCraftPipeline._scale_image_by_mask(
        image, mask, minimum_absolute_mean_ratio=.01,
        require_positive_mean=True, normalization_qc=qc)
    np.testing.assert_allclose(result.get_fdata().ravel(), [2, 2, -1])
    assert qc['Mean'] == pytest.approx(scale)
    assert qc['SignedMeanToMeanAbsoluteRatio'] == pytest.approx(.6)
    assert qc['VoxelCount'] == 3
    assert qc['Passed']


@pytest.mark.parametrize('scale', [1e-9, 1, 1e9])
def test_near_cancellation_fails_at_any_unit_scale(scale):
    image, mask = images(scale * np.array([1.001, -.999]))
    with pytest.raises(ValueError, match='Unstable normalization denominator'):
        StatCraftPipeline._scale_image_by_mask(image, mask, minimum_absolute_mean_ratio=.01)


@pytest.mark.parametrize('values', [[0, 0], [1, -1]])
def test_zero_mean_always_fails(values):
    image, mask = images(values)
    with pytest.raises(ValueError, match='zero mean'):
        StatCraftPipeline._scale_image_by_mask(image, mask)


def test_positive_mean_requirement_is_an_explicit_application_policy():
    image, mask = images([-2, -1])
    result = StatCraftPipeline._scale_image_by_mask(image, mask, minimum_absolute_mean_ratio=.01)
    np.testing.assert_allclose(result.get_fdata().ravel(), [4/3, 2/3])
    with pytest.raises(ValueError, match='positive mean'):
        StatCraftPipeline._scale_image_by_mask(image, mask, require_positive_mean=True)


def test_guard_uses_only_declared_support_not_outside_values():
    image, mask = images([2, 2, -1, -1e10])
    support = np.array([True, True, True, False]).reshape(image.shape)
    qc = {}
    result = StatCraftPipeline._scale_image_by_mask(
        image, mask, support=support, minimum_absolute_mean_ratio=.01,
        require_positive_mean=True, normalization_qc=qc)
    assert qc['SignedMeanToMeanAbsoluteRatio'] == pytest.approx(.6)
    np.testing.assert_allclose(result.get_fdata()[support], [2, 2, -1])


@pytest.mark.parametrize('ratio', [-.1, 1.1, np.nan, np.inf, True, '0.01'])
def test_invalid_ratio_fails(ratio):
    image, mask = images([1, 2])
    with pytest.raises(ValueError, match='minimum_absolute_mean_ratio'):
        StatCraftPipeline._scale_image_by_mask(image, mask, minimum_absolute_mean_ratio=ratio)


def test_invalid_sign_policy_fails():
    image, mask = images([1, 2])
    with pytest.raises(ValueError, match='require_positive_mean'):
        StatCraftPipeline._scale_image_by_mask(image, mask, require_positive_mean='false')
