"""Tests for the emsoft.indexing module."""

import numpy as np
import pytest
from emsoft.indexing import (
    normalize_patterns, circular_mask, match_patterns,
    confidence_index, adp_map, kam_map, convert_pattern_center,
)


class TestNormalizePatterns:
    def test_unit_norm(self):
        patterns = np.random.rand(10, 64)
        normed = normalize_patterns(patterns)
        norms = np.linalg.norm(normed.reshape(10, -1), axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-14)

    def test_zero_pattern(self):
        patterns = np.zeros((3, 16))
        normed = normalize_patterns(patterns)
        assert np.isfinite(normed).all()

    def test_shape_preserved(self):
        patterns = np.random.rand(5, 8, 8)
        normed = normalize_patterns(patterns)
        assert normed.shape == (5, 8, 8)


class TestCircularMask:
    def test_shape(self):
        mask = circular_mask(64, 48)
        assert mask.shape == (48, 64)

    def test_center_is_one(self):
        mask = circular_mask(100, 100)
        assert mask[50, 50] == 1.0

    def test_corner_is_zero(self):
        mask = circular_mask(100, 100)
        assert mask[0, 0] == 0.0


class TestMatchPatterns:
    def test_perfect_match(self):
        # Dictionary pattern should match itself perfectly
        dictionary = normalize_patterns(np.random.rand(5, 32))
        result = match_patterns(dictionary, dictionary, n_top=1)
        # Each pattern's best match should be itself
        np.testing.assert_array_equal(result['indices'][:, 0], np.arange(5))
        np.testing.assert_allclose(result['dot_products'][:, 0], 1.0, atol=1e-14)

    def test_n_top(self):
        dictionary = normalize_patterns(np.random.rand(20, 32))
        experimental = normalize_patterns(np.random.rand(5, 32))
        result = match_patterns(experimental, dictionary, n_top=3)
        assert result['indices'].shape == (5, 3)
        assert result['dot_products'].shape == (5, 3)

    def test_dot_products_sorted(self):
        dictionary = normalize_patterns(np.random.rand(10, 32))
        experimental = normalize_patterns(np.random.rand(3, 32))
        result = match_patterns(experimental, dictionary, n_top=5)
        # Top matches should be in descending order of dot product
        for i in range(3):
            dp = result['dot_products'][i]
            assert np.all(dp[:-1] >= dp[1:])

    def test_dot_products_in_range(self):
        dictionary = normalize_patterns(np.random.rand(10, 32))
        experimental = normalize_patterns(np.random.rand(5, 32))
        result = match_patterns(experimental, dictionary, n_top=3)
        assert np.all(result['dot_products'] >= -1.0)
        assert np.all(result['dot_products'] <= 1.0 + 1e-14)


class TestConfidenceIndex:
    def test_basic(self):
        dp = np.array([[0.95, 0.80], [0.90, 0.89]])
        ci = confidence_index(dp)
        np.testing.assert_allclose(ci, [0.15, 0.01])

    def test_needs_two(self):
        dp = np.array([[0.95]])
        with pytest.raises(ValueError):
            confidence_index(dp)


class TestADPMap:
    def test_shape(self):
        patterns = normalize_patterns(np.random.rand(12, 32))
        adp = adp_map(patterns, width=4, height=3)
        assert adp.shape == (3, 4)

    def test_identical_patterns(self):
        # All identical normalized patterns should give ADP = 1.0
        p = np.random.rand(1, 32)
        p = p / np.linalg.norm(p)
        patterns = np.tile(p, (6, 1))
        adp = adp_map(patterns, width=3, height=2)
        np.testing.assert_allclose(adp, 1.0, atol=1e-14)

    def test_wrong_size(self):
        with pytest.raises(ValueError):
            adp_map(np.zeros((10, 32)), width=3, height=2)


class TestKAMMap:
    def test_shape(self):
        eulers = np.random.rand(12, 3) * 0.1
        kam = kam_map(eulers, width=4, height=3)
        assert kam.shape == (3, 4)

    def test_identical_orientations(self):
        # All same orientation should give KAM = 0
        eulers = np.tile([0.1, 0.2, 0.3], (6, 1))
        kam = kam_map(eulers, width=3, height=2)
        np.testing.assert_allclose(kam, 0.0, atol=1e-10)

    def test_nonnegative(self):
        eulers = np.random.rand(12, 3) * 0.5
        kam = kam_map(eulers, width=4, height=3)
        assert np.all(kam >= 0)


class TestConvertPatternCenter:
    def test_emsoft_roundtrip(self):
        pc = [5.0, 10.0, 15000.0]
        result = convert_pattern_center(pc, 'EMsoft', 'EMsoft', 50.0, 640, 480)
        np.testing.assert_allclose(result, pc)

    def test_edax_roundtrip(self):
        pc = [0.5, 0.5, 0.6]
        emsoft = convert_pattern_center(pc, 'EDAX', 'EMsoft', 50.0, 640, 480)
        back = convert_pattern_center(emsoft, 'EMsoft', 'EDAX', 50.0, 640, 480)
        np.testing.assert_allclose(back, pc, atol=1e-12)

    def test_oxford_roundtrip(self):
        pc = [0.4, 0.6, 0.5]
        emsoft = convert_pattern_center(pc, 'Oxford', 'EMsoft', 50.0, 640, 480)
        back = convert_pattern_center(emsoft, 'EMsoft', 'Oxford', 50.0, 640, 480)
        np.testing.assert_allclose(back, pc, atol=1e-12)

    def test_unknown_vendor(self):
        with pytest.raises(ValueError):
            convert_pattern_center([0, 0, 0], 'Unknown', 'EMsoft', 50.0, 640, 480)

    def test_edax_center(self):
        # EDAX (0.5, 0.5, z*) should map to EMsoft (0, 0, L)
        pc = [0.5, 0.5, 0.5]
        emsoft = convert_pattern_center(pc, 'EDAX', 'EMsoft', 50.0, 640, 480)
        assert emsoft[0] == pytest.approx(0.0)  # xpc = 0
        assert emsoft[1] == pytest.approx(0.0)  # ypc = 0
