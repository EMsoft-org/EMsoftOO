"""Tests for the emsoft.ebsd module."""

import numpy as np
import pytest
from emsoft.ebsd import EBSDDetector, EBSDSimulator


def _make_test_master_pattern(npx=50):
    """Create a simple synthetic master pattern for testing.

    Returns a pattern with a cosine-like intensity: bright at the
    north pole (center of NH), dark at the equator.
    """
    n = 2 * npx + 1
    mLPNH = np.zeros((n, n), dtype=np.float64)
    mLPSH = np.zeros((n, n), dtype=np.float64)

    for i in range(n):
        for j in range(n):
            # Distance from center, normalized to [0, 1]
            x = (i - npx) / npx
            y = (j - npx) / npx
            r2 = x * x + y * y
            if r2 <= 1.0:
                mLPNH[i, j] = 1.0 - r2  # bright at center
                mLPSH[i, j] = r2         # bright at edges

    return mLPNH, mLPSH, npx


class TestEBSDDetector:
    def test_create(self):
        det = EBSDDetector(640, 480)
        assert det.numsx == 640
        assert det.numsy == 480

    def test_direction_cosines_shape(self):
        det = EBSDDetector(100, 80)
        assert det.rgx.shape == (100, 80)
        assert det.rgy.shape == (100, 80)
        assert det.rgz.shape == (100, 80)

    def test_direction_cosines_normalized(self):
        det = EBSDDetector(50, 40, L=15000.0, thetac=10.0)
        norms = np.sqrt(det.rgx**2 + det.rgy**2 + det.rgz**2)
        np.testing.assert_allclose(norms, 1.0, atol=1e-12)

    def test_repr(self):
        det = EBSDDetector(640, 480)
        assert '640x480' in repr(det)


class TestEBSDSimulator:
    def test_create(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern(20)
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        assert sim.npx == 20

    def test_create_infer_npx(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern(30)
        sim = EBSDSimulator(mLPNH, mLPSH)
        assert sim.npx == 30

    def test_wrong_shape(self):
        with pytest.raises(ValueError):
            EBSDSimulator(np.zeros((10, 10)), np.zeros((10, 10)), npx=20)

    def test_compute_single_pattern(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern()
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        det = EBSDDetector(64, 48, L=15000.0, thetac=10.0, delta=50.0)

        pattern = sim.compute_pattern(det, quaternion=[1, 0, 0, 0])
        assert pattern.shape == (48, 64)  # (numsy, numsx) image convention
        assert np.isfinite(pattern).all()

    def test_compute_pattern_euler(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern()
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        det = EBSDDetector(64, 48, L=15000.0, thetac=10.0, delta=50.0)

        pattern = sim.compute_pattern(det, euler=[30, 45, 60], degrees=True)
        assert pattern.shape == (48, 64)

    def test_different_orientations_differ(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern()
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        det = EBSDDetector(64, 48, L=15000.0, thetac=10.0, delta=50.0)

        p1 = sim.compute_pattern(det, quaternion=[1, 0, 0, 0])
        p2 = sim.compute_pattern(det, quaternion=[0.5, 0.5, 0.5, 0.5])
        # Different orientations should produce different patterns
        assert not np.allclose(p1, p2)

    def test_compute_batch_patterns(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern()
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        det = EBSDDetector(64, 48, L=15000.0, thetac=10.0, delta=50.0)

        quats = np.array([
            [1, 0, 0, 0],
            [0.5, 0.5, 0.5, 0.5],
            [0, 1, 0, 0],
        ], dtype=np.float64)
        patterns = sim.compute_patterns(det, quats)
        assert patterns.shape == (3, 48, 64)
        assert np.isfinite(patterns).all()

    def test_no_orientation_raises(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern(10)
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        det = EBSDDetector(32, 24)
        with pytest.raises(ValueError):
            sim.compute_pattern(det)

    def test_repr(self):
        mLPNH, mLPSH, npx = _make_test_master_pattern(25)
        sim = EBSDSimulator(mLPNH, mLPSH, npx)
        assert '25' in repr(sim)
