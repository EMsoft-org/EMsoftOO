"""Tests for the emsoft.crystallography module."""

import numpy as np
import pytest
from emsoft.crystallography import Crystal
from emsoft.symmetry import SpaceGroup


class TestCrystalCreation:
    def test_cubic(self):
        c = Crystal(0.35236, 0.35236, 0.35236, 90, 90, 90)
        lp = c.lattice_parameters
        np.testing.assert_allclose(lp[:3], [0.35236, 0.35236, 0.35236])
        np.testing.assert_allclose(lp[3:], [90, 90, 90])

    def test_hexagonal(self):
        c = Crystal(0.3209, 0.3209, 0.5211, 90, 90, 120)
        assert c.volume > 0

    def test_repr(self):
        c = Crystal(0.5, 0.5, 0.5, 90, 90, 90)
        assert 'Crystal' in repr(c)


class TestCrystalVolume:
    def test_cubic_volume(self):
        a = 0.35236
        c = Crystal(a, a, a, 90, 90, 90)
        expected = a ** 3
        assert c.volume == pytest.approx(expected, rel=1e-6)


class TestCrystalMetricTensors:
    def test_cubic_metric_tensor(self):
        a = 0.5
        c = Crystal(a, a, a, 90, 90, 90)
        dmt = c.direct_metric_tensor
        expected = np.diag([a**2, a**2, a**2])
        np.testing.assert_allclose(dmt, expected, atol=1e-14)

    def test_metric_tensor_symmetry(self):
        c = Crystal(0.3, 0.4, 0.5, 80, 85, 75)
        dmt = c.direct_metric_tensor
        np.testing.assert_allclose(dmt, dmt.T, atol=1e-14)

    def test_reciprocal_metric_tensor(self):
        a = 0.5
        c = Crystal(a, a, a, 90, 90, 90)
        rmt = c.reciprocal_metric_tensor
        expected = np.diag([1/a**2, 1/a**2, 1/a**2])
        np.testing.assert_allclose(rmt, expected, atol=1e-10)


class TestCrystalTransformations:
    def test_identity_transform(self):
        c = Crystal(0.5, 0.5, 0.5, 90, 90, 90)
        v = [1.0, 0.0, 0.0]
        # d->d should be identity
        result = c.transform(v, 'd', 'd')
        np.testing.assert_allclose(result, v, atol=1e-14)

    def test_cubic_d_to_c(self):
        a = 1.0
        c = Crystal(a, a, a, 90, 90, 90)
        # In cubic, direct = Cartesian (up to scale by a)
        v = [1.0, 0.0, 0.0]
        result = c.transform(v, 'd', 'c')
        np.testing.assert_allclose(result, [a, 0, 0], atol=1e-14)


class TestCrystalDotProduct:
    def test_cubic_dot(self):
        a = 1.0
        c = Crystal(a, a, a, 90, 90, 90)
        # In cubic, dot product in direct space: g_ij u^i v^j = a^2 * (u.v)
        d = c.dot([1, 0, 0], [1, 0, 0], space='d')
        assert d == pytest.approx(a**2)

    def test_orthogonal_dot(self):
        c = Crystal(1.0, 1.0, 1.0, 90, 90, 90)
        d = c.dot([1, 0, 0], [0, 1, 0], space='d')
        assert d == pytest.approx(0.0, abs=1e-14)


class TestCrystalLength:
    def test_cubic_length(self):
        a = 0.5
        c = Crystal(a, a, a, 90, 90, 90)
        L = c.length([1, 1, 1], space='d')
        assert L == pytest.approx(a * np.sqrt(3))


class TestCrystalAngle:
    def test_cubic_90deg(self):
        c = Crystal(1.0, 1.0, 1.0, 90, 90, 90)
        angle = c.angle([1, 0, 0], [0, 1, 0], space='d')
        assert angle == pytest.approx(np.pi / 2)

    def test_cubic_45deg(self):
        c = Crystal(1.0, 1.0, 1.0, 90, 90, 90)
        angle = c.angle([1, 0, 0], [1, 1, 0], space='d')
        assert angle == pytest.approx(np.pi / 4)


class TestCrystalNormalize:
    def test_normalize(self):
        c = Crystal(1.0, 1.0, 1.0, 90, 90, 90)
        v = c.normalize([3, 0, 0], space='d')
        L = c.length(v, space='d')
        assert L == pytest.approx(1.0)


class TestCrystalCross:
    def test_cubic_cross(self):
        c = Crystal(1.0, 1.0, 1.0, 90, 90, 90)
        r = c.cross([1, 0, 0], [0, 1, 0], in_space='d', out_space='d')
        # In cubic, cross product of [100] x [010] = [001]
        np.testing.assert_allclose(r / r[2], [0, 0, 1], atol=1e-14)


class TestSetupAtoms:
    def test_setup_ni(self):
        ni = Crystal(0.35236, 0.35236, 0.35236, 90, 90, 90)
        sg = SpaceGroup(225)
        # Should not raise
        ni.setup_atoms(sg, [(28, 0.0, 0.0, 0.0, 1.0, 0.003529)])

    def test_setup_multi_atom(self):
        # NaCl: two atom types
        nacl = Crystal(0.5640, 0.5640, 0.5640, 90, 90, 90)
        sg = SpaceGroup(225)
        nacl.setup_atoms(sg, [
            (11, 0.0, 0.0, 0.0, 1.0, 0.005),   # Na
            (17, 0.5, 0.0, 0.0, 1.0, 0.005),   # Cl
        ])


class TestInterplanarSpacing:
    def test_cubic_d_spacing(self):
        a = 0.35236  # nm
        c = Crystal(a, a, a, 90, 90, 90)
        d111 = c.interplanar_spacing([1, 1, 1])
        expected = a / np.sqrt(3)
        assert d111 == pytest.approx(expected, rel=1e-6)
