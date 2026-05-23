"""Tests for the emsoft.sampling module."""

import numpy as np
import pytest
from emsoft.sampling import FundamentalZone


class TestFundamentalZoneCreation:
    def test_create_cubic(self):
        fz = FundamentalZone(32)  # m-3m
        assert fz.point_group_number == 32

    def test_repr(self):
        fz = FundamentalZone(32)
        assert '32' in repr(fz)


class TestFZProperties:
    def test_cubic_type(self):
        fz = FundamentalZone(32)  # m-3m -> octahedral
        assert fz.fz_type == 4  # octahedral
        assert fz.fz_type_name == 'octahedral'

    def test_cubic_order(self):
        fz = FundamentalZone(32)
        # Octahedral FZ has order 0 (symmetry encoded in type, not order)
        assert fz.fz_order >= 0


class TestIsInsideFZ:
    def test_identity_inside(self):
        fz = FundamentalZone(32)
        # Identity rotation: Rodrigues = [0, 0, 1, 0]
        assert fz.is_inside([0, 0, 1, 0]) is True

    def test_small_rotation_inside(self):
        fz = FundamentalZone(32)
        # Small rotation around z: should be inside cubic FZ
        assert fz.is_inside([0, 0, 1, 0.1]) is True

    def test_large_rotation_outside(self):
        fz = FundamentalZone(32)
        # Very large Rodrigues parameter: should be outside cubic FZ
        assert fz.is_inside([0, 0, 1, 100.0]) is False


class TestMacKenzie:
    def test_mackenzie_shape(self):
        fz = FundamentalZone(32)
        angles, dist = fz.mackenzie_distribution(nsteps=100)
        assert angles.shape == (101,)
        assert dist.shape == (101,)

    def test_mackenzie_starts_at_zero(self):
        fz = FundamentalZone(32)
        angles, dist = fz.mackenzie_distribution(nsteps=100)
        assert angles[0] == pytest.approx(0.0)

    def test_mackenzie_nonnegative(self):
        fz = FundamentalZone(32)
        _, dist = fz.mackenzie_distribution(nsteps=100)
        assert np.all(dist >= 0)
