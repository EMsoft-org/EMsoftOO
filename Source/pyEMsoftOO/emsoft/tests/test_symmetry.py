"""Tests for the emsoft.symmetry module."""

import numpy as np
import pytest
from emsoft.symmetry import SpaceGroup


class TestSpaceGroupCreation:
    def test_create(self):
        sg = SpaceGroup(225)
        assert sg.number == 225

    def test_invalid_number(self):
        with pytest.raises(ValueError):
            SpaceGroup(0)
        with pytest.raises(ValueError):
            SpaceGroup(231)

    def test_repr(self):
        sg = SpaceGroup(225)
        assert '225' in repr(sg)
        assert 'Cubic' in repr(sg)


class TestSpaceGroupProperties:
    def test_cubic_system(self):
        sg = SpaceGroup(225)  # Fm-3m
        assert sg.crystal_system == 'Cubic'
        assert sg.crystal_system_number == 1

    def test_hexagonal_system(self):
        sg = SpaceGroup(194)  # P6_3/mmc
        assert sg.crystal_system == 'Hexagonal'

    def test_centrosymmetric(self):
        sg = SpaceGroup(225)  # Fm-3m is centrosymmetric
        assert sg.is_centrosymmetric is True

    def test_order(self):
        sg = SpaceGroup(225)
        assert sg.order > 0

    def test_matnum(self):
        sg = SpaceGroup(225)
        assert sg.n_matrices > 0

    def test_numpt(self):
        sg = SpaceGroup(225)
        assert sg.n_point_group_ops > 0


class TestReflectionAllowed:
    def test_fcc_allowed(self):
        sg = SpaceGroup(225)  # Fm-3m (FCC)
        # In FCC, all-odd or all-even Miller indices are allowed
        assert sg.is_reflection_allowed([1, 1, 1]) is True
        assert sg.is_reflection_allowed([2, 0, 0]) is True
        assert sg.is_reflection_allowed([2, 2, 0]) is True

    def test_fcc_forbidden(self):
        sg = SpaceGroup(225)  # Fm-3m (FCC)
        # Mixed odd/even are forbidden in FCC
        assert sg.is_reflection_allowed([1, 0, 0]) is False
        assert sg.is_reflection_allowed([1, 1, 0]) is False


class TestCalcOrbit:
    def test_general_position(self):
        sg = SpaceGroup(225)
        orbit = sg.calc_orbit([0.123, 0.456, 0.789])
        assert orbit.shape[1] == 3
        assert orbit.shape[0] > 1  # should have multiple equivalent positions

    def test_special_position(self):
        sg = SpaceGroup(225)
        # Origin is a special position
        orbit = sg.calc_orbit([0.0, 0.0, 0.0])
        assert orbit.shape[0] >= 1


class TestCalcStar:
    def test_cubic_star(self):
        sg = SpaceGroup(225)
        star = sg.calc_star([1, 0, 0], space='r')
        assert star.shape[1] == 3
        # {100} family in cubic has 6 members
        assert star.shape[0] == 6

    def test_111_star(self):
        sg = SpaceGroup(225)
        star = sg.calc_star([1, 1, 1], space='r')
        # {111} in cubic has 8 members
        assert star.shape[0] == 8


class TestCalcFamily:
    def test_cubic_family(self):
        sg = SpaceGroup(225)
        family = sg.calc_family([1, 0, 0], space='r')
        assert family.shape[1] == 3
        assert family.shape[0] == 6  # {100} has 6 members

    def test_111_family(self):
        sg = SpaceGroup(225)
        family = sg.calc_family([1, 1, 1], space='r')
        assert family.shape[0] == 8  # {111} has 8 members
