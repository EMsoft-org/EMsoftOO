"""Tests for the emsoft.diffraction module."""

import numpy as np
import pytest
from emsoft.diffraction import Diffraction
from emsoft.crystallography import Crystal
from emsoft.symmetry import SpaceGroup


def _make_ni_crystal():
    """Create a Ni crystal with atom positions for diffraction tests."""
    ni = Crystal(0.35236, 0.35236, 0.35236, 90, 90, 90)
    sg = SpaceGroup(225)  # Fm-3m
    # Ni: Z=28, at origin, full occupancy, DW=0.003529 nm^2
    ni.setup_atoms(sg, [(28, 0.0, 0.0, 0.0, 1.0, 0.003529)])
    return ni


class TestDiffractionCreation:
    def test_create(self):
        d = Diffraction(200.0, _make_ni_crystal())
        assert d.voltage == pytest.approx(200.0)

    def test_repr(self):
        d = Diffraction(300.0, _make_ni_crystal())
        assert '300.0' in repr(d)


class TestDiffractionProperties:
    def test_wavelength_200kv(self):
        d = Diffraction(200.0, _make_ni_crystal())
        # 200 keV electrons: lambda ~ 0.00251 nm
        assert 0.002 < d.wavelength < 0.003

    def test_wavelength_decreases_with_voltage(self):
        c = _make_ni_crystal()
        d100 = Diffraction(100.0, c)
        d200 = Diffraction(200.0, c)
        d300 = Diffraction(300.0, c)
        assert d100.wavelength > d200.wavelength > d300.wavelength

    def test_relativistic_correction(self):
        d = Diffraction(200.0, _make_ni_crystal())
        assert d.relativistic_correction > 1.0

    def test_sigma_positive(self):
        d = Diffraction(200.0, _make_ni_crystal())
        assert d.sigma > 0

    def test_psihat(self):
        d = Diffraction(200.0, _make_ni_crystal())
        assert d.psihat > 200000


class TestDiffractionMethod:
    def test_set_method(self):
        d = Diffraction(200.0, _make_ni_crystal())
        d.set_method('WK')
        d.set_method('DT')

    def test_invalid_method(self):
        d = Diffraction(200.0, _make_ni_crystal())
        with pytest.raises(ValueError):
            d.set_method('XX')


class TestStructureFactor:
    def test_calc_ucg_111(self):
        ni = _make_ni_crystal()
        d = Diffraction(200.0, ni)
        d.set_method('WK')
        result = d.calc_structure_factor(ni, [1, 1, 1])
        # 111 is allowed in FCC, should have finite extinction distance
        assert result['xg'] > 0

    def test_calc_ucg_forbidden(self):
        ni = _make_ni_crystal()
        d = Diffraction(200.0, ni)
        d.set_method('WK')
        result = d.calc_structure_factor(ni, [1, 0, 0])
        # 100 is forbidden in FCC, Ucg should be ~zero
        assert abs(result['Ucg_real']) < 1e-6
        assert abs(result['Ucg_imag']) < 1e-6
