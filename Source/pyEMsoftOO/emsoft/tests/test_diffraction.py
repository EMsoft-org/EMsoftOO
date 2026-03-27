"""Tests for the emsoft.diffraction module."""

import numpy as np
import pytest
from emsoft.diffraction import Diffraction


class TestDiffractionCreation:
    def test_create(self):
        d = Diffraction(200.0)
        assert d.voltage == pytest.approx(200.0)

    def test_repr(self):
        d = Diffraction(300.0)
        assert '300.0' in repr(d)


class TestDiffractionProperties:
    def test_wavelength_200kv(self):
        d = Diffraction(200.0)
        # 200 keV electrons: lambda ~ 0.00251 nm
        assert 0.002 < d.wavelength < 0.003

    def test_wavelength_decreases_with_voltage(self):
        d100 = Diffraction(100.0)
        d200 = Diffraction(200.0)
        d300 = Diffraction(300.0)
        assert d100.wavelength > d200.wavelength > d300.wavelength

    def test_relativistic_correction(self):
        d = Diffraction(200.0)
        # gamma > 1 for relativistic electrons
        assert d.relativistic_correction > 1.0

    def test_sigma_positive(self):
        d = Diffraction(200.0)
        assert d.sigma > 0

    def test_psihat(self):
        d = Diffraction(200.0)
        # Psihat should be larger than the accelerating voltage (in V)
        assert d.psihat > 200000


class TestDiffractionMethod:
    def test_set_method(self):
        d = Diffraction(200.0)
        d.set_method('WK')  # should not raise
        d.set_method('DT')  # should not raise

    def test_invalid_method(self):
        d = Diffraction(200.0)
        with pytest.raises(ValueError):
            d.set_method('XX')
