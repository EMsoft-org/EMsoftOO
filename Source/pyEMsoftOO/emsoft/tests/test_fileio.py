"""Tests for the emsoft.fileio module and Crystal.from_file()."""

import os
import numpy as np
import pytest

# Skip all tests if h5py is not installed
h5py = pytest.importorskip("h5py")

from emsoft.fileio import read_xtal
from emsoft.crystallography import Crystal
from emsoft.diffraction import Diffraction

# Locate test .xtal files relative to the repo root
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
_TRICLINIC_XTAL = os.path.join(_REPO_ROOT, 'resources', 'test_triclinic.xtal')
_NI_XTAL = os.path.join(_REPO_ROOT, '..', 'EMsoftData', 'DItutorial', 'Ni', 'Ni.xtal')


def _have_triclinic():
    return os.path.isfile(_TRICLINIC_XTAL)


def _have_ni():
    return os.path.isfile(_NI_XTAL)


class TestReadXtal:
    @pytest.mark.skipif(not _have_triclinic(), reason="test_triclinic.xtal not found")
    def test_read_triclinic(self):
        data = read_xtal(_TRICLINIC_XTAL)
        assert 'lattice_parameters' in data
        assert 'space_group_number' in data
        assert 'atom_types' in data
        assert 'atom_data' in data
        assert data['lattice_parameters'].shape == (6,)
        assert 1 <= data['space_group_number'] <= 230
        assert 1 <= data['crystal_system'] <= 7

    @pytest.mark.skipif(not _have_ni(), reason="Ni.xtal not found")
    def test_read_ni(self):
        data = read_xtal(_NI_XTAL)
        assert data['space_group_number'] == 225
        assert data['crystal_system_name'] == 'Cubic'
        assert data['n_atom_types'] == 1
        assert data['atom_types'][0] == 28  # Ni
        lp = data['lattice_parameters']
        assert lp[0] == pytest.approx(0.35236, rel=1e-3)  # a in nm


class TestCrystalFromFile:
    @pytest.mark.skipif(not _have_ni(), reason="Ni.xtal not found")
    def test_from_file_ni(self):
        ni = Crystal.from_file(_NI_XTAL)
        lp = ni.lattice_parameters
        assert lp[0] == pytest.approx(0.35236, rel=1e-3)
        assert ni.volume > 0

    @pytest.mark.skipif(not _have_ni(), reason="Ni.xtal not found")
    def test_from_file_diffraction(self):
        """Test that a crystal loaded from file works with Diffraction."""
        ni = Crystal.from_file(_NI_XTAL)
        d = Diffraction(200.0, ni)
        assert 0.002 < d.wavelength < 0.003

        # 111 should be allowed in Ni (FCC)
        result = d.calc_structure_factor(ni, [1, 1, 1])
        assert result['xg'] > 0

    @pytest.mark.skipif(not _have_triclinic(), reason="test_triclinic.xtal not found")
    def test_from_file_triclinic(self):
        c = Crystal.from_file(_TRICLINIC_XTAL)
        assert c.volume > 0
        lp = c.lattice_parameters
        assert lp.shape == (6,)
