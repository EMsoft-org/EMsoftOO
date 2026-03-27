"""Tests for the emsoft.diffraction module.

Note: The Diffraction class requires a fully initialized Crystal with
atom positions (loaded from a .xtal file via getCrystalData). Tests
using lattice-parameters-only Crystal objects will segfault because
the Diffraction_T constructor calls CalcUcg which needs atom data.

These tests are skipped until crystal file I/O is added to the Python
wrapper (Phase 4+).
"""

import pytest

pytestmark = pytest.mark.skip(
    reason="Diffraction requires Crystal with atom positions (needs .xtal file I/O)"
)
