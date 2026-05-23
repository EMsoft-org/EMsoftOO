"""
Electron diffraction calculations wrapping EMsoftOO mod_diffraction.

Provides the Diffraction class for computing relativistic electron
wavelengths, structure factors, extinction distances, and absorption lengths.
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_int = ctypes.c_int
c_char = ctypes.c_char
c_void_p = ctypes.c_void_p

# Physical constants (matching EMsoftOO mod_global)
_PLANCK = 6.62607015e-34      # J s
_REST_MASS = 9.1093837015e-31  # kg
_CHARGE = 1.602176634e-19      # C


def _setup_bindings():
    lib = get_lib()

    lib.emsoft_diff_create.argtypes = [c_double, c_void_p]
    lib.emsoft_diff_create.restype = c_void_p

    lib.emsoft_diff_destroy.argtypes = [c_void_p]
    lib.emsoft_diff_destroy.restype = None

    lib.emsoft_diff_calc_wavelength.argtypes = [c_void_p, c_void_p]
    lib.emsoft_diff_calc_wavelength.restype = None

    lib.emsoft_diff_get_voltage.argtypes = [c_void_p]
    lib.emsoft_diff_get_voltage.restype = c_double

    lib.emsoft_diff_get_wavelength.argtypes = [c_void_p]
    lib.emsoft_diff_get_wavelength.restype = c_double

    lib.emsoft_diff_get_relcor.argtypes = [c_void_p]
    lib.emsoft_diff_get_relcor.restype = c_double

    lib.emsoft_diff_get_sigma.argtypes = [c_void_p]
    lib.emsoft_diff_get_sigma.restype = c_double

    lib.emsoft_diff_get_psihat.argtypes = [c_void_p]
    lib.emsoft_diff_get_psihat.restype = c_double

    lib.emsoft_diff_set_method.argtypes = [c_void_p, c_char, c_char]
    lib.emsoft_diff_set_method.restype = None

    lib.emsoft_diff_calc_ucg.argtypes = [c_void_p, c_void_p, c_int * 3,
                                          ctypes.POINTER(c_double),
                                          ctypes.POINTER(c_double),
                                          ctypes.POINTER(c_double),
                                          ctypes.POINTER(c_double),
                                          ctypes.POINTER(c_double)]
    lib.emsoft_diff_calc_ucg.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


class Diffraction:
    """Electron diffraction parameters backed by the EMsoftOO Fortran library.

    Parameters
    ----------
    voltage : float
        Accelerating voltage in keV.
    crystal : emsoft.crystallography.Crystal
        Crystal unit cell (needed for wavelength calculation).

    Examples
    --------
    >>> from emsoft.crystallography import Crystal
    >>> ni = Crystal(0.35236, 0.35236, 0.35236, 90, 90, 90)
    >>> d = Diffraction(200.0, ni)
    >>> d.voltage
    200.0
    >>> d.wavelength  # in nm
    0.002508...
    """

    __slots__ = ('_handle', '_lib')

    def __init__(self, voltage, crystal):
        self._lib = _ensure_bindings()
        self._handle = self._lib.emsoft_diff_create(c_double(voltage),
                                                     crystal._handle)

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_diff_destroy(self._handle)
            self._handle = None

    def calc_wavelength(self, crystal):
        """Compute wavelength with mean inner potential correction from a Crystal.

        Parameters
        ----------
        crystal : emsoft.crystallography.Crystal
            Crystal with atom positions for V0 correction.
        """
        self._lib.emsoft_diff_calc_wavelength(self._handle, crystal._handle)

    @property
    def voltage(self):
        """Accelerating voltage in keV."""
        return self._lib.emsoft_diff_get_voltage(self._handle)

    @property
    def wavelength(self):
        """Relativistic electron wavelength in nm."""
        return self._lib.emsoft_diff_get_wavelength(self._handle)

    @property
    def relativistic_correction(self):
        """Relativistic correction factor (gamma)."""
        return self._lib.emsoft_diff_get_relcor(self._handle)

    @property
    def sigma(self):
        """Interaction constant (V^-1 nm^-1)."""
        return self._lib.emsoft_diff_get_sigma(self._handle)

    @property
    def psihat(self):
        """Relativistically corrected accelerating potential (V)."""
        return self._lib.emsoft_diff_get_psihat(self._handle)

    def set_method(self, method):
        """Set scattering factor calculation method.

        Parameters
        ----------
        method : str
            'WK' (Weickenmeier-Kohl), 'DT' (Doyle-Turner), or 'XR' (X-ray).
        """
        m = method.upper()
        if m not in ('WK', 'DT', 'XR'):
            raise ValueError(f"method must be 'WK', 'DT', or 'XR', got '{method}'")
        self._lib.emsoft_diff_set_method(self._handle,
                                         m[0].encode('ascii'),
                                         m[1].encode('ascii'))

    def calc_structure_factor(self, crystal, hkl):
        """Compute the structure factor for a reflection.

        Parameters
        ----------
        crystal : emsoft.crystallography.Crystal
            Crystal with atom positions calculated.
        hkl : array_like of 3 ints
            Miller indices.

        Returns
        -------
        dict with keys:
            'xg' : float — extinction distance (nm)
            'xgp' : float — absorption length (nm)
            'Ucg_real' : float — real part of Ucg (nm^-2)
            'Ucg_imag' : float — imaginary part of Ucg (nm^-2)
            'Vphase' : float — phase of Vg (radians)
        """
        g = (c_int * 3)(*[int(x) for x in hkl])
        xg = c_double()
        xgp = c_double()
        ucg_r = c_double()
        ucg_i = c_double()
        vphase = c_double()
        self._lib.emsoft_diff_calc_ucg(self._handle, crystal._handle, g,
                                        ctypes.byref(xg), ctypes.byref(xgp),
                                        ctypes.byref(ucg_r), ctypes.byref(ucg_i),
                                        ctypes.byref(vphase))
        return {
            'xg': xg.value,
            'xgp': xgp.value,
            'Ucg_real': ucg_r.value,
            'Ucg_imag': ucg_i.value,
            'Vphase': vphase.value,
        }

    def __repr__(self):
        return f'Diffraction(voltage={self.voltage:.1f} keV)'
