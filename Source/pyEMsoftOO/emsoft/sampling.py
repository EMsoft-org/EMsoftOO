"""
SO(3) sampling and fundamental zone operations wrapping EMsoftOO mod_so3.

Provides the FundamentalZone class for testing whether orientations lie
inside the Rodrigues fundamental zone and computing MacKenzie distributions.
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_int = ctypes.c_int
c_bool = ctypes.c_bool
c_void_p = ctypes.c_void_p

# FZ type descriptions
FZ_TYPES = {
    0: 'no symmetry',
    1: 'cyclic',
    2: 'dihedral',
    3: 'tetrahedral',
    4: 'octahedral',
}


def _setup_bindings():
    lib = get_lib()

    lib.emsoft_so3_create.argtypes = [c_int]
    lib.emsoft_so3_create.restype = c_void_p

    lib.emsoft_so3_destroy.argtypes = [c_void_p]
    lib.emsoft_so3_destroy.restype = None

    lib.emsoft_so3_get_fz_type_order.argtypes = [c_void_p,
                                                   ctypes.POINTER(c_int),
                                                   ctypes.POINTER(c_int)]
    lib.emsoft_so3_get_fz_type_order.restype = None

    lib.emsoft_so3_is_inside_fz.argtypes = [c_void_p, c_double * 4]
    lib.emsoft_so3_is_inside_fz.restype = c_bool

    lib.emsoft_so3_mackenzie.argtypes = [c_void_p, c_int,
                                          c_void_p, c_void_p]
    lib.emsoft_so3_mackenzie.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


class FundamentalZone:
    """Rodrigues fundamental zone for a crystallographic point group.

    Parameters
    ----------
    pgnum : int
        Point group number (1-32).

    Examples
    --------
    >>> fz = FundamentalZone(32)  # m-3m (cubic)
    >>> fz.fz_type_name
    'octahedral'
    >>> fz.is_inside([0, 0, 1, 0.1])  # Rodrigues vector
    True
    """

    __slots__ = ('_handle', '_lib', '_pgnum')

    def __init__(self, pgnum):
        self._lib = _ensure_bindings()
        self._pgnum = pgnum
        self._handle = self._lib.emsoft_so3_create(c_int(pgnum))

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_so3_destroy(self._handle)
            self._handle = None

    @property
    def point_group_number(self):
        """Point group number (1-32)."""
        return self._pgnum

    @property
    def fz_type(self):
        """Fundamental zone type number."""
        t = c_int()
        o = c_int()
        self._lib.emsoft_so3_get_fz_type_order(self._handle,
                                                ctypes.byref(t), ctypes.byref(o))
        return t.value

    @property
    def fz_order(self):
        """Fundamental zone order."""
        t = c_int()
        o = c_int()
        self._lib.emsoft_so3_get_fz_type_order(self._handle,
                                                ctypes.byref(t), ctypes.byref(o))
        return o.value

    @property
    def fz_type_name(self):
        """Human-readable FZ type name."""
        return FZ_TYPES.get(self.fz_type, 'unknown')

    def is_inside(self, rodrigues):
        """Test if a Rodrigues vector is inside the fundamental zone.

        Parameters
        ----------
        rodrigues : array_like of length 4
            Rodrigues vector [n1, n2, n3, tan(angle/2)].

        Returns
        -------
        bool
        """
        rod = (c_double * 4)(*rodrigues)
        return bool(self._lib.emsoft_so3_is_inside_fz(self._handle, rod))

    def mackenzie_distribution(self, nsteps=500):
        """Compute the theoretical MacKenzie misorientation distribution.

        Parameters
        ----------
        nsteps : int
            Number of angle bins from 0 to max misorientation angle.

        Returns
        -------
        angles : numpy.ndarray
            Misorientation angles in degrees.
        distribution : numpy.ndarray
            Probability density values.
        """
        misor = np.linspace(0, np.pi, nsteps + 1)
        mk = np.zeros(nsteps + 1, dtype=np.float64)
        self._lib.emsoft_so3_mackenzie(
            self._handle, c_int(nsteps),
            misor.ctypes.data_as(c_void_p),
            mk.ctypes.data_as(c_void_p))
        return np.degrees(misor), mk

    def __repr__(self):
        return f'FundamentalZone(pg={self._pgnum}, type={self.fz_type_name})'
