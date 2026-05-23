"""
Symmetry operations wrapping EMsoftOO mod_symmetry.

Provides the SpaceGroup class for crystallographic symmetry:
computing equivalent positions (orbits), equivalent reflections
(families), stars, and checking systematic absences.
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_int = ctypes.c_int
c_char = ctypes.c_char
c_bool = ctypes.c_bool
c_void_p = ctypes.c_void_p

# Crystal system names indexed by number (1-7)
CRYSTAL_SYSTEMS = {
    1: 'Cubic',
    2: 'Tetragonal',
    3: 'Orthorhombic',
    4: 'Hexagonal',
    5: 'Trigonal',
    6: 'Monoclinic',
    7: 'Triclinic',
}


def _setup_bindings():
    lib = get_lib()

    lib.emsoft_sg_create.argtypes = [c_int]
    lib.emsoft_sg_create.restype = c_void_p

    lib.emsoft_sg_destroy.argtypes = [c_void_p]
    lib.emsoft_sg_destroy.restype = None

    lib.emsoft_sg_get_number.argtypes = [c_void_p]
    lib.emsoft_sg_get_number.restype = c_int

    lib.emsoft_sg_get_order.argtypes = [c_void_p]
    lib.emsoft_sg_get_order.restype = c_int

    lib.emsoft_sg_get_matnum.argtypes = [c_void_p]
    lib.emsoft_sg_get_matnum.restype = c_int

    lib.emsoft_sg_get_numpt.argtypes = [c_void_p]
    lib.emsoft_sg_get_numpt.restype = c_int

    lib.emsoft_sg_get_xtal_system.argtypes = [c_void_p]
    lib.emsoft_sg_get_xtal_system.restype = c_int

    lib.emsoft_sg_get_centro.argtypes = [c_void_p]
    lib.emsoft_sg_get_centro.restype = c_bool

    lib.emsoft_sg_get_symmorphic.argtypes = [c_void_p]
    lib.emsoft_sg_get_symmorphic.restype = c_bool

    lib.emsoft_sg_is_g_allowed.argtypes = [c_void_p, c_int * 3]
    lib.emsoft_sg_is_g_allowed.restype = c_bool

    # CalcOrbit: handle, site(3), n_out, ctmp(maxn,3), maxn
    lib.emsoft_sg_calc_orbit.argtypes = [c_void_p, c_double * 3,
                                          ctypes.POINTER(c_int), c_void_p, c_int]
    lib.emsoft_sg_calc_orbit.restype = None

    # CalcStar: handle, kk(3), n_out, stmp(maxn,3), space, maxn
    lib.emsoft_sg_calc_star.argtypes = [c_void_p, c_double * 3,
                                         ctypes.POINTER(c_int), c_void_p, c_char, c_int]
    lib.emsoft_sg_calc_star.restype = None

    # CalcFamily: handle, ind(3), num_out, itmp(maxn,3), space, maxn
    lib.emsoft_sg_calc_family.argtypes = [c_void_p, c_int * 3,
                                           ctypes.POINTER(c_int), c_void_p, c_char, c_int]
    lib.emsoft_sg_calc_family.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


class SpaceGroup:
    """A crystallographic space group backed by the EMsoftOO Fortran library.

    Parameters
    ----------
    number : int
        Space group number (1-230).

    Examples
    --------
    >>> sg = SpaceGroup(225)  # Fm-3m (FCC)
    >>> sg.crystal_system
    'Cubic'
    >>> sg.is_centrosymmetric
    True

    >>> # Equivalent positions of (0.25, 0.25, 0.25)
    >>> orbit = sg.calc_orbit([0.25, 0.25, 0.25])

    >>> # Family of {111} planes
    >>> family = sg.calc_family([1, 1, 1], space='r')

    >>> # Check if (100) is allowed
    >>> sg.is_reflection_allowed([1, 0, 0])
    False
    """

    __slots__ = ('_handle', '_lib')

    def __init__(self, number):
        if not 1 <= number <= 230:
            raise ValueError(f'Space group number must be 1-230, got {number}')
        self._lib = _ensure_bindings()
        self._handle = self._lib.emsoft_sg_create(c_int(number))

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_sg_destroy(self._handle)
            self._handle = None

    # --- Properties ---

    @property
    def number(self):
        """Space group number (1-230)."""
        return self._lib.emsoft_sg_get_number(self._handle)

    @property
    def order(self):
        """Order of the space group (number of symmetry operations)."""
        return self._lib.emsoft_sg_get_order(self._handle)

    @property
    def n_matrices(self):
        """Number of symmetry matrices."""
        return self._lib.emsoft_sg_get_matnum(self._handle)

    @property
    def n_point_group_ops(self):
        """Number of point group operators."""
        return self._lib.emsoft_sg_get_numpt(self._handle)

    @property
    def crystal_system_number(self):
        """Crystal system number (1-7)."""
        return self._lib.emsoft_sg_get_xtal_system(self._handle)

    @property
    def crystal_system(self):
        """Crystal system name."""
        return CRYSTAL_SYSTEMS.get(self.crystal_system_number, 'Unknown')

    @property
    def is_centrosymmetric(self):
        """Whether the space group is centrosymmetric."""
        return bool(self._lib.emsoft_sg_get_centro(self._handle))

    @property
    def is_symmorphic(self):
        """Whether the space group is symmorphic."""
        return bool(self._lib.emsoft_sg_get_symmorphic(self._handle))

    # --- Symmetry operations ---

    def is_reflection_allowed(self, hkl):
        """Check if a reflection (h, k, l) is allowed (not systematically absent).

        Parameters
        ----------
        hkl : array_like of 3 ints
            Miller indices.

        Returns
        -------
        bool
        """
        g = (c_int * 3)(*[int(x) for x in hkl])
        return bool(self._lib.emsoft_sg_is_g_allowed(self._handle, g))

    def calc_orbit(self, site):
        """Compute the orbit of a fractional coordinate position.

        Parameters
        ----------
        site : array_like of length 3
            Fractional coordinates.

        Returns
        -------
        numpy.ndarray of shape (n, 3)
            Equivalent positions.
        """
        maxn = 192
        s = (c_double * 3)(*site)
        n = c_int(0)
        ctmp = np.zeros((maxn, 3), dtype=np.float64, order='F')
        self._lib.emsoft_sg_calc_orbit(self._handle, s, ctypes.byref(n),
                                        ctmp.ctypes.data_as(c_void_p), c_int(maxn))
        return ctmp[:n.value, :].copy()

    def calc_star(self, kk, space='r'):
        """Compute the star of a lattice vector.

        Parameters
        ----------
        kk : array_like of length 3
            Lattice vector indices.
        space : str
            'd' (direct) or 'r' (reciprocal).

        Returns
        -------
        numpy.ndarray of shape (n, 3)
            Star vectors.
        """
        maxn = 48
        k = (c_double * 3)(*[float(x) for x in kk])
        n = c_int(0)
        stmp = np.zeros((maxn, 3), dtype=np.float64, order='F')
        self._lib.emsoft_sg_calc_star(self._handle, k, ctypes.byref(n),
                                       stmp.ctypes.data_as(c_void_p),
                                       space[0].encode('ascii'), c_int(maxn))
        return stmp[:n.value, :].copy()

    def calc_family(self, hkl, space='r'):
        """Compute the family of symmetry-equivalent planes or directions.

        Parameters
        ----------
        hkl : array_like of 3 ints
            Miller indices.
        space : str
            'r' (reciprocal / planes) or 'd' (direct / directions).

        Returns
        -------
        numpy.ndarray of shape (n, 3), dtype int
            Equivalent indices.
        """
        maxn = 48
        ind = (c_int * 3)(*[int(x) for x in hkl])
        num = c_int(0)
        itmp = np.zeros((maxn, 3), dtype=np.int32, order='F')
        self._lib.emsoft_sg_calc_family(self._handle, ind, ctypes.byref(num),
                                         itmp.ctypes.data_as(c_void_p),
                                         space[0].encode('ascii'), c_int(maxn))
        return itmp[:num.value, :].copy()

    @property
    def multiplicity(self):
        """Multiplicity of the general position."""
        return self.order

    def __repr__(self):
        return f'SpaceGroup({self.number}, system={self.crystal_system})'
