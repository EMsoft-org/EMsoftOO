"""
Crystallography operations wrapping EMsoftOO mod_crystallography.

Provides the Crystal class for unit cell operations: coordinate
transformations, metric tensor calculations, dot products, angles,
lengths, and cross products in direct, reciprocal, and Cartesian spaces.
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_char = ctypes.c_char
c_void_p = ctypes.c_void_p


def _setup_bindings():
    lib = get_lib()

    lib.emsoft_cell_create.argtypes = [c_double * 6]
    lib.emsoft_cell_create.restype = c_void_p

    lib.emsoft_cell_destroy.argtypes = [c_void_p]
    lib.emsoft_cell_destroy.restype = None

    lib.emsoft_cell_get_latparm.argtypes = [c_void_p, c_double * 6]
    lib.emsoft_cell_get_latparm.restype = None

    lib.emsoft_cell_get_volume.argtypes = [c_void_p]
    lib.emsoft_cell_get_volume.restype = c_double

    for name in ['dmt', 'rmt', 'dsm', 'rsm']:
        fn = getattr(lib, f'emsoft_cell_get_{name}')
        fn.argtypes = [c_void_p, c_double * 9]
        fn.restype = None

    lib.emsoft_cell_trans_space.argtypes = [c_void_p, c_double * 3, c_double * 3,
                                            c_char, c_char]
    lib.emsoft_cell_trans_space.restype = None

    lib.emsoft_cell_calc_dot.argtypes = [c_void_p, c_double * 3, c_double * 3, c_char]
    lib.emsoft_cell_calc_dot.restype = c_double

    lib.emsoft_cell_calc_length.argtypes = [c_void_p, c_double * 3, c_char]
    lib.emsoft_cell_calc_length.restype = c_double

    lib.emsoft_cell_calc_angle.argtypes = [c_void_p, c_double * 3, c_double * 3, c_char]
    lib.emsoft_cell_calc_angle.restype = c_double

    lib.emsoft_cell_norm_vec.argtypes = [c_void_p, c_double * 3, c_char]
    lib.emsoft_cell_norm_vec.restype = None

    lib.emsoft_cell_calc_cross.argtypes = [c_void_p, c_double * 3, c_double * 3,
                                           c_double * 3, c_char, c_char]
    lib.emsoft_cell_calc_cross.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


def _space_char(space):
    """Convert a space string to a c_char."""
    s = space[0].lower()
    if s not in ('d', 'r', 'c'):
        raise ValueError(f"space must be 'd' (direct), 'r' (reciprocal), or 'c' (Cartesian), got '{space}'")
    return s.encode('ascii')


class Crystal:
    """A crystallographic unit cell backed by the EMsoftOO Fortran library.

    Parameters
    ----------
    a, b, c : float
        Lattice parameters in nanometers.
    alpha, beta, gamma : float
        Lattice angles in degrees.

    Examples
    --------
    >>> # Cubic Ni (a = 0.35236 nm)
    >>> ni = Crystal(0.35236, 0.35236, 0.35236, 90, 90, 90)
    >>> ni.volume
    0.04376...

    >>> # Transform [1,1,1] from direct to Cartesian
    >>> ni.transform([1, 1, 1], 'd', 'c')
    array([...])

    >>> # Interplanar angle
    >>> ni.angle([1,0,0], [1,1,0], space='r')
    0.7853...
    """

    __slots__ = ('_handle', '_lib')

    def __init__(self, a, b, c, alpha, beta, gamma):
        self._lib = _ensure_bindings()
        arr = (c_double * 6)(a, b, c, alpha, beta, gamma)
        self._handle = self._lib.emsoft_cell_create(arr)

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_cell_destroy(self._handle)
            self._handle = None

    # --- Properties ---

    @property
    def lattice_parameters(self):
        """Return [a, b, c, alpha, beta, gamma]."""
        out = (c_double * 6)()
        self._lib.emsoft_cell_get_latparm(self._handle, out)
        return np.array(out, dtype=np.float64)

    @property
    def volume(self):
        """Unit cell volume in nm^3."""
        return self._lib.emsoft_cell_get_volume(self._handle)

    @property
    def direct_metric_tensor(self):
        """Direct space metric tensor (3x3)."""
        out = (c_double * 9)()
        self._lib.emsoft_cell_get_dmt(self._handle, out)
        return np.array(out, dtype=np.float64).reshape(3, 3).T

    @property
    def reciprocal_metric_tensor(self):
        """Reciprocal space metric tensor (3x3)."""
        out = (c_double * 9)()
        self._lib.emsoft_cell_get_rmt(self._handle, out)
        return np.array(out, dtype=np.float64).reshape(3, 3).T

    @property
    def direct_structure_matrix(self):
        """Direct structure matrix (3x3)."""
        out = (c_double * 9)()
        self._lib.emsoft_cell_get_dsm(self._handle, out)
        return np.array(out, dtype=np.float64).reshape(3, 3).T

    @property
    def reciprocal_structure_matrix(self):
        """Reciprocal structure matrix (3x3)."""
        out = (c_double * 9)()
        self._lib.emsoft_cell_get_rsm(self._handle, out)
        return np.array(out, dtype=np.float64).reshape(3, 3).T

    # --- Computations ---

    def transform(self, v, from_space, to_space):
        """Transform a 3-vector between coordinate systems.

        Parameters
        ----------
        v : array_like of length 3
        from_space, to_space : str
            'd' (direct), 'r' (reciprocal), or 'c' (Cartesian).
        """
        vin = (c_double * 3)(*v)
        vout = (c_double * 3)()
        self._lib.emsoft_cell_trans_space(self._handle, vin, vout,
                                          _space_char(from_space),
                                          _space_char(to_space))
        return np.array(vout, dtype=np.float64)

    def dot(self, p, q, space='d'):
        """Dot product of two vectors in the given space."""
        pp = (c_double * 3)(*p)
        qq = (c_double * 3)(*q)
        return self._lib.emsoft_cell_calc_dot(self._handle, pp, qq,
                                              _space_char(space))

    def length(self, v, space='d'):
        """Length of a vector in the given space."""
        vv = (c_double * 3)(*v)
        return self._lib.emsoft_cell_calc_length(self._handle, vv,
                                                 _space_char(space))

    def angle(self, p, q, space='d'):
        """Angle (radians) between two vectors in the given space."""
        pp = (c_double * 3)(*p)
        qq = (c_double * 3)(*q)
        return self._lib.emsoft_cell_calc_angle(self._handle, pp, qq,
                                                _space_char(space))

    def normalize(self, v, space='d'):
        """Return a normalized copy of vector v in the given space."""
        vv = (c_double * 3)(*v)
        self._lib.emsoft_cell_norm_vec(self._handle, vv, _space_char(space))
        return np.array(vv, dtype=np.float64)

    def cross(self, p, q, in_space='d', out_space='d'):
        """Cross product of two vectors.

        Parameters
        ----------
        p, q : array_like of length 3
        in_space : str
            Space of input vectors.
        out_space : str
            Space of output vector.
        """
        pp = (c_double * 3)(*p)
        qq = (c_double * 3)(*q)
        rr = (c_double * 3)()
        self._lib.emsoft_cell_calc_cross(self._handle, pp, qq, rr,
                                         _space_char(in_space),
                                         _space_char(out_space))
        return np.array(rr, dtype=np.float64)

    def interplanar_spacing(self, hkl):
        """Compute d-spacing for reflection (h, k, l) in nm."""
        return 1.0 / self.length(hkl, space='r')

    def __repr__(self):
        lp = self.lattice_parameters
        return (f'Crystal(a={lp[0]:.5f}, b={lp[1]:.5f}, c={lp[2]:.5f}, '
                f'alpha={lp[3]:.2f}, beta={lp[4]:.2f}, gamma={lp[5]:.2f})')
