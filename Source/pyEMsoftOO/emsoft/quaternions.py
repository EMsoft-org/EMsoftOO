"""
Quaternion operations wrapping EMsoftOO mod_quaternions.

Provides the Quaternion class for single quaternion operations and
QuaternionArray for operating on arrays of quaternions efficiently.

All quaternions use the convention [w, x, y, z] where w is the scalar part.
"""

import ctypes
import numpy as np
from ._lib import get_lib

# C type aliases
c_double = ctypes.c_double
c_double_p = ctypes.POINTER(c_double)
c_int = ctypes.c_int
c_bool = ctypes.c_bool
c_void_p = ctypes.c_void_p


def _setup_bindings():
    """Set up ctypes function signatures for the C-interop layer."""
    lib = get_lib()

    # --- Quaternion_T ---

    # Constructors / destructors
    lib.emsoft_quat_create.argtypes = [c_double * 4]
    lib.emsoft_quat_create.restype = c_void_p

    lib.emsoft_quat_create_identity.argtypes = []
    lib.emsoft_quat_create_identity.restype = c_void_p

    lib.emsoft_quat_destroy.argtypes = [c_void_p]
    lib.emsoft_quat_destroy.restype = None

    # Getters / setters
    lib.emsoft_quat_get.argtypes = [c_void_p, c_double * 4]
    lib.emsoft_quat_get.restype = None

    lib.emsoft_quat_set.argtypes = [c_void_p, c_double * 4]
    lib.emsoft_quat_set.restype = None

    # Arithmetic (return new handle)
    lib.emsoft_quat_add.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_add.restype = c_void_p

    lib.emsoft_quat_subtract.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_subtract.restype = c_void_p

    lib.emsoft_quat_multiply.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_multiply.restype = c_void_p

    lib.emsoft_quat_divide.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_divide.restype = c_void_p

    lib.emsoft_quat_scale.argtypes = [c_void_p, c_double]
    lib.emsoft_quat_scale.restype = c_void_p

    # Conjugate, norm, normalization
    lib.emsoft_quat_conjugate.argtypes = [c_void_p]
    lib.emsoft_quat_conjugate.restype = c_void_p

    lib.emsoft_quat_norm.argtypes = [c_void_p]
    lib.emsoft_quat_norm.restype = c_double

    lib.emsoft_quat_normalize.argtypes = [c_void_p]
    lib.emsoft_quat_normalize.restype = None

    lib.emsoft_quat_flip.argtypes = [c_void_p]
    lib.emsoft_quat_flip.restype = None

    lib.emsoft_quat_pos.argtypes = [c_void_p]
    lib.emsoft_quat_pos.restype = None

    # Geometric operations
    lib.emsoft_quat_innerproduct.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_innerproduct.restype = c_double

    lib.emsoft_quat_angle.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_angle.restype = c_double

    lib.emsoft_quat_equal.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quat_equal.restype = c_bool

    # Vector rotation
    lib.emsoft_quat_rotate_vector.argtypes = [c_void_p, c_double * 3, c_double * 3]
    lib.emsoft_quat_rotate_vector.restype = None

    lib.emsoft_quat_rotate_vecarray.argtypes = [c_void_p, c_int,
                                                 ctypes.c_void_p, ctypes.c_void_p]
    lib.emsoft_quat_rotate_vecarray.restype = None

    # --- QuaternionArray_T ---

    lib.emsoft_quatarray_create.argtypes = [c_int, ctypes.c_void_p]
    lib.emsoft_quatarray_create.restype = c_void_p

    lib.emsoft_quatarray_create_empty.argtypes = [c_int]
    lib.emsoft_quatarray_create_empty.restype = c_void_p

    lib.emsoft_quatarray_destroy.argtypes = [c_void_p]
    lib.emsoft_quatarray_destroy.restype = None

    lib.emsoft_quatarray_size.argtypes = [c_void_p]
    lib.emsoft_quatarray_size.restype = c_int

    lib.emsoft_quatarray_get_element.argtypes = [c_void_p, c_int, c_double * 4]
    lib.emsoft_quatarray_get_element.restype = None

    lib.emsoft_quatarray_set_element.argtypes = [c_void_p, c_int, c_double * 4]
    lib.emsoft_quatarray_set_element.restype = None

    lib.emsoft_quatarray_multiply.argtypes = [c_void_p, c_void_p]
    lib.emsoft_quatarray_multiply.restype = c_void_p

    lib.emsoft_quatarray_normalize.argtypes = [c_void_p]
    lib.emsoft_quatarray_normalize.restype = None

    lib.emsoft_quatarray_rotate_vector.argtypes = [c_void_p, c_double * 3,
                                                    ctypes.c_void_p, c_int]
    lib.emsoft_quatarray_rotate_vector.restype = None

    return lib


# Module-level lazy initialization
_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


class Quaternion:
    """A quaternion [w, x, y, z] backed by the EMsoftOO Fortran library.

    Parameters
    ----------
    w : float
        Scalar part (default 1.0).
    x, y, z : float
        Vector part (default 0.0).
    components : array_like of length 4, optional
        If provided, overrides w/x/y/z with [w, x, y, z].

    Examples
    --------
    >>> q = Quaternion(1, 0, 0, 0)          # identity
    >>> q = Quaternion(components=[1, 0, 0, 0])
    >>> q1 * q2                              # Hamilton product
    >>> q.conjugate()                        # returns new quaternion
    >>> q.rotate([1.0, 0.0, 0.0])            # rotate a vector
    """

    __slots__ = ('_handle', '_lib')

    def __init__(self, w=1.0, x=0.0, y=0.0, z=0.0, *, components=None, _handle=None):
        self._lib = _ensure_bindings()
        if _handle is not None:
            self._handle = _handle
        elif components is not None:
            arr = (c_double * 4)(*components)
            self._handle = self._lib.emsoft_quat_create(arr)
        else:
            arr = (c_double * 4)(w, x, y, z)
            self._handle = self._lib.emsoft_quat_create(arr)

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_quat_destroy(self._handle)
            self._handle = None

    @classmethod
    def identity(cls):
        """Create the identity quaternion [1, 0, 0, 0]."""
        lib = _ensure_bindings()
        h = lib.emsoft_quat_create_identity()
        return cls(_handle=h)

    # --- Properties ---

    @property
    def components(self):
        """Return quaternion components as a numpy array [w, x, y, z]."""
        out = (c_double * 4)()
        self._lib.emsoft_quat_get(self._handle, out)
        return np.array(out, dtype=np.float64)

    @components.setter
    def components(self, value):
        arr = (c_double * 4)(*value)
        self._lib.emsoft_quat_set(self._handle, arr)

    @property
    def w(self):
        return self.components[0]

    @property
    def x(self):
        return self.components[1]

    @property
    def y(self):
        return self.components[2]

    @property
    def z(self):
        return self.components[3]

    # --- Arithmetic ---

    def __add__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        h = self._lib.emsoft_quat_add(self._handle, other._handle)
        return Quaternion(_handle=h)

    def __sub__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        h = self._lib.emsoft_quat_subtract(self._handle, other._handle)
        return Quaternion(_handle=h)

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            h = self._lib.emsoft_quat_multiply(self._handle, other._handle)
            return Quaternion(_handle=h)
        elif isinstance(other, (int, float)):
            h = self._lib.emsoft_quat_scale(self._handle, c_double(float(other)))
            return Quaternion(_handle=h)
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            h = self._lib.emsoft_quat_scale(self._handle, c_double(float(other)))
            return Quaternion(_handle=h)
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Quaternion):
            h = self._lib.emsoft_quat_divide(self._handle, other._handle)
            return Quaternion(_handle=h)
        return NotImplemented

    def __eq__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        return bool(self._lib.emsoft_quat_equal(self._handle, other._handle))

    def __neg__(self):
        q = Quaternion(_handle=self._lib.emsoft_quat_scale(
            self._handle, c_double(-1.0)))
        return q

    # --- Quaternion operations ---

    def conjugate(self):
        """Return the quaternion conjugate [w, -x, -y, -z]."""
        h = self._lib.emsoft_quat_conjugate(self._handle)
        return Quaternion(_handle=h)

    def norm(self):
        """Return the quaternion norm (magnitude)."""
        return self._lib.emsoft_quat_norm(self._handle)

    def normalize(self):
        """Normalize the quaternion in-place to unit length."""
        self._lib.emsoft_quat_normalize(self._handle)
        return self

    def flip(self):
        """Negate all components in-place."""
        self._lib.emsoft_quat_flip(self._handle)
        return self

    def positive(self):
        """Ensure the scalar part is positive in-place (q or -q convention)."""
        self._lib.emsoft_quat_pos(self._handle)
        return self

    def inner(self, other):
        """Compute the inner (dot) product with another quaternion."""
        return self._lib.emsoft_quat_innerproduct(self._handle, other._handle)

    def angle(self, other):
        """Compute the angle (in radians) between two unit quaternions."""
        return self._lib.emsoft_quat_angle(self._handle, other._handle)

    def rotate(self, v):
        """Rotate a 3-vector by this unit quaternion: v' = q v q*.

        Parameters
        ----------
        v : array_like of length 3
            The vector to rotate.

        Returns
        -------
        numpy.ndarray
            The rotated vector.
        """
        vin = (c_double * 3)(*v)
        vout = (c_double * 3)()
        self._lib.emsoft_quat_rotate_vector(self._handle, vin, vout)
        return np.array(vout, dtype=np.float64)

    def rotate_vectors(self, v):
        """Rotate an array of 3-vectors by this unit quaternion.

        Parameters
        ----------
        v : numpy.ndarray of shape (n, 3)
            The vectors to rotate.

        Returns
        -------
        numpy.ndarray of shape (n, 3)
            The rotated vectors.
        """
        v = np.asarray(v, dtype=np.float64)
        n = v.shape[0]
        # Fortran expects (3, n) column-major
        v_f = np.asfortranarray(v.T)
        vout_f = np.empty((3, n), dtype=np.float64, order='F')
        self._lib.emsoft_quat_rotate_vecarray(
            self._handle, c_int(n),
            v_f.ctypes.data_as(c_void_p),
            vout_f.ctypes.data_as(c_void_p))
        return vout_f.T.copy()

    # --- Representations ---

    def __repr__(self):
        c = self.components
        return f'Quaternion({c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}, {c[3]:.6f})'

    def __str__(self):
        c = self.components
        return f'[{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}, {c[3]:.6f}]'

    def to_array(self):
        """Return components as a numpy array [w, x, y, z]."""
        return self.components


class QuaternionArray:
    """An array of quaternions backed by the EMsoftOO Fortran library.

    Parameters
    ----------
    data : numpy.ndarray of shape (n, 4), optional
        Array of quaternions, each row [w, x, y, z].
    n : int, optional
        Create an empty array with n slots.

    Examples
    --------
    >>> qa = QuaternionArray(np.array([[1,0,0,0], [0,1,0,0]], dtype=np.float64))
    >>> qa[0]  # returns Quaternion
    >>> len(qa)
    2
    """

    __slots__ = ('_handle', '_lib')

    def __init__(self, data=None, *, n=None, _handle=None):
        self._lib = _ensure_bindings()
        if _handle is not None:
            self._handle = _handle
        elif data is not None:
            data = np.asarray(data, dtype=np.float64)
            nq = data.shape[0]
            # Fortran expects (4, n) column-major
            data_f = np.asfortranarray(data.T)
            self._handle = self._lib.emsoft_quatarray_create(
                c_int(nq), data_f.ctypes.data_as(c_void_p))
        elif n is not None:
            self._handle = self._lib.emsoft_quatarray_create_empty(c_int(n))
        else:
            raise ValueError("Provide either 'data' or 'n'")

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_quatarray_destroy(self._handle)
            self._handle = None

    def __len__(self):
        return self._lib.emsoft_quatarray_size(self._handle)

    def __getitem__(self, i):
        """Get the i-th quaternion (0-based index)."""
        if i < 0:
            i += len(self)
        if i < 0 or i >= len(self):
            raise IndexError(f'index {i} out of range for array of length {len(self)}')
        out = (c_double * 4)()
        self._lib.emsoft_quatarray_get_element(self._handle, c_int(i + 1), out)
        return Quaternion(components=list(out))

    def __setitem__(self, i, q):
        """Set the i-th quaternion (0-based index)."""
        if i < 0:
            i += len(self)
        if isinstance(q, Quaternion):
            c = q.components
        else:
            c = np.asarray(q, dtype=np.float64)
        arr = (c_double * 4)(*c)
        self._lib.emsoft_quatarray_set_element(self._handle, c_int(i + 1), arr)

    def __mul__(self, other):
        """Element-wise quaternion multiplication."""
        if not isinstance(other, QuaternionArray):
            return NotImplemented
        h = self._lib.emsoft_quatarray_multiply(self._handle, other._handle)
        return QuaternionArray(_handle=h)

    def normalize(self):
        """Normalize all quaternions in-place."""
        self._lib.emsoft_quatarray_normalize(self._handle)
        return self

    def rotate(self, v):
        """Rotate a single 3-vector by each quaternion in the array.

        Parameters
        ----------
        v : array_like of length 3
            The vector to rotate.

        Returns
        -------
        numpy.ndarray of shape (n, 3)
            Each row is the vector rotated by the corresponding quaternion.
        """
        n = len(self)
        vin = (c_double * 3)(*v)
        vout_f = np.empty((3, n), dtype=np.float64, order='F')
        self._lib.emsoft_quatarray_rotate_vector(
            self._handle, vin, vout_f.ctypes.data_as(c_void_p), c_int(n))
        return vout_f.T.copy()

    def to_array(self):
        """Return all quaternions as a numpy array of shape (n, 4)."""
        n = len(self)
        result = np.empty((n, 4), dtype=np.float64)
        out = (c_double * 4)()
        for i in range(n):
            self._lib.emsoft_quatarray_get_element(self._handle, c_int(i + 1), out)
            result[i, :] = out
        return result

    def __repr__(self):
        return f'QuaternionArray(n={len(self)})'
