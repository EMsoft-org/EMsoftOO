"""
Rotation representations wrapping EMsoftOO mod_rotations.

Provides the Rotation class that can be created from any of 9 rotation
representations and converted to any other:

- Euler angles (Bunge convention: phi1, Phi, phi2)
- Quaternion [w, x, y, z]
- Rotation matrix (3x3)
- Axis-angle [n1, n2, n3, angle]
- Rodrigues vector [n1, n2, n3, tan(angle/2)]
- Homochoric vector
- Cubochoric vector
- Stereographic projection vector
- Rotation vector

All conversions use the EMsoftOO Fortran library for numerical accuracy.
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_void_p = ctypes.c_void_p


def _setup_bindings():
    lib = get_lib()

    # Constructors (one per representation)
    for name in ['euler', 'quaternion', 'axisangle', 'rodrigues',
                  'homochoric', 'cubochoric', 'stereographic', 'rotvec']:
        fn = getattr(lib, f'emsoft_rot_from_{name}')
        fn.restype = c_void_p

    lib.emsoft_rot_from_euler.argtypes = [c_double * 3]
    lib.emsoft_rot_from_quaternion.argtypes = [c_double * 4]
    lib.emsoft_rot_from_matrix.argtypes = [c_double * 9]
    lib.emsoft_rot_from_matrix.restype = c_void_p
    lib.emsoft_rot_from_axisangle.argtypes = [c_double * 4]
    lib.emsoft_rot_from_rodrigues.argtypes = [c_double * 4]
    lib.emsoft_rot_from_homochoric.argtypes = [c_double * 3]
    lib.emsoft_rot_from_cubochoric.argtypes = [c_double * 3]
    lib.emsoft_rot_from_stereographic.argtypes = [c_double * 3]
    lib.emsoft_rot_from_rotvec.argtypes = [c_double * 3]

    lib.emsoft_rot_destroy.argtypes = [c_void_p]
    lib.emsoft_rot_destroy.restype = None

    # Extractors
    for name in ['euler', 'homochoric', 'cubochoric', 'stereographic', 'rotvec']:
        fn = getattr(lib, f'emsoft_rot_to_{name}')
        fn.argtypes = [c_void_p, c_double * 3]
        fn.restype = None

    for name in ['quaternion', 'axisangle', 'rodrigues']:
        fn = getattr(lib, f'emsoft_rot_to_{name}')
        fn.argtypes = [c_void_p, c_double * 4]
        fn.restype = None

    lib.emsoft_rot_to_matrix.argtypes = [c_void_p, c_double * 9]
    lib.emsoft_rot_to_matrix.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


class Rotation:
    """A rotation backed by the EMsoftOO Fortran library.

    Create from any representation, convert to any other. The Fortran
    Orientation_T object computes all representations simultaneously.

    Examples
    --------
    >>> r = Rotation.from_euler(30, 45, 60, degrees=True)
    >>> r.to_quaternion()
    array([...])
    >>> r.to_matrix()
    array([[...], [...], [...]])

    >>> r = Rotation.from_quaternion(0.5, 0.5, 0.5, 0.5)
    >>> r.to_euler(degrees=True)
    array([...])
    """

    __slots__ = ('_handle', '_lib')

    def __init__(self, _handle):
        self._lib = _ensure_bindings()
        self._handle = _handle

    def __del__(self):
        if hasattr(self, '_handle') and self._handle is not None:
            self._lib.emsoft_rot_destroy(self._handle)
            self._handle = None

    # --- Constructors ---

    @classmethod
    def from_euler(cls, phi1, Phi, phi2, degrees=False):
        """Create from Euler angles (Bunge convention).

        Parameters
        ----------
        phi1, Phi, phi2 : float
            Euler angles.
        degrees : bool
            If True, angles are in degrees (converted to radians internally).
        """
        lib = _ensure_bindings()
        if degrees:
            phi1, Phi, phi2 = np.radians(phi1), np.radians(Phi), np.radians(phi2)
        arr = (c_double * 3)(phi1, Phi, phi2)
        h = lib.emsoft_rot_from_euler(arr)
        return cls(h)

    @classmethod
    def from_quaternion(cls, w, x, y, z):
        """Create from unit quaternion [w, x, y, z]."""
        lib = _ensure_bindings()
        arr = (c_double * 4)(w, x, y, z)
        h = lib.emsoft_rot_from_quaternion(arr)
        return cls(h)

    @classmethod
    def from_matrix(cls, om):
        """Create from a 3x3 rotation matrix."""
        lib = _ensure_bindings()
        om = np.asarray(om, dtype=np.float64)
        # Pass as column-major for Fortran
        arr = (c_double * 9)(*om.T.flatten())
        h = lib.emsoft_rot_from_matrix(arr)
        return cls(h)

    @classmethod
    def from_axisangle(cls, axis, angle, degrees=False):
        """Create from axis-angle pair.

        Parameters
        ----------
        axis : array_like of length 3
            Rotation axis (will be normalized).
        angle : float
            Rotation angle.
        degrees : bool
            If True, angle is in degrees.
        """
        lib = _ensure_bindings()
        axis = np.asarray(axis, dtype=np.float64)
        axis = axis / np.linalg.norm(axis)
        if degrees:
            angle = np.radians(angle)
        arr = (c_double * 4)(axis[0], axis[1], axis[2], angle)
        h = lib.emsoft_rot_from_axisangle(arr)
        return cls(h)

    @classmethod
    def from_rodrigues(cls, ro):
        """Create from Rodrigues vector [n1, n2, n3, tan(angle/2)]."""
        lib = _ensure_bindings()
        arr = (c_double * 4)(*ro)
        h = lib.emsoft_rot_from_rodrigues(arr)
        return cls(h)

    @classmethod
    def from_homochoric(cls, ho):
        """Create from homochoric vector [h1, h2, h3]."""
        lib = _ensure_bindings()
        arr = (c_double * 3)(*ho)
        h = lib.emsoft_rot_from_homochoric(arr)
        return cls(h)

    @classmethod
    def from_cubochoric(cls, cu):
        """Create from cubochoric vector [c1, c2, c3]."""
        lib = _ensure_bindings()
        arr = (c_double * 3)(*cu)
        h = lib.emsoft_rot_from_cubochoric(arr)
        return cls(h)

    @classmethod
    def from_stereographic(cls, st):
        """Create from stereographic projection vector [s1, s2, s3]."""
        lib = _ensure_bindings()
        arr = (c_double * 3)(*st)
        h = lib.emsoft_rot_from_stereographic(arr)
        return cls(h)

    @classmethod
    def from_rotvec(cls, rv):
        """Create from rotation vector [v1, v2, v3] (direction is axis, magnitude is angle)."""
        lib = _ensure_bindings()
        arr = (c_double * 3)(*rv)
        h = lib.emsoft_rot_from_rotvec(arr)
        return cls(h)

    @classmethod
    def identity(cls):
        """Create the identity rotation."""
        return cls.from_quaternion(1.0, 0.0, 0.0, 0.0)

    # --- Extractors ---

    def to_euler(self, degrees=False):
        """Return Euler angles [phi1, Phi, phi2] in radians (or degrees)."""
        out = (c_double * 3)()
        self._lib.emsoft_rot_to_euler(self._handle, out)
        eu = np.array(out, dtype=np.float64)
        if degrees:
            eu = np.degrees(eu)
        return eu

    def to_quaternion(self):
        """Return unit quaternion [w, x, y, z]."""
        out = (c_double * 4)()
        self._lib.emsoft_rot_to_quaternion(self._handle, out)
        return np.array(out, dtype=np.float64)

    def to_matrix(self):
        """Return 3x3 rotation matrix."""
        out = (c_double * 9)()
        self._lib.emsoft_rot_to_matrix(self._handle, out)
        # Fortran returns column-major; reshape and transpose
        return np.array(out, dtype=np.float64).reshape(3, 3).T

    def to_axisangle(self, degrees=False):
        """Return axis-angle [n1, n2, n3, angle]."""
        out = (c_double * 4)()
        self._lib.emsoft_rot_to_axisangle(self._handle, out)
        ax = np.array(out, dtype=np.float64)
        if degrees:
            ax[3] = np.degrees(ax[3])
        return ax

    def to_rodrigues(self):
        """Return Rodrigues vector [n1, n2, n3, tan(angle/2)]."""
        out = (c_double * 4)()
        self._lib.emsoft_rot_to_rodrigues(self._handle, out)
        return np.array(out, dtype=np.float64)

    def to_homochoric(self):
        """Return homochoric vector [h1, h2, h3]."""
        out = (c_double * 3)()
        self._lib.emsoft_rot_to_homochoric(self._handle, out)
        return np.array(out, dtype=np.float64)

    def to_cubochoric(self):
        """Return cubochoric vector [c1, c2, c3]."""
        out = (c_double * 3)()
        self._lib.emsoft_rot_to_cubochoric(self._handle, out)
        return np.array(out, dtype=np.float64)

    def to_stereographic(self):
        """Return stereographic projection vector [s1, s2, s3]."""
        out = (c_double * 3)()
        self._lib.emsoft_rot_to_stereographic(self._handle, out)
        return np.array(out, dtype=np.float64)

    def to_rotvec(self):
        """Return rotation vector [v1, v2, v3]."""
        out = (c_double * 3)()
        self._lib.emsoft_rot_to_rotvec(self._handle, out)
        return np.array(out, dtype=np.float64)

    def __repr__(self):
        eu = self.to_euler(degrees=True)
        return f'Rotation(euler_deg=[{eu[0]:.2f}, {eu[1]:.2f}, {eu[2]:.2f}])'
