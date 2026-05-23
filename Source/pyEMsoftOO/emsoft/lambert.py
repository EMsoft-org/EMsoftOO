"""
Lambert projection utilities wrapping EMsoftOO mod_Lambert.

Provides functions for mapping between:
- Square grids and hemisphere (2D Lambert)
- Cubic grids and unit ball (3D Lambert)
- Stereographic projection (forward and inverse)
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_int = ctypes.c_int
c_void_p = ctypes.c_void_p


def _setup_bindings():
    lib = get_lib()

    for name in ['emsoft_lambert_square_to_sphere', 'emsoft_lambert_sphere_to_square']:
        fn = getattr(lib, name)
        fn.argtypes = [c_double * 2 if 'square_to' in name else c_double * 3,
                       c_double * 3 if 'square_to' in name else c_double * 2,
                       ctypes.POINTER(c_int)]
        fn.restype = None

    lib.emsoft_lambert_square_to_sphere.argtypes = [c_double * 2, c_double * 3,
                                                     ctypes.POINTER(c_int)]
    lib.emsoft_lambert_sphere_to_square.argtypes = [c_double * 3, c_double * 2,
                                                     ctypes.POINTER(c_int)]

    for name in ['emsoft_lambert_cube_to_ball', 'emsoft_lambert_ball_to_cube']:
        fn = getattr(lib, name)
        fn.argtypes = [c_double * 3, c_double * 3, ctypes.POINTER(c_int)]
        fn.restype = None

    lib.emsoft_lambert_stereo_forward.argtypes = [c_double * 3, c_double * 2,
                                                   ctypes.POINTER(c_int)]
    lib.emsoft_lambert_stereo_forward.restype = None

    lib.emsoft_lambert_stereo_inverse.argtypes = [c_double * 2, c_double * 3,
                                                   ctypes.POINTER(c_int)]
    lib.emsoft_lambert_stereo_inverse.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


def _check_err(ierr, operation):
    if ierr != 0:
        raise ValueError(f'{operation} failed with error code {ierr} '
                         f'(input out of valid range)')


def square_to_sphere(xy):
    """Map a 2D square point to a 3D hemisphere point (Lambert projection).

    Parameters
    ----------
    xy : array_like of length 2
        Point in [-1, 1]^2.

    Returns
    -------
    numpy.ndarray of length 3
        Point on the unit hemisphere.
    """
    lib = _ensure_bindings()
    inp = (c_double * 2)(*xy)
    out = (c_double * 3)()
    ierr = c_int()
    lib.emsoft_lambert_square_to_sphere(inp, out, ctypes.byref(ierr))
    _check_err(ierr.value, 'square_to_sphere')
    return np.array(out, dtype=np.float64)


def sphere_to_square(xyz):
    """Map a hemisphere point to a 2D square point (inverse Lambert).

    Parameters
    ----------
    xyz : array_like of length 3
        Point on the unit hemisphere (z >= 0).

    Returns
    -------
    numpy.ndarray of length 2
        Point in [-1, 1]^2.
    """
    lib = _ensure_bindings()
    inp = (c_double * 3)(*xyz)
    out = (c_double * 2)()
    ierr = c_int()
    lib.emsoft_lambert_sphere_to_square(inp, out, ctypes.byref(ierr))
    _check_err(ierr.value, 'sphere_to_square')
    return np.array(out, dtype=np.float64)


def cube_to_ball(cube):
    """Map a 3D cube point to a unit ball point (3D Lambert).

    Parameters
    ----------
    cube : array_like of length 3
        Point in the cube [-a, a]^3 where a = pi^(2/3)/2.

    Returns
    -------
    numpy.ndarray of length 3
        Point in the unit ball.
    """
    lib = _ensure_bindings()
    inp = (c_double * 3)(*cube)
    out = (c_double * 3)()
    ierr = c_int()
    lib.emsoft_lambert_cube_to_ball(inp, out, ctypes.byref(ierr))
    _check_err(ierr.value, 'cube_to_ball')
    return np.array(out, dtype=np.float64)


def ball_to_cube(ball):
    """Map a unit ball point to a 3D cube point (inverse 3D Lambert).

    Parameters
    ----------
    ball : array_like of length 3
        Point in the unit ball.

    Returns
    -------
    numpy.ndarray of length 3
        Point in the cube.
    """
    lib = _ensure_bindings()
    inp = (c_double * 3)(*ball)
    out = (c_double * 3)()
    ierr = c_int()
    lib.emsoft_lambert_ball_to_cube(inp, out, ctypes.byref(ierr))
    _check_err(ierr.value, 'ball_to_cube')
    return np.array(out, dtype=np.float64)


def stereo_forward(xyz):
    """Forward stereographic projection: hemisphere to plane.

    Parameters
    ----------
    xyz : array_like of length 3
        Point on the unit sphere.

    Returns
    -------
    numpy.ndarray of length 2
        Stereographic coordinates.
    """
    lib = _ensure_bindings()
    inp = (c_double * 3)(*xyz)
    out = (c_double * 2)()
    ierr = c_int()
    lib.emsoft_lambert_stereo_forward(inp, out, ctypes.byref(ierr))
    _check_err(ierr.value, 'stereo_forward')
    return np.array(out, dtype=np.float64)


def stereo_inverse(xy):
    """Inverse stereographic projection: plane to hemisphere.

    Parameters
    ----------
    xy : array_like of length 2
        Stereographic coordinates.

    Returns
    -------
    numpy.ndarray of length 3
        Point on the unit sphere.
    """
    lib = _ensure_bindings()
    inp = (c_double * 2)(*xy)
    out = (c_double * 3)()
    ierr = c_int()
    lib.emsoft_lambert_stereo_inverse(inp, out, ctypes.byref(ierr))
    _check_err(ierr.value, 'stereo_inverse')
    return np.array(out, dtype=np.float64)
