"""
EBSD pattern simulation wrapping EMsoftOO.

Provides the EBSDSimulator class for computing simulated EBSD patterns
from pre-computed master patterns. Given a detector geometry and crystal
orientations, the simulator interpolates the master pattern on a modified
Lambert projection to produce realistic EBSD patterns.

Typical workflow:
    1. Load a master pattern from an HDF5 file (via emsoft.fileio)
    2. Create an EBSDSimulator with the master pattern data
    3. Set detector geometry
    4. Compute patterns for one or more orientations
"""

import ctypes
import numpy as np
from ._lib import get_lib

c_double = ctypes.c_double
c_int = ctypes.c_int
c_void_p = ctypes.c_void_p


def _setup_bindings():
    lib = get_lib()

    lib.emsoft_ebsd_compute_detector.argtypes = [
        c_int, c_int,          # numsx, numsy
        c_double, c_double,    # xpc, ypc
        c_double, c_double,    # delta, thetac
        c_double, c_double,    # omega, L
        c_void_p, c_void_p, c_void_p  # rgx, rgy, rgz
    ]
    lib.emsoft_ebsd_compute_detector.restype = None

    lib.emsoft_ebsd_compute_pattern.argtypes = [
        c_int, c_int, c_int,   # numsx, numsy, npx
        c_void_p, c_void_p, c_void_p,  # rgx, rgy, rgz
        c_void_p, c_void_p,   # mLPNH, mLPSH
        c_void_p,              # quat(4)
        c_void_p               # pattern out
    ]
    lib.emsoft_ebsd_compute_pattern.restype = None

    lib.emsoft_ebsd_compute_patterns.argtypes = [
        c_int, c_int, c_int,   # numsx, numsy, npx
        c_void_p, c_void_p, c_void_p,  # rgx, rgy, rgz
        c_void_p, c_void_p,   # mLPNH, mLPSH
        c_void_p,              # quats(4, nquats)
        c_int,                 # nquats
        c_void_p               # patterns out
    ]
    lib.emsoft_ebsd_compute_patterns.restype = None

    return lib


_bindings_ready = False
_lib_ref = None


def _ensure_bindings():
    global _bindings_ready, _lib_ref
    if not _bindings_ready:
        _lib_ref = _setup_bindings()
        _bindings_ready = True
    return _lib_ref


class EBSDDetector:
    """EBSD detector geometry.

    Pre-computes direction cosines for each detector pixel, which are
    reused across all pattern computations.

    Parameters
    ----------
    numsx, numsy : int
        Detector dimensions in pixels.
    xpc, ypc : float
        Pattern center coordinates (in pixels, relative to detector center).
    delta : float
        Detector pixel size in microns.
    thetac : float
        Detector tilt angle in degrees.
    L : float
        Sample-to-detector distance in mm.
    omega : float, optional
        Sample tilt angle in degrees (default 0).
    """

    def __init__(self, numsx, numsy, xpc=0.0, ypc=0.0, delta=50.0,
                 thetac=10.0, L=15000.0, omega=0.0):
        self._lib = _ensure_bindings()
        self.numsx = numsx
        self.numsy = numsy

        # Compute direction cosines (Fortran column-major)
        self.rgx = np.empty((numsx, numsy), dtype=np.float64, order='F')
        self.rgy = np.empty((numsx, numsy), dtype=np.float64, order='F')
        self.rgz = np.empty((numsx, numsy), dtype=np.float64, order='F')

        self._lib.emsoft_ebsd_compute_detector(
            c_int(numsx), c_int(numsy),
            c_double(xpc), c_double(ypc),
            c_double(delta), c_double(thetac),
            c_double(omega), c_double(L),
            self.rgx.ctypes.data_as(c_void_p),
            self.rgy.ctypes.data_as(c_void_p),
            self.rgz.ctypes.data_as(c_void_p))

    def __repr__(self):
        return f'EBSDDetector({self.numsx}x{self.numsy})'


class EBSDSimulator:
    """Simulate EBSD patterns from a master pattern.

    Parameters
    ----------
    mLPNH : numpy.ndarray
        Northern hemisphere master pattern on modified Lambert grid.
        Shape: (2*npx+1, 2*npx+1) or higher-dimensional (energy/phase
        dimensions will be summed).
    mLPSH : numpy.ndarray
        Southern hemisphere master pattern, same shape as mLPNH.
    npx : int, optional
        Half-width of the master pattern grid. If not provided, inferred
        from the array shape as (shape[0] - 1) // 2.

    Examples
    --------
    >>> from emsoft.fileio import read_master_pattern
    >>> from emsoft.ebsd import EBSDSimulator, EBSDDetector
    >>>
    >>> mp = read_master_pattern('Ni_master.h5')
    >>> sim = EBSDSimulator(mp['mLPNH'], mp['mLPSH'])
    >>> det = EBSDDetector(640, 480, delta=50.0, thetac=10.0, L=15000.0)
    >>>
    >>> # Single pattern
    >>> pattern = sim.compute_pattern(det, quaternion=[1, 0, 0, 0])
    >>>
    >>> # Batch computation
    >>> quats = np.array([[1,0,0,0], [0.5,0.5,0.5,0.5]], dtype=np.float64)
    >>> patterns = sim.compute_patterns(det, quaternions=quats)
    """

    def __init__(self, mLPNH, mLPSH, npx=None):
        self._lib = _ensure_bindings()

        # Handle multi-dimensional master patterns (sum over extra dims)
        mLPNH = np.asarray(mLPNH, dtype=np.float64)
        mLPSH = np.asarray(mLPSH, dtype=np.float64)

        # If 3D or 4D (energy bins, phases), sum to get 2D
        while mLPNH.ndim > 2:
            mLPNH = mLPNH.sum(axis=-1)
            mLPSH = mLPSH.sum(axis=-1)

        if npx is None:
            npx = (mLPNH.shape[0] - 1) // 2

        self.npx = npx
        expected = 2 * npx + 1
        if mLPNH.shape != (expected, expected):
            raise ValueError(
                f'Master pattern shape {mLPNH.shape} does not match npx={npx} '
                f'(expected ({expected}, {expected}))')

        # Store as Fortran-contiguous for the C-interop layer
        self._mLPNH = np.asfortranarray(mLPNH)
        self._mLPSH = np.asfortranarray(mLPSH)

    def compute_pattern(self, detector, quaternion=None, euler=None, degrees=False):
        """Compute a single EBSD pattern.

        Parameters
        ----------
        detector : EBSDDetector
            Detector geometry.
        quaternion : array_like of length 4, optional
            Orientation as [w, x, y, z].
        euler : array_like of length 3, optional
            Orientation as Euler angles [phi1, Phi, phi2].
        degrees : bool
            If True, Euler angles are in degrees.

        Returns
        -------
        numpy.ndarray of shape (numsx, numsy)
            Simulated EBSD pattern.
        """
        quat = self._parse_orientation(quaternion, euler, degrees)

        pattern = np.empty((detector.numsx, detector.numsy),
                           dtype=np.float64, order='F')
        q = np.asfortranarray(quat)

        self._lib.emsoft_ebsd_compute_pattern(
            c_int(detector.numsx), c_int(detector.numsy), c_int(self.npx),
            detector.rgx.ctypes.data_as(c_void_p),
            detector.rgy.ctypes.data_as(c_void_p),
            detector.rgz.ctypes.data_as(c_void_p),
            self._mLPNH.ctypes.data_as(c_void_p),
            self._mLPSH.ctypes.data_as(c_void_p),
            q.ctypes.data_as(c_void_p),
            pattern.ctypes.data_as(c_void_p))

        return pattern.T  # Return as (numsy, numsx) for image convention

    def compute_patterns(self, detector, quaternions):
        """Compute EBSD patterns for multiple orientations.

        Parameters
        ----------
        detector : EBSDDetector
            Detector geometry.
        quaternions : numpy.ndarray of shape (n, 4)
            Array of unit quaternions [w, x, y, z].

        Returns
        -------
        numpy.ndarray of shape (n, numsy, numsx)
            Stack of simulated EBSD patterns.
        """
        quats = np.asarray(quaternions, dtype=np.float64)
        nq = quats.shape[0]

        # Fortran expects (4, nquats) column-major
        quats_f = np.asfortranarray(quats.T)
        patterns = np.empty((detector.numsx, detector.numsy, nq),
                            dtype=np.float64, order='F')

        self._lib.emsoft_ebsd_compute_patterns(
            c_int(detector.numsx), c_int(detector.numsy), c_int(self.npx),
            detector.rgx.ctypes.data_as(c_void_p),
            detector.rgy.ctypes.data_as(c_void_p),
            detector.rgz.ctypes.data_as(c_void_p),
            self._mLPNH.ctypes.data_as(c_void_p),
            self._mLPSH.ctypes.data_as(c_void_p),
            quats_f.ctypes.data_as(c_void_p),
            c_int(nq),
            patterns.ctypes.data_as(c_void_p))

        # Return as (nquats, numsy, numsx) for image convention
        return np.transpose(patterns, (2, 1, 0))

    def _parse_orientation(self, quaternion, euler, degrees):
        """Convert orientation input to a quaternion array."""
        if quaternion is not None:
            return np.asarray(quaternion, dtype=np.float64)
        elif euler is not None:
            from .rotations import Rotation
            eu = np.asarray(euler, dtype=np.float64)
            if degrees:
                eu = np.radians(eu)
            r = Rotation.from_euler(*eu)
            return r.to_quaternion()
        else:
            raise ValueError("Provide either 'quaternion' or 'euler'")

    def __repr__(self):
        return f'EBSDSimulator(npx={self.npx})'
