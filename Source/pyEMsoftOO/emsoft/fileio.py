"""
HDF5 file I/O for EMsoftOO data files.

Uses h5py to read .xtal crystal structure files and master pattern
HDF5 files without going through the Fortran HDF5 layer.

Requires: h5py (pip install h5py)
"""

import numpy as np

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

# Crystal system names (EMsoftOO convention: 1-7)
CRYSTAL_SYSTEM_NAMES = {
    1: 'Triclinic',
    2: 'Monoclinic',
    3: 'Orthorhombic',
    4: 'Tetragonal',
    5: 'Trigonal',
    6: 'Hexagonal',
    7: 'Cubic',
}

# Element symbols indexed by atomic number
ELEMENT_SYMBOLS = [
    '', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
]


def _require_h5py():
    if not HAS_H5PY:
        raise ImportError(
            "h5py is required for file I/O. Install with: pip install h5py")


def _read_string(dataset):
    """Read a string dataset, handling bytes/str differences."""
    val = dataset[()]
    if isinstance(val, bytes):
        return val.decode('utf-8', errors='replace')
    if isinstance(val, np.ndarray):
        v = val.flat[0]
        return v.decode('utf-8', errors='replace') if isinstance(v, bytes) else str(v)
    return str(val)


def read_xtal(filename):
    """Read an EMsoftOO .xtal crystal structure file.

    Parameters
    ----------
    filename : str
        Path to the .xtal file.

    Returns
    -------
    dict with keys:
        'lattice_parameters' : numpy.ndarray of shape (6,)
            [a, b, c, alpha, beta, gamma] in nm and degrees.
        'space_group_number' : int
            ITC space group number (1-230).
        'crystal_system' : int
            Crystal system number (1-7).
        'crystal_system_name' : str
            Crystal system name.
        'n_atom_types' : int
            Number of distinct atom types.
        'atom_types' : numpy.ndarray of ints
            Atomic numbers for each atom type.
        'atom_data' : numpy.ndarray of shape (n, 5)
            Each row: [x, y, z, occupancy, Debye-Waller].
        'source' : str
            Citation for the crystallographic data.
        'setting' : int
            Space group setting (1 or 2).
    """
    _require_h5py()

    with h5py.File(filename, 'r') as f:
        cd = f['CrystalData']

        data = {}
        data['lattice_parameters'] = cd['LatticeParameters'][()]
        data['space_group_number'] = int(cd['SpaceGroupNumber'][()])
        data['crystal_system'] = int(cd['CrystalSystem'][()])
        data['crystal_system_name'] = CRYSTAL_SYSTEM_NAMES.get(
            data['crystal_system'], 'Unknown')
        data['n_atom_types'] = int(cd['Natomtypes'][()])
        data['atom_types'] = cd['Atomtypes'][()]
        data['atom_data'] = cd['AtomData'][()]
        data['setting'] = int(cd['SpaceGroupSetting'][()])

        data['source'] = _read_string(cd['Source']) if 'Source' in cd else ''

    return data


def read_master_pattern(filename):
    """Read key data from an EMsoftOO EBSD master pattern HDF5 file.

    Parameters
    ----------
    filename : str
        Path to the master pattern .h5 file.

    Returns
    -------
    dict with keys:
        'mLPNH' : numpy.ndarray
            Northern hemisphere master pattern (modified Lambert projection).
        'mLPSH' : numpy.ndarray
            Southern hemisphere master pattern.
        'xtal_name' : str
            Crystal structure file used.
        'voltage' : float
            Accelerating voltage in keV.
        'npx' : int
            Half-width of the master pattern grid.
        'atom_types' : numpy.ndarray
            Atomic numbers.

    Notes
    -----
    The master pattern arrays are square grids of size (2*npx+1, 2*npx+1)
    in a modified Lambert projection. Use emsoft.lambert functions to
    map between grid coordinates and sphere directions.
    """
    _require_h5py()

    data = {}

    with h5py.File(filename, 'r') as f:
        # Navigate the EMsoft HDF5 structure
        # Typical path: EMData/EBSDmaster/mLPNH
        if 'EMData' in f:
            emd = f['EMData']
            # Find the master pattern group (EBSDmaster, ECPmaster, etc.)
            for key in ['EBSDmaster', 'ECPmaster', 'TKDmaster']:
                if key in emd:
                    mp = emd[key]
                    break
            else:
                raise KeyError("No master pattern group found in EMData")

            if 'mLPNH' in mp:
                data['mLPNH'] = mp['mLPNH'][()]
            if 'mLPSH' in mp:
                data['mLPSH'] = mp['mLPSH'][()]

        # Crystal data
        if 'CrystalData' in f:
            cd = f['CrystalData']
            if 'Atomtypes' in cd:
                data['atom_types'] = cd['Atomtypes'][()]

        # Namelist parameters
        if 'NMLparameters' in f:
            nml = f['NMLparameters']
            for key in ['EBSDMasterNameList', 'ECPMasterNameList', 'TKDMasterNameList']:
                if key in nml:
                    nl = nml[key]
                    if 'npx' in nl:
                        data['npx'] = int(nl['npx'][()])
                    if 'energyfile' in nl:
                        data['xtal_name'] = _read_string(nl['energyfile'])
                    break

        # Accelerating voltage (from MC data or namelist)
        if 'NMLparameters' in f:
            nml = f['NMLparameters']
            for key in ['MCCLNameList', 'MCOpenCLNameList']:
                if key in nml:
                    mcnl = nml[key]
                    if 'EkeV' in mcnl:
                        data['voltage'] = float(mcnl['EkeV'][()])
                    break

    return data
