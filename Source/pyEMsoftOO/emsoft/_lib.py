"""
Shared library loader for the EMsoftOO C-interop layer.
"""

import ctypes
import os
import sys

_lib = None


def _find_library():
    """Locate and load the libEMsoftOO_c shared library."""
    # Library name varies by platform
    if sys.platform == 'darwin':
        libname = 'libEMsoftOO_c.dylib'
    elif sys.platform == 'win32':
        libname = 'EMsoftOO_c.dll'
    else:
        libname = 'libEMsoftOO_c.so'

    # Search order:
    # 1. EMSOFTOO_LIB environment variable
    # 2. Same directory as this file
    # 3. ../lib relative to this file
    # 4. Standard library paths (let ctypes search)
    search_paths = []

    env_path = os.environ.get('EMSOFTOO_LIB')
    if env_path:
        search_paths.append(env_path)

    this_dir = os.path.dirname(os.path.abspath(__file__))
    search_paths.append(os.path.join(this_dir, libname))
    search_paths.append(os.path.join(this_dir, '..', 'lib', libname))
    search_paths.append(os.path.join(this_dir, '..', '..', 'lib', libname))

    for path in search_paths:
        if os.path.isfile(path):
            return ctypes.CDLL(path)

    # Fall back to system search
    return ctypes.CDLL(libname)


def get_lib():
    """Return the loaded shared library, loading it on first call."""
    global _lib
    if _lib is None:
        _lib = _find_library()
    return _lib
