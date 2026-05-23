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
    # 1. EMSOFTOO_LIB environment variable (file path or directory)
    # 2. Same directory as this file
    # 3. ../lib relative to this file
    # 4. Standard library paths (let ctypes search)
    search_paths = []

    env_path = os.environ.get('EMSOFTOO_LIB')
    if env_path:
        # Accept either a full file path or a directory
        if os.path.isfile(env_path):
            search_paths.append(env_path)
        elif os.path.isdir(env_path):
            search_paths.append(os.path.join(env_path, libname))

    this_dir = os.path.dirname(os.path.abspath(__file__))
    search_paths.append(os.path.join(this_dir, libname))
    search_paths.append(os.path.join(this_dir, '..', 'lib', libname))
    search_paths.append(os.path.join(this_dir, '..', '..', 'lib', libname))

    # Try each candidate path
    for path in search_paths:
        if os.path.isfile(path):
            try:
                return ctypes.CDLL(path)
            except OSError as e:
                raise OSError(
                    f"Found {path} but failed to load it:\n  {e}\n\n"
                    f"This usually means a dependent library (e.g. Fortran runtime) "
                    f"cannot be found.\nOn macOS, try:\n"
                    f"  export DYLD_LIBRARY_PATH=/path/to/EMsoftOO_SDK/lib:$DYLD_LIBRARY_PATH"
                ) from e

    # Fall back to system search
    try:
        return ctypes.CDLL(libname)
    except OSError:
        searched = '\n  '.join(search_paths) if search_paths else '(none)'
        raise OSError(
            f"Could not find {libname}. Set the EMSOFTOO_LIB environment variable\n"
            f"to the full path of the shared library or the directory containing it.\n"
            f"Example: setenv EMSOFTOO_LIB /path/to/EMsoftOOBuild/Release/lib/{libname}\n"
            f"Searched:\n  {searched}"
        )


def get_lib():
    """Return the loaded shared library, loading it on first call."""
    global _lib
    if _lib is None:
        _lib = _find_library()
    return _lib
