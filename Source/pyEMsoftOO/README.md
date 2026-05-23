# EMsoft Python Interface

Python wrappers for the EMsoftOO Fortran library.

## Requirements

- Python >= 3.9
- NumPy >= 1.20
- The `libEMsoftOO_c` shared library (built from `Source/EMsoftOOLib/c_interface/`)

## Installation

1. Build the C-interop shared library:

```bash
cd EMsoftOOBuild/Release
cmake -DCMAKE_BUILD_TYPE=Release -DEMsoftOO_SDK=/path/to/SDK ../../EMsoftOO
make EMsoftOO_c
```

2. Set the library path:

```bash
export EMSOFTOO_LIB=/path/to/EMsoftOOBuild/Release/lib/libEMsoftOO_c.dylib
```

3. Install the Python package:

```bash
cd Source/pyEMsoftOO
pip install -e ".[test]"
```

## Quick Start

```python
from emsoft.quaternions import Quaternion, QuaternionArray
import numpy as np

# Create quaternions
q1 = Quaternion(1, 0, 0, 0)        # identity
q2 = Quaternion(0.5, 0.5, 0.5, 0.5)  # 120 deg around [111]

# Arithmetic
q3 = q1 * q2       # Hamilton product
qc = q2.conjugate() # conjugate

# Rotate a vector
v = q2.rotate([1.0, 0.0, 0.0])

# Array operations
data = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0]], dtype=np.float64)
qa = QuaternionArray(data)
rotated = qa.rotate([1.0, 0.0, 0.0])  # rotate by each quaternion
```

## Running Tests

```bash
pytest emsoft/tests/
```
