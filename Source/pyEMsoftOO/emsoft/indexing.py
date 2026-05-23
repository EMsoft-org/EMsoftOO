"""
Dictionary indexing support for EBSD pattern matching.

Provides pattern preprocessing, normalized dot-product matching,
quality metrics (CI, ADP, KAM), and pattern center conversion
between vendor conventions.

All operations are implemented in numpy for simplicity and flexibility.
"""

import numpy as np


def normalize_patterns(patterns):
    """L2-normalize each pattern in a stack.

    Parameters
    ----------
    patterns : numpy.ndarray of shape (n, ...)
        Stack of n patterns (any shape per pattern).

    Returns
    -------
    numpy.ndarray of same shape
        Normalized patterns (each pattern has unit L2 norm).
    """
    patterns = np.asarray(patterns, dtype=np.float64)
    n = patterns.shape[0]
    flat = patterns.reshape(n, -1)
    norms = np.linalg.norm(flat, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-30)  # avoid division by zero
    flat = flat / norms
    return flat.reshape(patterns.shape)


def circular_mask(numsx, numsy):
    """Create a circular mask for EBSD patterns.

    Parameters
    ----------
    numsx, numsy : int
        Pattern dimensions.

    Returns
    -------
    numpy.ndarray of shape (numsy, numsx), dtype float64
        1.0 inside the circle, 0.0 outside.
    """
    cx, cy = numsx / 2.0, numsy / 2.0
    r = min(cx, cy)
    y, x = np.ogrid[:numsy, :numsx]
    mask = ((x - cx + 0.5)**2 + (y - cy + 0.5)**2) <= r**2
    return mask.astype(np.float64)


def match_patterns(experimental, dictionary, n_top=1):
    """Find best-matching dictionary entries for experimental patterns.

    Computes normalized dot products between all pairs of experimental
    and dictionary patterns and returns the top matches.

    Parameters
    ----------
    experimental : numpy.ndarray of shape (n_exp, npixels)
        Preprocessed and L2-normalized experimental patterns (flattened).
    dictionary : numpy.ndarray of shape (n_dict, npixels)
        Preprocessed and L2-normalized dictionary patterns (flattened).
    n_top : int
        Number of top matches to return per experimental pattern.

    Returns
    -------
    dict with keys:
        'indices' : numpy.ndarray of shape (n_exp, n_top)
            Indices into the dictionary of the top matches.
        'dot_products' : numpy.ndarray of shape (n_exp, n_top)
            Dot product values for the top matches (higher = better).
    """
    experimental = np.asarray(experimental, dtype=np.float64)
    dictionary = np.asarray(dictionary, dtype=np.float64)

    # Compute all dot products: (n_exp, n_dict)
    dp = experimental @ dictionary.T

    # Find top-N matches
    if n_top >= dp.shape[1]:
        indices = np.argsort(-dp, axis=1)
    else:
        indices = np.argpartition(-dp, n_top, axis=1)[:, :n_top]
        # Sort the top-N by dot product (descending)
        for i in range(len(indices)):
            order = np.argsort(-dp[i, indices[i]])
            indices[i] = indices[i][order]

    top_indices = indices[:, :n_top]
    top_dp = np.take_along_axis(dp, top_indices, axis=1)

    return {
        'indices': top_indices,
        'dot_products': top_dp,
    }


def confidence_index(dot_products):
    """Compute the confidence index from top dot products.

    CI = dp[0] - dp[1] (difference between best and second-best match).
    Higher values indicate more confident indexing.

    Parameters
    ----------
    dot_products : numpy.ndarray of shape (n, k) with k >= 2
        Top-k dot products for each pattern.

    Returns
    -------
    numpy.ndarray of shape (n,)
        Confidence index values.
    """
    if dot_products.shape[1] < 2:
        raise ValueError("Need at least 2 top matches to compute CI")
    return dot_products[:, 0] - dot_products[:, 1]


def adp_map(patterns, width, height):
    """Compute the Average Dot Product map.

    For each pattern, computes the average dot product with its
    spatial neighbors (up to 4: left, right, up, down).

    Parameters
    ----------
    patterns : numpy.ndarray of shape (n, npixels)
        L2-normalized patterns, ordered row by row.
    width, height : int
        Map dimensions (width * height must equal n).

    Returns
    -------
    numpy.ndarray of shape (height, width)
        ADP map values.
    """
    patterns = np.asarray(patterns, dtype=np.float64)
    n = patterns.shape[0]
    if n != width * height:
        raise ValueError(f"width*height={width*height} != n_patterns={n}")

    adp = np.zeros(n, dtype=np.float64)

    for idx in range(n):
        row = idx // width
        col = idx % width
        count = 0
        dp_sum = 0.0

        # Check 4-connected neighbors
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = row + dr, col + dc
            if 0 <= nr < height and 0 <= nc < width:
                nidx = nr * width + nc
                dp_sum += np.dot(patterns[idx], patterns[nidx])
                count += 1

        adp[idx] = dp_sum / max(count, 1)

    return adp.reshape(height, width)


def kam_map(eulers, width, height, symmetry_ops=None):
    """Compute the Kernel Average Misorientation map.

    For each point, computes the average misorientation angle with
    its spatial neighbors.

    Parameters
    ----------
    eulers : numpy.ndarray of shape (n, 3)
        Euler angles in radians, ordered row by row.
    width, height : int
        Map dimensions.
    symmetry_ops : numpy.ndarray of shape (nsym, 4), optional
        Symmetry quaternions for disorientation calculation.
        If None, raw misorientation angles are used.

    Returns
    -------
    numpy.ndarray of shape (height, width)
        KAM values in radians.
    """
    from .rotations import Rotation

    n = eulers.shape[0]
    if n != width * height:
        raise ValueError(f"width*height={width*height} != n_eulers={n}")

    # Convert all Euler angles to quaternions
    quats = np.empty((n, 4), dtype=np.float64)
    for i in range(n):
        r = Rotation.from_euler(*eulers[i])
        quats[i] = r.to_quaternion()

    kam = np.zeros(n, dtype=np.float64)

    for idx in range(n):
        row = idx // width
        col = idx % width
        count = 0
        angle_sum = 0.0

        q1 = quats[idx]

        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = row + dr, col + dc
            if 0 <= nr < height and 0 <= nc < width:
                nidx = nr * width + nc
                q2 = quats[nidx]
                # Misorientation: q_mis = q1* . q2
                q1c = np.array([q1[0], -q1[1], -q1[2], -q1[3]])
                # Quaternion product q1c * q2
                w = q1c[0]*q2[0] - q1c[1]*q2[1] - q1c[2]*q2[2] - q1c[3]*q2[3]
                angle = 2.0 * np.arccos(np.clip(abs(w), 0, 1))
                angle_sum += angle
                count += 1

        kam[idx] = angle_sum / max(count, 1)

    return kam.reshape(height, width)


def convert_pattern_center(pc, from_vendor, to_vendor, delta, numsx, numsy):
    """Convert pattern center between vendor conventions.

    Parameters
    ----------
    pc : array_like of length 3
        Pattern center [x, y, z] in the source convention.
    from_vendor : str
        Source convention: 'EMsoft', 'EDAX', 'Oxford', or 'Bruker'.
    to_vendor : str
        Target convention: 'EMsoft', 'EDAX', 'Oxford', or 'Bruker'.
    delta : float
        Detector pixel size in microns.
    numsx, numsy : int
        Detector dimensions in pixels.

    Returns
    -------
    numpy.ndarray of length 3
        Pattern center in the target convention.

    Notes
    -----
    EMsoft convention: [xpc, ypc, L] where xpc/ypc are in pixels from
    center, L is in microns.

    EDAX/TSL convention: [x*, y*, z*] as fractions of pattern width,
    origin at lower-left.

    Oxford convention: [x*, y*, z*] as fractions of pattern width/height,
    origin at lower-left.

    Bruker convention: [x*, y*, z*] as fractions of pattern width,
    origin at upper-left.
    """
    pc = np.asarray(pc, dtype=np.float64)

    # Convert to EMsoft first
    if from_vendor.upper() in ('EMSOFT',):
        xpc, ypc, L = pc[0], pc[1], pc[2]
    elif from_vendor.upper() in ('EDAX', 'TSL', 'EDAX/TSL'):
        xpc = numsx * (0.5 - pc[0])
        ypc = numsy * (0.5 - pc[1])
        L = numsx * delta * pc[2]
    elif from_vendor.upper() in ('OXFORD',):
        xpc = numsx * (0.5 - pc[0])
        ypc = numsy * (pc[1] - 0.5)
        L = numsx * delta * pc[2]
    elif from_vendor.upper() in ('BRUKER',):
        xpc = numsx * (0.5 - pc[0])
        ypc = numsy * (pc[1] - 0.5)
        L = numsx * delta * pc[2]
    else:
        raise ValueError(f"Unknown vendor: {from_vendor}")

    # Convert from EMsoft to target
    if to_vendor.upper() in ('EMSOFT',):
        return np.array([xpc, ypc, L])
    elif to_vendor.upper() in ('EDAX', 'TSL', 'EDAX/TSL'):
        return np.array([
            0.5 - xpc / numsx,
            0.5 - ypc / numsy,
            L / (numsx * delta)
        ])
    elif to_vendor.upper() in ('OXFORD',):
        return np.array([
            0.5 - xpc / numsx,
            0.5 + ypc / numsy,
            L / (numsx * delta)
        ])
    elif to_vendor.upper() in ('BRUKER',):
        return np.array([
            0.5 - xpc / numsx,
            0.5 + ypc / numsy,
            L / (numsx * delta)
        ])
    else:
        raise ValueError(f"Unknown vendor: {to_vendor}")
