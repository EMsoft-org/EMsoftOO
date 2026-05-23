"""Tests for the emsoft.lambert module."""

import numpy as np
import pytest
from emsoft import lambert


class TestSquareSphere:
    def test_origin(self):
        xyz = lambert.square_to_sphere([0.0, 0.0])
        # Origin maps to north pole
        np.testing.assert_allclose(xyz, [0, 0, 1], atol=1e-14)

    def test_roundtrip(self):
        xy_in = [0.3, 0.4]
        xyz = lambert.square_to_sphere(xy_in)
        xy_out = lambert.sphere_to_square(xyz)
        np.testing.assert_allclose(xy_out, xy_in, atol=1e-12)

    def test_unit_sphere(self):
        # Output should be on unit sphere
        xy = [0.5, 0.5]
        xyz = lambert.square_to_sphere(xy)
        assert np.linalg.norm(xyz) == pytest.approx(1.0, abs=1e-14)


class TestCubeBall:
    def test_origin(self):
        ball = lambert.cube_to_ball([0.0, 0.0, 0.0])
        np.testing.assert_allclose(ball, [0, 0, 0], atol=1e-14)

    def test_roundtrip(self):
        cube_in = [0.2, 0.3, 0.1]
        ball = lambert.cube_to_ball(cube_in)
        cube_out = lambert.ball_to_cube(ball)
        np.testing.assert_allclose(cube_out, cube_in, atol=1e-12)

    def test_inside_ball(self):
        # Output should be inside the unit ball
        cube = [0.5, 0.5, 0.5]
        ball = lambert.cube_to_ball(cube)
        assert np.linalg.norm(ball) <= 1.0 + 1e-14


class TestStereographic:
    def test_north_pole(self):
        xy = lambert.stereo_forward([0.0, 0.0, 1.0])
        np.testing.assert_allclose(xy, [0, 0], atol=1e-14)

    def test_roundtrip(self):
        # Start from a point on the sphere (northern hemisphere)
        xyz_in = np.array([0.3, 0.4, np.sqrt(1 - 0.3**2 - 0.4**2)])
        xy = lambert.stereo_forward(xyz_in)
        xyz_out = lambert.stereo_inverse(xy)
        np.testing.assert_allclose(xyz_out, xyz_in, atol=1e-12)

    def test_equator(self):
        # Point on equator: z=0
        xyz = [1.0, 0.0, 0.0]
        xy = lambert.stereo_forward(xyz)
        # Stereographic of equator point (1,0,0) -> (1,0)
        np.testing.assert_allclose(xy, [1, 0], atol=1e-14)
