"""Tests for the emsoft.rotations module."""

import numpy as np
import pytest
from emsoft.rotations import Rotation


class TestRotationCreation:
    def test_identity(self):
        r = Rotation.identity()
        q = r.to_quaternion()
        np.testing.assert_allclose(q, [1, 0, 0, 0], atol=1e-14)

    def test_from_euler_radians(self):
        r = Rotation.from_euler(0, 0, 0)
        q = r.to_quaternion()
        np.testing.assert_allclose(q, [1, 0, 0, 0], atol=1e-14)

    def test_from_euler_degrees(self):
        r = Rotation.from_euler(90, 0, 0, degrees=True)
        eu = r.to_euler(degrees=True)
        np.testing.assert_allclose(eu, [90, 0, 0], atol=1e-10)

    def test_from_quaternion(self):
        r = Rotation.from_quaternion(0.5, 0.5, 0.5, 0.5)
        q = r.to_quaternion()
        np.testing.assert_allclose(q, [0.5, 0.5, 0.5, 0.5], atol=1e-14)

    def test_from_matrix(self):
        # Identity matrix
        r = Rotation.from_matrix(np.eye(3))
        q = r.to_quaternion()
        np.testing.assert_allclose(abs(q[0]), 1.0, atol=1e-14)

    def test_from_axisangle(self):
        r = Rotation.from_axisangle([0, 0, 1], 90, degrees=True)
        ax = r.to_axisangle()
        np.testing.assert_allclose(ax[:3], [0, 0, 1], atol=1e-14)
        np.testing.assert_allclose(ax[3], np.pi / 2, atol=1e-14)

    def test_from_rodrigues(self):
        # Identity: [0, 0, 1, 0]
        r = Rotation.from_rodrigues([0, 0, 1, 0])
        q = r.to_quaternion()
        np.testing.assert_allclose(abs(q[0]), 1.0, atol=1e-14)


class TestRotationConversions:
    def test_euler_roundtrip(self):
        eu_in = [0.5, 1.0, 1.5]
        r = Rotation.from_euler(*eu_in)
        eu_out = r.to_euler()
        np.testing.assert_allclose(eu_out, eu_in, atol=1e-12)

    def test_quaternion_roundtrip(self):
        q_in = [0.5, 0.5, 0.5, 0.5]
        r = Rotation.from_quaternion(*q_in)
        q_out = r.to_quaternion()
        np.testing.assert_allclose(q_out, q_in, atol=1e-14)

    def test_matrix_roundtrip(self):
        # 90 degrees around z
        angle = np.pi / 2
        om_in = np.array([
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle),  np.cos(angle), 0],
            [0,              0,             1]
        ])
        r = Rotation.from_matrix(om_in)
        om_out = r.to_matrix()
        np.testing.assert_allclose(om_out, om_in, atol=1e-14)

    def test_euler_to_matrix_to_euler(self):
        r = Rotation.from_euler(30, 45, 60, degrees=True)
        om = r.to_matrix()
        r2 = Rotation.from_matrix(om)
        eu = r2.to_euler(degrees=True)
        np.testing.assert_allclose(eu, [30, 45, 60], atol=1e-10)

    def test_all_representations_consistent(self):
        r = Rotation.from_euler(10, 20, 30, degrees=True)
        # All representations should convert back to the same Euler angles
        for method_name in ['to_quaternion', 'to_matrix', 'to_axisangle',
                            'to_rodrigues', 'to_homochoric', 'to_cubochoric',
                            'to_stereographic', 'to_rotvec']:
            rep = getattr(r, method_name)()
            from_name = method_name.replace('to_', 'from_')
            if from_name == 'from_matrix':
                r2 = Rotation.from_matrix(rep)
            elif from_name == 'from_quaternion':
                r2 = Rotation.from_quaternion(*rep)
            elif from_name == 'from_axisangle':
                r2 = Rotation.from_axisangle(rep[:3], rep[3])
            elif from_name == 'from_rodrigues':
                r2 = Rotation.from_rodrigues(rep)
            else:
                r2 = getattr(Rotation, from_name)(rep)
            eu2 = r2.to_euler(degrees=True)
            np.testing.assert_allclose(eu2, [10, 20, 30], atol=1e-8,
                                       err_msg=f'Failed roundtrip via {method_name}')


class TestRotationRepr:
    def test_repr(self):
        r = Rotation.from_euler(30, 45, 60, degrees=True)
        s = repr(r)
        assert 'Rotation' in s
        assert '30.00' in s
