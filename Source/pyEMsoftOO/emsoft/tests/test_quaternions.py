"""
Tests for the emsoft.quaternions module.

Run with: pytest test_quaternions.py
"""

import numpy as np
import pytest
from emsoft.quaternions import Quaternion, QuaternionArray


class TestQuaternionCreation:
    def test_identity(self):
        q = Quaternion.identity()
        np.testing.assert_allclose(q.components, [1, 0, 0, 0])

    def test_from_components(self):
        q = Quaternion(0.5, 0.5, 0.5, 0.5)
        np.testing.assert_allclose(q.components, [0.5, 0.5, 0.5, 0.5])

    def test_from_array(self):
        q = Quaternion(components=[0.0, 1.0, 0.0, 0.0])
        assert q.w == pytest.approx(0.0)
        assert q.x == pytest.approx(1.0)

    def test_default_is_identity(self):
        q = Quaternion()
        np.testing.assert_allclose(q.components, [1, 0, 0, 0])


class TestQuaternionArithmetic:
    def test_add(self):
        q1 = Quaternion(1, 0, 0, 0)
        q2 = Quaternion(0, 1, 0, 0)
        q3 = q1 + q2
        np.testing.assert_allclose(q3.components, [1, 1, 0, 0])

    def test_subtract(self):
        q1 = Quaternion(1, 1, 0, 0)
        q2 = Quaternion(0, 1, 0, 0)
        q3 = q1 - q2
        np.testing.assert_allclose(q3.components, [1, 0, 0, 0])

    def test_multiply_quaternions(self):
        # i * j = k
        qi = Quaternion(0, 1, 0, 0)
        qj = Quaternion(0, 0, 1, 0)
        qk = qi * qj
        np.testing.assert_allclose(qk.components, [0, 0, 0, 1], atol=1e-15)

    def test_multiply_scalar(self):
        q = Quaternion(1, 0, 0, 0)
        q2 = q * 2.0
        np.testing.assert_allclose(q2.components, [2, 0, 0, 0])

    def test_rmul_scalar(self):
        q = Quaternion(1, 0, 0, 0)
        q2 = 3.0 * q
        np.testing.assert_allclose(q2.components, [3, 0, 0, 0])

    def test_divide_quaternions(self):
        q1 = Quaternion(1, 2, 3, 4)
        q2 = Quaternion(1, 2, 3, 4)
        q3 = q1 / q2
        # q / q should be identity
        np.testing.assert_allclose(q3.components, [1, 0, 0, 0], atol=1e-14)

    def test_negate(self):
        q = Quaternion(1, 2, 3, 4)
        qn = -q
        np.testing.assert_allclose(qn.components, [-1, -2, -3, -4])


class TestQuaternionOperations:
    def test_conjugate(self):
        q = Quaternion(1, 2, 3, 4)
        qc = q.conjugate()
        np.testing.assert_allclose(qc.components, [1, -2, -3, -4])

    def test_norm(self):
        q = Quaternion(1, 0, 0, 0)
        assert q.norm() == pytest.approx(1.0)

        q2 = Quaternion(1, 1, 1, 1)
        assert q2.norm() == pytest.approx(2.0)

    def test_normalize(self):
        q = Quaternion(2, 0, 0, 0)
        q.normalize()
        assert q.norm() == pytest.approx(1.0)
        np.testing.assert_allclose(q.components, [1, 0, 0, 0])

    def test_flip(self):
        q = Quaternion(1, 2, 3, 4)
        q.flip()
        np.testing.assert_allclose(q.components, [-1, -2, -3, -4])

    def test_positive(self):
        q = Quaternion(-1, 0, 0, 0)
        q.positive()
        assert q.w == pytest.approx(1.0)

    def test_inner_product(self):
        q1 = Quaternion(1, 0, 0, 0)
        q2 = Quaternion(1, 0, 0, 0)
        assert q1.inner(q2) == pytest.approx(1.0)

    def test_angle(self):
        q1 = Quaternion.identity()
        q2 = Quaternion.identity()
        assert q1.angle(q2) == pytest.approx(0.0)

    def test_equality(self):
        q1 = Quaternion(1, 0, 0, 0)
        q2 = Quaternion(1, 0, 0, 0)
        assert q1 == q2

    def test_inequality(self):
        q1 = Quaternion(1, 0, 0, 0)
        q2 = Quaternion(0, 1, 0, 0)
        assert not (q1 == q2)


class TestQuaternionRotation:
    def test_rotate_identity(self):
        q = Quaternion.identity()
        v = [1.0, 0.0, 0.0]
        vr = q.rotate(v)
        np.testing.assert_allclose(vr, [1, 0, 0], atol=1e-15)

    def test_rotate_90_around_z(self):
        # 90 degrees around z-axis: q = [cos(45), 0, 0, sin(45)]
        angle = np.pi / 2
        q = Quaternion(np.cos(angle / 2), 0, 0, np.sin(angle / 2))
        v = [1.0, 0.0, 0.0]
        vr = q.rotate(v)
        np.testing.assert_allclose(vr, [0, 1, 0], atol=1e-14)

    def test_rotate_vectors(self):
        angle = np.pi / 2
        q = Quaternion(np.cos(angle / 2), 0, 0, np.sin(angle / 2))
        vecs = np.array([[1, 0, 0], [0, 1, 0]], dtype=np.float64)
        vr = q.rotate_vectors(vecs)
        np.testing.assert_allclose(vr[0], [0, 1, 0], atol=1e-14)
        np.testing.assert_allclose(vr[1], [-1, 0, 0], atol=1e-14)


class TestQuaternionArray:
    def test_create(self):
        data = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64)
        qa = QuaternionArray(data)
        assert len(qa) == 2

    def test_get_element(self):
        data = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64)
        qa = QuaternionArray(data)
        q0 = qa[0]
        np.testing.assert_allclose(q0.components, [1, 0, 0, 0])
        q1 = qa[1]
        np.testing.assert_allclose(q1.components, [0, 1, 0, 0])

    def test_set_element(self):
        qa = QuaternionArray(n=2)
        qa[0] = Quaternion(1, 0, 0, 0)
        qa[1] = Quaternion(0, 0, 1, 0)
        q1 = qa[1]
        np.testing.assert_allclose(q1.components, [0, 0, 1, 0])

    def test_negative_index(self):
        data = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64)
        qa = QuaternionArray(data)
        q = qa[-1]
        np.testing.assert_allclose(q.components, [0, 1, 0, 0])

    def test_multiply(self):
        data1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64)
        data2 = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64)
        qa1 = QuaternionArray(data1)
        qa2 = QuaternionArray(data2)
        qa3 = qa1 * qa2
        assert len(qa3) == 2

    def test_normalize(self):
        data = np.array([[2, 0, 0, 0], [0, 0, 3, 0]], dtype=np.float64)
        qa = QuaternionArray(data)
        qa.normalize()
        q0 = qa[0]
        assert q0.norm() == pytest.approx(1.0)

    def test_rotate(self):
        # Two quaternions: identity and 90 degrees around z
        angle = np.pi / 2
        data = np.array([
            [1, 0, 0, 0],
            [np.cos(angle / 2), 0, 0, np.sin(angle / 2)]
        ], dtype=np.float64)
        qa = QuaternionArray(data)
        vr = qa.rotate([1.0, 0.0, 0.0])
        np.testing.assert_allclose(vr[0], [1, 0, 0], atol=1e-14)
        np.testing.assert_allclose(vr[1], [0, 1, 0], atol=1e-14)

    def test_to_array(self):
        data = np.array([[1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64)
        qa = QuaternionArray(data)
        arr = qa.to_array()
        np.testing.assert_allclose(arr, data)

    def test_repr(self):
        qa = QuaternionArray(n=5)
        assert repr(qa) == 'QuaternionArray(n=5)'
