"""Unit test of flux module

"""

import unittest as ut
from fluxer.eddycov import flux
import numpy as np
from numpy import testing as npt
import transforms3d.derivations.eulerangles as euler_deriv


class TestAngleVectorFuncs (ut.TestCase):
    def setUp(self):
        angles = np.linspace(0, 360, 10, dtype=int)
        vmags = np.linspace(0, 45, 10, dtype=int)
        # Expected returned angles grid
        E_angles, E_vmags = np.meshgrid(angles, vmags, indexing="ij")
        # Account for angle 0=360 and magnitude=0 means angle=0 conventions
        E_angles[0, :] = 360
        E_angles[:, 0] = 360
        # Expected returned coordinates grid
        E_x2d = np.array([[0, 5, 10, 15, 20, 25, 30, 35, 40, 45],
                          [0, 3.83022222, 7.66044443, 11.49066665,
                           15.32088886, 19.15111108, 22.98133329,
                           26.81155551, 30.64177772, 34.47199994],
                          [0, 0.86824089, 1.73648178, 2.60472267,
                           3.47296355, 4.34120444, 5.20944533, 6.07768622,
                           6.94592711, 7.814168],
                          [0, -2.5, -5, -7.5, -10, -12.5, -15, -17.5,
                           -20, -22.5],
                          [0, -4.6984631, -9.39692621, -14.09538931,
                           -18.79385242, -23.49231552, -28.19077862,
                           -32.88924173, -37.58770483, -42.28616794],
                          [0, -4.6984631, -9.39692621, -14.09538931,
                           -18.79385242, -23.49231552, -28.19077862,
                           -32.88924173, -37.58770483, -42.28616794],
                          [0, -2.5, -5, -7.5, -10, -12.5, -15, -17.5,
                           -20, -22.5],
                          [0, 0.86824089, 1.73648178, 2.60472267,
                           3.47296355, 4.34120444, 5.20944533, 6.07768622,
                           6.94592711, 7.814168],
                          [0, 3.83022222, 7.66044443, 11.49066665,
                           15.32088886, 19.15111108, 22.98133329,
                           26.81155551, 30.64177772, 34.47199994],
                          [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]])
        E_y2d = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 3.21393805, 6.42787610, 9.64181415,
                           1.28557522e+01, 1.60696902e+01, 1.92836283e+01,
                           2.24975663e+01, 2.57115044e+01, 2.89254424e+01],
                          [0, 4.92403877, 9.84807753, 1.47721163e+01,
                           1.96961551e+01, 2.46201938e+01, 2.95442326e+01,
                           3.44682714e+01, 3.93923101e+01, 4.43163489e+01],
                          [0, 4.33012702, 8.66025404, 1.29903811e+01,
                           1.73205081e+01, 2.16506351e+01, 2.59807621e+01,
                           3.03108891e+01, 3.46410162e+01, 3.89711432e+01],
                          [0, 1.7101007, 3.42020143, 5.13030215, 6.84040287,
                           8.55050358, 1.02606043e+01, 1.19707050e+01,
                           1.36808057e+01, 1.53909064e+01],
                          [0, -1.71010072, -3.42020143, -5.13030215,
                           -6.84040287, -8.55050358, -1.02606043e+01,
                           -1.19707050e+01, -1.36808057e+01, -1.53909064e+01],
                          [0, -4.33012702, -8.66025404, -1.29903811e+01,
                           -1.73205081e+01, -2.16506351e+01, -2.59807621e+01,
                           -3.03108891e+01, -3.46410162e+01, -3.89711432e+01],
                          [0, -4.92403877, -9.84807753, -1.47721163e+01,
                           -1.96961551e+01, -2.46201938e+01, -2.95442326e+01,
                           -3.44682714e+01, -3.93923101e+01, -4.43163489e+01],
                          [0, -3.21393805, -6.42787610, -9.64181415,
                           -1.28557522e+01, -1.60696902e+01, -1.92836283e+01,
                           -2.24975663e+01, -2.57115044e+01, -2.89254424e+01],
                          [0, -1.22464680e-15, -2.44929360e-15,
                           -3.67394040e-15, -4.89858720e-15, -6.12323400e-15,
                           -7.34788079e-15, -8.57252759e-15, -9.79717439e-15,
                           -1.10218212e-14]])
        self.angles = angles
        self.vmags = vmags
        self.E_angles = E_angles
        self.E_vmags = E_vmags
        self.E_x2d = E_x2d
        self.E_y2d = E_y2d

    def test_decompose_scalars(self):
        """Decomposing scalars should match expected"""
        msg = ("Calculated and expected {0} coordinates are not almost equal"
               " to 6 decimals")  # default
        for i, angle in enumerate(self.angles):
            for j, vmag in enumerate(self.vmags):
                xx, yy = flux.decompose(angle, vmag)
                E_x = self.E_x2d[i, j]
                E_y = self.E_y2d[i, j]
                npt.assert_almost_equal(E_x, xx, err_msg=msg.format("`x`"))
                npt.assert_almost_equal(E_y, yy, err_msg=msg.format("`y`"))

    def test_decompose_vectors(self):
        """Decomposing vectors should match expected"""
        xx, yy = flux.decompose(self.angles, self.vmags)
        E_x = np.diag(self.E_x2d)
        E_y = np.diag(self.E_y2d)
        msg = ("Calculated and expected {0} coordinates are not almost equal"
               " to 6 decimals")  # default
        npt.assert_array_almost_equal(E_x, xx, err_msg=msg.format("`x`"))
        npt.assert_array_almost_equal(E_y, yy, err_msg=msg.format("`y`"))

    def test_decompose_matrix(self):
        """Decomposing matrices should match expected"""
        xx, yy = flux.decompose(self.angles[:, np.newaxis],
                                self.vmags[np.newaxis, :])
        msg = ("Calculated and expected {0} coordinates are not almost equal"
               " to 6 decimals")  # default
        npt.assert_array_almost_equal(self.E_x2d, xx,
                                      err_msg=msg.format("`x`"))
        npt.assert_array_almost_equal(self.E_y2d, yy,
                                      err_msg=msg.format("`y`"))

    def test_recompose_scalars(self):
        """Recomposing scalars should match expected"""
        msg = ("Calculated and expected {0} are not almost equal"
               " to 6 decimals, at x[{1}], y[{2}]")  # default
        E_angles = self.E_angles
        E_vmags = self.E_vmags
        for i in range(len(self.angles)):
            for j in range(len(self.vmags)):
                angle_hat, vmag_hat = flux.recompose(self.E_x2d[i, j],
                                                     self.E_y2d[i, j])
                msg_angle = msg.format("angles", i, j)
                msg_vmag = msg.format("magnitudes", i, j)
                npt.assert_array_almost_equal(angle_hat, E_angles[i, j],
                                              err_msg=msg_angle)
                npt.assert_array_almost_equal(vmag_hat, E_vmags[i, j],
                                              err_msg=msg_vmag)

    def test_recompose_vectors(self):
        """Recomposing vectors should match expected"""
        E_angles = np.diag(self.E_angles)
        E_vmags = np.diag(self.E_vmags)
        angles_hat, vmags_hat = flux.recompose(np.diag(self.E_x2d),
                                               np.diag(self.E_y2d))
        msg = ("Calculated and expected {0} are not almost equal"
               " to 6 decimals")  # default
        npt.assert_array_almost_equal(angles_hat, E_angles,
                                      err_msg=msg.format("angles"))
        npt.assert_array_almost_equal(vmags_hat, E_vmags,
                                      err_msg=msg.format("magnitudes"))

    def test_recompose_matrix(self):
        """Recomposing matrices should match expected"""
        E_angles = self.E_angles
        E_vmags = self.E_vmags
        angles_hat, vmags_hat = flux.recompose(self.E_x2d, self.E_y2d)
        msg = ("Calculated and expected {0} are not almost equal"
               " to 6 decimals")  # default
        npt.assert_array_almost_equal(angles_hat, E_angles,
                                      err_msg=msg.format("angles"))
        npt.assert_array_almost_equal(vmags_hat, E_vmags,
                                      err_msg=msg.format("magnitudes"))


class TestRotations (ut.TestCase):
    def setUp(self):
        np.random.seed(123)
        self.uvw = np.random.randn(10, 3)
        self.euler_degs = np.array([[10, 20, 30]])
        self.degs = np.resize(self.euler_degs, self.uvw.shape)
        self.uvw_rot = flux.euler_rotate(self.uvw, self.degs)
        phi, theta, psi = euler_deriv.symbols("\\phi, \\theta, \\psi")
        self.phi, self.theta, self.psi = phi, theta, psi

        def euler_mat_fun(x):
            """Compute active rotation matrix from 3 angles in 1-D array x"""
            rotx = flux.rotation_matrix(x[0], 0)
            roty = flux.rotation_matrix(x[1], 1)
            rotz = flux.rotation_matrix(x[2], 2)
            return np.transpose(np.dot(rotz, np.dot(roty, rotx)))

        self.euler_mat_fun = euler_mat_fun

    def test_euler_rotations(self):
        """Compare rotations from fluxer and transforms3d"""

        # Derive the rotation matrix applied in euler_rotate

        # The transpose of a product of matrices is equal to the product of
        # their transposes in reverse order... so reverse the order of
        # multiplication, which obfuscates the meaning... The alternative
        # is to have the rotation matrix pre-multiply column vectors, which
        # is inconvenient.  Either way at least one transpose is required
        euler_mat0 = (euler_deriv.x_rotation(self.phi).T *
                      euler_deriv.y_rotation(self.theta).T *
                      euler_deriv.z_rotation(self.psi).T)
        subs = dict(zip([self.phi, self.theta, self.psi],
                        *np.radians(self.euler_degs)))
        M_zyx = np.array(euler_mat0.subs(subs)).astype(np.float)
        npt.assert_array_almost_equal(np.dot(self.uvw, M_zyx),
                                      self.uvw_rot)

        xrot_active = flux.rotation_matrix(self.euler_degs[0, 0], 0, True)
        yrot_active = flux.rotation_matrix(self.euler_degs[0, 1], 1, True)
        zrot_active = flux.rotation_matrix(self.euler_degs[0, 2], 2, True)
        euler_mat1 = np.dot(xrot_active, np.dot(yrot_active, zrot_active))
        npt.assert_array_almost_equal(np.dot(self.uvw, euler_mat1),
                                      self.uvw_rot)

        # Try with passive rotation and transpose the resulting rotation
        # matrix
        xrot_passive = flux.rotation_matrix(self.euler_degs[0, 0], 0)
        yrot_passive = flux.rotation_matrix(self.euler_degs[0, 1], 1)
        zrot_passive = flux.rotation_matrix(self.euler_degs[0, 2], 2)
        euler_mat2 = np.dot(zrot_passive,
                            np.dot(yrot_passive, xrot_passive))
        npt.assert_array_almost_equal(np.dot(self.uvw, euler_mat2.T),
                                      self.uvw_rot)

        # Now looping through rows, mimicking the real case of changing
        # euler angles.  First show derivation:
        uvw_rots = np.empty_like(self.uvw)
        for i, v in enumerate(self.degs):
            euler_mati = self.euler_mat_fun(v)
            uvw_roti = np.dot(self.uvw[i], euler_mati)
            uvw_rots[i] = uvw_roti
        npt.assert_array_almost_equal(uvw_rots, self.uvw_rot)
