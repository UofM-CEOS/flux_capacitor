"""Unit test of flux module

"""

import unittest as ut
from fluxer.eddycov import flux
import numpy as np
from numpy import testing as npt


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
