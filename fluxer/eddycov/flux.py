# pylint: disable=too-many-locals,invalid-name,no-member

"""Core functionality for the package

Library of main functions required for flux calculations.

"""

from collections import namedtuple
from itertools import groupby
import numpy as np
from scipy import interpolate as itpl
from scipy import signal
from scipy.stats import zscore
from astropy.convolution import convolve, Box1DKernel


__all__ = ["smooth_angle", "planarfit", "rotate_vectors",
           "wind3D_correct", "despike_VickersMahrt"]

# Valid 3-D wind rotation methods
_VECTOR_ROTATION_METHODS = {"DR", "TR", "PF"}


AngleCoordinates = namedtuple("AngleCoordinates", ["x", "y"])
Vector = namedtuple("Vector", ["angle", "magnitude"])
PlanarFitCoefs = namedtuple("PlanarFitCoefs", ["k_vct", "tilt_coefs"])
RotatedVectors = namedtuple("RotatedVectors", ["rotated", "phi_roll"])
CorrectedWind3D = namedtuple("CorrectedWind3D",
                             ["uvw_ship",
                              "euler_angles",
                              "euler_angles_angular_rates",
                              "euler_angles_accelerations",
                              "euler_angles_slow",
                              "euler_angles_fast",
                              "mount_offset_rotations",
                              "uvw_earth",
                              "uvw_angular",
                              "uvw_linear",
                              "ship_enu",
                              "uvw_enu",
                              "imu_enu"])
VickersMahrt = namedtuple("VickersMahrt",
                          ["x", "nspikes", "ntrends", "kclass"])


def decompose(angle, vmagnitude):
    """Decompose angle and magnitude into `x` and `y` coordinates

    Parameters
    ----------
    angle : array_like
        The angle(s) in degree units.
    vmagnitude : array_like
        The magnitude array associated with each angle.  It can be a scalar
        that is common to all `angle` elements.

    Returns
    -------
    namedtuple with ndarrays `x` and `y`, in that order.

    Examples
    --------
    >>> angles = np.arange(0, 360, 36)
    >>> vmags = np.arange(10)
    >>> decompose(angles, vmags)  # doctest: +NORMALIZE_WHITESPACE
    (array([ 0. , 0.80901699, 0.61803399, -0.92705098, -3.23606798,
            -5. , -4.85410197, -2.16311896, 2.47213595, 7.28115295]),
     array([  0.00000000e+00, 5.87785252e-01, 1.90211303e+00,
              2.85316955e+00, 2.35114101e+00, 6.12323400e-16,
             -3.52671151e+00, -6.65739561e+00, -7.60845213e+00,
             -5.29006727e+00]))

    """
    x = vmagnitude * np.cos(np.radians(angle))
    y = vmagnitude * np.sin(np.radians(angle))
    return AngleCoordinates(x, y)


def recompose(x, y):
    """Recompose angles and associated magnitudes from `x` and `y` vectors

    Parameters
    ----------
    x : array_like or scalar
        `x`-coordinates
    y : array_like or scalar
        `y`-coordinates

    Returns
    -------
    namedtuple with ndarrays `angle` and `magnitude`, in that order.

    """
    vmag = np.sqrt((x ** 2) + (y ** 2))
    ang = np.arctan2(y, x)
    try:                        # first assume array inputs
        ang[ang < 0] = ang[ang < 0] + (2 * np.pi)  # output range 0 - 2*pi
        ang[vmag == 0] = 0        # when magnitude is 0 the angle is also 0
        ang[ang == 0] = 2 * np.pi   # convention
    except (IndexError, TypeError):  # then check if scalar inputs
        if np.isscalar(ang) and ang < 0:
            ang = ang + (2 * np.pi)
        if np.isscalar(ang) and vmag == 0:
            ang = 0
        if np.isscalar(ang) and ang == 0:
            ang = 2 * np.pi

    return Vector(np.degrees(ang), vmag)


def smooth_angle(angle, vmagnitude=1, kernel_width=21):
    """Smooth angles by decomposing them, applying a boxcar average

    Smoothing is done by using a 1D kernel smoothing filter of a given
    width.

    Parameters
    ----------
    angle : numpy.ndarray
        The angle(s) in degree units.
    vmagnitude : numpy.ndarray, optional
        The magnitude array associated with each angle.  It can be a scalar
        that is common to all `angle` elements.
    kernel_width : int
        The width of the filter kernel.

    Returns
    -------
    namedtuple with ndarrays `angle` and `magnitude`, in that order.
    """
    x, y = decompose(angle, vmagnitude)
    x_smooth = convolve(x, Box1DKernel(kernel_width), boundary="extend")
    y_smooth = convolve(y, Box1DKernel(kernel_width), boundary="extend")
    return recompose(x_smooth, y_smooth)


def level3D_motion(accel, ang_rate, roll_range, pitch_range):
    """Level 3D acceleration and angular rate measured by motion sensor"""
    pass                        # IMPLEMENT THIS?


def level3D_anemometer(wind_speed, roll, pitch):
    """Level 3D anemometer measurements, given mean roll and pitch"""
    pass                        # IMPLEMENT THIS?


def planarfit(vectors):
    """Calculate planar fit coefficients for coordinate rotation

    See Handbook of Micrometeorology (Lee et al. 2004).  Ported from
    getPlanarCoeffs.m from Patric Sturm <pasturm@ethz.ch>.

    Parameters
    ----------
    vectors : numpy.ndarray
        A 2-D (Nx3) array with x, y, and z vectors, expressed in a
        right-handed coordinate system.  These vectors may correspond to u,
        v, and w wind speed vectors, or inertial acceleration components.

    Returns
    -------
    namedtuple with (index and name in brackets):
    numpy.ndarray [0, 'k_vct']
        1-D array (1x3) unit vector parallel to the new z-axis.
    numpy.ndarray [1, 'tilt_coefs']
        1-D array (1x3) Tilt coefficients.

    """
    vct_u = vectors[:, 0]
    vct_v = vectors[:, 1]
    vct_w = vectors[:, 2]
    vct_nrows = vectors.shape[0]
    sum_u = sum(vct_u)
    sum_v = sum(vct_v)
    sum_w = sum(vct_w)
    sum_uv = np.dot(vct_u, vct_v)
    sum_uw = np.dot(vct_u, vct_w)
    sum_vw = np.dot(vct_v, vct_w)
    sum_u2 = np.dot(vct_u, vct_u)
    sum_v2 = np.dot(vct_v, vct_v)
    H_arr = np.array([[vct_nrows, sum_u, sum_v],
                      [sum_u, sum_u2, sum_uv],
                      [sum_v, sum_uv, sum_v2]])
    g_arr = np.array([sum_w, sum_uw, sum_vw])
    tilt_coef = np.linalg.solve(H_arr, g_arr)
    # Determine unit vector parellel to new z-axis
    k_2 = 1 / np.sqrt(1 + tilt_coef[1] ** 2 + tilt_coef[2] ** 2)
    k_0 = -tilt_coef[1] * k_2
    k_1 = -tilt_coef[2] * k_2
    k_vct = np.array([k_0, k_1, k_2])
    return PlanarFitCoefs(k_vct, tilt_coef)


def rotate_vectors(vectors, method="PF", **kwargs):
    """Transform vectors to reference mean streamline coordinate system

    Use double rotation, triple rotation, or planar fit methods (Wilczak et
    al. 2001; Handbook of Micrometeorology).

    This is a general coordinate rotation tool, so can handle inputs such
    as wind speed and acceleration from inertial measurement units.

    Parameters
    ----------
    vectors : numpy.ndarray
        A 2-D (Nx3) array with x, y, and z vector components, expressed in
        a right-handed coordinate system.  These may represent u, v, and w
        wind speed vectors, or inertial acceleration.
    method : str
        One of: "DR", "TR", "PF" for double rotation, triple rotation, or
        planar fit.
    k_vector : numpy.ndarray (optional)
        1-D array (1x3) unit vector parallel to the new z-axis, when
        "method" is "PF" (planar fit).  If not supplied, then it is
        calculated.

    Returns
    -------
    namedtuple with (index, name in brackets):
    numpy.ndarray [0, 'rotated']
        2-D array (Nx3) Array with rotated vectors
    numpy.ndarray [1, 'phi_theta']
        1-D array (1x2) Phi (roll) and Theta (pitch) rotation angles

    """
    if method not in _VECTOR_ROTATION_METHODS:
        msg = "method must be one of "
        raise ValueError(msg + ', '.join("\"{}\"".format(m) for m in
                                         _VECTOR_ROTATION_METHODS))

    if method == "PF":
        if "k_vector" in kwargs:
            k_vct = kwargs.get("k_vector")
        else:
            k_vct, tilt_coef = planarfit(vectors)
        j_vct = np.cross(k_vct, np.mean(vectors, 0))
        j_vct = j_vct / np.sqrt(np.sum(j_vct ** 2))
        i_vct = np.cross(j_vct, k_vct)
        vcts_mat = np.column_stack((i_vct, j_vct, k_vct))
        vcts_new = np.dot(vectors, vcts_mat)
        phi = np.arccos(np.dot(k_vct, np.array([0, 0, 1])))
        theta = np.arctan2(np.mean(-vectors[:, 1], 0),
                           np.mean(vectors[:, 0], 0))
    else:
        # First rotation to set mean v to 0
        theta = np.arctan2(np.mean(vectors[:, 1]),
                           np.mean(vectors[:, 0]))
        rot1 = np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0],
                         [0, 0, 1]])
        vcts1 = np.dot(vectors, rot1)
        # Second rotation to set mean w to 0
        phi = np.arctan2(np.mean(vcts1[:, 2]),
                         np.mean(vcts1[:, 0]))
        rot2 = np.array([[np.cos(phi), 0, -np.sin(phi)],
                         [0, 1, 0],
                         [np.sin(phi), 0, np.cos(phi)]])
        vcts_new = np.dot(vcts1, rot2)
        # Third rotation to set mean vw to 0
        if method == "TR":
            mean_vw = np.mean(vcts_new[:, 1] * vcts_new[:, 2])
            psi = 0.5 * np.arctan2(2 * mean_vw,
                                   np.mean(vcts_new[:, 1] ** 2) -
                                   np.mean(vcts_new[:, 2] ** 2))
            rot3 = np.array([[1, 0, 0],
                             [0, np.cos(psi), -np.sin(psi)],
                             [0, np.sin(psi), np.cos(psi)]])
            vcts_new = np.dot(vcts_new, rot3)

    return RotatedVectors(vcts_new, np.array([phi, theta]))


def euler_rotate(X, euler):
    """Rotate vector matrix given Euler transformation matrix"""
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    phi, theta, psi = euler[:, 0], euler[:, 1], euler[:, 2]
    x_new = (x * np.cos(theta) * np.cos(psi) +
             y * (np.sin(phi) * np.sin(theta) * np.cos(psi) -
                  np.cos(phi) * np.sin(psi)) +
             z * (np.cos(phi) * np.sin(theta) * np.cos(psi) +
                  np.sin(phi) * np.sin(psi)))
    y_new = (x * np.cos(theta) * np.sin(psi) +
             y * (np.sin(phi) * np.sin(theta) * np.sin(psi) +
                  np.cos(phi) * np.cos(psi)) +
             z * (np.cos(phi) * np.sin(theta) * np.sin(psi) -
                  np.sin(phi) * np.cos(psi)))
    z_new = (x * (-np.sin(theta)) + y * (np.cos(theta) * np.sin(phi)) +
             z * (np.cos(theta) * np.cos(phi)))
    return np.column_stack((x_new, y_new, z_new))


def wind3D_correct(wind_speed, acceleration, angle_rate, heading, speed,
                   anemometer_pos, sample_freq, Tcf, Ta,
                   tilt_motion=np.array([0.0, 0.0]),
                   tilt_anemometer=np.array([0.0, 0.0])):
    """Correct wind vector measurements from a moving platform

    This is a port of Scott Miller's `motion` Matlab function, which
    implements Miller et al. (2008) approach.  Coordinate frame is assumed
    to be right-handed:

    x - positive towards the bow.
    y - positive towards port side.
    z - positive upwards.

    Parameters
    ----------
    wind_speed : numpy.ndarray
        A 2-D (Nx3) array with u, v, and w wind speed (m/s) vectors.
    acceleration : numpy.ndarray
        A 2-D (Nx3) array with x, y, and z acceleration (m/s/s) vectors.
    angle_rate : numpy.ndarray
        A 2-D (Nx3) array with angular rates (radians/s) vectors.
    heading : array_like
        A vector with ship heading measurements (deg) in right-hand
        coordinate system.
    speed : array_like
        A vector with ship speed (m/s).
    anemometer_pos : array_like
        [x, y, z] position vector of anemometer, relative to motion sensor
        (m).
    sample_freq : float
        Sampling frequency (Hz).
    Tcf : float
        Complimentary filter period (s).
    Ta : float
        High-pass filter cutoff for accelerations (s).  This can be a
        scalar if the same cutoff is used for the three components, or a 3
        element vector to indicate cutoff period for the `x`, `y`, `z`
        components.
    tilt_motion : array_like
        [roll, pitch] vector with mean tilt (radians) relative to a
        horizontal plane.  `Roll` angle offset is positive is port side up,
        and `pitch` is positive for bow down (see Miller 2008 for info, set
        to [0 0] if unknown).
    tilt_anemometer : array_like
        [roll, pitch] vector with mean tilt (radians) relative to a
        horizontal plane.  `Roll` angle offset is positive is port side up,
        and `pitch` is positive for bow down (see Miller 2008 for info, set
        to [0 0] if unknown).

    Returns
    -------
    namedtuple with (index, name in brackets):
    numpy.ndarray [0, 'uvw_ship']
        2-D array (Nx3) with corrected wind vectors in ship-referenced
        frame with z-axis parallel to gravity.
    numpy.ndarray [1, 'euler_angles']
        2-D array (Nx3) with Euler angles.
    numpy.ndarray [2, 'euler_angles_angular_rates']
        2-D array (NX3) with Euler angles from rate sensors (unfiltered).
    numpy.ndarray [3, 'euler_angles_accelerations']
        2-D array (NX3) with Euler angles from accelerometers (unfiltered).
    numpy.ndarray [4, 'euler_angles_slow']
        2-D array (NX3) with slow Euler angles (low pass filtered).
    numpy.ndarray [5, 'euler_angles_fast']
        2-D array (NX3) with fast Euler angles (high pass filtered).
    numpy.ndarray [6, 'mount_offset_rotations']
        2-D array (3X3) with mounting offset rotation matrix (see Miller 2008).
    numpy.ndarray [7, 'uvw_earth']
        2-D array (NX3) with measured velocity, rotated to the earth frame.
    numpy.ndarray [8, 'uvw_angular']
        2-D array (NX3) with velocity induced by angular motion.
    numpy.ndarray [9, 'uvw_linear]
        2-D array (NX3) with velocity induced by platform linear motion.
    numpy.ndarray [10, 'ship_enu']
        2-D array (NX3) with ship velocity in eastward, northward, up
        frame.  (low pass filtered).
    numpy.ndarray [11, 'uvw_enu']
        2-D array (NX3) with corrected wind in eastward, northward, up
        frame.
    numpy.ndarray [12, 'imu_enu']
        2-D (NX3) array with displacement of the motion sensor.

    Reference
    ---------
    Miller,S.D., Hristov,T.S., Edson,J.B., and C.A. Friehe, 2008:
        Platform Motion Effects on Measurements of Turbulence and Air-Sea
        Exchange Over the Open Ocean, J. Atmo. Ocean. Tech. 25(9),
        1683-1694, DOI: 10.1175/2008JTECHO547.1.

    """
    if len([Ta]) == 1:
        Ta = Ta * np.ones((3, 1))

    Astop, Apass = 10.0, 0.5     # stopband attenuation; passband Ripple (dB)
    # Stop, passband cutoffs
    wp = 1.0 / (2.0 * Tcf) / (sample_freq / 2.0)
    ws = 1.0 / Tcf / (sample_freq / 2.0)
    N, Wn = signal.buttord(wp, ws, Apass, Astop)
    bc, ac = signal.butter(N, Wn, "high")
    # WATCH THIS: we need to make padlen the same as in Matlab
    pdl = 3 * (max(len(ac), len(bc)) - 1)

    # EULER ANGLES: (see Edson et al., 1998)
    # Low frequency tilt using accelerometers
    g = np.sqrt(np.sum(np.mean(acceleration, 0) ** 2))
    ax, ay = acceleration[:, 0], acceleration[:, 1]
    phi_lf = (np.arctan2(ay, g) -
              signal.filtfilt(bc, ac, np.arctan2(ay, g), padlen=pdl))
    # High frequency angles using angle rates
    rm = signal.detrend(angle_rate, 0)
    rx, ry, rz = rm[:, 0], rm[:, 1], rm[:, 2]  # don't do this in ipdb!
    phi_hf = signal.filtfilt(bc, ac,
                             ((np.cumsum(rx) - 0.5 * rx - 0.5 * rx[0]) /
                              sample_freq), padlen=pdl)
    phi = phi_lf + phi_hf
    theta_lf = (np.arctan2(-ax * np.cos(phi), g) -
                signal.filtfilt(bc, ac, np.arctan2(-ax * np.cos(phi), g),
                                padlen=pdl))

    gyro = np.unwrap(np.radians(heading))  # Put heading in radians
    gyro0 = gyro[0]               # Initial heading
    gyro = gyro - gyro0           # Remove initial heading
    # Low-pass filter to retain low-frequency content, but not the offset
    psi_lf = gyro - signal.filtfilt(bc, ac, gyro, padlen=pdl)

    theta_hf = signal.filtfilt(bc, ac,
                               ((np.cumsum(ry) - 0.5 * ry - 0.5 * ry[0]) /
                                sample_freq), padlen=pdl)
    psi_hf = signal.filtfilt(bc, ac,
                             ((np.cumsum(rz) - 0.5 * rz - 0.5 * rz[0]) /
                              sample_freq), padlen=pdl)
    theta = theta_lf + theta_hf
    # s = psi_lf + psi_hf             # SPL: not used...

    # NONLINEAR ROTATION OF ANGULAR RATES FROM MOTION SENSOR-FRAME TO
    # EARTH-FRAME This next step puts the measured angular rates into the
    # earth reference frame.  Motion sensor frame-based anglular rates are
    # first rolled based on the roll angle (phim) estimated above; they are
    # then pitched *in a rolled frame* based upon the pitch angle (thetam)
    # estimated above; finally, they are yawed *in a rolled and pitched
    # frame* based upon both the pitch and roll (thetam and phim)
    # estimates. These sequential rotations are exact for small angle RATES
    # (see Goldstein) - I've found this step doesn't have much effect on
    # the corrected winds.
    F12 = np.tan(theta) * np.sin(phi)
    F13 = np.tan(theta) * np.cos(phi)
    F22 = np.cos(phi)
    F23 = -np.sin(phi)
    F32 = np.sin(phi) / np.cos(theta)
    F33 = np.cos(phi) / np.cos(theta)

    phi_dot = rx + ry * F12 + rz * F13
    theta_dot = ry * F22 - rz * F23
    psi_dot = ry * F32 + rz * F33

    # The angular rates phi_dot, theta_dot, and psi_dot are now in the
    # earth based frame, which is what we want for creating the
    # rotation-transformation matrix. Re-integrate and high pass filter
    # these angular rates to get the fast angles.
    phi_hf = signal.filtfilt(bc, ac,
                             (np.cumsum(phi_dot) -
                              0.5 * phi_dot -
                              0.5 * phi_dot[0]) / sample_freq, padlen=pdl)
    theta_hf = signal.filtfilt(bc, ac,
                               (np.cumsum(theta_dot) -
                                0.5 * theta_dot -
                                0.5 * theta_dot[0]) / sample_freq,
                               padlen=pdl)
    psi_hf = signal.filtfilt(bc, ac,
                             (np.cumsum(psi_dot) -
                              0.5 * psi_dot -
                              0.5 * psi_dot[0]) / sample_freq, padlen=pdl)

    # combine high- and low-frequency angles; make 2D matrix
    EA = np.column_stack((phi_lf + phi_hf,
                          theta_lf + theta_hf,
                          psi_lf + psi_hf))

    # SHIP AND ANEMOMETER OFFSET ANGLES
    # Pitch and roll tilt of platform w.r.t. earth (see Wilczak et al. 2000)
    phi_em, theta_em = tilt_motion[0], tilt_motion[1]
    PHI_em = np.array([[1, 0, 0], [0, np.cos(phi_em), -np.sin(phi_em)],
                       [0, np.sin(phi_em), np.cos(phi_em)]]).T
    THETA_em = np.array([[np.cos(theta_em), 0, np.sin(theta_em)],
                         [0, 1, 0],
                         [-np.sin(theta_em), 0, np.cos(theta_em)]]).T

    # Pitch and roll tilt of anemometer w.r.t. earth
    phi_ea, theta_ea = tilt_anemometer[0], tilt_anemometer[1]
    PHI_ea = np.array([[1, 0, 0], [0, np.cos(phi_ea), -np.sin(phi_ea)],
                       [0, np.sin(phi_ea), np.cos(phi_ea)]]).T
    THETA_ea = np.array([[np.cos(theta_ea), 0, np.sin(theta_ea)],
                         [0, 1, 0],
                         [-np.sin(theta_ea), 0, np.cos(theta_ea)]]).T

    # Anemometer WRT motion sensor
    M_em = np.dot(THETA_em, PHI_em)
    M_ea = np.dot(THETA_ea, PHI_ea)
    M_ma = np.dot(np.transpose(M_em), M_ea)

    # WIND VECTOR ROTATION
    # Using coordinate rotations
    Ur = euler_rotate(np.dot(M_ma, wind_speed.T).T, EA)

    # PLATFORM ANGULAR VELOCITY
    uam = np.column_stack((anemometer_pos[2] * ry -
                           anemometer_pos[1] * rz,
                           anemometer_pos[0] * rz -
                           anemometer_pos[2] * rx,
                           anemometer_pos[1] * rx -
                           anemometer_pos[0] * ry))
    Ua = euler_rotate(uam, EA)  # rotate to earth frame

    # PLATFORM LINEAR VELOCITY
    # High-pass filters for accelerometers
    wp = 1.0 / (2.0 * Ta[0]) / (sample_freq / 2.0)
    ws = 1.0 / Ta[0] / (sample_freq / 2.0)   # stop, passband cutoffs
    N, Wn = signal.buttord(wp, ws, Apass, Astop)
    bax, aax = signal.butter(N, Wn, "high")
    # Again, watch padlen
    pdlx = 3 * (max(len(aax), len(bax)) - 1)

    wp = 1.0 / (2.0 * Ta[1]) / (sample_freq / 2.0)
    ws = 1.0 / Ta[1] / (sample_freq / 2.0)   # stop, passband cutoffs
    N, Wn = signal.buttord(wp, ws, Apass, Astop)
    bay, aay = signal.butter(N, Wn, "high")
    pdly = 3 * (max(len(aay), len(bay)) - 1)

    wp = 1.0 / (2.0 * Ta[2]) / (sample_freq / 2.0)
    ws = 1.0 / Ta[2] / (sample_freq / 2.0)   # stop, passband cutoffs
    N, Wn = signal.buttord(wp, ws, Apass, Astop)
    baz, aaz = signal.butter(N, Wn, "high")
    pdlz = 3 * (max(len(aaz), len(baz)) - 1)

    # Rotate accelerations
    ae = euler_rotate(acceleration, EA)
    ae[:, 2] = ae[:, 2] - g     # subtract gravity
    n = np.size(ae, 0)
    up = (np.cumsum(ae, 0) -
          0.5 * ae - 0.5 * (np.ones((n, 1)) * ae[0, :])) / sample_freq
    Up = np.column_stack((signal.filtfilt(bax, aax, up[:, 0], padlen=pdlx),
                          signal.filtfilt(bay, aay, up[:, 1], padlen=pdly),
                          signal.filtfilt(baz, aaz, up[:, 2], padlen=pdlz)))

    # SHIP SPEED: this is in the positive x-direction (assuming the bow is
    # positive x direction).  We low-pass filter this signal so there isnt
    # double counting of linear motion from the integrated accelerometers.
    # Also, the y- and z-components are zero (ie, captured by the
    # integrated accelerometers).
    n = len(speed)
    Us = np.column_stack((speed -
                          signal.filtfilt(bax, aax, speed, padlen=pdlx),
                          np.zeros((n, 1)), np.zeros((n, 1))))

    UVW = Ur + Ua + Up + Us     # corrected wind vector

    # Organize outputs
    EA_rate = np.column_stack(((np.cumsum(rx) - 0.5 * rx - 0.5 * rx[0]) /
                               sample_freq,
                               (np.cumsum(ry) - 0.5 * ry - 0.5 * ry[0]) /
                               sample_freq,
                               (np.cumsum(rz) - 0.5 * rz - 0.5 * rz[0]) /
                               sample_freq))
    EA_acc = np.column_stack((np.arctan2(ay, g),
                              np.arctan2(-ax, g),
                              gyro))
    EA_slow = np.column_stack((phi_lf, theta_hf, psi_hf))
    EA_fast = np.column_stack((phi_hf, theta_hf, psi_hf))
    u_ea = np.column_stack((np.zeros((n, 1)),
                            np.zeros((n, 1)),
                            (gyro0 + np.pi / 2) * np.ones((n, 1))))
    U_ship = euler_rotate(Us, u_ea)
    U_earth = euler_rotate(UVW, u_ea)
    # Platform displacement
    xp = (np.cumsum(Up, 0) - 0.5 * Up -
          0.5 * (np.ones((n, 1)) * Up[0, :])) / sample_freq
    Xp = np.column_stack((signal.filtfilt(bax, aax, xp[:, 0], padlen=pdlx),
                          signal.filtfilt(bay, aay, xp[:, 1], padlen=pdly),
                          signal.filtfilt(baz, aaz, xp[:, 2], padlen=pdlz)))

    return CorrectedWind3D(UVW, EA, EA_rate, EA_acc, EA_slow, EA_fast,
                           M_ma, Ur, Ua, Up, U_ship, U_earth, Xp)


def window_indices(idxs, width, step=None):
    """List of sliding window indices across an index vector

    Parameters
    ----------
    idx : numpy.ndarray
        A 1-D index vector.
    width : int
        Window width.
    step : int
        Step size for sliding windows.

    Returns
    -------
    List of tuples, each with the indices for a window.

    """
    return zip(*(idxs[i::step] for i in range(width)))


def get_VickersMahrt(x, zscore_thr, nrep_thr):
    """Vickers Mahrt computations in a window

    Parameters
    ----------
    x : numpy.ndarray
        A 1-D signal vectors to be despiked.
    zscore_thr : float
        The zscore beyond which an observation is considered to be an
        outlier.
    nrep_thr : int
        The maximum number of consecutive outliers that should occur for a
        spike to be detected.

    Returns
    -------
    Tuple with (index, name in brackets):
    numpy.ndarray [0, 'x']
        1-D array with interpolated input.
    numpy.int [1, 'nspikes']
        Number of spikes detected.
    numpy.int [2, 'ntrends']
        Number of outlier trends detected.
    numpy.ndarray [3, 'kclass']
        1-D array of the same size as input, indicating the classification
        `k` for each measurement. k=0: measurement within plausibility
        range, k=[-1 or 1]: measurement outside plausibility range, abs(k)
        > 1: measurement is part of an outlier trend.

    """
    z = zscore(x)
    # Discern between outliers above and below the threshold
    isout_hi, isout_lo = (z > zscore_thr), (z < -zscore_thr)
    n_outs = sum(isout_hi) + sum(isout_lo)
    # Set categorical x: 0 (ok), 1 (upper outlier), -1 (lower outlier)
    xcat = np.zeros(x.shape, dtype=np.int)
    xcat[isout_hi] = 1
    xcat[isout_lo] = -1
    if n_outs > 0:
        # Create tuples for each sequence indicating whether it's outliers
        # and its length.
        grps = [(val, len(list(seq))) for val, seq in groupby(xcat)]
        vals = np.array([k[0] for k in grps])
        lens = np.array([k[1] for k in grps])
        is_spike = (vals != 0) & (lens <= nrep_thr)
        nspikes = sum(is_spike)
        # We tally trends as well
        is_trend = (vals != 0) & (lens > nrep_thr)
        ntrends = sum(is_trend)
        # If we have trends, loop through each one, knowing the length of
        # the spike and where we are along the input series.
        if ntrends > 0:
            trends = zip(vals[is_trend], lens[is_trend],
                         np.cumsum(lens)[is_trend])
            for i in trends:
                # Double the categorical value for trends and consider
                # these OK for interpolation. So abs(xcat) > 1 are trends.
                xcat[(i[2] - i[1]):i[2]] = i[0] * 2
        # Now we are left with true outliers to interpolate
        x_new = x.copy()
        xidx = np.arange(len(x))             # simple index along x
        isok = (xcat == 0) | (abs(xcat) > 1)  # ok if 0 or trend
        s = itpl.InterpolatedUnivariateSpline(xidx[isok],
                                              x[isok], k=1)
        x_itpl = s(xidx[~ isok])
        x_new[~ isok] = x_itpl
        return VickersMahrt(x_new, nspikes, ntrends, xcat)
    else:
        return VickersMahrt(x, 0, 0, xcat)


def despike_VickersMahrt(x, width, zscore_thr, nreps, step=None,
                         nrep_thr=None, interp_nan=True):
    """Vickers and Mahrt (1997) signal despiking procedure

    The interpolating function is created by the
    InterpolatedUnivariateSpline function from the scipy package, and uses
    a single knot to approximate a simple linear interpolation, so as to
    keep the original signal as untouched as possible.

    Parameters
    ----------
    x : numpy.ndarray
        A 1-D signal vectors to be despiked.
    width : int
        Window width.
    step : int
        Step size for sliding windows.  Default is one-half window width.
    zscore_thr : float
        The zscore beyond which an observation is considered to be an
        outlier.
    nrep_thr : int
        The maximum number of consecutive outliers that should occur for a
        spike to be detected.  Default: 3.
    nreps: int
        How many times to run the procedure.
    interp.nan : bool
        Whether missing values should be interpolated.  Interpolated values
        are computed after despiking.

    Returns
    -------
    Tuple with (index, name in brackets):
    numpy.ndarray [0, 'x']
        1-D array with despiked input.
    numpy.int [1, 'nspikes']
        Number of spikes detected.
    numpy.int [2, 'ntrends']
        Number of outlier trends detected.
    numpy.int [3, 'kclass']
        Number of iterations performed.

    """
    if step is None:            # set default step as
        step = width / 2        # one-half window size
    if nrep_thr is None:
        nrep_thr = 3
    nspikes, ntrends = 0, 0
    xout = x.copy()
    # Following EddyUH implementation, fill missing values with nearest
    # value for the purpose of spike and trend detection.  This ensures
    # that we can always calculate zscores.
    is_missing = np.isnan(np.array(xout))  # need to coerce to np array
    xidx = np.arange(len(xout))            # simple index along x
    f_itpl = itpl.interp1d(xidx[~ is_missing], xout[~ is_missing],
                           kind="nearest", fill_value="extrapolate")
    x_nonan = f_itpl(xidx)
    # Get a series of tuples with indices for each window
    idxl = window_indices(range(len(x)), width, step)
    nloops = 0
    while nloops < nreps:
        nspikes_loop = 0
        for w in idxl:
            winidx = [i for i in w]  # indices of current window
            xwin = x_nonan[winidx]   # values for the current window
            xnew, nsp, ntr, xmask = get_VickersMahrt(xwin, zscore_thr,
                                                     nrep_thr)
            nspikes_loop += nsp
            ntrends += ntr
            x_nonan[winidx] = xnew
        nloops += 1
        # Increase zscore_thr by 0.3, instead of 0.1 as in V&M (1997),
        # following EddyUH implementation
        zscore_thr += 0.3
        if nspikes_loop > 0:
            nspikes += nspikes_loop
        else:                   # Stop if we haven't found any new spikes
            break
    # Interpolate through missing values, if requested (default).
    nmissing = np.count_nonzero(is_missing)
    if (nmissing > 0) and interp_nan:
        s = itpl.InterpolatedUnivariateSpline(xidx[~ is_missing],
                                              xout[~ is_missing], k=1)
        x_itpl = s(xidx[is_missing])
        xout[is_missing] = x_itpl

    return VickersMahrt(xout, nspikes, ntrends, nloops)
