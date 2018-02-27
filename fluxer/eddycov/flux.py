# pylint: disable=too-many-locals,invalid-name,no-member

"""Core functionality for the package

Library of main functions required for flux calculations.

"""

from collections import namedtuple
from itertools import groupby
import numpy as np
from scipy import interpolate as itpl
from scipy import signal, integrate
from scipy.stats import zscore
from astropy.convolution import convolve, Box1DKernel


__all__ = ["smooth_angle", "planarfit", "rotate_coordinates",
           "rotate_wind3d", "wind3D_correct", "despike_VickersMahrt"]

# Valid 3-D wind rotation methods
_VECTOR_ROTATION_METHODS = {"DR", "TR", "PF"}


AngleCoordinates = namedtuple("AngleCoordinates", ["x", "y"])
Vector = namedtuple("Vector", ["angle", "magnitude"])
PlanarFitCoefs = namedtuple("PlanarFitCoefs",
                            ["k_vct", "tilt_coefs", "phi", "theta"])
RotatedVectors = namedtuple("RotatedVectors", ["rotated", "phi_theta"])
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
        that is common to all ``angle`` elements.

    Returns
    -------
    AngleCoordinates : namedtuple
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
    """Recompose angles and associated magnitudes from ``x`` and ``y`` vectors

    Parameters
    ----------
    x : array_like or scalar
        ``x``-coordinates
    y : array_like or scalar
        ``y``-coordinates

    Returns
    -------
    Vector : namedtuple
        namedtuple with ndarrays ``angle`` and ``magnitude``, in that
        order.

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
        that is common to all ``angle`` elements.
    kernel_width : int, optional
        The width of the filter kernel.

    Returns
    -------
    namedtuple
        namedtuple with ndarrays ``angle`` and ``magnitude``, in that
        order.

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
    getPlanarFitCoeffs.m from Patric Sturm <pasturm@ethz.ch>.

    Parameters
    ----------
    vectors : numpy.ndarray
        A 2-D (Nx3) array with `x`, `y`, and `z` vectors, expressed in a
        right-handed coordinate system.  These vectors may correspond to
        `u`, `v`, and `w` wind speed vectors, or inertial acceleration
        components.

    Returns
    -------
    PlanarFitCoefs : namedtuple
        namedtuple with (index and name in brackets):

        numpy.ndarray [0, `k_vct`]
            1-D array (1x3) unit vector parallel to the new z-axis.
        numpy.ndarray [1, `tilt_coefs`]
            1-D array (1x3) Tilt coefficients `b0`, `b1`, `b2`.
        numpy.float [2, `phi`]
            Scalar representing roll angle :math:`\\phi`.
        numpy.float [3, `theta`]
            Scalar representing pitch angle :math:`\\theta`.

    """
    vct_u = vectors[:, 0]
    vct_v = vectors[:, 1]
    vct_w = vectors[:, 2]
    vct_nrows = vectors.shape[0]
    sum_u = sum(vct_u)
    sum_v = sum(vct_v)
    sum_w = sum(vct_w)
    dot_uv = np.dot(vct_u, vct_v)
    dot_uw = np.dot(vct_u, vct_w)
    dot_vw = np.dot(vct_v, vct_w)
    dot_u2 = np.dot(vct_u, vct_u)
    dot_v2 = np.dot(vct_v, vct_v)
    H_arr = np.array([[vct_nrows, sum_u, sum_v],
                      [sum_u, dot_u2, dot_uv],
                      [sum_v, dot_uv, dot_v2]])
    g_arr = np.array([sum_w, dot_uw, dot_vw])
    tilt_coef = np.linalg.solve(H_arr, g_arr)
    # Estimated \phi (roll) and \theta (pitch) tilt angles
    phi_denom = np.sqrt(1 + (tilt_coef[2] ** 2))
    phi_sin = tilt_coef[2] / phi_denom
    phi_cos = 1 / phi_denom
    phi = np.arctan2(phi_sin, phi_cos)
    theta_denom = np.sqrt(1 + (tilt_coef[1] ** 2) + (tilt_coef[2] ** 2))
    theta_sin = -tilt_coef[1] / theta_denom
    theta_cos = np.sqrt((tilt_coef[2] ** 2) + 1) / theta_denom
    theta = np.arctan2(theta_sin, theta_cos)
    # Determine unit vector parallel to new z-axis
    k_2 = 1 / np.sqrt(1 + tilt_coef[1] ** 2 + tilt_coef[2] ** 2)
    k_0 = -tilt_coef[1] * k_2
    k_1 = -tilt_coef[2] * k_2
    k_vct = np.array([k_0, k_1, k_2])
    return PlanarFitCoefs(k_vct, tilt_coef, phi, theta)


def rotation_matrix(theta, axis, active=False):
    """Generate rotation matrix for a given axis

    The default (``active=False``) output matrix corresponds to passive
    ("alias") rotations, if used to post-multiply row vectors.

    Parameters
    ----------
    active: bool, optional
        Whether to return active transformation matrix.

    Note
    ----
    See ``rotate_coordinates`` for rest of parameters.

    Returns
    -------
    R_theta : numpy.ndarray
        3x3 rotation matrix

    """
    theta = np.radians(theta)
    if axis == 0:
        R_theta = np.array([[1, 0, 0],
                            [0, np.cos(theta), -np.sin(theta)],
                            [0, np.sin(theta), np.cos(theta)]])
    elif axis == 1:
        R_theta = np.array([[np.cos(theta), 0, np.sin(theta)],
                            [0, 1, 0],
                            [-np.sin(theta), 0, np.cos(theta)]])
    else:
        R_theta = np.array([[np.cos(theta), -np.sin(theta), 0],
                            [np.sin(theta), np.cos(theta), 0],
                            [0, 0, 1]])
    if active:
        R_theta = np.transpose(R_theta)
    return R_theta


def rotate_coordinates(vectors, theta=0, axis=2, rotate_vectors=False):
    """Rotate vector coordinate system or vectors themselves around an axis

    A right-handed coordinate system is assumed, where positive rotations
    are clockwise when looking along the rotation axis from the origin
    towards the positive direction.  With the default
    ``rotate_vectors=False``, each vector :math:`i` in input coordinate
    system :math:`0` is rotated around the given axis by angle
    :math:`\\theta`, so that it can be expressed in coordinate system
    :math:`1` as follows:

    .. math:: \\vec{v}_{i,1} = \\vec{v}_{i,0} R_{\\theta}

    where the input and output coordinate systems are different.  If the
    input and output coordinate systems are the same,
    i.e. ``rotate_vectors=True``, the vectors are transformed according to:

    .. math:: \\vec{v}_{i,1} = \\vec{v}_{i,0} R_{\\theta}^\\intercal

    Parameters
    ----------
    vectors : array_like
        An nx3 array of vectors with their `x`, `y`, `z` components
    theta : numeric, optional
        The angle (degrees) by which to perform the rotation.  Default is
        0, which means return the coordinates of the vector in the rotated
        coordinate system, when rotate_vectors=False.
    axis : int, optional
        Axis around which to perform the rotation (`x` = 0; `y` = 1; `z` =
        2)
    rotate_vectors : bool, optional
        Whether to return the coordinates of each vector in the rotated
        coordinate system or to rotate the vectors themselves onto the same
        coordinate system.  Default is to perform a passive transformation,
        i.e. rotation of the coordinate system.

    Returns
    -------
    array_like
    The vector array with the same shape as input with rotated components.

    """
    R_theta = rotation_matrix(theta, axis, active=rotate_vectors)
    # Multiplying a vector by R_theta, as defined above, transforms the
    # vector coordinates from the original coordinate system to the new,
    # rotated coordinate frame.  Note that the rotation matrix
    # post-multiplies the vector matrix.  If rotating the vectors
    # themselves, then the *transposed* rotation matrix post-multiplies the
    # matrix of vectors.
    return np.dot(vectors, R_theta)


def rotate_wind3d(wind3D, method="PF", **kwargs):
    """Transform 3D wind vectors to reference mean streamline coordinate system

    Use double rotation, triple rotation, or planar fit methods (Wilczak et
    al. 2001; Handbook of Micrometeorology).

    This is a general coordinate rotation tool, so can handle inputs such
    as wind speed and acceleration from inertial measurement units.

    Parameters
    ----------
    wind3D : numpy.ndarray
        A 2-D (Nx3) array with `x`, `y`, and `z` vector components,
        expressed in a right-handed coordinate system.  These may represent
        `u`, `v`, and `w` wind speed vectors, or inertial acceleration.
    method : {"DR", "TR", "PF"}, optional
        One of: "DR", "TR", "PF" for double rotation, triple rotation, or
        planar fit.
    k_vector : numpy.ndarray, optional
        1-D array (1x3) unit vector parallel to the new z-axis, when
        "method" is "PF" (planar fit).  If not supplied, then it is
        calculated.

    Returns
    -------
    RotatedVectors : namedtuple
        namedtuple with (index, name in brackets):

        numpy.ndarray [0, `rotated`]
            2-D array (Nx3) Array with rotated vectors
        numpy.ndarray [1, `phi_theta`]
            1-D array (1x2) :math:`\\phi` and :math:`\\theta` rotation
            angles.  The former is the estimated angle between the vertical
            unit coordinate vector in the rotated frame and the vertical
            unit vector in the measured uv plane, while the latter is wind
            direction in the measured uv plane.  Note these are *not* roll
            and pitch angles of the measurement coordinate frame relative
            to the reference frame.

    """
    if method not in _VECTOR_ROTATION_METHODS:
        msg = "method must be one of "
        raise ValueError(msg + ', '.join("\"{}\"".format(m) for m in
                                         _VECTOR_ROTATION_METHODS))

    if method == "PF":
        if "k_vector" in kwargs:
            k_vct = kwargs.get("k_vector")
        else:
            pfit = planarfit(wind3D)
            k_vct = pfit.k_vct
        j_vct = np.cross(k_vct, np.mean(wind3D, 0))
        j_vct = j_vct / np.sqrt(np.sum(j_vct ** 2))
        i_vct = np.cross(j_vct, k_vct)
        vcts_mat = np.column_stack((i_vct, j_vct, k_vct))
        vcts_new = np.dot(wind3D, vcts_mat)
        phi = np.arccos(np.dot(k_vct, np.array([0, 0, 1])))
        theta = np.arctan2(np.mean(wind3D[:, 1], 0),
                           np.mean(wind3D[:, 0], 0))
    else:
        # First rotation to set mean v to 0
        theta = np.arctan2(np.mean(wind3D[:, 1]),
                           np.mean(wind3D[:, 0]))
        rot1 = np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0],
                         [0, 0, 1]])
        vcts1 = np.dot(wind3D, rot1)
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


def _rot_xyz(phi, theta, psi):
    """Compute passive rotation matrix from 3 angles in 1-D array x

    The matrix is meant to post-multiply row vectors.

    Parameters
    ----------
    phi : float
        Angle (degrees) for the rotation around the `x`-axis.
    theta : float
        Angle (degrees) for the rotation around the `y`-axis.
    psi : float
        Angle (degrees) for the rotation around the `z`-axis.

    Returns
    -------
    numpy.ndarray
        2-D array (3x3) with rotation matrix corresponding to the rotation
        sequence :math:`\\phi`, :math:`\\theta`, :math:`\\psi`.

    """
    rotx = rotation_matrix(phi, 0)
    roty = rotation_matrix(theta, 1)
    rotz = rotation_matrix(psi, 2)
    return np.transpose(np.dot(rotz, np.dot(roty, rotx)))


def euler_rotate(xyz, xyz_angles):
    """Rotate vectors using Euler angles

    Parameters
    ----------
    xyz : array_like
        2-D array (Nx3) with row vectors to be rotated.
    xyz_angles : array_like
        2-D array (Nx3) with row vectors of angles (degrees) to be used for
        rotation of ``xyz``, where the first column specifies the angle for
        the first rotation around the `x`-axis, the second column the angle
        for the second rotation around the `y`-axis, and the third column
        the angle for the last rotation around the `z`-axis.

    Returns
    -------
    numpy.ndarray
        Array of the same shape as ``xyz`` with rotated vectors.

    """
    xyz_rots = np.empty_like(xyz)
    for i, v in enumerate(xyz_angles):
        phi, theta, psi = v
        euler_mati = _rot_xyz(phi, theta, psi)
        xyz_roti = np.dot(xyz[i], euler_mati)
        xyz_rots[i] = xyz_roti
    return xyz_rots


def euler_rate_rotate(euler_angles, omega):
    """Rotate angular rates, given Euler angles

    The Euler angular rate column vector :math:`\\Omega` consists of the
    angular rates :math:`\\dot \\phi`, :math:`\\dot \\theta`, and
    :math:`\\dot \\psi` around the axes in the inertial frame.  It is
    related to the angular rate column vector :math:`\\omega`, consisting
    of the angular rates :math:`\\omega_x`, :math:`\\omega_y`, and
    :math:`\\omega_z` in the body frame via a non-orthogonal transformation
    matrix :math:`R`:

    .. math:: \\Omega = R \\omega

    This function performs this operation by using the transpose of the
    input body-frame angular rate row vectors, and transposing the result
    to return row vectors also.

    Parameters
    ----------
    euler_angles : array_like
        2-D array (Nx2) with Euler :math:`\\phi` (roll) and :math:`\\theta`
        (pitch) angles (radians) around the `x`- and `y`-axes in columns 1
        and 2, respectively, representing the orientation of the sensor
        orientation relative to the output frame.
    omega : float 2-D array
        (Nx3) with angular rates :math:`\\omega_x`, :math:`\\omega_y`, and
        :math:`\\omega_z` around the body-frame `x`-, `y`-, and `z`-axes in
        columns 1, 2, and 3, respectively.

    Returns
    -------
    array_like
        Array with same shape as `omega` with rotated angular rates.

    """
    def rot_mat(phi, theta):
        """Transformation matrix pre-multiplying column vectors"""
        return np.array([[1, np.sin(phi) * np.tan(theta),
                          np.cos(phi) * np.tan(theta)],
                         [0, np.cos(phi), -np.sin(phi)],
                         [0, np.sin(phi) / np.cos(theta),
                          np.cos(phi) / np.cos(theta)]])
    omega_rots = np.empty_like(omega)
    for i, v in enumerate(euler_angles):
        phi, theta = v
        # Multiply and transpose to original shape as row vectors
        omega_rots[i] = np.dot(rot_mat(phi, theta), omega[i].T).T
    return omega_rots


def _cumtrapz(x, sample_freq):
    """Cumulatively integrate signal given sampling frequency

    Parameters
    ----------
    x : array_like
        A 1- or 2-D array with signal such as angular rate (deg/s) or
        linear rate measurements.
    sample_freq : int, float
        The sampling frequency in units required for integration.

    Returns
    -------
    array_like
        Array with same shape as `x`

    """
    y = integrate.cumtrapz(x, dx=(1.0 / sample_freq), axis=0, initial=0)
    return y


def _butterworth_coefs(cutoff_period, sample_freq, Astop=10.0, Apass=0.5):
    """Compute Butterworth filter coefficients

    Parameters
    ----------
    cutoff_period : int, float
        Cutoff period (s) defining the cutoff frequency for designing the
        filter.
    sample_freq : int, float
        The sampling frequency in units required for integration.
    Apass : float, optional
        Stopband attenuation
    Astop : float, optional
        Passband ripple (dB)

    Returns
    -------
    tuple
        Tuple with (index, name in brackets):

        numpy.ndarray [0]
            Numerator (b) polynomial of the filter.
        numpy.ndarray [1]
            Denominator (a) polynomial of the filter.
        float [2]
            Padding length for the filter.

    """
    # Passband and stopband cutoffs, normalized to Nyquist frequency
    wp = 1.0 / (2.0 * cutoff_period) / (sample_freq / 2.0)
    ws = 1.0 / cutoff_period / (sample_freq / 2.0)
    N, Wn = signal.buttord(wp, ws, Apass, Astop)
    bc, ac = signal.butter(N, Wn, "high")
    # WATCH THIS: we need to make padlen the same as in Matlab
    pdl = 3 * (max(len(ac), len(bc)) - 1)
    return bc, ac, pdl


def wind3D_correct(wind_speed, acceleration, angle_rate, heading, speed,
                   anemometer_pos, sample_freq, Tcf, Ta,
                   tilt_motion=np.array([0.0, 0.0]),
                   tilt_anemometer=np.array([0.0, 0.0])):
    """Correct wind vector measurements from a moving platform

    This is a port of S. Miller's ``motion`` Matlab function, which
    implements Miller et al. (2008) approach.  Coordinate frame is assumed
    to be right-handed:

    * `x` - positive towards the bow.
    * `y` - positive towards port side.
    * `z` - positive upwards.

    Parameters
    ----------
    wind_speed : numpy.ndarray
        A 2-D (Nx3) array with `u`, `v`, and `w` wind speed (m/s) vectors.
    acceleration : numpy.ndarray
        A 2-D (Nx3) array with `x`, `y`, and `z` acceleration (m/s/s)
        vectors.
    angle_rate : numpy.ndarray
        A 2-D (Nx3) array with angular rates (radians/s) vectors.
    heading : array_like
        A vector with ship heading measurements (deg) in right-hand
        coordinate system.
    speed : array_like
        A vector with ship speed (m/s).
    anemometer_pos : array_like
        [`x`, `y`, `z`] position vector of anemometer, relative to motion
        sensor (m).
    sample_freq : float
        Sampling frequency (Hz).
    Tcf : float
        Complimentary filter period (s).
    Ta : float
        High-pass filter cutoff for accelerations (s).  This can be a
        scalar if the same cutoff is used for the three components, or a 3
        element vector to indicate cutoff period for the `x`, `y`, `z`
        components.
    tilt_motion : array_like, optional
        [roll, pitch] vector with mean tilt (radians) relative to a
        horizontal plane.  Roll angle offset is positive is port side up,
        and pitch is positive for bow down (see Miller 2008 for info, set
        to [0 0] if unknown).
    tilt_anemometer : array_like, optional
        [roll, pitch] vector with mean tilt (radians) relative to a
        horizontal plane.  Roll angle offset is positive is port side up,
        and pitch is positive for bow down (see Miller 2008 for info, set
        to [0 0] if unknown).

    Returns
    -------
    CorrectedWind3D : namedtuple
        namedtuple with (index, name in brackets):

        numpy.ndarray [0, `uvw_ship`]
            2-D array (Nx3) with corrected wind vectors in ship-referenced
            frame with z-axis parallel to gravity.
        numpy.ndarray [1, `euler_angles`]
            2-D array (Nx3) with Euler angles.
        numpy.ndarray [2, `euler_angles_angular_rates`]
            2-D array (NX3) with Euler angles from rate sensors
            (unfiltered).
        numpy.ndarray [3, `euler_angles_accelerations`]
            2-D array (NX3) with Euler angles from accelerometers
            (unfiltered).
        numpy.ndarray [4, `euler_angles_slow`]
            2-D array (NX3) with slow Euler angles (low pass filtered).
        numpy.ndarray [5, `euler_angles_fast`]
            2-D array (NX3) with fast Euler angles (high pass filtered).
        numpy.ndarray [6, `mount_offset_rotations`]
            2-D array (3X3) with mounting offset rotation matrix (see
            Miller 2008).
        numpy.ndarray [7, `uvw_earth`]
            2-D array (NX3) with measured velocity, rotated to the earth
            frame.
        numpy.ndarray [8, `uvw_angular`]
            2-D array (NX3) with velocity induced by angular motion.
        numpy.ndarray [9, `uvw_linear`]
            2-D array (NX3) with velocity induced by platform linear
            motion.
        numpy.ndarray [10, `ship_enu`]
            2-D array (NX3) with ship velocity in eastward, northward, up
            frame.  (low pass filtered).
        numpy.ndarray [11, `uvw_enu`]
            2-D array (NX3) with corrected wind in eastward, northward, up
            frame.
        numpy.ndarray [12, `imu_enu`]
            2-D (NX3) array with displacement of the motion sensor.

    References
    ----------
    Miller,S.D., Hristov,T.S., Edson,J.B., and C.A. Friehe, 2008:
        Platform Motion Effects on Measurements of Turbulence and Air-Sea
        Exchange Over the Open Ocean, J. Atmo. Ocean. Tech. 25(9),
        1683-1694, DOI: 10.1175/2008JTECHO547.1.

    """
    Astop, Apass = 10.0, 0.5     # stopband attenuation; passband Ripple (dB)
    # Stop, passband cutoffs
    bc, ac, pdl = _butterworth_coefs(Tcf, sample_freq,
                                     Astop=Astop, Apass=Apass)

    # EULER ANGLES: (see Edson et al., 1998)
    # Low frequency tilt using accelerometers and gyro
    g = np.sqrt(np.sum(np.mean(acceleration, 0) ** 2))
    gyro = np.unwrap(-np.radians(heading))  # Put heading in radians and RHS
    gyro0 = gyro[0]               # Initial heading
    gyro = gyro - gyro0           # Remove initial heading
    # Uncorrected Euler angles from accelerometers and gyro
    EA_acc = np.column_stack((np.arcsin(acceleration[:, 1] / g),  # phi
                              np.arcsin(-acceleration[:, 0] / g),  # theta
                              gyro))                               # psi
    # High frequency angles using angular rates
    rm = signal.detrend(angle_rate, 0, "constant")
    EA_rate = _cumtrapz(rm, sample_freq)  # Estimated Euler angles
    EA_rate_hf = signal.filtfilt(bc, ac, EA_rate, padlen=pdl, axis=0)
    theta_lf = (EA_acc[:, 1] -
                signal.filtfilt(bc, ac, EA_acc[:, 1], padlen=pdl))
    theta = theta_lf + EA_rate_hf[:, 1]
    phi_pitched = EA_acc[:, 0] / np.cos(theta_lf)
    phi_lf = phi_pitched - signal.filtfilt(bc, ac, phi_pitched, padlen=pdl)
    phi = phi_lf + EA_rate_hf[:, 0]  # complementary-filtered phi
    # Low-pass filter to retain low-frequency content, but not the offset
    psi_lf = gyro - signal.filtfilt(bc, ac, gyro, padlen=pdl)
    axg = np.arctan2(-acceleration[:, 0] * np.cos(phi), g)
    theta_lf = axg - signal.filtfilt(bc, ac, axg, padlen=pdl)
    theta = theta_lf + EA_rate_hf[:, 1]
    EA_slow = np.column_stack((phi_lf, theta_lf, psi_lf))

    # NONLINEAR ROTATION OF ANGULAR RATES FROM MOTION SENSOR-FRAME TO
    # EARTH-FRAME This next step puts the measured angular rates into the
    # earth reference frame.  Motion sensor frame-based angular rates are
    # first rolled based on the roll angle (phim) estimated above; they are
    # then pitched *in a rolled frame* based upon the pitch angle (thetam)
    # estimated above; finally, they are yawed *in a rolled and pitched
    # frame* based upon both the pitch and roll (thetam and phim)
    # estimates. These sequential rotations are exact for small angle RATES
    # (see Goldstein) - I've found this step doesn't have much effect on
    # the corrected winds.
    omega = euler_rate_rotate(np.column_stack((phi, theta)), rm)
    # The angular rates omega are now in the earth based frame, which is
    # what we want for creating the rotation-transformation
    # matrix. Re-integrate and high pass filter these angular rates to get
    # the fast angles.
    EA_fast = signal.filtfilt(bc, ac,
                              _cumtrapz(omega, sample_freq),
                              padlen=pdl, axis=0)

    # combine high- and low-frequency angles; make 2D transformation matrix
    EA = EA_slow + EA_fast

    # SHIP AND ANEMOMETER OFFSET ANGLES
    # Pitch and roll tilt of platform w.r.t. earth (see Wilczak et al. 2000)
    phi_em, theta_em = tilt_motion[0], tilt_motion[1]
    PHI_em = rotation_matrix(np.degrees(phi_em), axis=0, active=True)
    THETA_em = rotation_matrix(np.degrees(theta_em), axis=1, active=True)
    # Pitch and roll tilt of anemometer w.r.t. earth
    phi_ea, theta_ea = tilt_anemometer[0], tilt_anemometer[1]
    PHI_ea = rotation_matrix(np.degrees(phi_ea), axis=0, active=True)
    THETA_ea = rotation_matrix(np.degrees(theta_ea), axis=1, active=True)

    # Anemometer WRT motion sensor
    M_em = np.dot(THETA_em, PHI_em)
    M_ea = np.dot(THETA_ea, PHI_ea)
    M_ma = np.dot(np.transpose(M_em), M_ea)

    # WIND VECTOR ROTATION
    # Using coordinate rotations
    EA_degs = np.degrees(EA)    # only for getting rotation matrix
    Ur = euler_rotate(np.dot(M_ma, wind_speed.T).T, EA_degs)

    # PLATFORM ANGULAR VELOCITY
    uam = np.column_stack((anemometer_pos[2] * rm[:, 1] -
                           anemometer_pos[1] * rm[:, 2],
                           anemometer_pos[0] * rm[:, 2] -
                           anemometer_pos[2] * rm[:, 0],
                           anemometer_pos[1] * rm[:, 0] -
                           anemometer_pos[0] * rm[:, 1]))
    Ua = euler_rotate(uam, EA_degs)  # rotate to earth frame

    # PLATFORM LINEAR VELOCITY
    # High-pass filters for accelerometers
    if np.isscalar(Ta):
        Ta = np.repeat(Ta, 3)
    bax, aax, pdlx = _butterworth_coefs(Ta[0], sample_freq,
                                        Astop=Astop, Apass=Apass)
    bay, aay, pdly = _butterworth_coefs(Ta[1], sample_freq,
                                        Astop=Astop, Apass=Apass)
    baz, aaz, pdlz = _butterworth_coefs(Ta[2], sample_freq,
                                        Astop=Astop, Apass=Apass)
    # Rotate accelerations
    ae = euler_rotate(acceleration, EA_degs)
    ae[:, 2] = ae[:, 2] - g     # subtract gravity
    up = _cumtrapz(ae, sample_freq)
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
    u_ea = np.column_stack((np.zeros((n, 1)),
                            np.zeros((n, 1)),
                            (gyro0 + np.pi / 2) * np.ones((n, 1))))
    u_ea_degs = np.degrees(u_ea)
    U_ship = euler_rotate(Us, u_ea_degs)
    U_earth = euler_rotate(UVW, u_ea_degs)
    # Platform displacement
    xp = _cumtrapz(Up, sample_freq)
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
    step : int, optional
        Step size for sliding windows.

    Returns
    -------
    list
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
    VickersMahrt : tuple
        Tuple with (index, name in brackets):

        numpy.ndarray [0, `x`]
            1-D array with interpolated input.
        numpy.int [1, `nspikes`]
            Number of spikes detected.
        numpy.int [2, `ntrends`]
            Number of outlier trends detected.
        numpy.ndarray [3, `kclass`]
            1-D array of the same size as input, indicating the
            classification `k` for each measurement. k=0: measurement
            within plausibility range, k=[-1 or 1]: measurement outside
            plausibility range, abs(k) > 1: measurement is part of an
            outlier trend.

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
    step : int, optional
        Step size for sliding windows.  Default is one-half window width.
    zscore_thr : float
        The zscore beyond which an observation is considered to be an
        outlier.  Default is zero.
    nrep_thr : int, optional
        The maximum number of consecutive outliers that should occur for a
        spike to be detected.  Default: 3.
    nreps: int, optional
        How many times to run the procedure.  Default is zero.
    interp.nan : bool, optional
        Whether missing values should be interpolated.  Interpolated values
        are computed after despiking.  Default is True.

    Returns
    -------
    VickersMahrt : tuple
        Tuple with (index, name in brackets):

        numpy.ndarray [0, `x`]
            1-D array with despiked input.
        numpy.int [1, `nspikes`]
            Number of spikes detected.
        numpy.int [2, `ntrends`]
            Number of outlier trends detected.
        numpy.int [3, `kclass`]
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
