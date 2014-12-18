# $Id$

import numpy as np
from scipy import interpolate as itpl
from scipy import signal
from astropy.convolution import convolve, Box1DKernel


def shot_filter(x, sigma_thr=3):
    """Perform spline interpolation for extreme values in Series `x`.
    
    Extreme values are those that are larger than `sigma_thr` * `sigma`,
    where `sigma` is the sample standard deviation of Pandas Series `x`.

    Interpolation uses the time index of `x`.  The interpolating function
    is created by the InterpolatedUnivariateSpline function from the scipy
    package, and uses a single knot to approximate a simple linear
    interpolation, so as to keep the original signal as untouched as
    possible.

    Returns:
      A Pandas Series as `x` with interpolated values at index location of
      extremes.

    """
    x_new = x.copy().astype('d')
    is_ok = abs(x_new - np.mean(x_new)) < (sigma_thr * np.std(x_new))
    x_ok = x_new[is_ok]
    t_ok = x.index[is_ok]
    t_bad = x.index[~is_ok]
    # # Simple linear interpolation; extrapolation impossible
    # f_itpl = itpl.interp1d(t_ok.values.astype('d'), x_ok)
    # x_itpl = f_itpl(t_bad.values.astype('d'))
    # Trying a 1D B-spline for more realistic intra- and extrapolation.
    # This needs more work, and may not be worth it...
    s = itpl.InterpolatedUnivariateSpline(t_ok.values.astype('d'),
                                          x_ok, k=1)
    x_itpl = s(t_bad.values.astype('d'))
    x_new[t_bad] = x_itpl
    return x_new


def decompose(angle, vmagnitude):
    """Decompose angle and magnitude into `x` and `y` vector(s).

    Parameters
    ----------
    angle : array_like
            The angle(s) in degree units.
    vmagnitude : array_like
                 The magnitude array associated with each angle.  It can be
                 a scalar that is common to all `angle` elements.

    Returns
    -------
    Tuple with ndarrays `x` and `y`, in that order.

    """
    x = vmagnitude * np.cos(np.radians(angle))
    y = vmagnitude * np.sin(np.radians(angle))
    return x, y


def recompose(x, y):
    """Recompose angles and associated magnitudes from `x` and `y` vectors.

    Parameters
    ----------
    x : array_like
        `x`-coordinates
    y : array_like
        `y`-coordinates

    Returns
    -------
    Tuple with ndarrays `angle` and `vmagnitude`, in that order.

    """
    vmag = np.sqrt((x ** 2) + (y ** 2))
    ang = np.arctan2(y, x)
    ang[ang < 0] = ang[ang < 0] + (2 * np.pi) # output range 0 - 2*pi
    ang[vmag == 0] = 0          # when magnitude is 0 the angle is also 0
    ang[ang == 0] = 2 * np.pi   # convention
    return np.degrees(ang), vmag


def smooth_angle(angle, vmagnitude=1, kernel_width=21):
    """Smooth angles by decomposing them, applying a boxcar average.

    Smoothing is done by using a 1D kernel smoothing filter of a given
    width.

    Parameters
    ----------
    angle : numpy.ndarray
            The angle(s) in degree units.
    vmagnitude : numpy.ndarray, optional
                 The magnitude array associated with each angle.  It can be
                 a scalar that is common to all `angle` elements.
    kernel_width : int
                   The width of the filter kernel.

    Returns
    -------
    Tuple with ndarrays `angle` and `vmagnitude`, in that order.

    """
    x, y = decompose(angle, vmagnitude)
    x_smooth = convolve(x, Box1DKernel(kernel_width), boundary="extend")
    y_smooth = convolve(y, Box1DKernel(kernel_width), boundary="extend")
    return recompose(x_smooth, y_smooth)


def level3D_motion(accel, ang_rate, roll_range, pitch_range):
    """Level 3D acceleration and angular rate measured by motion sensor.

    """
    pass                        # IMPLEMENT THIS?


def level3D_anemometer(wind_speed, roll, pitch):
    """Level 3D anemometer measurements, given mean roll and pitch.

    """
    pass                        # IMPLEMENT THIS?


def euler_translate(X, euler):
    """Translate vector matrix using a translation matrix."""
    x, y, z = X[:, 0], X[:, 1], X[:, 2]
    pitch, tilt, s = euler[:, 0], euler[:, 1], euler[:, 2]
    x_new =  (x * np.cos(tilt) * np.cos(s) + 
              y * ( np.sin(pitch) * np.sin(tilt) * np.cos(s) -
                    np.cos(pitch) * np.sin(s)) +
              z * ( np.cos(pitch) * np.sin(tilt) * np.cos(s) +
                    np.sin(pitch) * np.sin(s)))
    y_new =  (x * np.cos(tilt) * np.sin(s) +
              y * ( np.sin(pitch) * np.sin(tilt) * np.sin(s) +
                    np.cos(pitch) * np.cos(s)) +
              z * ( np.cos(pitch) * np.sin(tilt) * np.sin(s) -
                    np.sin(pitch) * np.cos(s)))
    z_new =  (x * (-np.sin(tilt)) + y * ( np.cos(tilt) * np.sin(pitch)) +
              z * ( np.cos(tilt) * np.cos(pitch)))
    return np.column_stack((x_new, y_new, z_new))
    

def wind3D_correct(wind_speed, acceleration, angle_rate, heading, speed,
                   anemometer_pos, sample_freq, Tcf, Ta,
                   tilt_motion, tilt_anemometer):
    """Correct wind vector measurements from a moving platform.

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
        A vector with ship heading measurements (deg).
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
    Tuple with (index in brackets):
    numpy.ndarray [0]
        2-D array (Nx3) with corrected wind vectors in ship-referenced
        frame with z-axis parallel to gravity.
    numpy.ndarray [1]
        2-D array (Nx3) with Euler angles.
    numpy.ndarray [2]
        2-D array (NX3) with Euler angles from rate sensors (unfiltered).
    numpy.ndarray [3]
        2-D array (NX3) with Euler angles from accelerometers (unfiltered).
    numpy.ndarray [4]
        2-D array (NX3) with slow Euler angles (low pass filtered).
    numpy.ndarray [5]
        2-D array (NX3) with fast Euler angles (high pass filtered).
    numpy.ndarray [6]
        2-D array (3X3) with mounting offset rotation matrix (see Miller 2008).
    numpy.ndarray [7]
        2-D array (NX3) with measured velocity, rotated to the earth frame.
    numpy.ndarray [8]
        2-D array (NX3) with velocity induced by angular motion.
    numpy.ndarray [9]
        2-D array (NX3) with velocity induced by platform linear motion.
    numpy.ndarray [10]
        2-D array (NX3) with ship velocity in eastward, northward, up
        frame.  (low pass filtered).
    numpy.ndarray [11]
        2-D array (NX3) with corrected wind in eastward, northward, up
        frame.
    numpy.ndarray [12]
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
    p_lf = np.arctan2(ay, g) - \
           signal.filtfilt(bc, ac, np.arctan2(ay, g), padlen=pdl)
    # High frequency angles using angle rates
    rm = signal.detrend(angle_rate, 0)
    rx, ry, rz = rm[:, 0], rm[:, 1], rm[:, 2] # don't do this in ipdb!
    p_hf = signal.filtfilt(bc, ac,
                           ((np.cumsum(rx) - 0.5 * rx - 0.5 * rx[0]) /
                            sample_freq), padlen=pdl)
    pitch = p_lf + p_hf
    t_lf = np.arctan2(-ax * np.cos(pitch), g) - \
           signal.filtfilt(bc, ac, np.arctan2(-ax * np.cos(pitch), g),
                           padlen=pdl)

    gyro = np.unwrap(-np.radians(heading)) # Put heading in radians and RHS
    heading_rhs = gyro[0]               # Initial heading
    gyro = gyro - heading_rhs           # Remove the intial heading
    # Low-pass filter to retain low-frequency content, but not the offset
    s_lf = gyro - signal.filtfilt(bc, ac, gyro, padlen=pdl)

    t_hf = signal.filtfilt(bc, ac,
                           ((np.cumsum(ry) - 0.5 * ry - 0.5 * ry[0]) /
                            sample_freq), padlen=pdl)
    s_hf = signal.filtfilt(bc, ac,
                           ((np.cumsum(rz) - 0.5 * rz - 0.5 * rz[0]) /
                            sample_freq), padlen=pdl)
    tilt = t_lf + t_hf
    # s = s_lf + s_hf             # SPL: not used...

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
    F12 = np.tan(tilt) * np.sin(pitch)
    F13 = np.tan(tilt) * np.cos(pitch)
    F22 = np.cos(pitch)
    F23 = -np.sin(pitch)
    F32 = np.sin(pitch) / np.cos(tilt)
    F33 = np.cos(pitch) / np.cos(tilt)

    pdot = rx + ry * F12 + rz * F13
    tdot =      ry * F22 - rz * F23
    sdot =      ry * F32 + rz * F33

    # The angular rates pdot, tdot, and sdot are now in the earth based
    # frame, which is what we want for creating the rotation-transformation
    # matrix. Re-integrate and high pass filter these angular rates to get
    # the fast angles.
    p_hf = signal.filtfilt(bc, ac,
                           (np.cumsum(pdot) - 0.5 * pdot - 0.5 * pdot[0]) /
                           sample_freq, padlen=pdl)
    t_hf = signal.filtfilt(bc, ac,
                           (np.cumsum(tdot) - 0.5 * tdot - 0.5 * tdot[0]) /
                           sample_freq, padlen=pdl)
    s_hf = signal.filtfilt(bc, ac,
                           (np.cumsum(sdot) - 0.5 * sdot - 0.5 * sdot[0]) /
                           sample_freq, padlen=pdl)

    # combine high- and low-frequency angles; make 2D matrix
    EA = np.column_stack((p_lf+p_hf, t_lf+t_hf, s_lf+s_hf))

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
    Ur = euler_translate(np.dot(M_ma, wind_speed.T).T, EA)

    ## PLATFORM ANGULAR VELOCITY
    uam = np.column_stack((anemometer_pos[2] * ry -
                           anemometer_pos[1] * rz,
                           anemometer_pos[0] * rz -
                           anemometer_pos[2] * rx,
                           anemometer_pos[1] * rx -
                           anemometer_pos[0] * ry))
    Ua = euler_translate(uam, EA) # rotate to earth frame

    ## PLATFORM LINEAR VELOCITY
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
    ae = euler_translate(acceleration, EA)
    ae[:, 2] = ae[:, 2] - g     # subtract gravity
    n = np.size(ae, 0)
    up = (np.cumsum(ae, 0) -
          0.5 * ae - 0.5 * (np.ones((n, 1)) * ae[0, :])) / sample_freq
    Up = np.column_stack((signal.filtfilt(bax, aax, up[:, 0], padlen=pdlx),
                          signal.filtfilt(bay, aay, up[:, 1], padlen=pdly),
                          signal.filtfilt(baz, aaz, up[:, 2], padlen=pdlz)))

    ## SHIP SPEED: this is in the positive x-direction (assuming the bow is
    ## positive x direction).  We low-pass filter this signal so there isnt
    ## double counting of linear motion from the integrated accelerometers.
    ## Also, the y- and z-components are zero (ie, captured by the
    ## integrated accelerometers).
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
    EA_slow = np.column_stack((p_lf, t_hf, s_hf))
    EA_fast = np.column_stack((p_hf, t_hf, s_hf))
    u_ea = np.column_stack((np.zeros((n, 1)),
                            np.zeros((n, 1)),
                            (heading_rhs + np.pi / 2) * np.ones((n, 1))))
    U_ship = euler_translate(Us, u_ea)
    U_earth = euler_translate(UVW, u_ea)
    # Platform displacement
    xp = (np.cumsum(Up, 0) - 0.5 * Up -
          0.5 * (np.ones((n, 1)) * Up[0, :])) / sample_freq
    Xp = np.column_stack((signal.filtfilt(bax, aax, xp[:, 0], padlen=pdlx),
                          signal.filtfilt(bay, aay, xp[:, 1], padlen=pdly),
                          signal.filtfilt(baz, aaz, xp[:, 2], padlen=pdlz)))

    return (UVW, EA, EA_rate, EA_acc, EA_slow, EA_fast, M_ma, Ur, Ua, Up,
            U_ship, U_earth, Xp)