
function [u, varargout] = motion_octave(um, am, rm, gyro, speed, r, sf, Tcf,
					Ta, tilt_ms, tilt_anem, OutputVars)
# =========================================================================
# MOTION.M - corrects wind vector measurements from a moving platform. The
#     method is based on a 3-dimensional wind sensor, and a 6 degree of
#     freedom motion sensor that measures 3 orthogonal accelerations and 3
#     orthogonal angle rates.
#
#    Copyright (C) Scott Miller 1995-2008, 2015.
#    Version: 2.0
#    Last Revised: Dec 4, 2008
#
# REFERENCE:
#
# Miller,S.D., Hristov,T.S., Edson,J.B., and C.A. Friehe, 2008: Platform
#     Motion Effects on Measurements of Turbulence and Air-Sea Exchange
#     Over the Open Ocean, J. Atmo. Ocean. Tech. 25(9), 1683-1694,
#     DOI: 10.1175/2008JTECHO547.1.
#
# AUTHOR: Scott Miller (smiller@albany.edu)
#          Atmospheric Sciences Research Center
#          State University of New York at Albany
#
# DISCLAIMER: This software is provided "as is" and is intended for 
#     research purposes without warranty of any kind.
#
# DETAILS:
#
#   This scripts makes use of the filter design functions in the Matlab
#   Signal Processing Toolbox, so will not work without that toolbox unless
#   the filters are specified explicitly in this code.
#
#   Coordinate Frame: (right handed system)
#     x - positive to the bow
#     y - positive to port side
#     z - positive up
#
# INPUTS:
# -------
#   uvwm    - NX3 measured (u,v,w) velocity vectors (m/s)
#   am      - NX3 measured (ax,ay,az) accelerations (m/s/s)
#   rm      - NX3 measured (rx,ry,rz) angle rates (radians/sec)
#   gyro    - ship heading (heading, deg)
#   speed   - ship speed (m/s)
#   r       - position vector of sonic w.r.t. motion sensor (m)
#   sf      - sampling frequency (Hz)
#   Tcf     - complimentary filter period (seconds)
#   Ta      - high-pass filter cutoff for accelerations
#             (seconds). This can be a scalar if the same cutoff is
#             used for the three components, or a 3 element vector
#             to indicate cutoff period for the [x y z] components
#   tilt_ms - [roll pitch]: mean motion sensor tilt in radians with
#             respect to a horizontal plane. 'roll' angle offset is
#             positive port side up, 'pitch' is positive for bow down
#             (see Miller 2008 for info, set to [0 0] if unknown)
#   tilt_anem - [roll pitch]: mean anemometer tilt in radians with
#             respect to a horizontal plane. 'roll' angle offset
#             is positive port side up, 'pitch' is positive for bow down
#             (see Miller 2008 for info, set to [0 0] if unknown)
#
# OUTPUTS:
# --------
# uvw - corrected wind vector in ship-referenced frame with z-axis parallel
#       to gravity
#
#   Optional outputs according to the optional input list OutputVars
#   The same number and order of output arguments should be on the
#   left hand sign of the function call e.g., OutputVars={'EA';'UEARTH'}
#   Options are:
#     'ea'      - (NX3) the euler angle array
#     'ea_rate' - (NX3) euler angles from rate sensors (unfiltered)
#     'ea_acc'  - (NX3) euler angles from accelerometers (unfiltered)
#     'ea_slow' - (NX3) slow euler angles (low pass filtered)
#     'ea_fast' - (NX3) fast euler angles (high pass filtered)
#     'M'       - (3X3) mounting offset rotation matrix (see Miller 2008)
#     'urotation' - (NX3) the measured velocity rotated to the earth frame
#     'uangular'  - (NX3) velocity induced by angular motion
#     'ulinear'   - (NX3) velocity induced by platform linear motion
#     'uship'     - (NX3) ship velocity in eastward, northward, up frame (low pass filtered)
#     'uearth'    - (NX3) corrected wind in eastward, northward, up frame.
#     'xyz'       - (NX3) displacement of the motion sensor
#
# Here is an example of the call to this function from one of my scripts:
#
# [uvw,xyz,uvwe,uship] = motion_octave(uvwm,acc,rate,heading,speed,[1 -.5 0.25],...
#                         10,10,20,[0 0],[0 0],{'xyz';'uearth';'uship'});
#
  
  if length(Ta) == 1
    Ta = Ta * ones(3, 1);
  end

  ## Define complimentary filter
  Astop = 10; Apass = 0.5;    # stopband attenuation; passband Ripple (dB)
  wp = 1 / (2 * Tcf) / (sf / 2); ws = 1 / Tcf / (sf / 2);  # stop, passband cutoffs
  [N, Wn] = buttord(wp, ws, Apass, Astop);

  ## N=5; Wn=0.0079; # commented out for testing mar 29, 06
  [bc, ac]=butter(N, Wn, 'high');

  ## EULER ANGLES: (see Edson et al., 1998)
  ## low frequency tilt using accelerometers
  g = sqrt(sum(mean(am) .^ 2));
  ax = am(:, 1); ay = am(:, 2); az = am(:, 3);
  p_lf = atan2(ay, g) - filtfilt(bc, ac, atan2(ay, g));
  ## high frequency angles using angle rates.
  rm = detrend(rm, 0);
  rx = rm(:, 1); ry = rm(:, 2); rz =rm(:, 3);
  p_hf = filtfilt(bc, ac, ((cumsum(rx) - 0.5 * rx - 0.5 * rx(1)) / sf));
  p = reshape(p_lf, length(p_lf), 1) + reshape(p_hf, length(p_hf), 1);
  t_lf = atan2(-ax .* cos(p), g) - ...
	 filtfilt(bc, ac, atan2(-ax .* cos(p), g));

  gyro=unwrap(deg2rad(-gyro)); # put gyro into radians and rhs
  headingrhs=gyro(1);        # get initial heading (changed from mean
			     # heading, Oct 2006)
  gyro=gyro-headingrhs;      # remove the intial heading
  # low-pass filter to retain low-frequency content, but not the offset
  s_lf=(gyro-filtfilt(bc, ac, gyro));

  t_hf=filtfilt(bc, ac, ((cumsum(ry) - 0.5 * ry - 0.5 * ry(1)) / sf));
  s_hf=filtfilt(bc, ac, ((cumsum(rz) - 0.5 * rz - 0.5 * rz(1)) / sf));
  t=reshape(t_lf, length(t_lf), 1) + reshape(t_hf, length(t_hf), 1);
  s=reshape(s_lf, length(s_lf), 1) + reshape(s_hf, length(s_hf), 1);

  ## NONLINEAR ROTATION OF ANGULAR RATES FROM MOTION SENSOR-FRAME TO
  ## EARTH-FRAME This next step puts the measured angular rates into the
  ## earth reference frame.  Motion sensor frame-based anglular rates are
  ## first rolled based on the roll angle (phim) estimated above; they are
  ## then pitched *in a rolled frame* based upon the pitch angle (thetam)
  ## estimated above; finally, they are yawed *in a rolled and pitched
  ## frame* based upon both the pitch and roll (thetam and phim)
  ## estimates. These sequential rotations are exact for small angle RATES
  ## (see Goldstein) - I've found this step doesn't have much effect on the
  ## corrected winds.

  F12 = tan(t) .* sin(p);
  F13 = tan(t) .* cos(p);
  F22 = cos(p);
  F23 = -sin(p);
  F32 = sin(p) ./ cos(t);
  F33 = cos(p) ./ cos(t);
  
  pdot = rx + ry .* F12 + rz .* F13;
  tdot =      ry .* F22 - rz .* F23;
  sdot =      ry .* F32 + rz .* F33;
  
  ## The angular rates pdot, tdot, and sdot are now in the earth based
  ## frame, which is what we want for creating the rotation-transformation
  ## matrix. Re-integrate and high pass filter these angular rates to get
  ## the fast angles.

  p_hf = filtfilt(bc, ac,
		  (cumsum(pdot) - 0.5 * pdot - 0.5 * pdot(1)) / sf);
  t_hf = filtfilt(bc, ac,
		  (cumsum(tdot) - 0.5 * tdot - 0.5 * tdot(1)) / sf);
  s_hf = filtfilt(bc, ac,
		  (cumsum(sdot) - 0.5 * sdot - 0.5 * sdot(1)) / sf);
  
  ## combine high- and low-frequency angles
  p_freqs = reshape(p_lf, length(p_lf), 1) + reshape(p_hf, length(p_hf), 1);
  t_freqs = reshape(t_lf, length(t_lf), 1) + reshape(t_hf, length(t_hf), 1);
  s_freqs = reshape(s_lf, length(s_lf), 1) + reshape(s_hf, length(s_hf), 1);
  EA = [p_freqs t_freqs s_freqs];

  ## Ship AND ANEMOMETER OFFSET ANGLES pitch and roll tilt of platform
  ## w.r.t. earth (see Wilczak et al. 2000)

  phi_em = tilt_ms(1);
  theta_em = tilt_ms(2);
  PHI_em = [1 0 0; 0 cos(phi_em) -sin(phi_em); 0 sin(phi_em) cos(phi_em)];
  THETA_em = [cos(theta_em) 0 sin(theta_em); 0 1 0; ...
	      -sin(theta_em) 0 cos(theta_em)];

  ## pitch and roll tilt of anemometer w.r.t. earth
  phi_ea = tilt_anem(1);
  theta_ea = tilt_anem(2);
  PHI_ea = [1 0 0; 0 cos(phi_ea) -sin(phi_ea); 0 sin(phi_ea) cos(phi_ea)];
  THETA_ea = [cos(theta_ea) 0 sin(theta_ea); 0 1 0; ...
	      -sin(theta_ea) 0 cos(theta_ea)];

  ## anemometer w.r.t. motion sensor
  M_em = THETA_em * PHI_em;
  M_ea = THETA_ea * PHI_ea;
  M_ma = transpose(M_em) * M_ea;

  ## WIND VECTOR ROTATION
  ## using coordinate rotations
  [ur] = trans([M_ma * um']', EA);

  ## PLATFORM ANGULAR VELOCITY
  uam(:, 1) = r(3) * rm(:, 2) - r(2) * rm(:, 3);
  uam(:, 2) = r(1) * rm(:, 3) - r(3) * rm(:, 1);
  uam(:, 3) = r(2) * rm(:, 1) - r(1) * rm(:, 2);
  [ua] = trans(uam, EA);  # rotate to earth frame

  ## PLATFORM LINEAR VELOCITY
  ## high-pass filters for accelerometers
  Astop = 10; Apass = 0.5;    # stopband attenuation; passband Ripple (dB)

  wp = 1 / (2 * Ta(1)) / (sf / 2); ws = 1 / Ta(1) / (sf / 2);  # stop, passband cutoffs
  [n, wn] = buttord(wp, ws, Apass, Astop);
  [bax, aax] = butter(n, wn, 'high');

  wp = 1 / (2 * Ta(2)) / (sf / 2); ws = 1 / Ta(2) / (sf / 2);  # stop, passband cutoffs
  [n, wn] = buttord(wp, ws, Apass, Astop);
  [bay, aay] = butter(n, wn, 'high');

  wp = 1 / (2 * Ta(3)) / (sf / 2); ws = 1 / Ta(3) / (sf / 2);  # stop, passband cutoffs
  [n, wn] = buttord(wp, ws, Apass, Astop);
  [baz, aaz] = butter(n, wn, 'high');

  ##rotate accelerations
  [ae] = trans(am, EA);
  ae(:, 3) = ae(:, 3) - g;   # subtract gravity
  n = size(ae, 1);
  up = (cumsum(ae) - 0.5 * ae - 0.5 * (ones(n, 1) * ae(1, :))) / sf;
  up(:, 1) = filtfilt(bax, aax, up(:, 1));  # hp filter
  up(:, 2) = filtfilt(bay, aay, up(:, 2));
  up(:, 3) = filtfilt(baz, aaz, up(:, 3));

  ## SHIP SPEED: this is in the positive x-direction (assuming the bow is
  ## positive x direction).  We low-pass filter this signal so there isnt
  ## double counting of linear motion from the integrated
  ## accelerometers. Also, the y- and z-components are zero (ie, captured
  ## by the integrated accelerometers).
  n = length(speed);
  us = [speed - filtfilt(bax, aax, speed) zeros(n, 1) zeros(n, 1)];

  u = [ur + ua + up + us];   # corrected wind vector

  ##process optional outputs

  if exist('OutputVars','var')
    for i=1:size(OutputVars, 1)
      switch upper((OutputVars{i}))
        case {'EA'}
          varargout(i) = {EA};
        case {'M'}
          varargout(i) = {M_ma};
        case {'EA_SLOW'}
          varargout(i) = {[p_lf t_lf s_lf]};
        case {'EA_FAST'}
          varargout(i) = {[p_hf t_hf s_hf]};
        case {'EA_ACC'}
          varargout(i) = {[atan2(ay, g) atan2(-ax, g) gyro]};
        case {'EA_RATE'}
	  x_col = (cumsum(rx) - 0.5 * rx - 0.5 * rx(1)) / sf;
	  y_col = (cumsum(ry) - 0.5 * ry - 0.5 * ry(1)) / sf;
	  z_col = (cumsum(rz) - 0.5 * rz - 0.5 * rz(1)) / sf;
          varargout(i)={[x_col y_col z_col]};
        case {'UROTATION'}
          varargout(i) = {ur};
        case {'UANGULAR'}
          varargout(i) = {ua};
        case {'ULINEAR'}
          varargout(i) = {up};
        case {'USHIP'}
	  hdgrhs = (headingrhs + pi / 2) * ones(n, 1);
          us_e = trans(us, [zeros(n, 1) zeros(n,1) hdgrhs]);
          varargout(i)={us_e};
        case {'UEARTH'}
	  ## rotate to earth frame, where x is east, y is north
	  hdgrhs = (headingrhs + pi / 2) * ones(n, 1);
          ue = trans(u, [zeros(n, 1) zeros(n, 1) hdgrhs]);
          varargout(i) = {ue};
        case {'XYZ'}
	  ## PLATFORM DISPLACEMENT
          xp = (cumsum(up) - 0.5 * up -
		0.5 * (ones(n, 1) * up(1, :))) / sf;
          xp(:, 1) = filtfilt(bax, aax, xp(:, 1));  # hp filter
          xp(:, 2) = filtfilt(bay, aay, xp(:, 2));
          xp(:, 3) = filtfilt(baz, aaz, xp(:, 3));
          varargout(i) = {xp};
        otherwise
          disp('Unknown output requested (motion.m)');
      end
    end
  end
  return
endfunction

function [xout] = trans(xin, euler)
  x = xin(:, 1); y = xin(:, 2); z = xin(:, 3);
  p = euler(:, 1); t = euler(:, 2); s = euler(:, 3);
  xout = x .* cos(t) .* cos(s) + ...
	 y .* (sin(p) .* sin(t) .* cos(s) - cos(p) .* sin(s)) + ...
	 z .* (cos(p) .* sin(t) .* cos(s) + sin(p) .* sin(s));
  yout =  x .* cos(t) .* sin(s) + ...
	  y .* (sin(p) .* sin(t) .* sin(s) + cos(p) .* cos(s)) + ...
	  z .* (cos(p) .* sin(t) .* sin(s) - sin(p) .* cos(s));
  zout = x .* (-sin(t)) + ...
	 y .* (cos(t) .* sin(p)) + ...
	 z .* (cos(t) .* cos(p));
  xout = [xout yout zout];
  return
endfunction
