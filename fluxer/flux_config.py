
"""Utility for configuring user variables for flux analyses."""

import configparser as cfg
from collections import OrderedDict
from os import getcwd
import os.path as osp
import glob
import re


__all__ = ["parse_config"]

# We need to control names.  We leave the date time columns (always the
# first two in input) up to pandas to determine, as they're always safe.
_incols_all = dict(longitude=float, latitude=float, speed_over_ground=float,
                   course_over_ground=float, heading=float, pitch=float,
                   roll=float, heave=float, atmospheric_pressure=float,
                   air_temperature=float, relative_humidity=float,
                   surface_temperature=float, wind_speed=float,
                   wind_direction=float, true_wind_speed=float,
                   true_wind_direction=float, PAR=float, K_down=float,
                   temperature_thermopile=float, temperature_case=float,
                   temperature_dome=float, LW_down=float,
                   UV_sensor_temperature=float, UV_a=float, UV_b=float,
                   UV_broad=float, acceleration_x=float, acceleration_y=float,
                   acceleration_z=float, rate_x=float, rate_y=float,
                   rate_z=float, wind_speed_u=float, wind_speed_v=float,
                   wind_speed_w=float, air_temperature_sonic=float,
                   sound_speed=float, anemometer_status=float,
                   op_CO2_fraction=float, op_CO2_density=float,
                   op_CO2_absorptance=float, op_H2O_fraction=float,
                   op_H2O_density=float, op_H2O_absorptance=float,
                   op_pressure=float, op_temperature=float,
                   op_temperature_base=float, op_temperature_spar=float,
                   op_temperature_bulb=float, op_cooler_voltage=float,
                   op_bandwidth=float, op_delay_interval=float,
                   op_bad_chopper_wheel_temperature_flag=bool,
                   op_bad_detector_temperature_flag=bool,
                   op_bad_optical_wheel_rate_flag=bool,
                   op_bad_sync_flag=bool, op_CO2_signal_strength=float,
                   op_AGC=float, op_analyzer_status=float,
                   cp_analyzer_status=float, cp_CO2_fraction=float,
                   cp_CO2_density=float, cp_CO2_dry_fraction=float,
                   cp_CO2_absorptance=float, cp_H2O_fraction=float,
                   cp_H2O_density=float, cp_H2O_dry_fraction=float,
                   cp_H2O_absorptance=float, cp_pressure=float,
                   cp_temperature=float, cp_temperature_in=float,
                   cp_temperature_out=float, cp_temperature_block=float,
                   cp_temperature_cell=float, cp_CO2_signal_strength=float,
                   cp_H2O_signal_strength=float, uw_CO2_cellA=float,
                   uw_CO2_cellB=float, uw_CO2_fraction=float,
                   uw_H2O_cellA=float, uw_H2O_cellB=float,
                   uw_H2O_fraction=float, uw_temperature_analyzer=float,
                   uw_pressure_analyzer=float, equ_temperature=float,
                   equ_pressure=float, H2O_flow=float,
                   air_flow_analyzer=float, equ_speed_pump=float,
                   ventilation_flow=float, condensation_atm=float,
                   condensation_equ=float, drip_1=float, drip_2=float,
                   condenser_temperature=float, ctd_pressure=float,
                   ctd_temperature=float, ctd_conductivity=float,
                   ctd_O2_saturation=float, ctd_O2_concentration=float,
                   uw_pH=float, uw_redox_potential=float,
                   temperature_external=float, temperature_in=float,
                   bad_wind_direction_flag=bool,
                   very_bad_wind_direction_flag=bool, bad_ice_flag=bool,
                   bad_atmospheric_pressure_flag=bool, bad_ctd_flag=bool,
                   bad_CO2_flag=bool, bad_H2O_flag=bool,
                   bad_H2O_flow_flag=bool, bad_pressure_analyzer_flag=bool,
                   bad_temperature_analyzer_flag=bool,
                   bad_equ_temperature_flag=bool,
                   bad_temperature_external_flag=bool)

# Simple dictionary to list our defaults
_dflts = {
    'ec_input_directory': getcwd(),
    'ec_file_pattern': "*[0-9].csv",
    'ec_colnames': ("time_20min,time_study,longitude,latitude," +
                    "speed_over_ground,course_over_ground,heading," +
                    "pitch,roll,heave,atmospheric_pressure," +
                    "air_temperature,relative_humidity," +
                    "surface_temperature,wind_speed,wind_direction," +
                    "true_wind_speed,true_wind_direction,PAR,K_down," +
                    "temperature_thermopile,temperature_case," +
                    "temperature_dome,LW_down,UV_sensor_temperature," +
                    "UV_b,UV_broad,acceleration_x,acceleration_y," +
                    "acceleration_z,rate_x,rate_y,rate_z,wind_speed_u," +
                    "wind_speed_v,wind_speed_w,air_temperature_sonic," +
                    "sound_speed,anemometer_status,op_CO2_fraction," +
                    "op_CO2_density,op_CO2_absorptance,op_H2O_fraction," +
                    "op_H2O_density,op_H2O_absorptance,op_pressure," +
                    "op_temperature,op_temperature_base," +
                    "op_temperature_spar,op_temperature_bulb," +
                    "op_cooler_voltage,op_bandwidth,op_delay_interval," +
                    "op_bad_chopper_wheel_temperature_flag," +
                    "op_bad_detector_temperature_flag," +
                    "op_bad_optical_wheel_rate_flag,op_bad_sync_flag," +
                    "op_CO2_signal_strength,cp_analyzer_status," +
                    "cp_CO2_fraction,cp_CO2_density,cp_CO2_dry_fraction," +
                    "cp_CO2_absorptance,cp_H2O_fraction,cp_H2O_density," +
                    "cp_H2O_dry_fraction,cp_H2O_absorptance,cp_pressure," +
                    "cp_temperature,cp_temperature_in,cp_temperature_out," +
                    "cp_temperature_block,cp_temperature_cell," +
                    "cp_CO2_signal_strength,cp_H2O_signal_strength"),
    'ec_sample_frequency': "10.0",
    'ec_summary_file': "fluxes.csv",
    'ec_despike_win_width': "3000",
    'ec_despike_step': "1500",
    'ec_despike_nreps': "20",
    'motion2anemometer_pos': "0.0, 0.0, 0.0",
    'ec_complementary_filter_period': "10.0",
    'ec_accel_highpass_cutoff': "20.0",
    'uw_input_directory': getcwd(),
    'uw_file_pattern': "*[0-9].csv",
    'uw_summary_file': "underway_pCO2.csv",
    'uw_colnames': ("time_30min,time_study,longitude,latitude," +
                    "speed_over_ground,course_over_ground,heading" +
                    "atmospheric_pressure,air_temperature," +
                    "relative_humidity,surface_temperature,wind_speed," +
                    "wind_direction,true_wind_speed,true_wind_direction," +
                    "PAR,K_down,LW_down,UV_b,UV_a,UV_broad," +
                    "air_temperature_sonic,cp_CO2_fraction," +
                    "cp_H2O_fraction,cp_pressure,cp_temperature," +
                    "op_CO2_density,op_H2O_density,op_pressure," +
                    "op_temperature,equ_temperature,equ_pressure," +
                    "equ_speed_pump,uw_CO2_fraction,uw_H2O_fraction," +
                    "uw_temperature_analyzer,uw_pressure_analyzer," +
                    "H2O_flow,air_flow_analyzer,ctd_pressure," +
                    "ctd_temperature,ctd_conductivity,ctd_O2_saturation," +
                    "ctd_O2_concentration,uw_pH,uw_redox_potential," +
                    "temperature_external,temperature_in," +
                    "bad_wind_direction_flag," +
                    "very_bad_wind_direction_flag,bad_ice_flag," +
                    "bad_atmospheric_pressure_flag,bad_ctd_flag," +
                    "bad_CO2_flag,bad_H2O_flag,bad_H2O_flow_flag," +
                    "bad_pressure_analyzer_flag," +
                    "bad_temperature_analyzer_flag," +
                    "bad_equ_temperature_flag," +
                    "bad_temperature_external_flag"),
    'uw_intake_depth': "5.0",
    'anemometer2d_height': "16.0"}
# Scalar option names
_scalar_opts = ["ec_sample_frequency", "ec_despike_win_width",
                "ec_despike_step", "ec_despike_nreps",
                "ec_complementary_filter_period", "ec_accel_highpass_cutoff",
                "anemometer2d_height", "uw_intake_depth"]
_vector_opts = ["motion2anemometer_pos"]

def parse_config(cfg_file):
    """Parse configuration file for essential variables during flux analysis.

    Parameters
    ----------
    cfg_file : string
        Path to configuration file.

    Returns
    -------
    Ordered dictionary with variables for each section.  A list of input
    files found, given the input directory and file pattern, is also
    appended to the 'Inputs' section.

    """
    # Ordered dictionary based on dflts to give to the parser
    dflt_dict = OrderedDict((
        ("EC Inputs",
         OrderedDict((
             ("input_directory", _dflts["ec_input_directory"]),
             ("file_pattern", _dflts["ec_file_pattern"]),
             ("colnames", _dflts["ec_colnames"]),
             ("sample_frequency", _dflts["ec_sample_frequency"]),
         ))
     ),
        ("EC Outputs",
         OrderedDict((
             ("summary_file", _dflts["ec_summary_file"]),
         ))
     ),
        ("EC Despiking",
         OrderedDict((
             ("despike_win_width", _dflts["ec_despike_win_width"]),
             ("despike_step", _dflts["ec_despike_step"]),
             ("despike_nreps", _dflts["ec_despike_nreps"]),
         ))
     ),
        ("EC Motion Correction",
         OrderedDict((
             ("motion2anemometer_pos",
              _dflts["motion2anemometer_pos"]),
             ("complementary_filter_period",
              _dflts["ec_complementary_filter_period"]),
             ("accel_highpass_cutoff",
              _dflts["ec_accel_highpass_cutoff"]),
         ))
     ),
        ("UW Inputs",
         OrderedDict((
             ("input_directory", _dflts["uw_input_directory"]),
             ("file_pattern", _dflts["uw_file_pattern"]),
             ("colnames", _dflts["uw_colnames"]),
             ("uw_intake_depth", _dflts["uw_intake_depth"]),
             ("anemometer2d_height", _dflts["anemometer2d_height"]),
         ))
     ),
        ("UW Outputs",
         OrderedDict((
             ("summary_file", _dflts["uw_summary_file"]),
         ))
     ),))

    # Set up the parser to interpolate across sections
    config = cfg.ConfigParser(interpolation=cfg.ExtendedInterpolation())
    config.read_dict(dflt_dict)     # set up our specific defaults
    config.read_file(open(cfg_file)) # replace defaults with what
                                     # we're given
    # Copy where we'll replace strings with other types
    py_dict = dflt_dict.copy()
    # Loop through all items and clean them up to generate our variables as
    # lists of strings
    for sec in config.sections():
        for opt in config.options(sec):
            opt_value = config.get(sec, opt)
            # Just replace and skip if we have the same as defaults
            if (dflt_dict[sec][opt] == opt_value):
                py_dict[sec][opt] = opt_value.split(",")
            else:
                # Otherwise move on and remove double quotes, newlines, and
                # spaces
                clean_opt = re.sub('["\n ]+', "", opt_value)
                config.set(sec, opt, clean_opt)
                # Replace values with lists, splitting on comma character,
                # on our local dictionary
                py_dict[sec][opt] = clean_opt.split(",")
            # Extract single elements, and convert to floats and arrays where
            # appropriate
            if (len(py_dict[sec][opt]) == 1):
                py_dict[sec][opt] = py_dict[sec][opt][0]
            if (opt in _scalar_opts):
                py_dict[sec][opt] = float(py_dict[sec][opt])
            elif (opt in _vector_opts):
                py_dict[sec][opt] = [float(x) for x in py_dict[sec][opt]]

    # Check input directories exist
    if (not osp.exists(py_dict["EC Inputs"]["input_directory"])):
        raise Exception("Input directory doesn't exist")
    if (not osp.exists(py_dict["UW Inputs"]["input_directory"])):
        raise Exception("Input directory doesn't exist")

    # Check if we have a valid sampling frequency
    if (py_dict["EC Inputs"]["sample_frequency"] <= 0):
        raise Exception("Sampling frequency must be greater than zero")

    # Check if we have a valid despiking parameters
    if (py_dict["EC Despiking"]["despike_win_width"] <= 0):
        raise Exception("Despiking window width must be greater than zero")
    if (py_dict["EC Despiking"]["despike_step"] <= 0):
        raise Exception("Despiking step size must be greater than zero")
    if (py_dict["EC Despiking"]["despike_nreps"] < 0):
        raise Exception("The number of despiking iterations cannot be negative")

    # Sort input file lists
    ec_input_files = glob.glob(osp.join(config["EC Inputs"]["input_directory"],
                                        config["EC Inputs"]["file_pattern"]))
    ec_input_files.sort()
    py_dict["EC Inputs"]["input_files"] = ec_input_files
    uw_input_files = glob.glob(osp.join(config["UW Inputs"]["input_directory"],
                                        config["UW Inputs"]["file_pattern"]))
    uw_input_files.sort()
    py_dict["UW Inputs"]["input_files"] = uw_input_files

    # Check if we have all legal names for header of input files (don't
    # care about first 2 time columns).
    illegal_names = ((set(py_dict["EC Inputs"]["colnames"][2:]) |
                      set(py_dict["UW Inputs"]["colnames"][2:])) -
                     set(_incols_all.keys()))
    if len(illegal_names) > 0:
        raise Exception("There are illegal column names in config file:")
        print illegal_names
    else:
        ec_dtypes = {key: _incols_all[key] for key in
                     py_dict["EC Inputs"]["colnames"][2:]}
        py_dict["EC Inputs"]["dtypes"] = ec_dtypes
        uw_dtypes = {key: _incols_all[key] for key in
                     py_dict["UW Inputs"]["colnames"][2:]}
        py_dict["UW Inputs"]["dtypes"] = uw_dtypes

    return py_dict
