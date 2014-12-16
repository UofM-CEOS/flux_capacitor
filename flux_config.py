#! /usr/bin/env python
# $Id$

"""Utility for configuring user variables for flux analyses."""

import configparser as cfg
from collections import OrderedDict
from os import getcwd
import os.path as osp
import re

# Simple dictionary to list our defaults
dflts = {
    "input_directory": getcwd(),
    "file_pattern": '*[0-9].csv',
    "colnames": 'time_20min,time_study,longitude,latitude,speed_over_ground,course_over_ground,heading,pitch,roll,heave,atmospheric_pressure,air_temperature,relative_humidity,surface_temperature,wind_speed,wind_direction,true_wind_speed,true_wind_direction,PAR,K_down,temperature_thermopile,temperature_case,temperature_dome,LW_down,UV_sensor_temperature,UV_b,UV_broad,acceleration_x,acceleration_y,acceleration_z,rate_x,rate_y,rate_z,wind_speed_u,wind_speed_v,wind_speed_w,air_temperature_sonic,sound_speed,anemometer_status,op_CO2_fraction,op_CO2_density,op_CO2_absorptance,op_H2O_fraction,op_H2O_density,op_H2O_absorptance,op_pressure,op_temperature,op_temperature_base,op_temperature_spar,op_temperature_bulb,op_cooler_voltage,op_bandwidth,op_delay_interval,bad_chopper_wheel_temperature_flag,bad_detector_temperature_flag,bad_optical_wheel_rotation_rate_flag,bad_sync_flag,op_AGC,cp_analyzer_status,cp_CO2_fraction,cp_CO2_density,cp_CO2_dry_fraction,cp_CO2_absorptance,cp_H2O_fraction,cp_H2O_density,cp_H2O_dry_fraction,cp_H2O_absorptance,cp_pressure,cp_temperature,cp_temperature_in,cp_temperature_out,cp_temperature_block,cp_temperature_cell,cp_CO2_signal_strength,cp_H2O_signal_strength',
    "sample_frequency": '10.0',
    "summary_file": 'fluxes.csv',
    "motion2anemometer_pos": '0.0, 0.0, 0.0',
    "complementary_filter_period": '10.0',
    "accel_highpass_cutoff": '20.0'}
# Scalar option names
scalar_opts = ['sample_frequency', 'complementary_filter_period',
               'accel_highpass_cutoff']
vector_opts = ['motion2anemometer_pos']

def parse_config(cfg_file):
    """Parse configuration file for essential variables for flux analysis.

    Parameters
    ----------
    cfg_file : string
        Path to configuration file.

    Returns
    -------
    Ordered dictionary with variables for each section.

    """
    # Ordered dictionary based on dflts to give to the parser
    dflt_dict = OrderedDict((
        ('Inputs',
         OrderedDict((
             ('input_directory', dflts['input_directory']),
             ('file_pattern', dflts['file_pattern']),
             ('colnames', dflts['colnames']),
             ('sample_frequency', dflts['sample_frequency']),
         ))
     ),
        ('Outputs',
         OrderedDict((
             ('summary_file', dflts['summary_file']),
         ))
     ),
        ('Motion Correction',
         OrderedDict((
             ('motion2anemometer_pos',
              dflts['motion2anemometer_pos']),
             ('complementary_filter_period',
              dflts['complementary_filter_period']),
             ('accel_highpass_cutoff',
              dflts['accel_highpass_cutoff']),
         ))
     ),))

    # Set up the parser to interpolate across sections
    config = cfg.ConfigParser(interpolation=cfg.ExtendedInterpolation())
    config.read_dict(dflt_dict)     # set up our specific defaults
    config.read_file(open(cfg_file)) # replace defaults with what
                                     # we're given
    # Copy where we'll replace strings with other types
    py_dict = dflt_dict.copy()
    # Loop through all items and clean them to generate our variables as
    # lists of strings
    for sec in config.sections():
        for opt in config.options(sec):
            opt_value = config.get(sec, opt)
            # Just replace and skip if we have the same as defaults
            if (dflt_dict[sec][opt] == opt_value):
                py_dict[sec][opt] = opt_value.split(',')
                continue
            # Otherwise move on and remove double quotes, newlines, and
            # spaces
            clean_opt = re.sub('["\n ]+', '', opt_value)
            config.set(sec, opt, clean_opt)
            # Replace values with lists, splitting on comma character, on
            # our local dictionary
            py_dict[sec][opt] = clean_opt.split(',')

    # Loop again to extract single elements, and convert to floats and
    # arrays where needed
    for sec in config.sections():
        for opt in config.options(sec):
            # If we have a single element list, remove the container list
            if (len(py_dict[sec][opt]) == 1):
                py_dict[sec][opt] = py_dict[sec][opt][0]
            if (opt in scalar_opts):
                py_dict[sec][opt] = float(py_dict[sec][opt])
            elif (opt in vector_opts):
                py_dict[sec][opt] = [float(x) for x in py_dict[sec][opt]]

    # Check input directory exists
    if (not osp.exists(py_dict['Inputs']['input_directory'])):
        raise Exception("Input directory doesn't exist")

    # Check if we have a valid sampling frequency
    if (py_dict['Inputs']['sample_frequency'] <= 0):
        raise Exception("Sampling frequency must be greater than zero")

    # # Check if directory for summary file exists
    # if (not osp.exists(osp.dirname(py_dict['Outputs']['summary_file']))):
    #     raise Exception("Directory for summary file doesn't exist")

    return py_dict
