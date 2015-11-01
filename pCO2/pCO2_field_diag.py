"""Diagnostic plots for underway data in the field.

We take the output of the AWK pre-processing scripts that gather all files
together, and then here we examine variables.

We have the following file structure (2015) for pCO2:

[1]  time [YYYY-MM-DD HH:MM:SS]
[2]  record_type [string]
[3]  uw_diag [integer]
[4]  equ temperature [D+]
[5]  std value [D+]
[6]  "uw CO2 cell A (raw value)" [D+]
[7]  "uw CO2 cell B (raw value)" [D+]
[8]  "uw CO2 fraction (um/m)" [D+]
[9]  "uw H2O cell A (raw value)" [D+]
[10] "uw H2O cell B (raw value)" [D+]
[11] "uw H2O fraction (mm/m)" [D+]
[12] uw temperature analyzer [D+]
[13] uw pressure analyzer [D+]
[14] equ pressure [D+]
[15] "H2O flow" [D+]
[16] air flow analyzer [D+]
[17] equ speed pump [D+]
[18] ventilation flow [D+]
[19] condensation_atm [D+]
[20] condensation_equ [D+]
[21] drip 1 [D+]
[22] drip 2 [D+]
[23] condenser temperature [D+]
[24] temperature dry box [D+]
[25] ctd pressure [D+]
[26] ctd temperature [D+]
[27] ctd conductivity [D+]
[28] "ctd O2 saturation" [D+]
[29] "ctd O2 concentration" [D+]
[30] "uw pH" [D+]
[31] uw redox potential [D+]
[32] temperature external [D+]

And for external temperature data files:

[1] time [YYYY-MM-DD HH:MM:SS]
[2] battery voltage [D+]
[3] logger_temperature [D+]
[4] water_temperature [D+]
[5] battery_voltage_sd [D+]
[6] logger_temperature_sd [D+]
[7] water_temperature_sd [D+]

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('ggplot')

pCO2_names = ["time", "record_type", "uw_diag", "equ_temperature",
              "std_value", "uw_CO2_cellA", "uw_CO2_cellB",
              "uw_CO2_fraction", "uw_H2O_cellA", "uw_H2O_cellB",
              "uw_H2O_fraction", "uw_temperature_analyzer",
              "uw_pressure_analyzer", "equ_pressure", "H2O_flow",
              "air_flow_analyzer", "equ_speed_pump", "ventilation_flow",
              "condensation_atm", "condensation_equ", "drip_1", "drip_2",
              "condenser_temperature", "temperature_dry_box",
              "ctd_pressure", "ctd_temperature", "ctd_conductivity",
              "ctd_O2_saturation", "ctd_O2_concentration", "uw_pH",
              "uw_redox_potential", "temperature_external"]
tsw_names = ["time", "battery_voltage", "logger_temperature",
             "water_temperature", "battery_voltage_sd",
             "logger_temperature_sd", "water_temperature_sd"]

# CSV file
pCO2 = pd.read_csv("/media/sluque/Data_2015/Data/ArcticNet/2015/UW_pCO2/"
                   "UW_pCO2.csv", header=None, names=pCO2_names,
                   parse_dates=["time"], index_col="time")
# pCO2["record_type"] = pCO2["record_type"].astype("category")

tsw = pd.read_csv("/media/sluque/Data_2015/Data/ArcticNet/2015/UW_pCO2/"
                  "UW_water_temperature.csv",
                  parse_dates=["time"], index_col="time")
# Alternative from TSG
tsg = pd.read_csv("/media/sluque/Data_2015/Data/ArcticNet/2015/UW_pCO2/"
                  "TSG.csv",
                  parse_dates=["time"], index_col="time")

pCO2_equ = pCO2[pCO2.record_type == "EQU"]
pCO2_std2 = pCO2[pCO2.record_type == "STD2"]
pCO2_std3 = pCO2[pCO2.record_type == "STD3"]
pCO2_std4 = pCO2[pCO2.record_type == "STD4"]

# Y-limits
CO2_lims = (100, 1000)
CO2zero_lims = (-1, 2)
CO2std3_lims = (445, 465)
CO2std4_lims = (600, 605)
H2O_lims = (2, 18)
H2Ozero_lims = (-0.1, 0.3)
analyzer_temp_lims = (5, 40)
pressure_lims = (970, 1040)
gasflow_lims = (50, 120)

equ_temp_lims = (-5, 15)
equ_pressure_lims = (-5, 5)
equ_cond_lims = (0, 20)

H2O_temp_lims = (-5, 30)

# LICOR EQU plots
fig, axs = plt.subplots(5, 1, sharex=True)
fig.set_size_inches((11.5, 12.5))
pCO2_equ[["uw_CO2_fraction"]].plot(ax=axs[0], title="EQU samples",
                                   ylim=CO2_lims, legend=False)
axs[0].set_ylabel('CO2 fraction\n($\mu$mol/mol)')
axs[0].set_xlabel('')
pCO2_equ[["uw_H2O_fraction"]].plot(ax=axs[1], ylim=H2O_lims, legend=False)
axs[1].set_ylabel('H2O fraction\n(mmol/mol)')
axs[1].set_xlabel('')
pCO2_equ[["uw_temperature_analyzer"]].plot(ax=axs[2],
                                           ylim=analyzer_temp_lims,
                                           legend=False)
axs[2].set_ylabel('Temperature analyzer\n(C$^\circ$)')
axs[2].set_xlabel('')
pCO2_equ[["uw_pressure_analyzer"]].plot(ax=axs[3], ylim=pressure_lims,
                                        legend=False)
axs[3].set_ylabel('Pressure analyzer\n(mbar)')
axs[3].set_xlabel('')
pCO2_equ[["air_flow_analyzer"]].plot(ax=axs[4], rot=0,
                                     ylim=gasflow_lims, legend=False)
axs[4].set_ylabel('Flow rate analyzer\n(ml/min)')
axs[4].set_xlabel('')
plt.savefig("licor_equ.png", bbox_inches="tight"); plt.close()
# LICOR STD2 plots
fig, axs = plt.subplots(5, 1, sharex=True)
fig.set_size_inches((11.5, 12.5))
pCO2_std2[["uw_CO2_fraction"]].plot(ax=axs[0], title="STD2 samples",
                                    ylim=CO2zero_lims, style=".-",
                                    legend=False)
axs[0].set_ylabel('CO2 fraction\n($\mu$mol/mol)')
axs[0].set_xlabel('')
pCO2_std2[["uw_H2O_fraction"]].plot(ax=axs[1], ylim=H2Ozero_lims,
                                    style=".-", legend=False)
axs[1].set_ylabel('H2O fraction\n(mmol/mol)')
axs[1].set_xlabel('')
pCO2_std2[["uw_temperature_analyzer"]].plot(ax=axs[2],
                                            ylim=analyzer_temp_lims,
                                            style=".-", legend=False)
axs[2].set_ylabel('Temperature analyzer\n(C$^\circ$)')
axs[2].set_xlabel('')
pCO2_std2[["uw_pressure_analyzer"]].plot(ax=axs[3], ylim=pressure_lims,
                                         style=".-", legend=False)
axs[3].set_ylabel('Pressure analyzer\n(mbar)')
axs[3].set_xlabel('')
pCO2_std2[["air_flow_analyzer"]].plot(ax=axs[4], rot=0, ylim=gasflow_lims,
                                      style=".-", legend=False)
axs[4].set_ylabel('Flow rate analyzer\n(ml/min)')
axs[4].set_xlabel('')
plt.savefig("licor_std2.png", bbox_inches="tight"); plt.close()
# LICOR STD3 plots
fig, axs = plt.subplots(5, 1, sharex=True)
fig.set_size_inches((11.5, 12.5))
pCO2_std3[["uw_CO2_fraction"]].plot(ax=axs[0], title="STD3 samples",
                                    ylim=CO2std3_lims, style=".-",
                                    legend=False)
axs[0].set_ylabel('CO2 fraction\n($\mu$mol/mol)')
axs[0].set_xlabel('')
pCO2_std3[["uw_H2O_fraction"]].plot(ax=axs[1], style=".-",
                                    ylim=H2Ozero_lims, legend=False)
axs[1].set_ylabel('H2O fraction\n(mmol/mol)')
axs[1].set_xlabel('')
pCO2_std3[["uw_temperature_analyzer"]].plot(ax=axs[2], style=".-",
                                            ylim=analyzer_temp_lims,
                                            legend=False)
axs[2].set_ylabel('Temperature analyzer\n(C$^\circ$)')
axs[2].set_xlabel('')
pCO2_std3[["uw_pressure_analyzer"]].plot(ax=axs[3], style=".-",
                                         ylim=pressure_lims, legend=False)
axs[3].set_ylabel('Pressure analyzer\n(mbar)')
axs[3].set_xlabel('')
pCO2_std3[["air_flow_analyzer"]].plot(ax=axs[4], rot=0, style=".-",
                                      ylim=gasflow_lims, legend=False)
axs[4].set_ylabel('Flow rate analyzer\n(ml/min)')
axs[4].set_xlabel('')
plt.savefig("licor_std3.png", bbox_inches="tight"); plt.close()
# LICOR STD4 plots
fig, axs = plt.subplots(5, 1, sharex=True)
fig.set_size_inches((11.5, 12.5))
pCO2_std4[["uw_CO2_fraction"]].plot(ax=axs[0], title="STD4 samples",
                                    ylim=CO2std4_lims, style=".-",
                                    legend=False)
axs[0].set_ylabel('CO2 fraction\n($\mu$mol/mol)')
axs[0].set_xlabel('')
pCO2_std4[["uw_H2O_fraction"]].plot(ax=axs[1], style=".-",
                                    ylim=H2Ozero_lims, legend=False)
axs[1].set_ylabel('H2O fraction\n(mmol/mol)')
axs[1].set_xlabel('')
pCO2_std4[["uw_temperature_analyzer"]].plot(ax=axs[2], style=".-",
                                            ylim=analyzer_temp_lims,
                                            legend=False)
axs[2].set_ylabel('Temperature analyzer\n(C$^\circ$)')
axs[2].set_xlabel('')
pCO2_std4[["uw_pressure_analyzer"]].plot(ax=axs[3], style=".-",
                                         ylim=pressure_lims, legend=False)
axs[3].set_ylabel('Pressure analyzer\n(mbar)')
axs[3].set_xlabel('')
pCO2_std4[["air_flow_analyzer"]].plot(ax=axs[4], rot=0, style=".-",
                                      ylim=gasflow_lims, legend=False)
axs[4].set_ylabel('Flow rate analyzer\n(ml/min)')
axs[4].set_xlabel('')
plt.savefig("licor_std4.png", bbox_inches="tight"); plt.close()

# Equilibrator plots
fig, axs = plt.subplots(3, 1, sharex=True)
fig.set_size_inches((11, 9))
pCO2_equ[["equ_temperature"]].plot(ax=axs[0], ylim=equ_temp_lims,
                                   legend=False)
axs[0].set_ylabel('Equilibrator temperature\n(C$^\circ$)')
axs[0].set_xlabel('')
pCO2_equ[["equ_pressure"]].plot(ax=axs[1], ylim=equ_pressure_lims,
                                legend=False)
axs[1].set_ylabel('Equilibrator pressure\n(mbar)')
axs[1].set_xlabel('')
pCO2_equ[["condensation_equ"]].plot(ax=axs[2], rot=0,
                                    ylim=equ_cond_lims, legend=False)
axs[2].set_ylabel('Equilibrator condensation')
axs[2].set_xlabel('')
plt.savefig("equilibrator_equ.png", bbox_inches="tight"); plt.close()

# CTD plots
fig, axs = plt.subplots(3, 1, sharex=True)
fig.set_size_inches((11, 9))
pCO2_equ[["ctd_temperature"]].plot(ax=axs[0], style=".-",
                                   ylim=equ_temp_lims,
                                   legend=False)
axs[0].set_ylabel('CTD temperature\n(C$^\circ$)')
axs[0].set_xlabel('')
pCO2_equ[["equ_temperature"]].plot(ax=axs[0], ylim=equ_temp_lims,
                                   legend=True)
pCO2_equ[["ctd_pressure"]].plot(ax=axs[1], legend=False)
axs[1].set_ylabel('CTD pressure\n(dbar)')
axs[1].set_xlabel('')
pCO2_equ[["ctd_conductivity"]].plot(ax=axs[2], rot=0, legend=False)
axs[2].set_ylabel('CTD conductivity\n(ms/cm)')
axs[2].set_xlabel('')
plt.savefig("ctd_equ.png", bbox_inches="tight"); plt.close()

# External temperature plots
fig, axs = plt.subplots(2, 1, sharex=True)
fig.set_size_inches((11, 7))
tsw[["water_temperature"]].plot(ax=axs[0], style=".-", ylim=H2O_temp_lims,
                                legend=False)
tsg[["water_temperature"]].plot(ax=axs[0], ylim=H2O_temp_lims, legend=False)
axs[0].legend(axs[0].get_lines(), ["CR23X", "TSG"])
axs[0].set_ylabel('Water temperature\n(C$^\circ$)')
axs[0].set_xlabel('')
tsw[["water_temperature_sd"]].plot(ax=axs[1], rot=0,
                                   ylim=(0, 3), legend=False)
axs[1].set_ylabel('Water temperature SD\n(C$^\circ$)')
axs[1].set_xlabel('')
plt.savefig("external_temperature.png", bbox_inches="tight"); plt.close()

# # And here's the reason: temperature_external constant until mid-september
# # in the database, whereas Tim's data vary; why?
# plt.switch_backend("Agg")
# plt.figure(figsize=(13, 6))
# pCO2_tim.temp2.plot(rot=0)
# pCO2_db.temperature_external.plot(rot=0)
# plt.ylabel("External Temperature (C)")
# leg = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), frameon=False,
#                  borderaxespad=0, ncol=2)
# leg.get_texts()[0].set_text("Tim")
# leg.get_texts()[1].set_text("Database")
# plt.tight_layout()
# plt.savefig("external_temperature_2010.png", bbox_extra_artists=(leg,),
#             bbox_inches='tight')
# plt.close()
