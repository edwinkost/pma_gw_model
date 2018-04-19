# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pcraster as pcr
import numpy as np

# set the clone map
pcr.setclone("input/clone.map")
clone_map = pcr.boolean(1.0)

# get the cell dimension from the clone map (length and width, in meter)
cell_length = pcr.clone().cellSize()   
cell_width  = pcr.clone().cellSize()
# - cell area (m2)
cell_area   = cell_length * cell_width

# landmask
landmask = pcr.readmap("input/clone.map")

# original dem
input_dem = pcr.readmap("input/DEM_150929_110004.map")

# interpolate using inverse distance function
interpolated_dem = pcr.inversedistance(landmask, input_dem, 2, 0, 10) 

# save the interpolated dem to a file
pcr.report(interpolated_dem, "output_transient/interpolated_dem.map")
# and display it
pcr.aguila("output_transient/interpolated_dem.map")

# initial head - estimated from the steady state model
initial_head_steady_state = pcr.readmap("output_steady_state/groundwater_head.map")

# timestep
for timestep in range(0, 900):
    
    # initialize modflow object 
    modflow_object = pcr.initialise(pcr.clone())
    
    # defining the layer
    thickness = 50.0
    bottom_elevation = interpolated_dem - thickness
    top_elevation = interpolated_dem
    modflow_object.createBottomLayer(bottom_elevation, top_elevation)
    
    # set the DIS parameter, see http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/dis.html
    #  - time and spatial units
    ITMUNI = 4 # indicating that the time unit is "days"
    LENUNI = 2 # indicating that the spatial unit is "meters"
    # - PERLEN: duration of stress period (days)
    PERLEN = 1./24. # hourly stress period
    # - NSTP: number of sub time steps within the PERLEN
    NSTP = 1  
    # - TSMULT # always 1 by default
    TSMULT = 1
    # - SSTR: transient (0) or steady state (1)
    SSTR = 0
    modflow_object.setDISParameter(ITMUNI, LENUNI, PERLEN, NSTP, TSMULT, SSTR)
    
    # set the IBOND of the BAS package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bas.html
    # - in the ocean region (x < -75 m), assume the heads will follow the tides 
    ibound = pcr.ifthenelse(pcr.xcoordinate(clone_map) < -75., pcr.nominal(-1.), pcr.nominal(1.))
    pcr.report(ibound, "output_transient/ibound.map")
    # pcr.aguila("output_transient/ibound.map")
    modflow_object.setBoundary(ibound, 1)
    
    timestep_in_day  = timestep * PERLEN   # day
    
    # tide water level (m, relative to MSL)
    tide_amplitude       = 1.0        # meter
    tide_periode_in_hour = 12.4       # hour
    tide_periode_in_day  = 12.4 / 24. # day
    pi = 3.14159265359
    #tide_water_level = tide_amplitude * pcr.sin( (2.0 * pi * timestep_in_day / (tide_periode_in_day )) * 0.0174532925 )
    tide_water_level = tide_amplitude * np.sin( (2.0 * pi * timestep_in_day / (tide_periode_in_day )) )
    print(tide_water_level)
    
    # set the initial head of the BAS package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bas.html
    # - for the first test, use the DEM as initial head
    if timestep == 0:
        initial_head = initial_head_steady_state
    else:
        initial_head = groundwater_head
    # - in the ocean (ibound = -1), groundwater head is equal to the tide
    initial_head = pcr.ifthenelse(pcr.scalar(ibound) > 0, initial_head, tide_water_level)
    modflow_object.setInitialHead(initial_head, 1)
    
    # set the conductivities for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
    # - sand conductivity in m.day-1
    sand_conductivity = pcr.spatial(pcr.scalar(10.))
    # - horizontal and vertical conductivity
    hConductivity = sand_conductivity 
    vConductivity = sand_conductivity # for 1 layer case, this is just dummy and never used
    # layer type, we use LAYTYPE = 0 (harmonic mean) and LAYCON = 0 (confined, constant transmissivities and storage coefficients)
    modflow_object.setConductivity(00, hConductivity, vConductivity, 1)
    
    # set the storage coefficients for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
    # - sand porosity (m3.m-3)
    sand_porosity = pcr.spatial(pcr.scalar(0.25))
    # - prtimary and secondary storage coefficients 
    primary_storage_coefficient   = sand_porosity
    secondary_storage_coefficient = primary_storage_coefficient  # for LAYCON = 0 (and 1), this is just dummy and never used
    modflow_object.setStorage(primary_storage_coefficient, secondary_storage_coefficient, 1)
    
    # set the RIVER package
    bottom_morphology                  = interpolated_dem
    bed_thickness                      = 0.001
    tide_water_level_entering_the_land = pcr.ifthenelse(bottom_morphology < tide_water_level, tide_water_level, bottom_morphology)
    # conductance for the RIVER package (m2.day-1)
    bed_conductance = sand_porosity * cell_area / bed_thickness
    # - inactive the RIV package for dry areas
    bed_conductance = pcr.ifthenelse(bottom_morphology < tide_water_level, bed_conductance, 0.0)
    # - set the RIV package
    modflow_object.setRiver(tide_water_level_entering_the_land, bottom_morphology, bed_conductance, 1)
    
    # set the SOLVER package 
    MXITER = 50                 # maximum number of outer iterations           # Deltares use 50
    ITERI  = 30                 # number of inner iterations                   # Deltares use 30
    NPCOND = 1                  # 1 - Modified Incomplete Cholesky, 2 - Polynomial matrix conditioning method
    HCLOSE = 0.001              # HCLOSE (unit: m) 
    RCLOSE = 0.001              # RCLOSE (unit: m3)
    RELAX  = 1.00               # relaxation parameter used with NPCOND = 1
    NBPOL  = 2                  # indicates whether the estimate of the upper bound on the maximum eigenvalue is 2.0 (but we don ot use it, since NPCOND = 1) 
    DAMP   = 1                  # no damping (DAMP introduced in MODFLOW 2000)
    modflow_object.setPCG(MXITER, ITERI, NPCOND, HCLOSE, RCLOSE, RELAX, NBPOL, DAMP)
    
    # execute the MODFLOW run and write all modflow temporary files to a certain folder
    temporary_folder = "temp"
    modflow_object.run(temporary_folder)
    
    # get the output
    # - groundwater head (m)
    groundwater_head = modflow_object.getHeads(1)
    groundwater_head_filename = "output_transient/groundwater_head" + "_" + str(timestep) + ".map"
    pcr.report(groundwater_head, groundwater_head_filename)
    # pcr.aguila(groundwater_head_filename)
    # - groundwater head (m) in time series
    groundwater_tss_file_name = "output_transient/h0000000" + ".00" + str(timestep)
    if timestep >= 10: groundwater_tss_file_name = "output_transient/h0000000" + ".0" + str(timestep)
    if timestep >= 100: groundwater_tss_file_name = "output_transient/h0000000" + "." + str(timestep)
    if timestep >= 1000: groundwater_tss_file_name = "output_transient/h0000001" + ".00" + str(timestep - 1000)
    if timestep >= 1010: groundwater_tss_file_name = "output_transient/h0000001" + ".0" + str(timestep - 1000)
    pcr.report(groundwater_head, groundwater_tss_file_name)
