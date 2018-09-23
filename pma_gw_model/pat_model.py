#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pcraster as pcr

#~ from pcraster import *
import pcraster as pcr

from pcraster.framework import *

class PantaiAirTanahModel(DynamicModel, MonteCarloModel):

    def __init__(self):

        # In this part ("__init__"), we initiate pcraster frameworks and set the clone map and set the output folder. 
        
        # initiate pcraster dynamic and monte carlo frameworks   
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        
        # set the clone map based on DEM
        self.clone_map = "input_files/DEM_150929_110004_correct.map"
        pcr.setclone(self.clone_map)
        #~ # - landmask - needed if we want to mask out some areas/cells
        #~ self.landmask = pcr.readmap(self.clone_map)
        
        # set and create the output folder 
        self.output_folder = "output/test/"
        
        cleaning_previous_output_folder = False
        try: 
            os.makedirs()
        except:
            if cleaning_previous_output_folder: 
			    cmd = 'rm -r ' + self.output_folder
			    os.system(cmd)
        
    def premcloop(self):

        # In this part (premcloop), we initiate parameters/variables/objects that are the same (constant) values throughout all monte carlo (MC) samples. 

        # get the cell dimension from the clone map (length and width, in meter)
        self.cell_length = pcr.clone().cellSize()   
        self.cell_width  = pcr.clone().cellSize()
        # - cell area (m2)
        self.cell_area   = self.cell_length * self.cell_width

        # digital elevation model (m)
        input_dem = pcr.readmap("input/DEM_150929_110004_correct.map")
    
    def initial(self):

        # In this part (premcloop), we initiate parameters/variables/objects that are changing throughout all monte carlo samples. 

        # initialize modflow object - this object is unique for each sample
        self.modflow_object = pcr.initialise(pcr.clone())
        
        # defining the layer (one layer model), thickness (m), top and bottom elevations 
        self.thickness = 15.0
        # - thickness value is suggested by Kim Cohen (neede)
        self.top_elevation    = self.input_dem
        self.bottom_elevation = self.top_elevation - self.thickness
        # - set one modflow layer model
        self.modflow_object.createBottomLayer(bottom_elevation, top_elevation)
        
        # set the DIS parameter, see http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/dis.html
        #  - time and spatial units
        ITMUNI = 4 # indicating that the time unit is "days"
        LENUNI = 2 # indicating that the spatial unit is "meters"
        # - PERLEN: duration of stress period (days)
        # PERLEN = 1./24. # hourly stress period
        PERLEN = 600. / (24.*60.*60.)   # 600 seconds stress period 
        # - NSTP: number of sub time steps within the PERLEN
        NSTP = 1  
        # - TSMULT # always 1 by default
        TSMULT = 1
        # - SSTR: transient (0) or steady state (1)
        SSTR = 0
        modflow_object.setDISParameter(ITMUNI, LENUNI, PERLEN, NSTP, TSMULT, SSTR)
        
        # set the IBOND of the BAS package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bas.html
        # - in the ocean region (x < -75 m), assume the heads will follow the tides 
        # ibound = pcr.ifthenelse(pcr.xcoordinate(clone_map) < -75., pcr.nominal(-1.), pcr.nominal(1.))
        ibound = pcr.ifthenelse(pcr.xcoordinate(clone_map) < -200., pcr.nominal(-1.), pcr.nominal(1.))
        pcr.report(ibound, "output_transient/ibound.map")
        if timestep == 0: pcr.aguila("output_transient/ibound.map")
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
        sand_conductivity = pcr.spatial(pcr.scalar(10.))  # get Sebastian value
        # - horizontal and vertical conductivity
        hConductivity = sand_conductivity 
        vConductivity = sand_conductivity # for 1 layer case, this is just dummy and never used
        # layer type, we use LAYTYPE = 0 (harmonic mean) and LAYCON = 0 (confined, constant transmissivities and storage coefficients)
        modflow_object.setConductivity(00, hConductivity, vConductivity, 1)
        
        # set the storage coefficients for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
        # - sand porosity (m3.m-3)
        sand_porosity = pcr.spatial(pcr.scalar(0.25))   # get Sebastian value
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

        self.snow = scalar(0)
        self.temperatureLapseRate = 0.005 + (mapnormal() * 0.001)
        self.report(self.temperatureLapseRate, "lapse")
        self.temperatureCorrection = self.elevationAboveMeteoStation\
             * self.temperatureLapseRate
    
    def dynamic(self):

        # run modflow and report
        
        temperatureObserved = self.readDeterministic("tavgo")
        precipitationObserved = self.readDeterministic("pr")
        precipitation = max(0, precipitationObserved * (mapnormal() * 0.2 + 1.0))
        temperature = temperatureObserved - self.temperatureCorrection
        snowFall = ifthenelse(temperature < 0, precipitation, 0)
        self.snow = self.snow + snowFall
        potentialMelt = ifthenelse(temperature > 0, temperature\
             * self.degreeDayFactor, 0)
        actualMelt = min(self.snow, potentialMelt)
        self.snow = max(0, self.snow - actualMelt)
        rain = ifthenelse(temperature >= 0, precipitation, 0)
        discharge = accuflux(self.ldd, actualMelt + rain)
        self.report(self.snow, "s")
        self.report(discharge, "q")
        
    def postmcloop(self):
        
        
        names = ["s", "q"]
        mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        percentiles = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())

myModel = SnowModel()
dynamicModel = DynamicFramework(myModel, lastTimeStep=180, firstTimestep=1)
mcModel = MonteCarloFramework(dynamicModel, nrSamples=10)
mcModel.run()
