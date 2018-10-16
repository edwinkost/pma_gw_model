#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

#~ from pcraster import *
from pcraster.framework import *

import pcraster as pcr
import numpy as np

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
        self.output_folder = "/scratch-shared/edwinhs/output_yvonne/test/"
        # - create output folder
        cleaning_previous_output_folder = True
        try: 
            os.makedirs(self.output_folder)
        except:
            if cleaning_previous_output_folder: 
			    cmd = 'rm -r ' + self.output_folder + "/*"
			    os.system(cmd)
        
    def premcloop(self):

        # In this part (premcloop), we initiate parameters/variables/objects that are the same (constant) values throughout all monte carlo (MC) samples. 

        # get the cell dimension from the clone map (length and width, in meter)
        self.cell_length = pcr.clone().cellSize()   
        self.cell_width  = pcr.clone().cellSize()
        # - cell area (m2)
        self.cell_area   = self.cell_length * self.cell_width

        # digital elevation model (m)
        self.input_dem = pcr.readmap("input_files/DEM_150929_110004_correct.map")
        
        # read/generate some lists for perturbing model parameter values
        # TODO
        
        # read tide
        file_tide = open("input_files/tide_setup_modflowtime_egmond.txt", "r")
        self.time_and_tide = file_tide.readlines()
        file_tide.close()
        
        
    def initial(self):

        # In this part (premcloop), we initiate parameters/variables/objects that are changing throughout all monte carlo samples. 

        
        # defining the layer (one layer model), thickness (m), top and bottom elevations 
        self.thickness = 15.0
        # - thickness value is suggested by Kim Cohen (put reference here)
        self.top_elevation    = self.input_dem
        self.bottom_elevation = self.top_elevation - self.thickness

        # DIS parameters, see http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/dis.html
        #  - time and spatial units
        self.ITMUNI = 4 # indicating that the time unit is "days"
        self.LENUNI = 2 # indicating that the spatial unit is "meters"
        # - PERLEN: duration of stress period (days)
        # -- 10 minute stress period = 600 seconds stress period 
        self.length_of_stress_period = 600. / (24.*60.*60.)  
        self.PERLEN = self.length_of_stress_period 
        # - NSTP: number of sub time steps within the PERLEN
        self.NSTP = 1  
        # - TSMULT # always 1 by default
        self.TSMULT = 1
        # - SSTR: transient (0) or steady state (1)
        self.SSTR = 0
        
        # values for the IBOND of the BAS package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bas.html
        # - Alternative 1: all cells are active 
        self.ibound = pcr.spatial(pcr.nominal(1.))
        #~ # - Alternative 2: in the ocean region (x < -75 m), assume the heads will follow the tides 
        #~ self.ibound = pcr.ifthenelse(pcr.xcoordinate(clone_map) < -75., pcr.nominal(-1.), pcr.nominal(1.))
        #~ pcr.aguila(self.ibound)
        
        # the initial head for the BAS package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bas.html
        # - for the first test, use the DEM as initial head
        self.initial_head = self.input_dem
        
        # conductivities for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
        # - sand conductivity in m.day-1 # TODO: Find the value from Sebastian paper. 
        self.sand_conductivity = pcr.spatial(pcr.scalar(10.))
        #
        if self.currentSampleNumber() == 1: self.sand_conductivity = pcr.spatial(pcr.scalar(10.))
        if self.currentSampleNumber() == 2: self.sand_conductivity = pcr.spatial(pcr.scalar(2.0))
        if self.currentSampleNumber() == 3: self.sand_conductivity = pcr.spatial(pcr.scalar(1.0))
        if self.currentSampleNumber() == 4: self.sand_conductivity = pcr.spatial(pcr.scalar(0.1))
        if self.currentSampleNumber() == 5: self.sand_conductivity = pcr.spatial(pcr.scalar(.01))
        if self.currentSampleNumber() == 6: self.sand_conductivity = pcr.spatial(pcr.scalar(20.))
        if self.currentSampleNumber() == 7: self.sand_conductivity = pcr.spatial(pcr.scalar(15.))
        if self.currentSampleNumber() == 8: self.sand_conductivity = pcr.spatial(pcr.scalar(5.0))
        #
        # - horizontal and vertical conductivity
        self.hConductivity = self.sand_conductivity 
        self.vConductivity = self.hConductivity          
        # - for one layer model, vConductivity is just dummy and never used
        # - layer type, we use LAYTYPE = 0 (harmonic mean) and LAYCON = 0 (confined, constant transmissivities and storage coefficients)
        
        # storage coefficients for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
        # - sand porosity (m3.m-3)       # TODO: Find the value from Sebastian paper.
        self.sand_porosity = pcr.spatial(pcr.scalar(0.25))   # get Sebastian value
        # - primary and secondary storage coefficients 
        self.primary_storage_coefficient   = self.sand_porosity
        self.secondary_storage_coefficient = self.primary_storage_coefficient  
        # - for LAYCON = 0 (and 1), secondary_storage_coefficient is just dummy and never used      

        # parameter values for the SOLVER package 
        self.MXITER = 50                 # maximum number of outer iterations           # Deltares use 50
        self.ITERI  = 30                 # number of inner iterations                   # Deltares use 30
        self.NPCOND = 1                  # 1 - Modified Incomplete Cholesky, 2 - Polynomial matrix conditioning method
        self.HCLOSE = 0.001              # HCLOSE (unit: m) 
        self.RCLOSE = 0.001              # RCLOSE (unit: m3)
        self.RELAX  = 1.00               # relaxation parameter used with NPCOND = 1
        self.NBPOL  = 2                  # indicates whether the estimate of the upper bound on the maximum eigenvalue is 2.0 (but we don ot use it, since NPCOND = 1) 
        self.DAMP   = 1                  # no damping (DAMP introduced in MODFLOW 2000)
        
        # starting time step index
        self.time_step_index = 0

        
    def dynamic(self):

        # In this part (dynamic), the parameters/variables/objects are changing over time.  

        # timestep in day unit
        # - a stress period contains a time step (10 minute length)
        self.timestep_in_day = self.currentTimeStep() * self.length_of_stress_period
        
        # time step index
        self.time_step_index = self.time_step_index + 1
        

        # initialize modflow object - this object is unique for each sample and also changing over time
        self.modflow_object = None
        self.modflow_object = pcr.initialise(pcr.clone())

        # set one modflow layer model
        self.modflow_object.createBottomLayer(self.bottom_elevation, self.top_elevation)

        # set ibound 
        self.modflow_object.setBoundary(self.ibound, 1)

        # set initial head
        # - for the first time step, this will be taken from the 'initial' part
        self.modflow_object.setInitialHead(self.initial_head, 1)

        # set conductivity values
        # - layer type, we use LAYTYPE = 0 (harmonic mean) and LAYCON = 0 (confined, constant transmissivities and storage coefficients)
        self.modflow_object.setConductivity(00, self.hConductivity, self.vConductivity, 1)

        # set storage coefficients
        # - for LAYCON = 0 (and 1), secondary_storage_coefficient is just dummy and never used
        self.modflow_object.setStorage(self.primary_storage_coefficient, self.secondary_storage_coefficient, 1)

        # set DIS parameters
        self.modflow_object.setDISParameter(self.ITMUNI, self.LENUNI, self.PERLEN, self.NSTP, self.TSMULT, self.SSTR)

        # set SOLVER package 
        self.modflow_object.setPCG(self.MXITER, self.ITERI, self.NPCOND, self.HCLOSE, self.RCLOSE, self.RELAX, self.NBPOL, self.DAMP)

        # tide water level (m, relative to MSL) - assume a simple sinusoidal function
        tide_amplitude       = 1.0        # meter
        tide_periode_in_hour = 12.4       # hour
        tide_periode_in_day  = 12.4 / 24. # day
        #~ self.tide_water_level = tide_amplitude * np.sin( (2.0 * np.pi * self.timestep_in_day / (tide_periode_in_day )) )
        #
        
        #~ print self.time_and_tide[self.time_step_index-1].split()[1]
        #~ print self.time_and_tide[self.time_step_index].split()[1]
        
        # TODO: Read this from the file
        self.tide_water_level = 0.5 * (float(self.time_and_tide[self.time_step_index-1].split()[1]) + float(self.time_and_tide[self.time_step_index].split()[1]))
        print(self.tide_water_level)

        #~ # - far in the ocean (ibound = -1), groundwater head is equal to the tide - NOT NEEDED
        #~ self.initial_head = pcr.ifthenelse(pcr.scalar(self.ibound) > 0, self.initial_head, self.tide_water_level)


        # set the RIVER package
        # - coastal morphology (elevation)
        self.bottom_morphology = self.top_elevation
        # - conductance for the RIVER package (m2.day-1)
        # -- reset bed conductance for every time step
        self.bed_conductance = None
        # -- assume very thin river bed thickness (m)
        bed_thickness = 0.001
        self.bed_conductance = self.sand_porosity * self.cell_area / bed_thickness
        # - conductance for the RIVER package (m2.day-1)
        self.tide_water_level_entering_the_land = pcr.ifthenelse(self.bottom_morphology < self.tide_water_level, self.tide_water_level, self.bottom_morphology)
        # - inactive the RIV package for dry areas
        self.bed_conductance = pcr.ifthenelse(self.bottom_morphology < self.tide_water_level, self.bed_conductance, 0.0)
        # - set the RIV package
        self.modflow_object.setRiver(self.tide_water_level_entering_the_land, self.bottom_morphology, self.bed_conductance, 1)

        
        # run modflow
        # - go to the output folder before executing MODFLOW
        os.chdir(self.output_folder)
        # - execute the MODFLOW run and write all modflow temporary files to a certain folder
        temporary_folder = str(self.currentSampleNumber())
        try:
            os.makedirs(temporary_folder)
        except:
            pass
        self.modflow_object.run(temporary_folder)
        
        # get the output
        # - groundwater head (m)
        self.groundwater_head = self.modflow_object.getHeads(1)
        self.report(self.groundwater_head, "h")
        
        # set the calculate head for the next time step
        self.initial_head = self.groundwater_head


    def postmcloop(self):
        
        # - go to the output folder before doing some statistics
        os.chdir(self.output_folder)

        names = ["h"]
        mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
        percentiles = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())

# 

myModel = PantaiAirTanahModel()
dynamicModel = DynamicFramework(myModel, lastTimeStep=5500, firstTimestep=1)
mcModel = MonteCarloFramework(dynamicModel, nrSamples=16)
mcModel.setForkSamples(fork = True, nrCPUs=16)
mcModel.run()