#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import time
import datetime

#~ from pcraster import *
from pcraster.framework import *

import pcraster as pcr
import numpy as np

import output_netcdf_writer

class PantaiMukaAirTanahModel(DynamicModel, MonteCarloModel):

    def __init__(self, model_setup):

        # In this part ("__init__"), we initiate pcraster frameworks and set the clone map.
        
        # initiate pcraster dynamic and monte carlo frameworks   
        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)
        
        # make model_setup available for the entire method
        self.model_setup = model_setup
        
        # set the clone map based on DEM
        self.clone_map = self.model_setup['dem_file_name']
        pcr.setclone(self.clone_map)
        # - landmask - needed if we want to mask out some areas/cells
        self.landmask = pcr.defined(pcr.readmap(self.clone_map))

        # output folder
        self.output_folder = self.model_setup['output_folder']
        
        # initiate a netcdf writer
        self.netcdf_writer = output_netcdf_writer.OutputNetCDF(pcr.clone())
        
    def premcloop(self):

        # In this part (premcloop), we initiate parameters/variables/objects that are the same (constant) values throughout all monte carlo (MC) samples. 
        # - Note that you can also write this part in "__init__" (see above). 

        # get the map properties of the clone map (length and width, in meter)
        self.cell_length = pcr.clone().cellSize()   
        self.cell_width  = pcr.clone().cellSize()
        # - cell area (m2)
        self.cell_area   = self.cell_length * self.cell_width

        # digital elevation model (m)
        self.input_dem = pcr.readmap(self.model_setup['dem_file_name'])
        
        # before entering the initial part, save the current working direcory
        self.script_directory = os.getcwd()
        
        
    def initial(self):

        # make sure that you start from the script directory
        os.chdir(self.script_directory)
        
        # In this part (premcloop), we initiate parameters/variables/objects that are changing throughout all monte carlo samples. 

        msg  = "\n" 
        msg += "Sample number: " + str(self.currentSampleNumber()) 
        msg += "\n" 
        print(msg)
        
        # conductivities for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
        inp_soil_conductivity = self.model_setup['soil_conductivity'][int(str(self.currentSampleNumber())) - 1]
        self.soil_conductivity = pcr.spatial(pcr.scalar(inp_soil_conductivity))

        #
        # - horizontal and vertical conductivity
        self.hConductivity = self.soil_conductivity 
        self.vConductivity = self.hConductivity          
        # - for one layer model, vConductivity is just dummy and never used
        # - layer type, we use LAYTYPE = 0 (harmonic mean) and LAYCON = 0 (confined, constant transmissivities and storage coefficients)
        
        # storage coefficients for the BCF package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bcf.html
        # - sand porosity (m3.m-3)                                                                                         # TODO: Find the value from Sebastian paper.
        self.sand_porosity = pcr.spatial(pcr.scalar(0.25))
        #
        # - primary and secondary storage coefficients 
        self.primary_storage_coefficient   = self.sand_porosity
        self.secondary_storage_coefficient = self.primary_storage_coefficient  
        # - for LAYCON = 0 (and 1), secondary_storage_coefficient is just dummy and never used      


        
        # The following are CURRENTLY just the same for all samples. 
        ############################################################################
        
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
        
        # ibound with regional groundwater head values
        if self.model_setup['regional_groundwater_head']['activation']:
            self.ibound = pcr.ifthenelse(pcr.xcoordinate(pcr.boolean(1.0)) >= self.model_setup['regional_groundwater_head']['starting_x'], pcr.nominal(-1.), self.ibound)

        
        # parameter values for the SOLVER package 
        self.MXITER = 50                 # maximum number of outer iterations           # Deltares use 50
        self.ITERI  = 30                 # number of inner iterations                   # Deltares use 30
        self.NPCOND = 1                  # 1 - Modified Incomplete Cholesky, 2 - Polynomial matrix conditioning method
        self.HCLOSE = 0.001              # HCLOSE (unit: m) 
        self.RCLOSE = 0.001              # RCLOSE (unit: m3)
        self.RELAX  = 1.00               # relaxation parameter used with NPCOND = 1
        self.NBPOL  = 2                  # indicates whether the estimate of the upper bound on the maximum eigenvalue is 2.0 (but we don ot use it, since NPCOND = 1) 
        self.DAMP   = 1                  # no damping (DAMP introduced in MODFLOW 2000)
        


        # the initial head for the BAS package, see: http://pcraster.geo.uu.nl/pcraster/4.1.0/doc/modflow/bas.html
        # - Gerben recommends to start using 0.8 m
        self.initial_head = pcr.spatial(pcr.scalar(0.8))
        
        # regional groundwater head
        if self.model_setup['regional_groundwater_head']['activation']:
            self.initial_head = pcr.ifthenelse(pcr.xcoordinate(pcr.boolean(1.0)) >= self.model_setup['regional_groundwater_head']['starting_x'], self.model_setup['regional_groundwater_head']['value'], self.initial_head)


        # initialise timeoutput object for reporting time series in txt files
        # - groundwater head
        self.head_obs_point  = pcr.readmap("input_files/groundwater_well_coordinates.map")
        self.reportGwHeadTss = TimeoutputTimeseries("groundwater_head", self, self.head_obs_point, noHeader=False)


        # Save model parameter values (and other information) to a txt file.
        # - output directory (that will contain result)
        output_directory = self.output_folder + "/" + str(self.currentSampleNumber())
        try:
            os.makedirs(output_directory)
        except:
            pass
        # - file name for this information
        information_file = output_directory + "/" + "info.txt"
        file_info = open(information_file, 'w')
        write_line  = "" 
        # - DEM
        write_line += "DEM (m): " + str(self.model_setup['dem_file_name'])
        write_line += "\n"
        # - tide
        write_line += "Tide input file name: " + str(self.model_setup['tide_file_name'])
        write_line += "\n"
        # - starting date
        write_line += "Starting date and time: " + str(self.model_setup['start_datetime'])
        write_line += "\n"
        # - soil conductivity
        write_line += "Soil conductivity (m.day-1): " + str(inp_soil_conductivity)
        write_line += "\n"
        # - regional groundwater head (optional):
        if self.model_setup['regional_groundwater_head']['activation']:
            write_line += "Regional groundwater head is fixed to : " + str(self.model_setup['regional_groundwater_head']['value'])
            write_line += "\n"
            write_line += "for all cells in  with x coordinate >=: " + str(self.model_setup['regional_groundwater_head']['starting_x'])
            write_line += "\n"
        else:
            write_line += "NO regional groundwater head. "
            write_line += "\n"
        # - sample number
        write_line += "PCRaster sample number: " + str(self.currentSampleNumber())
        write_line += "\n"
        file_info.write(write_line)
        # - close the file
        file_info.close()
        
        
        # initiate netcdf files for groundwater head
        netcdf_variable_name  = "groundwater_head"
        self.netcdf_file_name = self.output_folder + "/" + str(self.currentSampleNumber()) + "/" + netcdf_variable_name + ".nc"
        self.model_setup['netcdf_attributes']['notes'] = write_line
        self.netcdf_writer.create_netcdf_file(self.netcdf_file_name, self.model_setup['netcdf_attributes'])
        # - create variable
        netcdf_variable_unit  = "m"
        self.netcdf_writer.create_variable(self.netcdf_file_name, netcdf_variable_name, netcdf_variable_unit)   


        
    def dynamic(self):

        # In this part (dynamic), the parameters/variables/objects are changing over time.  

        msg  = "\n" 
        msg += "Sample number: " + str(self.currentSampleNumber())  + " ; " + "time step number: " + str(self.currentTimeStep())
        msg += "\n" 
        print(msg)


        # time step index (from pcraster)
        self.time_step_index = self.currentTimeStep()

        # timestep in day unit
        # - a stress period contains a time step (10 minute length)
        self.timestep_in_day = self.currentTimeStep() * self.length_of_stress_period
        

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

        #~ # tide water level (m, relative to MSL) - assume a simple sinusoidal function
        #~ tide_amplitude       = 1.0        # meter
        #~ tide_periode_in_hour = 12.4       # hour
        #~ tide_periode_in_day  = 12.4 / 24. # day
        #~ self.tide_water_level = tide_amplitude * np.sin( (2.0 * np.pi * self.timestep_in_day / (tide_periode_in_day )) )
        #
        
        # tide water level from the file (m, relative to MSL???)
        # - average from two measurements
        self.tide_water_level = 0.5 * (float(self.model_setup['tide_series'][self.time_step_index-1].split()[1]) + float(self.model_setup['tide_series'][self.time_step_index].split()[1]))


        #~ # - far in the ocean (ibound = -1), groundwater head is equal to the tide - NOT NEEDED (all cells are active)
        #~ self.initial_head = pcr.ifthenelse(pcr.scalar(self.ibound) > 0, self.initial_head, self.tide_water_level)


        # set the RIVER package
        # - coastal morphology (elevation)
        self.bottom_morphology = self.top_elevation
        # - conductance for the RIVER package (m2.day-1)
        # -- reset bed conductance for every time step
        self.bed_conductance = None
        # -- assume very thin river bed thickness (m)
        bed_thickness = 0.001
        self.bed_conductance = self.soil_conductivity * self.cell_area / bed_thickness
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
        try:
            temporary_folder = str(self.currentSampleNumber())
            os.makedirs(temporary_folder)
        except:
            pass
        try:
            # new pcraster
            self.modflow_object.run(temporary_folder)
        except:
            # old pcraster
            os.chdir(temporary_folder)
            self.modflow_object.run()
            pass

        # get the output
        # - groundwater head (m)
        self.groundwater_head = self.modflow_object.getHeads(1)
        try:
            # new pcraster
            self.report(self.groundwater_head, "h")
        except:
            # old pcraster
            os.chdir(self.output_folder)
            self.report(self.groundwater_head, "h")
            pass
        

        # set the calculate head for the next time step
        self.initial_head = self.groundwater_head
        
        
        # sampling timeseries for given locations
        self.reportGwHeadTss.sample(self.groundwater_head)
        
        # write netcdf files for groundwater head
        netcdf_variable_name  = "groundwater_head"
        field_value = pcr.pcr2numpy(self.groundwater_head, 1e20)
        # - time bounds for netcdf files
        lowerTimeBound = self.model_setup['start_datetime'] + (self.time_step_index - 1) * datetime.timedelta(days = self.length_of_stress_period)
        upperTimeBound = lowerTimeBound + datetime.timedelta(days = self.length_of_stress_period)
        time_bounds = [lowerTimeBound, upperTimeBound]
        self.netcdf_writer.data_to_netcdf(self.netcdf_file_name, netcdf_variable_name, field_value, time_bounds)   

        
        # make sure that you return to the output folder
        os.chdir(self.output_folder)


    def postmcloop(self):
        
        calculate_statistics = False
        
        if len(self.sampleNumbers()) > 1 and calculate_statistics:
            
            print("Get some statistics.")
            
            # - go to the output folder before doing some statistics
            os.chdir(self.output_folder)
            
            names = ["h"]
            mcaveragevariance(names, self.sampleNumbers(), self.timeSteps())
            percentiles = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
            mcpercentiles(names, percentiles, self.sampleNumbers(), self.timeSteps())



def main():
    
    # initiate a dictionary that will contain input values (and some other information)  
    model_setup = {}
    

    # SET YOUR OUTPUT FOLDER HERE 
    # - cartesius
    model_setup['output_folder']    = "/scratch-shared/edwinsut/test_output_yvonne/test_with_regional_head/"
    model_setup['output_folder']    = "/scratch-shared/edwinsut/test_output_yvonne/test_without_regional_head/"
    #~ # - speedy
    #~ model_setup['output_folder'] = "/scratch/edwin/test_output_yvonne/test_using_old_pcraster/"
    #~ # - WINDOWS
    #~ model_setup['output_folder'] = "C:/test/"
    
    # create output folder
    cleaning_previous_output_folder = True
    try: 
        os.makedirs(model_setup['output_folder'])
    except:
        if cleaning_previous_output_folder: 
          cmd = 'rm -r ' + model_setup['output_folder'] + "/*"
          os.system(cmd)
    
    ####################################################################################################
    
    

    # INPUT SECTION
    # - all input file names are RELATIVE to the directory where THIS SCRIPT is stored. 
    
    # tide 
    # - the file should includes some extra data for spinning up purpose.
    model_setup['tide_file_name'] = "input_files/tide_setup_modflowtime_egmond_version_2018-10-05.txt"
    # - read the tide file
    file_tide = open(model_setup['tide_file_name'], "r")
    # - tide water levels (m) for every 10 minutes (600 seconds)
    model_setup['tide_series'] = file_tide.readlines()
    file_tide.close()

    # the starting date of the model run based on tide file (Yvonne: 7 September 2015)
    model_setup['start_datetime'] =  datetime.datetime(int("2015"), int("09"), int("07"), int("00"), int("00"))     
    
    # the number of timesteps based on the length of tide file
    model_setup['number_of_time_steps'] = len(model_setup['tide_series']) - 1
    #~ # - for testing
    #~ model_setup['number_of_time_steps'] = 3
    
    
    # elevation (DEM, in meter)
    model_setup['dem_file_name'] = "input_files/DEM_150929_110004_Jetski.map"
    
    
    # conductivity values (m.day-1) in a list: value for every sample
    model_setup['soil_conductivity'] = [100.0,
                                        50.0,
                                        45.0,
                                        40.0,
                                        35.0,
                                        30.0,
                                        25.0,
                                        20.0,
                                        15.0,
                                        10.0,
                                        5.0,
                                        2.0,
                                        1.0]
    
    # number of samples
    # - based on number of conductivity values
    number_of_samples = len(model_setup['soil_conductivity'])
    

    # regional groundwater head (m from MSL) boundary condition (inland side)
    model_setup['regional_groundwater_head'] = {}
    # - value (m, from MSL)
    model_setup['regional_groundwater_head']['value'] = 1.00
    # - position, all cells in the 'right' side of this will have fixed heads as defined above
    model_setup['regional_groundwater_head']['starting_x'] = 26.00 
    model_setup['regional_groundwater_head']['activation'] = True 
    model_setup['regional_groundwater_head']['activation'] = False
    
    
    # PS: There are also more input files/values that are harcoded in the other parts.   


    # attribute information for netcdf files
    model_setup['netcdf_attributes'] = {}
    model_setup['netcdf_attributes']['institution'] = "Utrecht University, Dept. of Physical Geography"
    model_setup['netcdf_attributes']['title'      ] = "Groundwater simulation (Egmond aan Zee)"
    model_setup['netcdf_attributes']['source'     ] = "Model scripts were written by Edwin H. Sutanudjaja (E.H.Sutanudjaja@uu.nl)."
    model_setup['netcdf_attributes']['references' ] = "Smit et al. (in prep.)"
    model_setup['netcdf_attributes']['description'] = "Groundwater simulation (Egmond aan Zee)"
    model_setup['netcdf_attributes']['created by' ] = "The model run was performed by Yvonne Smit (Y.Smit@uu.nl). Model scripts were written by Edwin H. Sutanudjaja (E.H.Sutanudjaja@uu.nl)."
    model_setup['netcdf_attributes']['notes'      ] = ""

    

    
    ####################################################################################################


    # STARTING THE MODEL
    
    myModel = PantaiMukaAirTanahModel(model_setup)
    dynamicModel = DynamicFramework(myModel, lastTimeStep = model_setup['number_of_time_steps'], firstTimestep = 1)
    
    # define the number of samples here
    mcModel = MonteCarloFramework(dynamicModel, nrSamples = number_of_samples)
    #~ mcModel = MonteCarloFramework(dynamicModel, nrSamples = 1)
    
    # - forking only work for linux
    mcModel.setForkSamples(fork = True, nrCPUs = 20)
    #~ mcModel.setForkSamples(fork = True, nrCPUs = 5)
    
    # run
    mcModel.run()

        
if __name__ == '__main__':
    sys.exit(main())

