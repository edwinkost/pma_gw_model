#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import sys
import datetime
import time
import re
import glob
import subprocess
import netCDF4 as nc
import numpy as np
import virtualOS as vos

class OutputNetCDF():
    
    def __init__(self, pcraster_clone, projection_info = None, netcdf_y_orientation_follow_cf_convention = True, netcdf_format = 'NETCDF4', zlib = True):
        		
        # pcraster clone (containing the information of pcraster.clone().nrRows(), pcraster.clone().nrCols(), pcraster.clone().cellSize(), pcraster.clone().west(), pcraster.clone().north())
        self.x_min       = pcraster_clone.west()
        self.y_max       = pcraster_clone.north()
        self.num_of_cols = pcraster_clone.nrCols()
        self.num_of_rows = pcraster_clone.nrRows()
        self.cell_width  = pcraster_clone.cellSize()     
        self.cell_length = pcraster_clone.cellSize()
        self.x_max       = self.x_min + self.num_of_cols * self.cell_width
        self.y_min       = self.y_max - self.num_of_rows * self.cell_length

        # projection info 
        # - TODO: I need the following info from Yvonne:
        self.projection_long_name = 'unknown'
        self.projection_EPSG_code = 'unknown'
        self.projection_proj4_params = 'unknown'
        self.projection_grid_mapping_name = 'unknown'
        if projection_info == "longlat":
        # - an example for latlon
            self.projection_long_name = 'wgs84'
            self.projection_EPSG_code = 'EPSG:4326'
            self.projection_proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
            self.projection_grid_mapping_name = 'latitude_longitude'

        # let users decide what their preference regarding latitude order 
        self.netcdf_y_orientation_follow_cf_convention = netcdf_y_orientation_follow_cf_convention
        
        # cf convention
        self.cf_convention = ''
        if self.netcdf_y_orientation_follow_cf_convention: self.cf_convention = 'CF-1.4'
        # - Edwin assume the convention still unknown.
        self.cf_convention = 'unknown' 

        # netcdf format and zlib setup 
        self.netcdf_format = netcdf_format
        self.zlib          = zlib
        
    def set_netcdf_attributes(self, netcdf_setup_dictionary = None):

        attributeDictionary                  = {}
        if netcdf_setup_dictionary != None:
            attributeDictionary['institution' ]  = netcdf_setup_dictionary['institution']
            attributeDictionary['title'       ]  = netcdf_setup_dictionary['title'      ]
            attributeDictionary['source'      ]  = netcdf_setup_dictionary['source'     ]
            attributeDictionary['references'  ]  = netcdf_setup_dictionary['references' ]
            attributeDictionary['description' ]  = netcdf_setup_dictionary['description']
            attributeDictionary['created by'  ]  = netcdf_setup_dictionary['created by' ]
            attributeDictionary['notes'       ]  = netcdf_setup_dictionary['notes'      ]
        else:
            attributeDictionary['institution' ]  = "Utrecht University, Dept. of Physical Geography"
            attributeDictionary['title'       ]  = "Groundwater simulation (Egmond aan Zee)"
            attributeDictionary['source'      ]  = "Model scripts were written by Edwin H. Sutanudjaja (E.H.Sutanudjaja@uu.nl)."
            attributeDictionary['references'  ]  = "Smit et al. (in prep.)"
            attributeDictionary['description' ]  = "Groundwater simulation (Egmond aan Zee)"
            attributeDictionary['created by'  ]  = "The model run was performed by Yvonne Smit (Y.Smit@uu.nl). Model scripts were written by Edwin H. Sutanudjaja (E.H.Sutanudjaja@uu.nl)."
            attributeDictionary['notes'       ]  = ""
                                             
        attributeDictionary["history"     ]  = 'created on ' + datetime.datetime.today().isoformat(' ')
        attributeDictionary["date_created"]  = datetime.datetime.today().isoformat(' ')
                                             
        attributeDictionary["Conventions" ]  = self.cf_convention
        
        return attributeDictionary

    def create_netcdf_file(self, netcdf_file_name, netcdf_setup_dictionary):

        #~ self.x_min       = pcraster_clone.west()
        #~ self.y_max       = pcraster_clone.north()
        #~ self.num_of_cols = pcraster_clone.nrCols()
        #~ self.num_of_rows = pcraster_clone.nrRows()
        #~ self.cell_width  = pcraster_clone.cellSize()     
        #~ self.cell_length = pcraster_clone.cellSize()
        #~ self.x_max       = self.x_min + self.num_of_cols * self.cell_width
        #~ self.y_min       = self.y_max - self.num_of_rows * self.cell_length

        # cell centres coordinates
        x_coordinates = np.linspace(self.x_min + 0.5*self.cell_width , self.x_max + 0.5*self.cell_width , self.num_of_cols)
        y_coordinates = np.linspace(self.y_max - 0.5*self.cell_length, self.y_min + 0.5*self.cell_length, self.num_of_rows) 
        if self.netcdf_y_orientation_follow_cf_convention: y_coordinates = y_coordinates[::-1]

        # prepare the file
        ncFileName = netcdf_file_name
        rootgrp = nc.Dataset(ncFileName, 'w', format = self.netcdf_format)

        # create dimensions - time is unlimited, y and x are fixed
        rootgrp.createDimension('time', None)
        rootgrp.createDimension('y', len(y_coordinates) )
        rootgrp.createDimension('x', len(x_coordinates))
        # - dimension for time bounds
        rootgrp.createDimension("nv", 2) # 
        
        # time
        date_time = rootgrp.createVariable('time', 'f8', ('time',))
        date_time.units = 'Days since 1901-01-01' 
        date_time.calendar = 'gregorian'
        date_time.standard_name = 'time'
        date_time.long_name = 'Days since 1900-01-01 00:00:00'
        date_time.bounds = 'time_bounds'
        
        # time bounds
        time_bounds = rootgrp.createVariable('time_bounds', 'f8', ('time', 'nv',))
        time_bounds.units = 'Days since 1900-01-01 00:00:00'
        
        # y ('latitude', vertical diection)
        y = rootgrp.createVariable('y', 'f4', ('y',))
        y.long_name = 'y'
        y.units = 'metres'
        y.standard_name = 'y'

        # x ('longitude', horizontal direction)
        x = rootgrp.createVariable('x', 'f4', ('x',))
        x.long_name = 'x'
        x.units = 'metres'
        x.standard_name = 'x'

        # set latitude and and longitude values
        y[:] = y_coordinates
        x[:] = x_coordinates

        # projection info 
        projection = rootgrp.createVariable('projection', 'c')
        projection.long_name = self.projection_long_name
        projection.EPSG_code = self.projection_EPSG_code
        projection.proj4_params = self.projection_proj4_params
        projection.grid_mapping_name = self.projection_grid_mapping_name
        
        # set netcdf attribute information
        attributeDictionary = self.set_netcdf_attributes(netcdf_setup_dictionary)
        for k, v in attributeDictionary.items(): setattr(rootgrp,k,v)

        # sync and close the file
        rootgrp.sync()
        rootgrp.close()

    def create_variable(self, ncFileName, varName, varUnit, fill_value = 1e20, longName = None, comment = None):

        rootgrp = nc.Dataset(ncFileName,'a')

        # short and long variable names
        shortVarName = varName
        longVarName  = longName
        if longVarName == None: longVarName = shortVarName
        # - comment
        if comment == None: comment = ''

        # creating the variable
        var = rootgrp.createVariable(shortVarName, 'f4', ('time', 'y', 'x',), fill_value, zlib = self.zlib)
        var.standard_name = shortVarName
        var.long_name = longVarName
        var.comment = comment
        
        # variable unit
        var.units = varUnit

        # sync and close the file
        rootgrp.sync()
        rootgrp.close()

    def data_to_netcdf(self, ncFileName, shortVarName, varField, timeBounds, timeStamp = None, posCnt = None):

        rootgrp = nc.Dataset(ncFileName, 'a')
        
        lowerTimeBound = timeBounds[0]
        upperTimeBound = timeBounds[1]
        if timeStamp == None: timeStamp = lowerTimeBound + (upperTimeBound - lowerTimeBound) / 2

        # time
        date_time = rootgrp.variables['time']
        if posCnt == None: posCnt = len(date_time)
        date_time[posCnt] = nc.date2num(timeStamp, date_time.units, date_time.calendar)
        
        # time bounds
        time_bounds = rootgrp.variables['time_bounds']
        time_bounds[posCnt, 0] = nc.date2num(lowerTimeBound, date_time.units, date_time.calendar)
        time_bounds[posCnt, 1] = nc.date2num(upperTimeBound, date_time.units, date_time.calendar)
        
        # flip variable if necessary (to follow cf_convention)
        if self.netcdf_y_orientation_follow_cf_convention: varField = np.flipud(varField)
        
        # the variable
        rootgrp.variables[shortVarName][posCnt,:,:] = varField

        rootgrp.sync()
        rootgrp.close()

