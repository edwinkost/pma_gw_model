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
    
    def __init__(self, global_map = True, netcdf_y_orientation_follow_cf_convention = True, netcdf_format = 'NETCDF4', zlib = True):
        		
        # corner cordinates (lat/lon system)
        if global_map == True:
            self.x_min = -180.
            self.y_min =  -90.
            self.x_max =  180. 
            self.y_max =   90. 
        # TODO: Make the option for a non global map.     
        
        # projection info 
        self.projection_long_name = 'wgs84'
        self.projection_EPSG_code = 'EPSG:4326'
        self.projection_proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        self.projection_grid_mapping_name = 'latitude_longitude'

        # let users decide what their preference regarding latitude order 
        self.netcdf_y_orientation_follow_cf_convention = netcdf_y_orientation_follow_cf_convention
        
        # cf convention
        self.cf_convention = ''
        if self.netcdf_y_orientation_follow_cf_convention: self.cf_convention = 'CF-1.4'

        # netcdf format and zlib setup 
        self.netcdf_format = netcdf_format
        self.zlib          = zlib
        
    def set_netcdf_attributes(self, netcdf_setup_dictionary):

        attributeDictionary                  = {}
        attributeDictionary['institution' ]  = netcdf_setup_dictionary['institution']
        attributeDictionary['title'       ]  = netcdf_setup_dictionary['title'      ]
        attributeDictionary['source'      ]  = netcdf_setup_dictionary['source'     ]
        attributeDictionary['references'   ] = netcdf_setup_dictionary['references' ]
        attributeDictionary['description' ]  = netcdf_setup_dictionary['description']
        attributeDictionary['created by'  ]  = netcdf_setup_dictionary['created by' ]
                                             
        attributeDictionary["history"     ]  = 'created on ' + datetime.datetime.today().isoformat(' ')
        attributeDictionary["date_created"]  = datetime.datetime.today().isoformat(' ')
                                             
        attributeDictionary["Conventions" ]  = self.cf_convention
        
        return attributeDictionary

    def create_netcdf_file(self, netcdf_setup_dictionary):

        # cell centres coordinates (lat/lon - arc degree)
        deltaLon = netcdf_setup_dictionary['resolution_arcmin'] / 60.0
        deltaLat = deltaLon
        nrCols   = int((self.x_max - self.x_min) / deltaLon)
        nrRows   = int((self.y_max - self.y_min) / deltaLat)
        longitudes = np.linspace(self.x_min + 0.5*deltaLon, self.x_max + 0.5*deltaLon, nrCols)
        latitudes  = np.linspace(self.y_max - 0.5*deltaLat, self.y_min + 0.5*deltaLat, nrRows) 
        if self.netcdf_y_orientation_follow_cf_convention: latitudes = latitudes[::-1]

        # prepare the file
        ncFileName = netcdf_setup_dictionary['file_name']
        rootgrp = nc.Dataset(ncFileName, 'w', format = self.netcdf_format)

        # create dimensions - time is unlimited, lat and lon are fixed
        rootgrp.createDimension('time', None)
        rootgrp.createDimension('lat', len(latitudes) )
        rootgrp.createDimension('lon', len(longitudes))
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
        
        # latitude
        lat = rootgrp.createVariable('lat', 'f4', ('lat',))
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lat.standard_name = 'latitude'

        # longitude
        lon = rootgrp.createVariable('lon', 'f4', ('lon',))
        lon.standard_name = 'longitude'
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'

        # set latitude and and longitude values
        lat[:] = latitudes
        lon[:] = longitudes

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

    def create_variable(self, ncFileName, varName, varUnit, longName = None, comment = None):

        rootgrp = nc.Dataset(ncFileName,'a')

        # short and long variable names
        shortVarName = varName
        longVarName  = longName
        if longVarName == None: longVarName = shortVarName
        # - comment
        if comment == None: comment = ''

        # creating the variable
        var = rootgrp.createVariable(shortVarName, 'f4', ('time', 'lat', 'lon',), fill_value = vos.MV, zlib = self.zlib)
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

    def dictionary_of_data_to_netcdf(self, ncFileName, dataDictionary, timeBounds, timeStamp = None, posCnt = None):

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

        shortVarNameList = dataDictionary.keys()
        for shortVarName in shortVarNameList:
            
            varField = dataDictionary[shortVarName]
            # flip variable if necessary (to follow cf_convention)
            if self.netcdf_y_orientation_follow_cf_convention: varField = np.flipud(varField)
            
            # the variable
            rootgrp.variables[shortVarName][posCnt,:,:] = varField

        rootgrp.sync()
        rootgrp.close()
