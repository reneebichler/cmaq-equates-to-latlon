## ------------------------------------------------------------------------------------
## Description
## ------------------------------------------------------------------------------------

## Author: Renee Bichler, 2025
## Find me on GitHub: reneebichler

## Within this code, we're going to download EPA EQUATES data based on the CMAQ model.
## https://www.epa.gov/cmaq/equates
## YOU CAN DO IT!! ╰( ^o^)╮╰( ^o^)╮

print("##################################################################################")
print("##")
print("##                            Download EPA EQUATES Data")
print("##                                 ╰( ^o^)╮╰( ^o^)╮")
print("##")
print("##################################################################################")

## ------------------------------------------------------------------------------------
## Libraries
## ------------------------------------------------------------------------------------

import pyproj
import pycno
import pyrsig
import pandas as pd

## ------------------------------------------------------------------------------------
## Main
## ------------------------------------------------------------------------------------

## Name of the location
locname = 'conus'

## Set up a bounding box for CONUS
bbox = (-125.0, 24.5, -66.5, 49.5)

## Define the CMAQ key for NO2 concentrations
cmaqkey = 'cmaq.equates.conus.conc.NO2'

## Create a date sequence for the year 2019
dates = pd.date_range(start='2019-09-01', end='2019-12-31', freq='D').date

## Create a for loop to process each date
for bdate in dates:

    cur_date = pd.to_datetime(bdate)
    bdate = cur_date.strftime('%Y%m%d')

    print(f'Processing: 12US1_CMAQ_hourly_3D_{locname}_equates_conc_{bdate}')

    try:

        ## Get a pyrsig API object
        api = pyrsig.RsigApi(bdate = bdate, bbox = bbox, workdir = locname, gridfit = True)

        ## Get a column form CMAQ
        cds = api.to_ioapi(cmaqkey, bdate = bdate)

        ## Export cds to netcdf
        cds.to_netcdf(f'/DATA/EQUATES/12US1_CMAQ_hourly_3D_{locname}_equates_conc_{bdate}.nc')
    
    except Exception as e:
        print(f'Error processing {bdate}: {e}')
        continue

print("##################################################################################")
print("##")
print("##                                      Success!")
print("##                                  Job completed!")
print("##                                     ( ┘^o^)┘")
print("##")
print("##################################################################################")
