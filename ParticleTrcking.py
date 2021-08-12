# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import Idw, Contour, RadiusVariable, Kriging, Spline
import geopandas as gpd
import pandas as pd
import numpy as np
import flopy.modflow as mf
import flopy.utils.binaryfile as bf
import os

arcpy.env.overwriteOutput = True
date1 = pd.to_datetime("2011-01-01)
date2 = pd.to_datetime("2015-01-01)
#%% Read in model heads
base_hds_path = "S:/LAX/UVRGA.M001.SRVCS/Flow Model/Model/Transient/UVRGA_Transient_Alt2.hds"
base_hdobj = bf.HeadFile(base_hds_path)


sfq = flopy.utils.SfrFile('S:\\LAX\\UVRGA.M001.SRVCS\\Flow Model\\GWV\\SFR\\postprocessing\\UVRGA_Transient_streamflow.dat')
sfr_data=sfq.get_dataframe()
sfr_data=sfr_data.reset_index()
sp_set=np.zeros((len(sfr_data.index)))
for i in sfr_data.index:
    sp_set[i]=sfr_data.kstpkper[i][1]
    
sfr_data.insert(0,'SP',sp_set.astype(int)+1)

sfr_data=sfr_data.drop(columns=['width','Cond','gradient','k','i','j', 'Qprecip','Qet','Qovr'])


sp_map=pd.read_csv('S:/LAX/UVRGA.M001.SRVCS/Flow Model/Model/Dates_full_stress_period_map.csv', index_col=1)
sp_map.index=pd.to_datetime(sp_map.index)

#%% 