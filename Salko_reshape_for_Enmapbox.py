# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:44:33 2025

@author: mabo6137
"""

import pandas as pd
import geopandas as gpd
import json
import os
from matplotlib import pyplot as plt

i=0
# %%
i+=1
sl_path = "C:\\Users\\mabo6137\\Documents\\hyperspectral\\Spectral_Library_Salko"
sl_fn = "Reflectance_spectra_of_peatland_vegetation_Finland_Estonia_smoothed.csv"
output_fn = "SpectralLibrary_Salko_test{}.csv".format(i)
# choose only arctic and sub-arctic sites
sites_list = ["Kilpisjarvi", "Lompolonjankka", "Halssiaapa", "Viiankiaapa"]
# read salko et al. spectral library, skip first three rows
sl_data = pd.read_csv(os.path.join(sl_path, sl_fn), header=3)

#%% write class to handle spectral library
#%%
class SpectralLibrary:
    ''' A spectral Library. of a very specific kind: table with wavelengths'''
    def __init__(self, ds, nanfilter=False, swirfilter = False, sites_list=[]):
        if nanfilter == True:
            ds = ds.dropna(axis='columns', how='all')
        if swirfilter == True:
            # class 1: noisy data at 2025-284 nm -> replace with NA
            ds.loc[:, 'wl2025' : 'wl2048'] = ds.loc[:, 'wl2025' : 'wl2048'].where(ds['SWIR_class'] !=1)
            # class 2: noisy data at 2025–2310 nm 
            ds.loc[:, 'wl2025' : 'wl2310'] = ds.loc[:, 'wl2025' : 'wl2310'].where(ds['SWIR_class'] !=2)
        # implement a filter for sites
        if sites_list: 
            ds = ds.loc[ds.Site.isin(sites_list)]
            ds.reset_index(drop=True)
            print(ds.Site)
        # assign data
        self.data = ds
        # make spectra accessible
        spectra = ds.filter(regex='wl')
        spectra = spectra.rename(columns = lambda x: x.removeprefix("wl"))
        self.spectra_subset = spectra
        # add some useful information
        self.data = self.data.rename(columns = lambda x: x.removeprefix("wl"))
        self.band_nr = self.spectra_subset.shape[1]
        self.band_wavelengths = [int(x) for x in list(self.spectra_subset.columns)]
        self.obs_nr =self.spectra_subset.shape[0]
        self.spectra_names = self.data.Plot_ID
        #self.spectra_subset.index 
        self.pft_subset = ds.filter(regex= 'PFT')
        self.pft_max = self.pft_subset.idxmax(axis='columns')
    def return_bandnames(self):
        return([str(wl) + " nm" for wl in self.band_wavelengths])
    
# %% open library file and create information     

sl = SpectralLibrary(sl_data, nanfilter = True, swirfilter = True, sites_list = sites_list)

# %% write ascii files for ENVI
csv_params = {'sep':"\t", 'encoding': 'ascii', 'na_rep':0}
# NA value 0 leads to "match" status.

ascii_outfile = os.path.join(sl_path, "SalkoSubArctic_spectra.txt")
pft_table_outfile = os.path.join(sl_path, "SalkoArctic_maxPFT.txt")
wavelengths_outfile = os.path.join(sl_path, "SalkoArctic_wavelengths.txt")

# save spectra
out = sl.spectra_subset.transpose()
out["wavelengths"] = sl.band_wavelengths
out.to_csv(ascii_outfile, **csv_params) 
# index has to be True!!

# save spectra names (each in new line)
names = (sl.pft_max + "_" + sl.spectra_names).str[4:]
names.transpose().to_csv(pft_table_outfile,
                              header=False, index=False,
                              **csv_params)

# save list of wavelengths
with open(wavelengths_outfile, "w") as file:
    for line in sl.band_wavelengths:
        file.write("%s\n" % line)

# TODO add spectra names as columns names also
# THIS WORKS!!! --> necessary input: wavelength column, columns = observations

#%% save as GeoJSON for QGIS EnMAPBox
output_fn = "SpectralLibrarySalko.geojson"

band_wavelengths_list = sl.band_wavelengths 
bad_bands_list = [1]*len(sl.band_wavelengths) # all bands are valid at this point


# create list of features
features_list=[]
for site in sl.data.index:
    measured_reflectance_list = list(sl.spectra_subset.loc[site])
    band_wavelengths_list = sl.band_wavelengths 
    PlotID = sl.spectra_names.loc[site]
    coordinates = [float(sl.data.Coordinate_x.loc[site]),
                   float(sl.data.Coordinate_y.loc[site])]
    _feature={
        "type": "Feature",
        "properties": {
 		"name": PlotID,
 		"profiles":{
         		"y": measured_reflectance_list, 
         		"x": band_wavelengths_list,
         		"xUnit": "Nanometers", 
         		"bbl": bad_bands_list
         		}
 		},
 	"geometry": {
 		"type": "Point",
 		"coordinates": coordinates
 		}
     }
    
    features_list.append(_feature)


_collection = {
"type": "FeatureCollection",
     "name": "SpectralLibrary_Salko2024_arctic",
     "description": "Peatland Vegetation in northern Finland",
     "features": features_list
 }



with open(os.path.join(sl_path, output_fn), 'w') as file:
    json.dump(_collection, file, ensure_ascii=False)

# %% cut actual spectral measurements
output_fn = "SpectralLibrary_Salko_forQGIS.csv"
# test: only use shortest wavelenghts

spectra = sl.filter(regex='wl')
spectra= spectra.rename(columns = lambda x: x.removeprefix("wl"))
cs= [int(item) for item in list(spectra.columns)] # to get columns as ints

# transform them to json tables
# spectrajson = spectra.to_json(orient='split')
# spectrajson = spectra.apply(lambda x: x.to_json(orient='split'), axis=1)

spectrajson = spectra.apply(lambda x: {"x":cs, "y":list(x.values), "xUnit":"Nanometers"}, axis=1)

# re-attach to table
output = sl.filter(items = ("Plot_ID", "Coordinate_x", "Coordinate_y"))
output['profiles'] = spectrajson
# save as csv again to use in QGIS

output.to_csv(os.path.join(sl_path, output_fn))

# %% try to create ENVI file
# header
envi_fn = "SpectralLibrary_Salko_envi.sli"
envi_header_fn = "SpectralLibrary_Salko_envi.hdr"
# get properties
hdr_entries = [
    "ENVI", 
    "samples = {}".format(sl.band_nr, nanfilter = True),
    "description = {Spectra of Finnish and Estonian peatlands, provided by Salko et al., https://doi.org/10.1016/j.ecoinf.2024.102772}",
    "lines = {}".format(sl.obs_nr),
    "bands = 1",
    "header_offset = 3", 
    "file type = ENVI Spectral Library",
    "data type = 4",
    "interleave = bsq",
    "byte order = 0",
    "wavelength units = Nanometers",
    "reflectance scale factor = 1.000000",
    "band_names = {}".format(set(sl.return_bandnames())),
    "wavelength = {}". format(set(sl.band_wavelengths)),
    "spectra_names = {}".format(set(sl.spectra_names)),
    "z plot range = {0,1}"
]
# TODO: improve spectra names using Warp(Band 1: Name) (some envi shenanigans)


''' more possibilities:
    "z plot titles = {Wavelength, Reflectance}
    z plot range = {0.00, 1.00}
'''

with open(os.path.join(sl_path, envi_header_fn), "w") as hdr:
    for line in hdr_entries:
        hdr.write("%s\n" % line)

