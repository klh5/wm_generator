import datacube
import sys
import xarray as xr
import numpy as np
import geopandas as gpd
from datacube.virtual import construct_from_yaml
from datacube.storage.masking import mask_invalid_data
from osgeo import gdal, osr

site = sys.argv[1]

grid = gpd.read_file('/scratch/a.klh5/mangrove_data/shapefiles/{}.shp'.format(site))

bounds = grid.total_bounds

xmin = bounds[0]
xmax = bounds[2]
ymin = bounds[3]
ymax = bounds[1]
crs = int(grid.crs['init'].split(':')[1])

srs = osr.SpatialReference()
srs.ImportFromEPSG(crs)

combined_ls_sref = construct_from_yaml("""
    collate:
        - product: ls4_arcsi_sref_global_mangroves
          measurements: [green, NIR, red]
        - product: ls5_arcsi_sref_global_mangroves   
          measurements: [green, NIR, red]
        - product: ls7_arcsi_sref_global_mangroves   
          measurements: [green, NIR, red]
        - product: ls8_arcsi_sref_global_mangroves   
          measurements: [green, NIR, red]
    """)
    
def getDataset(crs, xmin, xmax, ymin, ymax):
    
    "Fetch all data for the given area."
    
    print("Fetching data...")
    
    fetch_ds = combined_ls_sref.query(dc, x=(xmin, xmax), y=(ymin, ymax), crs='EPSG:{}'.format(crs), time=('2009-01-01', '2011-12-31'))

    grouped_ds = combined_ls_sref.group(fetch_ds, resolution=(-30, 30), output_crs='EPSG:{}'.format(crs))
            
    ds = combined_ls_sref.fetch(grouped_ds)    
    
    ds = mask_invalid_data(ds)
    
    print("Done.")

    return(ds)
    
def getNDWI(ds):
    
    print("Generating NDWI...")
    
    ds['NDWI'] = (ds.green - ds.NIR) / (ds.green + ds.NIR)
    
    avg = ds.NDWI.mean('time', skipna=True)
    
    print("Producing mask...")
    
    wm = xr.where(avg >= -0.3, 1, 0)
    
    avg = avg.fillna(255)
        
    wm = wm.astype('uint8')
        
    wm = wm.sortby("y", ascending=False) 
        
    print("Done.")
    
    return(wm)

def outputToFile(output):
    
    outfile = '{}.kea'.format(site)
    
    # Output to KEA file
    x_size = len(output.x.values)
    y_size = len(output.y.values)
    x_min = np.amin(output.x.values)
    y_max = np.amax(output.y.values)
    
    geo_transform = (x_min, 30, 0.0, y_max, 0.0, -30)
    
    driver = gdal.GetDriverByName('KEA')
    output_raster = driver.Create(outfile, x_size, y_size, 1, 1) # Only one band, byte data type since there are only 2 values
    
    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.SetGeoTransform(geo_transform)
    
    raster_band = output_raster.GetRasterBand(1)
    raster_band.SetNoDataValue(255)
    raster_band.SetDescription("mask")
    
    raster_band.WriteArray(output.values)
        
    output_raster = None

dc = datacube.Datacube()

ds = getDataset(crs, xmin, xmax, ymin, ymax)

wm = getNDWI(ds)

outputToFile(wm)
















