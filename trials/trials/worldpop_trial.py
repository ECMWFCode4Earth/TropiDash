# (failed) trials on downloading tiff files by API form worldpop database
# trials using json
import requests
import json

pop = requests.get("https://www.worldpop.org/rest/data/pop")
jpop = json.loads(pop.content)
jpop["data"]

# pop = requests.get("https://www.worldpop.org/rest/data/pop/wpgp")
# jpop = json.loads(pop.content)
# jpop["data"][0]
requests.get("https://www.worldpop.org/rest/data/pop_density/pd_ic_1km").json()

requests.get("https://www.worldpop.org/rest/data/pop").json()

requests.get("https://www.worldpop.org/rest/data/pop/wpgp1km/data").json()


#to try
requests.get("https://www.worldpop.org/rest/data/pop/wpgp1km?id=24777").json()


#normalize raster values
with rasterio.open(os.path.join(cwd,"TropiDash\data\impacts\ppp_2020_1km_Aggregated_resampled_10km_sum_clipped_3402na_fix.tif"), "r+") as src:
    with rasterio.open("data/impacts/normalized_pop_v2.tif", 'w',  **src.profile) as dst:
                band = src.read(1)
                band[band == src.nodata] = -1000
                # band = (band - band.ravel().min()) / (band.ravel()[band.ravel()<src.nodata].max() - band.ravel().min())
                bmax = band.ravel()[band.ravel()<r.nodata].max()
                bmin = band.ravel()[band.ravel()>0].min()
                band = (band - bmin) / (bmax - bmin)
                band[band<0] = src.nodata
                dst.write(band,1)