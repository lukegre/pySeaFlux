# SeaFlux

A package to calculate all things air-sea CO<sub>2</sub> fluxes. 

## Installation 
I recommend the pip installation of the package
    
    pip install seaflux


## Basic usage

```python
import seaflux as sf
import xarray as xr

xds = xr.open_dataset('some_file_path.nc')
area = xds.pco2.area()
ice = xds.ice_frac.fillna(0)

url = "url from noaa mbl page (the actual download link)"
apco2 = sf.atm.nooa_mbl_to_pCO2(url, xds.pres_hPa, xds.tempC, xds.salt)
spco2 = sf.fCO2_to_pCO2(xds.fCO2, xds.tempC, xds.pres_hPa)
dpco2 = xds.spco2 - xds.apco2

# note that Wanninkhof (2014) is scales the CCMPv2 wind to 16cm/hr
kw = xds.gas_transfer_velocity.k_W14(xds.wind_ccmp_ms, xds.tempC)
sol = sf.solubility.solubility_weiss1974(xds.salt, xds.tempC - 273.15)

flux = kw * sol * dpco2 * (1 - ice)
# flux is in molC/m2/day
# to get to global (gC/yr) we multiply by days/yr, area (m2), g/molC
global_flux = (flux * area * 365 * 12.011).sum(['lat, 'lon'])

```
