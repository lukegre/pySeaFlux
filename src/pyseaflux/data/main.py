"""creates the SeaFlux data set."""
from . import aux_vars, pco2atm, spco2



def get_era5_data(era5_config: dict):
    
    def make_era5_urls(yyyy_mm: str, era5_url:dict)->xr.Dataset:
        from pyseaflux.data.download import make_urls
        
        variables = (
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "surface_pressure",
        )
    
        dt = pd.DateOffset(months=1, days=-1)
        t0 = pd.Timestamp(yyyy_mm)
        t1 = t0 + dt
        times = dict(start=t0, end=t1, freq='1D')
    
        urls = []
        for v in variables:
            urls += make_urls(era5_url, **times, variable=v)
        urls = np.sort(urls).tolist()
        return urls

    


if __name__ == "__main__":
    main()
