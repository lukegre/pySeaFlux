import xarray as xr
from loguru import logger


def get_data(url_list, **kwargs):

    ds_raw = get_netcdfs(
        url_list, 
        kwargs.get('storage_options', {}),
        kwargs.get('open_kwargs', {}),
        kwargs.get('downloader', None),
    )
    
    ds_processed = process_data(ds_raw, **kwargs)

    return ds_processed


def process_data(ds, **kwargs):
    from .processors import run_processors, rename
    
    ds_rename = rename(ds, **kwargs.get('variables', {}))
    ds_processed = run_processors(ds_rename, kwargs.get('processors', []))

    return ds_processed


def _fsspec_open_local(url, storage_options)->list:
    import fsspec

    flist = fsspec.open_local(url, **storage_options)

    return flist


def _fsspec_open_files(urls, storage_options)->list:
    import fsspec

    if isinstance(urls, str):
        urls = [urls]

    olist = []
    for url in urls:
        flist = fsspec.open_files(url, **storage_options)
        olist += [f.open() for f in flist]

    return olist


def get_netcdfs(urls, storage_options, open_kwargs={}, downloader=None):

    if isinstance(urls, str):
        urls = [urls]

    has_zip_protocol = any(['zip:' in u for u in urls])
    
    if downloader is None:
        if has_zip_protocol:
            downloader = _fsspec_open_files
        else:
            downloader = _fsspec_open_local

    logger.info(f"Opening {len(urls)} urls")
    flist = downloader(urls, storage_options)

    ds = xr.open_mfdataset(flist, **open_kwargs)
    
    return ds


def make_urls(url_fmt, start, end, freq, **kwargs):
    import pandas as pd

    dates = pd.date_range(start, end, freq=freq)
    urls = [url_fmt.format(t=t, **kwargs) for t in dates]

    return urls
