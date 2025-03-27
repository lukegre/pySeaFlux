import pathlib
import pandas as pd
import xarray as xr
from loguru import logger


def get_data(data_config_entry:dict):
    from copy import deepcopy

    # make sure that nothing changes in the loaded config
    entry = deepcopy(data_config_entry)

    # make groups of urls - based on entry.urls[0].time.output_freq
    output_freq = entry.get('output_options', {}).get('output_freq', None)
    url_groups = make_urls(entry.urls, output_freq)

    ds_groups = ()
    for url_list in url_groups:
        ds_groups += get_data_from_url_list(url_list, entry),

    logger.info(f'Combing {len(ds_groups)} datasets into one group with time freq of {output_freq}')
    ds = xr.combine_by_coords(ds_groups, combine_attrs='override')

    return ds


def get_data_from_url_list(url_list:tuple[str], entry:dict):
    from dask.diagnostics import ProgressBar
    
    # path for the processed data 
    output_path = entry.get('output_options', {}).get('output_storage')
    if output_path is None:
        raise KeyError('you must provide `["output_options"]["output_storage"]`')
    
    flist = download_netcdfs(
        url_list, 
        entry.get('fsspec_options'),
        downloader=entry.get('downloader', None))
    
    ds_raw = xr.open_mfdataset(flist, **entry.get('open_kwargs', {}))
    ds_processed = process_dataset(ds_raw, **entry)

    sname = make_output_name(ds_processed, entry)
    spath = pathlib.Path(output_path) / sname
    spath.parent.mkdir(parents=True, exist_ok=True)

    if not spath.exists():
        logger.info(f"Computing and saving to: {spath.name}")
        with ProgressBar():
            ds = ds_processed.compute()
            ds.to_netcdf(spath)
    else:
        logger.info(f"File exists: {spath.name}")
        ds = xr.open_dataset(spath, chunks={})

    if entry.get('output_options', {}).get('delete_raw_files'):
        if not isinstance(flist[0], str):
            raise ValueError("You can only delete raw files when you're not working with zip files")
        logger.warning(
            f'The downloaded raw files (n={len(flist)}) will be removed. '
            'Note that directories are not removed.')
        flist = [pathlib.Path(f).unlink() for f in flist]


def process_dataset(ds, **kwargs):
    from .processors import run_processors, rename
    
    ds_rename = rename(ds, **kwargs.get('variables', {}))
    ds_processed = run_processors(ds_rename, kwargs.get('processors', []))

    return ds_processed


def download_netcdfs(urls, storage_options, downloader=None):
    
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
    
    return flist


def make_output_name(ds: xr.Dataset, item_config: dict, template:str="{name}-{var}-{time}.nc"):

    kwargs = {
        'name': item_config['name'],
    }
    
    time = ds.time.to_index()
    time_start = time.min().strftime('%Y%m%d')
    time_end = time.max().strftime('%Y%m%d')
    if time_start == time_end:
        kwargs['time'] = time_start
    else:
        kwargs['time'] = f"{time_start}_{time_end}"

    vars = list(ds.data_vars)
    if len(vars) <= 2:
        kwargs['var'] = '_'.join(vars)
    elif len(vars) > 2:
        kwargs['var'] = f"{vars[0]}_and_others"

    fname = template.format(**kwargs)
    return fname


def make_urls(urls_specs: list, output_freq=None)->tuple:

    def _make_urls_for_url_entry(entry, output_freq=None):
        from itertools import product
        from copy import deepcopy
        
        entry = deepcopy(entry)
    
        # pop the two manditory values from our dictionary
        url = entry.pop('url')
        times = _process_url_time(entry.pop('time', None), output_freq=output_freq)
    
        # remaining entries are converted to list for product
        keys = entry.keys()
        values = entry.values()
        # make sure that there aren't any nasty surprises
        _check_all_entries_lists(values)
    
        # convert the entries into a list of kwargs that we can iterate over
        kwarg_list = [dict(zip(keys, v)) for v in product(*values)]
    
        urls:tuple = ()
        for time_group in times:
            url_group = ()
            for t in time_group:
                for kw in kwarg_list:
                    url_group += url.format(t=t, **kw),
            urls += url_group,
        
        return urls
    
    def _process_url_time(time_props:dict, output_freq:str=None)->list[pd.DatetimeIndex]:
        from copy import deepcopy
        
        time_props = deepcopy(time_props)
        if time_props is None:
            return [['']]
        
        valid_start = time_props.pop('start')
        valid_end = time_props.pop('end')
        freq_out = (valid_end - valid_start) if output_freq is None else output_freq
        dt = time_props.pop('file_freq')
    
        times = []
        date_ranges = pd.date_range(valid_start, valid_end, freq=freq_out)
        date_ranges = zip(date_ranges[:-1], date_ranges[1:])
        for t0, t1 in date_ranges:
            times += pd.date_range(t0, t1, freq=dt, inclusive='left'),
    
        return times
    
    def _check_all_entries_lists(values:list)->None:
        
        all_list_type = all([isinstance(item, (list, tuple)) for item in values]), 
        if not all_list_type:
            message = (
                'Entries in a urls item must follow the following types:\n'
                '{ url: str,  time: {start, end, freq},  *other: list}.\n'
                f'You have key `{key}` that is `{type(value)}` that must be a list.')
            raise TypeError(message)
    
    urls = ()
    for url_spec in urls_specs:
        urls += _make_urls_for_url_entry(url_spec, output_freq=output_freq)
    return urls

