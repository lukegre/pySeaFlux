"""
Utilities for data set homogenisation, metadata, and plotting
"""

import munch
from loguru import logger


def load_yaml_config(fname: str, **kwargs) -> munch.Munch:
    """
    Load and parse a YAML configuration file with Jinja2 template support.

    This function reads a YAML file, renders it as a Jinja2 template using the provided
    keyword arguments, and then parses the rendered content. If no keyword arguments
    are provided, the function will perform a second pass, using values from the YAML
    file itself as template variables, which enables self-referential configurations.

    Parameters
    ----------
    fname : str
        Path to the YAML configuration file.
    **kwargs : Dict[str, Any]
        Keyword arguments to use for rendering the Jinja2 template.
        These values will replace corresponding template variables in the YAML file.

    Returns
    -------
    config_dict : dict
        A dictionary-like object with attribute-style access to the configuration properties

    Examples
    --------
    >>> config = load_yaml_config('config.yaml')
    >>> print(config.some_property)

    >>> config = load_yaml_config('config.yaml', custom_var='value')
    >>> print(config.some_property)
    """
    import jinja2
    import yaml

    # Load the file and render it as a Jinja2 template
    with open(fname, "r") as f:
        template = jinja2.Template(f.read())
        rendered = template.render(**kwargs)

    # Parse the rendered YAML content
    config_dict = yaml.safe_load(rendered)

    # If no keyword arguments were provided, perform a second pass
    # using the values from the YAML file itself as template variables
    if kwargs == {}:
        config_dict = load_yaml_config(fname, **config_dict)
        config_dict = validate_data_config(config_dict)

    # Convert to Munch object for attribute-style access
    return munch.munchify(config_dict)


def is_package_func(func_name):
    from voluptuous import Invalid
    from . import processors
    from . import custom_funcs

    source_libs = (processors, custom_funcs)
    for library in source_libs:
        func = getattr(library, func_name, None)
        if func is not None:
            break

    if func is None:
        raise Invalid(
            f"Could not find the function `{func_name}` in "
            f"{str([sl.__name__ for sl in source_libs])}. "
            "Please edit `custom_funcs` to add the function"
        )
    return func


def validate_data_config(config: dict):
    from voluptuous import Schema, Required, Optional, All, Any
    import datetime

    data_config_schema = Schema(
        {
            "release": object,
            "atm_co2": dict,
            str: {
                Required("name"): str,
                Optional("metadata"): dict,
                Required("urls"): [
                    {
                        Required("url"): str,
                        Optional("time"): {
                            Required("start"): datetime.date,
                            Required("end"): datetime.date,
                            Required("file_freq"): str,
                        },
                        Optional(str): All([str]),
                    }
                ],
                Required("fsspec_options"): {
                    Required("cache_storage"): str,
                    Optional("same_names"): bool,
                    Optional("cache_mapper"): All(str, is_package_func),
                    Optional(str): Any(str, dict),
                },
                Required("output_options"): {
                    Required("output_storage"): str,
                    Optional("output_freq"): str,
                    Optional("delete_raw_files"): bool,
                },
                Required("variables"): {str: str},
                Optional("processors"): [is_package_func],
            },
        }
    )

    return data_config_schema(config)


def check_delete_raw_files(catalog):
    from voluptuous import Schema, Invalid, Optional

    def warn_delete_raw_files(value):
        if value:
            raise Invalid(
                'Downloaded raw files for "{0}" will be deleted from '
                '"{1}" after each batch of final output is saved to {2}'
            )

    schema_catch_delete_raw_files = Schema(
        {
            str: object,
            str: {
                Optional("output_options"): {
                    str: object,
                    Optional("delete_raw_files"): warn_delete_raw_files,
                },
                str: object,
            },
        }
    )

    try:
        schema_catch_delete_raw_files(catalog)
    except Invalid as e:
        lvl0_path = catalog[e.path[0]]["fsspec_options"]["cache_storage"]
        lvl1_path = catalog[e.path[0]]["output_options"]["output_storage"]
        logger.warning(e.msg.format(e.path[0], lvl0_path, lvl1_path))
