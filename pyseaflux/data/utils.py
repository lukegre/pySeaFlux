"""
Utilities for data set homogenisation, metadata, and plotting
"""
import munch
from typing import Union
from .processesors_custom import *


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

    # Evaluate lambda functions
    config_dict = eval_lambdas(config_dict)

    # Convert to Munch object for attribute-style access
    return munch.munchify(config_dict)


def eval_lambdas(config):
    """
    Recursively evaluate lambda expressions in a configuration dictionary.

    Searches through a dictionary for string values that start with "lambda"
    and evaluates them into actual callable functions. This allows defining
    functions in YAML configuration files.

    Parameters
    ----------
    config : dict
        The configuration dictionary to process

    Returns
    -------
    dict
        The processed configuration dictionary with lambda strings converted to functions

    Examples
    --------
    >>> config = {'func': 'lambda x: x * 2', 'nested': {'func2': 'lambda x: x + 10'}}
    >>> processed = eval_lambdas(config)
    >>> processed['func'](5)
    10
    >>> processed['nested']['func2'](5)
    15

    Notes
    -----
    This function uses eval() which can pose security risks if the configuration
    comes from untrusted sources.
    """
    
    for key, value in config.items():
        if isinstance(value, dict):
            config[key] = eval_lambdas(value)
        elif isinstance(value, str) and value.startswith("lambda"):
            config[key] = eval(value)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                if isinstance(item, dict):
                    value[i] = eval_lambdas(item)
                elif isinstance(item, str) and item.startswith("lambda"):
                    value[i] = eval(item)
    return config


def munch_to_dict(munch_obj: Union[munch.Munch, dict]) -> dict:
    """
    Converts a Munch object to a dictionary

    Parameters
    ----------
    munch_obj : Union[munch.Munch, dict]
        A Munch object or dictionary to convert

    Returns
    -------
    dict
        A dictionary representation of the input object

    Examples
    --------
    >>> from munch import Munch
    >>> m = Munch(a=1, b=Munch(c=2))
    >>> d = munch_to_dict(m)
    >>> isinstance(d, dict)
    True
    >>> d
    {'a': 1, 'b': {'c': 2}}
    """

    if isinstance(munch_obj, munch.Munch):
        return munch_obj.toDict()
    else:
        return munch_obj
