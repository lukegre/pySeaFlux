import numpy as np
from functools import wraps


class UnitError(Exception):
    pass


def unit_checker(arr, func, lim, vname):
    """
    Raise an error if the unit_checker are outside the given limit
    """
    from numpy import array, any

    arr = array(arr, ndmin=1)
    outside = func(arr, lim)
    if any(outside):
        fname = func.__name__
        nouts = outside.sum()
        prep = "are" if nouts > 1 else "is"
        plural = "s" if nouts > 1 else ""
        msg = (
            f"{vname} has {nouts} measurement{plural} that "
            f"{prep} {fname} than the limit ({lim}). \nUse "
            "use pd.Series.where/clip to remove/limit this data."
        )
        raise UnitError(msg)


def temp_K(temp_K):
    from numpy import less, greater

    varname = "temperature (K)"
    unit_checker(temp_K, less, 271.15, varname)
    unit_checker(temp_K, greater, 318.15, varname)


def pres_atm(pres_atm):
    from numpy import greater, less

    varname = "Pressure (atm)"
    unit_checker(pres_atm, greater, 1.5, varname)
    unit_checker(pres_atm, less, 0.5, varname)


def CO2_mol(CO2_mol):
    from numpy import greater, less

    varname = "CO2 mole fraction (ppm)"
    unit_checker(CO2_mol, greater, 0.01, varname)
    unit_checker(CO2_mol, less, 50e-6, varname)


def salt(salt):
    from numpy import greater, less

    varname = "Salinity"
    unit_checker(salt, less, 5, varname)
    unit_checker(salt, greater, 50, varname)


def wind_ms(wind_ms):
    from numpy import greater, less

    varname = "Wind speed (m/s)"
    unit_checker(wind_ms, less, 0, varname)
    unit_checker(wind_ms, greater, 40, varname)


class MetaArray(np.ma.MaskedArray):
    """
    Builds on np.ma.masked_array but adds a meta data column.
    Unique meta data is added when two MetaArrays are operated on.
    """

    def __new__(cls, input_array, mask=None, meta=None, **options):
        if isinstance(input_array, MetaArray):
            if meta is not None:
                metaA = np.array(input_array.meta)
                metaB = np.array(meta)
                meta = metaA + metaB
            if mask is not None:
                mask = input_array.mask | mask
        obj = np.ma.masked_array(input_array, mask=mask, **options).view(cls)
        obj._optinfo['meta'] = meta
        return obj

    @property
    def meta(self):
        return np.array(self._optinfo.get('meta', None))

    @property
    def ignore_mask(self):
        out = self._optinfo.get('ignore_mask', False)
        print(out)
        return out

    def to_dict(self):
        return dict(
            data=self.data,
            mask=self.mask,
            meta=self.meta)

    def __repr__(self):
        """
        Literal string representation.
        """
        import builtins
        if self._baseclass is np.ndarray:
            name = 'array'
        else:
            name = self._baseclass.__name__

        prefix = 'meta_{}('.format(name)

        dtype_needed = (
            not np.core.arrayprint.dtype_is_implied(self.dtype) or
            np.all(self.mask) or
            self.size == 0
        )

        # determine which keyword args need to be shown
        keys = ['data', 'mask', 'meta', 'fill_value']
        if dtype_needed:
            keys.append('dtype')

        # array has only one row (non-column)
        is_one_row = builtins.all(dim == 1 for dim in self.shape[:-1])

        # choose what to indent each keyword with
        min_indent = 2
        if is_one_row:
            # first key on the same line as the type, remaining keys
            # aligned by equals
            indents = {}
            indents[keys[0]] = prefix
            for k in keys[1:]:
                n = builtins.max(min_indent, len(prefix + keys[0]) - len(k))
                indents[k] = ' ' * n
            prefix = ''  # absorbed into the first indent
        else:
            # each key on its own line, indented by two spaces
            indents = {k: ' ' * min_indent for k in keys}
            prefix = prefix + '\n'  # first key on the next line

        # format the field values
        reprs = {}
        reprs['data'] = np.array2string(
            self._insert_masked_print(),
            separator=", ",
            prefix=indents['data'] + 'data=',
            suffix=',')
        reprs['mask'] = np.array2string(
            self._mask,
            separator=", ",
            prefix=indents['mask'] + 'mask=',
            suffix=',')
        reprs['meta'] = np.array2string(
            self.meta,
            separator=", ",
            prefix=indents['meta'] + 'meta=',
            suffix=',')
        reprs['fill_value'] = repr(self.fill_value)
        if dtype_needed:
            reprs['dtype'] = np.core.arrayprint.dtype_short_repr(self.dtype)

        # join keys with values and indentations
        result = ',\n'.join(
            '{}{}={}'.format(indents[k], k, reprs[k])
            for k in keys
        )
        return prefix + result + ')'

    @staticmethod
    def _ignore_mask(func):
        def inner(*args):
            output = args[0].copy()
            empty = np.ndarray(output.shape, dtype=str)
            meta = [getattr(a, 'meta', empty) for a in args]
            meta = np.array([np.array(a, ndmin=1) for a in meta]).T
            meta = [''.join(set(a)) for a in meta]
            args = [np.array(a, ndmin=1) for a in args]
            result = func(*args)
            output.data[:] = result
            output._optinfo['meta'] = meta
            return output
        return inner

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        data = [np.ma.getdata(a) for a in inputs]
        empty = np.ndarray(self.shape, dtype=str)
        meta = np.array([getattr(a, 'meta', empty) for a in inputs]).T
        meta = [''.join(set(a)) for a in meta]
        result = ufunc(*data)
        output = self.copy()
        output.data[:] = result
        output._optinfo['meta'] = meta
        return output

    def __add__(self, other):
        return self._ignore_mask(np.ma.add)(self, other)

    def __sub__(self, other):
        return self._ignore_mask(np.ma.subtract)(self, other)

    def __div__(self, other):
        return self._ignore_mask(np.ma.divide)(self, other)

    def __mul__(self, other):
        return self._ignore_mask(np.ma.multiply)(self, other)

    def __mod__(self, other):
        return self._ignore_mask(np.ma.mod)(self, other)

    def __pow__(self, other):
        return self._ignore_mask(np.ma.power)(self, other)

    def __floordiv__(self, other):
        return self._ignore_mask(np.ma.floor_divide)(self, other)


def check_limits(*funcs, return_data=True, **limits):
    def is_outside(arr, limit, func_name='', var_name='', return_info=True):
        arr = np.array(arr, ndmin=1)

        mask_lo = arr < limit[0]
        mask_up = arr > limit[1]
        mask = mask_lo | mask_up

        if return_info:
            info = np.ndarray(arr.shape, dtype='O')
            info[:] = ''
            info[mask_lo] = f'{func_name}: {var_name} < {limit[0]}; '
            info[mask_up] = f'{func_name}: {var_name} > {limit[1]}; '
            return mask, info
        else:
            return mask

    def wrapper(f):
        @wraps(f)
        def inner(*args, return_data=return_data, **kwargs):
            from inspect import getfullargspec

            spec = getfullargspec(f)
            # extra_keys = [k for k in limits if k not in spec.args]
            output = []
            for i, key in enumerate(spec.args):
                if key not in limits:
                    continue
                if i >= len(args):
                    continue
                val = args[i]
                lim = limits[key]
                assert len(lim) == 2, f'input of limit ({key}) must be len 2'
                assert lim[0] < lim[1], f'limits ({key}) must be [min, max]'

                output += is_outside(val, lim, f.__name__, key),

            output = np.array(output)
            mask = output[:, 0].any(0)
            meta = output[:, 1].sum(0)

            if return_data:
                args = [np.array(a) for a in args]
                data = f(*args, **kwargs)
                return MetaArray(data, mask=mask, meta=meta)
            else:
                return mask, meta

        return inner

    if len(funcs) == 1:
        f = funcs[0]
        return wrapper(f)
    else:
        return wrapper
