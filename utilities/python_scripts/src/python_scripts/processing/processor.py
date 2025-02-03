#! /usr/bin/env python3
import python_scripts as ps
from python_scripts import utils
from python_scripts.processing import dataio

import logging
import pandas as pd
import numpy as np
import gvar as gv
import gvar.dataset as ds
import typing as t
from python_scripts.a2a import contract as a2a

ACTION_ORDER = ['build_high', 'average', 'sum', 'time_average', 'real', 'permkey_average', 'permkey_normalize',
                'normalize', 'index', 'drop', 'gvar']


def stdjackknife(series: pd.Series) -> pd.Series:
    """Builds an array of standard deviations,
    which you can take the mean and standard deviation
    from for the error and the error on the error
    """
    marray = np.ma.array(series.to_numpy(), mask=False)
    array_out = np.empty_like(marray)

    marray.mask[0] = True
    for i in range(len(array_out)):
        marray.mask[:] = False
        marray.mask[i] = True
        array_out[i] = marray.std()

    return pd.Series(array_out, index=series.index)


def group_apply(df: pd.DataFrame, func: t.Callable, data_col: str,
                ungrouped_cols: t.List, invert=False) -> pd.DataFrame:

    all_cols = list(df.index.names) + list(df.columns)

    if not invert:
        ungrouped = ungrouped_cols + [data_col]
        grouped = [x for x in all_cols if x not in ungrouped]
    else:
        grouped = ungrouped_cols
        ungrouped = [x for x in all_cols if x not in grouped] + [data_col]

    df_out = df.reset_index().groupby(by=grouped)[
        ungrouped
    ].apply(func)

    while None in df_out.index.names:
        df_out = df_out.droplevel(df_out.index.names.index(None))

    return df_out


def gvar(df: pd.DataFrame) -> pd.DataFrame:

    cfgs = len(df)
    j_knife = group_apply(df, stdjackknife, ['series', 'cfg'])

    def ds_avg(series: pd.Series, scale: float = 1.0) -> pd.Series:
        array = ds.avg_data(series.to_numpy()) / scale
        return pd.Series(
            array,
            index=series.index.droplevel(['series', 'cfg']).drop_duplicates()
        )

    df_out = group_apply(j_knife, lambda x: ds_avg(
        x, float(cfgs - 1)), ['series', 'cfg'])
    df_out['noise'] = df_out.pop('corr')
    df_out['signal'] = group_apply(df, ds_avg, ['series', 'cfg'])['corr']
    df_out['nts'] = np.divide(df_out['noise'], df_out['signal'])
    return df_out


def buffer(df: pd.DataFrame, data_col: str, key_index: str) -> gv.BufferDict:

    buff = gv.BufferDict()

    nt = df.index.get_level_values('dt').nunique()

    labels_dt_last = sorted(df.index.names,
                            key=lambda x: 0 if x == 'dt' else -1)

    if key_index in df.columns:
        group_param = {'by': key_index}
    else:
        group_param = {'level': key_index}

    for key, xs in df.groupby(**group_param):
        buff[key] = xs.reorder_levels(labels_dt_last)[data_col] \
            .to_numpy().reshape((-1, nt))

    return buff

    # Shaun example code for dicts:
    # dset = gv.BufferDict()
    # dset['local'] = localArray
    # dset['onelink'] = onelinkArray
    # dsetGvar = ds.avg_data(dset)
    # localMinusOnelink = dsetGvar['local'] - dsetGvar['onelink']


# def build_high(df: pd.DataFrame, data_col) -> pd.DataFrame:

#     high = df.xs('ama', level='dset').sort_index()[data_col] \
#         - df.xs('ranLL', level='dset').sort_index()[data_col]
#     high = high.to_frame(data_col)
#     high['dset'] = 'high'
#     high.set_index('dset', append=True, inplace=True)
#     high = high.reorder_levels(df.index.names)

#     return pd.concat([df, high])


def drop(df, data_col, *args):
    for key in args:
        assert isinstance(key, str)

        if key in df.index.names:
            df.reset_index(key, drop=True, inplace=True)
        elif key in df.columns:
            _ = df.pop(key)
        else:
            raise ValueError(f'Drop Failed - No index or column `{key}` found.')
    return df


def index(df, data_col, *args):

    indices = [i for i in args]
    assert all([isinstance(i, str) for i in indices])

    if indices:
        if 'series.cfg' in indices:

            series: pd.DataFrame
            cfg: pd.DataFrame
            for key in ['series', 'cfg']:
                if key in df.index.names:
                    df.reset_index(key, inplace=True)

            series = df.pop('series')
            cfg = df.pop('cfg')

            df['series.cfg'] = series + '.' + cfg

            if 'series.cfg' in df.index.names:
                df.reset_index('series.cfg', drop=True, inplace=True)

        df.reset_index(inplace=True)
        df.set_index(indices, inplace=True)
        df.sort_index(inplace=True)
    return df


def real(df, data_col, apply_real: bool = True):
    if apply_real:
        df[data_col] = df[data_col].apply(np.real)
    return df


def normalize(df, data_col, divisor):
    return df['corr'].apply(lambda x: x / float(divisor)).to_frame()


def sum(df: pd.DataFrame, data_col, *sum_indices) -> pd.DataFrame:
    """Sums `data_col` column in `df` over columns or indices specified in `avg_indices`
        """
    return group_apply(df,lambda x: x[data_col].mean(), data_col,list(sum_indices))


def average(df: pd.DataFrame, data_col, *avg_indices) -> pd.DataFrame:
    """Averages `data_col` column in `df` over columns or indices specified in `avg_indices`,
    one at a time.
    """
    df_out = df
    for col in avg_indices:
        df_out = group_apply(df_out,lambda x: x[data_col].mean(), data_col,[col]).to_frame(data_col)

    return df_out

def permkey_split(df: pd.DataFrame, data_col, permkey_col: str = 'permkey') -> pd.DataFrame:
    if permkey_col in df.index.names:
        df.reset_index(permkey_col, inplace=True)

    df[permkey_col] = df[permkey_col].str.replace('e', '')
    df[permkey_col] = df[permkey_col].str.replace('v[0-9]+', ',', regex=True)
    df[permkey_col] = df[permkey_col].str.replace('w', '')
    df[permkey_col] = df[permkey_col].str.rstrip(',')
    df[permkey_col] = df[permkey_col].str.lstrip(',')
    key_len = df.iloc[0][permkey_col].count(',')
    assert all(df[permkey_col].str.count(',') == key_len)
    n_high = int(key_len + 1)

    df[[f'{permkey_col}{i}' for i in range(n_high)]] = df[permkey_col].str.split(',', expand=True)
    df.rename({f'{permkey_col}{n_high-1}':'noise_max'})
    df.drop(permkey_col, inplace=True, axis='columns')
    return df

def permkey_normalize(df: pd.DataFrame, data_col, permkey_col: str = 'permkey') -> pd.DataFrame:
    df_out = df
    if f'{permkey_col}0' not in df.columns:
        df_out = permkey_split(df_out, data_col, permkey_col)

    perm_cols = [x for x in df_out.columns if permkey_col in x]

    assert 'noise_max' in df.columns
    n_high_modes = df['noise_max'].max()
    n_unique_comb = df[perm_cols+['noise_max']].drop_duplicates().count()
    index0_max = (n_high_modes - len(perm_cols))
    df[data_col] = df[data_col]*n_unique_comb/index0_max

    def norm_func(x):
        x[data_col] = x[data_col] / x[data_col].count()
        return x

    for p in perm_cols:
        p_inv = [x for x in perm_cols if x != p]
        df_out = group_apply(df_out,norm_func, data_col,p_inv)

    return df_out


def permkey_average(df: pd.DataFrame, data_col, permkey_col: str = 'permkey') -> pd.DataFrame:
    df_out = permkey_split(df, data_col, permkey_col)
    perm_cols = [x for x in df_out.columns if permkey_col in x]

    return average(df_out, data_col, *perm_cols)


def time_average(df: pd.DataFrame, data_col: str, *avg_indices) -> pd.DataFrame:
    """Averages `data_col` column in `df` over columns or indices specified in `avg_indices`,
    one at a time.
    """
    assert len(avg_indices) == 2

    def apply_func(x):
        nt = int(np.sqrt(len(x)))
        assert nt ** 2 == len(x)
        corr = x[data_col].to_numpy().reshape((nt, nt))
        return pd.DataFrame({data_col: a2a.time_average(corr)}, index=pd.Index(range(nt), name='dt'))


    return group_apply(df,apply_func,data_col,list(avg_indices))

# def fold(df: pd.DataFrame, apply_fold: bool = True) -> pd.DataFrame:
#
#     if not apply_fold:
#         return df
#
#     assert len(df.columns) == 2
#
#     data_col = df.columns[-1]
#
#     array = df.sort_values('dt')[data_col].to_numpy()
#     nt = len(array)
#     folded_len = nt // 2 + 1
#     array[1:nt // 2] = (array[1:nt // 2] + array[-1:nt // 2:-1]) / 2.0
#
#     return pd.DataFrame(
#         array[:folded_len],
#         index=pd.Index(range(folded_len), name='dt'),
#         columns=[data_col]
#     )
#
#

def call(df, func_name, data_col, *args, **kwargs):

    func = globals().get(func_name, None)
    if callable(func):
        return func(df, data_col, *args, **kwargs)
    else:
        raise AttributeError(
            f"Function '{func_name}' not found or is not callable.")


def execute(df: pd.DataFrame, actions: t.Dict) -> pd.DataFrame:

    df_out = df
    data_col = actions.pop('data_col', 'corr')

    for key in sorted(actions.keys(), key=ACTION_ORDER.index):
        assert key in ACTION_ORDER
        param = actions[key]
        if isinstance(param, t.Dict):
            df_out = call(df_out, key, data_col, **param)
        elif isinstance(param, t.List):
            df_out = call(df_out, key, data_col, *param)
        else:
            if param:
                df_out = call(df_out, key, data_col, param)
            else:
                df_out = call(df_out, key, data_col)

    return df_out


def main(**kwargs):
    ps.setup()
    logging_level: str
    if kwargs:
        logging_level = kwargs.pop('logging_level', 'INFO')
        params = kwargs
    else:
        try:
            params = utils.load_param('params.yaml')['process_files']
        except KeyError:
            raise ValueError("Expecting `process_files` key in params.yaml file.")

        logging_level = params.pop('logging_level', 'INFO')

    logging.getLogger().setLevel(logging_level)

    ps.set_parallel_load(False)

    result = {}
    for key in params['run']:
        run_params = params[key]

        result[key] = dataio.main(**run_params)
        actions = run_params.get('actions', {})
        out_files = run_params.get('out_files', {})
        index = out_files.pop('index', None)

        if index:
            actions.update({'index': index})

        if 'actions' in run_params:
            result[key] = execute(result[key], run_params['actions'])

        if out_files:
            out_type = out_files['type']
            if out_type == 'dictionary':
                filestem = out_files['filestem']
                depth = int(out_files['depth'])
                dataio.write_dict(result[key], filestem, depth)
            elif out_type == 'dataframe':
                filestem = out_files['filestem']
                dataio.write_frame(result[key], filestem)
            else:
                raise NotImplementedError(
                    f"No support for out file type {out_type}."
                )

    return result


if __name__ == '__main__':
    result = main()
