#! /usr/bin/env python3
import logging
import sys
import typing as t

import gvar as gv
import gvar.dataset as ds
import numpy as np
import pandas as pd

import python_scripts as ps
from python_scripts import utils
from python_scripts.a2a import contract as a2a
from python_scripts.nanny import config
from python_scripts.processing import dataio

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
                ungrouped_cols: t.List, invert: bool=False) -> pd.DataFrame:

    """Applies `func` to `data_col` in `df` grouped by `ungrouped_cols`.

    Parameters
    ----------
    df : pd.DataFrame
        
    func : t.Callable
        
    data_col : str
        
    ungrouped_cols : t.List
        
    invert : bool, optional
        

    Returns
    -------
    pd.DataFrame
        

    """
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


def gvar(df: pd.DataFrame, data_col: str, key_index: t.Optional[str] = None) -> gv.BufferDict:

    tvar = 't' if 't' in df.index.names else 'dt'

    nt = df.index.get_level_values(tvar).nunique()

    labels_dt_last = sorted(df.index.names,
                            key=lambda x: 0 if x == tvar else -1)

    result = {}
    if not key_index:
        result[data_col] = ds.avg_data(df.reorder_levels(labels_dt_last)[data_col] \
                                       .to_numpy().reshape((-1, nt)))
    else:
        if key_index in df.columns:
            group_param = {'by': key_index}
        else:
            group_param = {'level': key_index}

        for key, xs in df.groupby(**group_param):
            result[key] = ds.avg_data(xs.reorder_levels(labels_dt_last)[data_col] \
                .to_numpy().reshape((-1, nt)))

    return pd.DataFrame(result,index=pd.Index(range(nt),name=tvar))


def buffer(df: pd.DataFrame, data_col: str, key_index: str) -> gv.BufferDict:
    tvar = 't' if 't' in df.index.names else 'dt'

    buff = gv.BufferDict()

    nt = df.index.get_level_values(tvar).nunique()

    labels_dt_last = sorted(df.index.names,
                            key=lambda x: 0 if x == tvar else -1)

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

def permkey_split_old(df: pd.DataFrame, data_col, permkey_col: str = 'permkey') -> pd.DataFrame:
    df[permkey_col] = df[permkey_col].str.replace('e', '')
    df[permkey_col] = df[permkey_col].str.replace('v[0-9]+', ',', regex=True)
    df[permkey_col] = df[permkey_col].str.replace('w', '')
    df[permkey_col] = df[permkey_col].str.rstrip(',')
    df[permkey_col] = df[permkey_col].str.lstrip(',')
    key_len = df.iloc[0][permkey_col].count(',')
    assert all(df[permkey_col].str.count(',') == key_len)
    n_high = int(key_len + 1)

    df[[f'{permkey_col}{i}' for i in range(n_high)]] = df[permkey_col].str.split(',', expand=True)
    df.drop(permkey_col, inplace=True, axis='columns')
    return df

def permkey_split(df: pd.DataFrame, data_col, permkey_col: str = 'permkey') -> pd.DataFrame:
    if permkey_col in df.index.names:
        df.reset_index(permkey_col, inplace=True)

    if '_' not in df.iloc[0][permkey_col]:
        return permkey_split_old(df, data_col, permkey_col)

    df[permkey_col] = df[permkey_col].str.replace('(e_|_e)', '', regex=True)
    key_len = df.iloc[0][permkey_col].count('_')
    assert all(df[permkey_col].str.count('_') == key_len)
    n_high = int(key_len + 1)//2

    df[[f'{permkey_col}{i}' for i in range(n_high)]] = df[permkey_col].str.split('_', expand=True)[list(range(n_high))]
    df.drop(permkey_col, inplace=True, axis='columns')
    return df

def permkey_normalize(df: pd.DataFrame, data_col, permkey_col: str = 'permkey') -> pd.DataFrame:
    df_out = df
    if f'{permkey_col}0' not in df_out.columns:
        df_out = permkey_split(df_out, data_col, permkey_col)

    perm_cols = [x for x in df_out.columns if permkey_col in x]

    n_high_modes = df_out[f'{permkey_col}{len(perm_cols) - 1}'].astype(int).max() + 1
    n_unique_comb = df_out[perm_cols].drop_duplicates()[f'{permkey_col}0'].count()
    n_index_modes = n_high_modes - (len(perm_cols) - 1)
    df_out[data_col] = df_out[data_col] * n_unique_comb / n_index_modes
    for p in perm_cols[:-1]:
        n_index_modes += 1
        df_out[data_col] = df_out[data_col] / (n_index_modes - df_out[p].astype(int) - 1)

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
    tvar = 't' if 't' in df.index.names else 'dt'

    def apply_func(x):
        nt = int(np.sqrt(len(x)))
        assert nt ** 2 == len(x)
        corr = x[data_col].to_numpy().reshape((nt, nt))
        return pd.DataFrame({data_col: a2a.time_average(corr)}, index=pd.Index(range(nt), name=tvar))


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


def main(*args, **kwargs):
    ps.setup()
    logging_level: str

    if kwargs:
        logging_level = kwargs.pop('logging_level', 'INFO')
        proc_params = kwargs
    else:
        params = utils.load_param('params.yaml')
        if len(args) == 1 and isinstance(args[0],str):
            step = args[0]
            job_config = config.get_job_config(params, step)
            submit_config = config.get_submit_config(params, job_config)
            outfile_config = config.get_outfile_config(params)
            proc_params = config.processing_params(job_config, submit_config, outfile_config)
        else:
            try:
                proc_params = params['process_files']
            except KeyError:
                raise ValueError("Expecting `process_files` key in params.yaml file.")

        logging_level = proc_params.pop('logging_level', 'INFO')

    logging.getLogger().setLevel(logging_level)

    result = {}
    for key in proc_params['run']:
        run_params = proc_params[key]

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
    if len(sys.argv) == 2:
        step = sys.argv[1]
        result = main(step)
    else:
        result = main()

