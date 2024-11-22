#! /usr/bin/env python3
import logging
import pandas as pd
import numpy as np
import gvar.dataset as ds
import typing as t

ACTION_ORDER = ['build_high', 'fold', 'stdjackknife', 'average',
                'normalize', 'gvar']

GROUPED_ACTIONS = {
    'fold': ['dt'],
    'stdjackknife': ['series', 'cfg']
}


def gvar(df: pd.DataFrame):

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

    # Shaun example code for dicts:
    # dset = gv.BufferDict()
    # dset['local'] = localArray
    # dset['onelink'] = onelinkArray
    # dsetGvar = ds.avg_data(dset)
    # localMinusOnelink = dsetGvar['local'] - dsetGvar['onelink']


def build_high(df: pd.DataFrame) -> pd.DataFrame:

    high = df.xs('ama', level='dset').sort_index()['corr'] \
        - df.xs('ranLL', level='dset').sort_index()['corr']
    high = high.to_frame('corr')
    high['dset'] = 'high'
    high.set_index('dset', append=True, inplace=True)
    high = high.reorder_levels(df.index.names)

    return pd.concat([df, high])


def normalize(df, divisor, *args, **kwargs):
    return df['corr'].apply(lambda x: x.real / float(divisor)).to_frame()


def sum(df: pd.DataFrame, average=False, *avg_indices) -> pd.DataFrame:
    """Averages `df` attribute over columns specified in `indices`
    """
    logging.debug(df.index.names)
    logging.debug(avg_indices)
    assert all([isinstance(x, str) for x in avg_indices])
    assert all([i in df.index.names for i in avg_indices])
    assert len(df.columns) == 1

    df_group_keys = [k for k in df.index.names if k not in avg_indices]

    if average:
        df_out = df.groupby(level=df_group_keys).mean()
    else:
        df_out = df.groupby(level=df_group_keys).sum()

    return df_out


def average(df: pd.DataFrame, *avg_indices) -> pd.DataFrame:
    return sum(df, True, *avg_indices)


def fold(series: pd.Series) -> pd.Series:
    array = series.to_numpy()
    nt = len(array)

    array[1:nt // 2] = (array[1:nt // 2] + array[-1:nt // 2:-1]) / 2.0

    time_indices = series.index.get_level_values('dt') < (nt // 2 + 1)
    return pd.Series(
        array[:nt // 2 + 1],
        index=series[time_indices].index
    )


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


def group_apply(df: pd.DataFrame, func: t.Callable, ungrouped_cols: t.List) \
        -> pd.DataFrame:

    assert 'dt' in df.index.names
    assert 'corr' in df
    index_names = [x for x in df.index.names if x not in ungrouped_cols]

    return df['corr'].groupby(level=index_names, group_keys=False) \
                     .apply(func).to_frame()


def call(df, func_name, *args, **kwargs):
    func = globals().get(func_name, None)
    if callable(func):
        if func_name in GROUPED_ACTIONS:
            return group_apply(
                df, func, GROUPED_ACTIONS[func_name], *args, **kwargs)
        else:
            return func(df, *args, **kwargs)
    else:
        raise AttributeError(
            f"Function '{func_name}' not found or is not callable.")


def execute(df: pd.DataFrame, actions: t.Dict) -> pd.DataFrame:

    df_out = df
    for key in sorted(actions.keys(), key=ACTION_ORDER.index):
        assert key in ACTION_ORDER
        param = actions[key]
        if isinstance(param, t.Dict):
            df_out = call(df_out, key, **param)
        elif isinstance(param, t.List):
            df_out = call(df_out, key, *param)
        else:
            if param:
                df_out = call(df_out, key, param)
            else:
                df_out = call(df_out, key)

    return df_out
