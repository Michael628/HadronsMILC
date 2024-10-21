#! /usr/bin/env python3
import logging
import pandas as pd
import numpy as np
import gvar as gv
import gvar.dataset as ds
import typing as t

ACTION_ORDER = ['remove', 'fold', 'average',
                'build_high', 'build_lmi', 'normalize', 'gvar']


def remove(df: pd.DataFrame, columns: t.List[str]) -> None:
    if isinstance(columns, str):
        columns = [columns]

    df = df.drop(columns=columns)


def fold(df: pd.DataFrame, time_col: str, data_col: str) -> None:

    data = pd.concat([df.pop(time_col), df.pop(data_col)])

    df.drop_duplicates(inplace=True)

    data = np.concatenate(data.to_numpy()).reshape((size, -1))

    halftime = int(data.shape[-1]/2)

    corr_out = np.zeros((size, halftime+1,), dtype=data.dtype)

    for i in range(1, halftime):
        corr_out[..., i] = (data[..., i]+data[..., -i])/2
    corr_out[..., 0] = data[..., 0]
    corr_out[..., halftime] = data[..., halftime]

    df['data'] = [corr_out[index] for index in range(size)]


def stdjackknife(df, corr):
    # This generates an array of standard deviations,
    # which you can take the mean and standard deviation
    # from for the error and the error on the error
    corr_jackknife_std = []
    Ncf = len(corr)
    for i in range(Ncf):
        b = np.delete(corr, i, axis=0)
        jknife_std = np.nanstd(b, axis=0)/np.sqrt(Ncf-1)
        corr_jackknife_std.append(jknife_std)
    return np.array(corr_jackknife_std)


def gvar(df, scale):

    def to_gvar(data):
        corr = np.concatenate(data).reshape((-1, len(data[0])))

        Ncfgs = float(len(corr))

        j_knife = self.stdjackknife(corr)*np.sqrt(scale)
        noise = ds.avg_data(j_knife)/(Ncfgs-1)
        signal = ds.avg_data(corr)  # , norm=Ncfgs*scaleVar)

        temp = {
            "gvar_type": pd.Series(['signal', 'noise', 'nts'],
                                   index=[0, 1, 2]),
            "data": pd.Series([
                signal,
                noise,
                np.divide(noise, signal)
            ], index=[0, 1, 2])
        }
        temp = pd.DataFrame(
            temp,
            columns=list(df.columns) +
            ['gvar_type']
        )
        return temp

    # Shaun example code for dicts:
    # dset = gv.BufferDict()
    # dset['local'] = localArray
    # dset['onelink'] = onelinkArray
    # dsetGvar = ds.avg_data(dset)
    # localMinusOnelink = dsetGvar['local'] - dsetGvar['onelink']

    try:
        scale = float(scale)
    except TypeError as e:
        scale = 1.0

    group = [
        x for x in df
        if x not in ['series', 'cfg', 'data']
    ]

    # if len(group) != 0:
    # df = \
    # df.groupby(
    # group,
    # dropna=False
    # )

    # df = pd.DataFrame(df['data'].apply(np.stack))
    # df = df['data'].apply(to_gvar)
    return
    # df = df.reset_index()[group+['data', 'gvar_type']]

    # result = None
    # for index, row in df.iterrows():
    #     new_data = to_gvar(row['data'])
    #     if result is None:
    #         result = new_data
    #     else:
    #         result = pd.concat([result, new_data], ignore_index=True)
    #     result = result.fillna(
    #         row[[x for x in row.index if x != 'data']])

    # # result.drop(columns=['series', 'cfg'])
    # df = result


def build_lmi(df: pd.DataFrame) -> None:
    df = df.sort_values(by=[
        x for x in df.columns if x != 'data'])
    result = df[df['dset'] == 'high'].copy()
    sum_data = result['data'].to_numpy() + \
        df[df['dset'] == 'a2aLL']['data'].to_numpy()

    result['data'] = sum_data
    result['dset'] = 'lmi'

    df = pd.concat([df, result],
                   ignore_index=True)


def build_high(df: pd.DataFrame) -> None:
    df = df.sort_values(by=[
        x for x in df.columns if x != 'data'])
    result = df[df['dset'] == 'ama'].copy()
    diff_data = result['data'].to_numpy() - \
        df[df['dset'] == 'ranLL']['data'].to_numpy()

    result['data'] = diff_data
    result['dset'] = 'high'

    df = pd.concat([df, result],
                   ignore_index=True)


def normalize(df, divisor, *args, **kwargs):
    df['data'] = df['data'].apply(
        lambda x: x.real/float(divisor))


def average(df: pd.DataFrame, avg_indices: t.List[str]) -> pd.DataFrame:
    """Averages `df` attribute over columns specified in `indices`
    """
    assert all([i in df.index.names for i in avg_indices])
    assert len(df.columns) == 1

    df_group_keys = [k for k in df.index.names if k not in avg_indices]

    df_out = df.reset_index().set_index(df_group_keys)['corr'].groupby(
        df_group_keys).transform('mean').to_frame('corr').drop_duplicates()

    return df_out


def call(df, func_name, *args, **kwargs):
    func = globals().get(func_name, None)
    if callable(func):
        return func(df, *args, **kwargs)
    else:
        raise AttributeError(
            f"Function '{func_name}' not found or is not callable.")


def execute(df, actions: t.Dict):

    df_out = df
    for key in sorted(actions.keys(), key=ACTION_ORDER.index):
        assert key in ACTION_ORDER
        param = actions[key]
        if isinstance(param, t.Dict):
            df_out = call(df_out, key, **param)
        elif isinstance(param, t.List):
            df_out = call(df_out, key, *param)
        else:
            df_out = call(df_out, key)

    return df_out
