#! /usr/bin/env python3
import sys
import logging
import copy
from dataclasses import dataclass, field
import pandas as pd
import numpy as np
import gvar as gv
import gvar.dataset as ds

ACTION_ORDER = ['remove', 'fold', 'average',
                'build_high', 'build_lmi', 'normalize', 'gvar']


@dataclass
class DataProcess:

    _actions: dict = field(init=False, repr=False)

    def __post_init__(self):
        self._df = None

    @property
    def actions(self) -> dict:
        return self._actions

    @actions.setter
    def actions(self, actions: dict) -> None:
        if isinstance(actions, dict):
            self._actions = actions
        else:
            raise ValueError("actions must be assigned a dictionary.")

    def remove(self, columns, *args, **kwargs):
        if isinstance(columns, str):
            columns = [columns]
        try:
            self._df = self._df.drop(columns=columns)
        except KeyError as e:
            logging.info(e)

    def fold(self, *args, **kwargs):

        data = self._df.pop('data')
        size = len(data)

        data = np.concatenate(data.to_numpy()).reshape((size, -1))

        halftime = int(data.shape[-1]/2)

        corr_out = np.zeros((size, halftime+1,), dtype=data.dtype)

        for i in range(1, halftime):
            corr_out[..., i] = (data[..., i]+data[..., -i])/2
        corr_out[..., 0] = data[..., 0]
        corr_out[..., halftime] = data[..., halftime]

        self._df['data'] = [corr_out[index] for index in range(size)]

    def stdjackknife(self, corr):
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

    def gvar(self, scaleVar):

        def to_gvar(data):
            corr = np.concatenate(data).reshape((-1, len(data[0])))

            Ncfgs = float(len(corr))

            j_knife = self.stdjackknife(corr)*np.sqrt(scaleVar)
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
                columns=list(self._df.columns) +
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
            scaleVar = float(scaleVar)
        except TypeError as e:
            scaleVar = 1.0

        group = [
            x for x in self._df
            if x not in ['series', 'cfg', 'data']
        ]

        if len(group) != 0:
            self._df = \
                self._df.groupby(
                    group,
                    dropna=False
                )

        self._df = self._df['data'].apply(np.stack)

        return 

        result = None
        for index, row in self._df.iterrows():
            new_data = to_gvar(row['data'])
            if result is None:
                result = new_data
            else:
                result = pd.concat([result, new_data], ignore_index=True)
            result = result.fillna(
                row[[x for x in row.index if x != 'data']])

        # result.drop(columns=['series', 'cfg'])
        self._df = result

    def build_lmi(self, *args, **kwargs):
        self._df = self._df.sort_values(by=[
            x for x in self._df.columns if x != 'data'])
        result = self._df[self._df['dset'] == 'high'].copy()
        sum_data = result['data'].to_numpy() + \
            self._df[self._df['dset'] == 'a2aLL']['data'].to_numpy()

        result['data'] = sum_data
        result['dset'] = 'lmi'

        self._df = pd.concat([self._df, result],
                             ignore_index=True)

    def build_high(self, *args, **kwargs):
        self._df = self._df.sort_values(by=[
            x for x in self._df.columns if x != 'data'])
        result = self._df[self._df['dset'] == 'ama'].copy()
        diff_data = result['data'].to_numpy() - \
            self._df[self._df['dset'] == 'ranLL']['data'].to_numpy()

        result['data'] = diff_data
        result['dset'] = 'high'

        self._df = pd.concat([self._df, result],
                             ignore_index=True)

    def normalize(self, divisor, *args, **kwargs):
        self._df['data'] = self._df['data'].apply(
            lambda x: x.real/float(divisor))

    def average(self, indices: list[str]) -> None:
        """Averages `df` attribute over columns specified in `indices`
        """
        if len(indices) != 0:
            self._df = self._df.groupby(  # group rows for averaging
                [
                    x
                    for x in self._df
                    if x not in indices+['data']
                ],
                dropna=False
            )['data'].apply(  # average over 'data' column
                lambda x: np.mean(np.stack(x), axis=0)
            ).reset_index()

    def call(self, method_name, *args, **kwargs):
        method = getattr(self, method_name, None)
        if callable(method):
            return method(*args, **kwargs)
        else:
            raise AttributeError(
                f"Method '{method_name}' not found or is not callable.")

    def execute(self, df) -> dict:

        self._df = df.copy()
        for key in sorted(self.actions.keys(), key=ACTION_ORDER.index):
            param = self.actions[key]
            self.call(key, param)
        return self._df
