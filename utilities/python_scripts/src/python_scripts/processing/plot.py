import pandas as pd
from dataclasses import dataclass
import matplotlib.pyplot as plt
import gvar as gv


def plot_cols(df: pd.DataFrame) -> None:
    plt.figure(dpi=100, figsize=(15, 10))

    for col in df.columns:
        y = gv.mean(df[col].to_numpy())
        yerr = gv.sdev(df[col].to_numpy())
        plt.errorbar(df.index, y, yerr=yerr, marker="o", linestyle='None', capsize=7, markersize=7, label=col)
    # plt.yscale("log")
    plt.legend()
    plt.show()

@dataclass
class DataPlot:

    _plot_info = {
        "signal": "C",
        "noise": "\\sigma",
        "nts": "N/S"
    }

    @staticmethod
    def tfm(spacing, tsteps):
        return [i*spacing for i in range(tsteps)]

    @staticmethod
    def tovera(tsteps):
        return list(range(tsteps))

    @classmethod
    def plot(cls, df: pd.DataFrame, params=None):

        if params is None:
            params = {}

        spacing = params.get('spacing', '')
        label_func = params.get('label', None)

        df = df.sort_values(by='gvar_type')

        for plot_type, df_group in df.groupby('gvar_type'):
            fig = plt.figure(dpi=100, figsize=(15, 10))
            trange = None
            for _, row in df_group.iterrows():
                data = row['data']
                if trange is None:
                    if spacing:
                        trange = cls.tfm(float(spacing), len(data))
                    else:
                        trange = cls.tovera(len(data))

                if label_func is not None:
                    label = label_func(row)
                else:
                    label = row['dset']+row['setkey']
                y = gv.mean(data)
                yerr = gv.sdev(data)
                ebs = plt.errorbar(trange, y, yerr=yerr, marker="o",
                                   linestyle='None', capsize=7, markersize=7,
                                   label=f"{label}")
                ebs[-1][0].set_linewidth(0.3)
                # [eb.set_markeredgewidth(2) for eb in ebs[-2]]

                if plot_type == "nts":
                    plt.plot(trange, [1.0]*len(trange), linestyle="--")
                    plt.ylim([1e-7, 4])
            plt.yscale("log")
            plt.legend(fontsize=20)
            plt.title((
                f"{spacing+' fm' if spacing else ''} Correlator Comparison"
                f" ({plot_type})"), fontsize=20)
            plt.ylabel(f"${cls._plot_info[plot_type]}(t)$", fontsize=20)
            plt.xlabel("t (fm)", fontsize=20)
            # plt.xlim([-0.1,3.3])
            # plt.xlim([1,1.5])
            fig.set_facecolor('white')
            plt.show()
        # plt.savefig(f"plots/{gamma}_04fm_{dset}")
