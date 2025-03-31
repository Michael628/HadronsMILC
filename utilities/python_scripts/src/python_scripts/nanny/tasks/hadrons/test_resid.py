import itertools
import os.path
import re
import typing as t
from dataclasses import fields

import pandas as pd
from pydantic.dataclasses import dataclass

from python_scripts import Gamma, utils
from python_scripts.nanny import TaskBase
from python_scripts.nanny.config import OutfileList
from python_scripts.nanny.tasks.hadrons import SubmitHadronsConfig, templates
from python_scripts.nanny.tasks.hadrons.lmi import LMITask


# ============TestResid Task Configuration===========
@dataclass
class TestResidTask(TaskBase):
    # ============Epack===========

    epack: t.Optional[LMITask.EpackTask] = None
    high_modes: t.Optional[LMITask.HighModes] = None
    resid: t.Optional[t.List[float]] = None

    def __post_init__(self):
        assert isinstance(self.resid, t.List) and len(self.resid) > 0
        assert all([isinstance(r, float) for r in self.resid])

        assert not self.high_modes.skip_cg

        self.resid = sorted(self.resid, reverse=True)

    @classmethod
    def from_dict(cls, kwargs) -> "TestResidTask":
        """Creates a new instance of TestResidTaskConfig from a dictionary."""
        params = utils.deep_copy_dict(kwargs)
        for f in fields(cls):
            field_name = f.name
            field_type = f.type.__args__[0]
            if field_name in kwargs:
                params[field_name] = field_type.from_dict(kwargs[field_name])

        return cls(**params)

    @property
    def mass(self):
        """Returns list of labels for masses required by task components."""
        res = []

        if self.epack and not self.epack.load:
            res.append("zero")
        if self.high_modes:
            res += self.high_modes.mass

        return list(set(res))


# ============Functions for building params and checking outfiles===========
def input_params(
    tasks: TestResidTask,
    submit_config: SubmitHadronsConfig,
    outfile_config_list: OutfileList,
) -> t.Tuple[t.List[t.Dict], t.Optional[t.List[str]]]:
    submit_conf_dict = submit_config.string_dict()

    run_tsources = list(map(str, submit_config.tsource_range))

    gauge_filepath = outfile_config_list.gauge_links.filestem.format(**submit_conf_dict)
    gauge_fat_filepath = outfile_config_list.fat_links.filestem.format(
        **submit_conf_dict
    )
    gauge_long_filepath = outfile_config_list.long_links.filestem.format(
        **submit_conf_dict
    )

    modules = [
        templates.load_gauge("gauge", gauge_filepath),
        templates.load_gauge("gauge_fat", gauge_fat_filepath),
        templates.load_gauge("gauge_long", gauge_long_filepath),
        templates.cast_gauge("gauge_fatf", "gauge_fat"),
        templates.cast_gauge("gauge_longf", "gauge_long"),
    ]

    for mass_label in tasks.mass:
        name = f"stag_mass_{mass_label}"
        mass = str(submit_config.mass[mass_label])
        modules.append(
            templates.action(
                name=name, mass=mass, gauge_fat="gauge_fat", gauge_long="gauge_long"
            )
        )

    if tasks.high_modes and tasks.high_modes.solver == "mpcg":
        for mass_label in tasks.high_modes.mass:
            name = f"istag_mass_{mass_label}"
            mass = str(submit_config.mass[mass_label])
            modules.append(
                templates.action_float(
                    name=name,
                    mass=mass,
                    gauge_fat="gauge_fatf",
                    gauge_long="gauge_longf",
                )
            )

    if tasks.epack:
        epack_path = ""
        multifile = str(tasks.epack.multifile).lower()
        if tasks.epack.load or tasks.epack.save_eigs:
            epack_path = outfile_config_list.eig.filestem.format(**submit_conf_dict)

        # Load or generate eigenvectors
        if tasks.epack.load:
            modules.append(
                templates.epack_load(
                    name="epack",
                    filestem=epack_path,
                    size=submit_conf_dict["eigs"],
                    multifile=multifile,
                )
            )
        else:
            modules.append(templates.op("stag_op", "stag_mass_zero"))
            modules.append(
                templates.irl(
                    name="epack",
                    op="stag_op_schur",
                    alpha=submit_conf_dict["alpha"],
                    beta=submit_conf_dict["beta"],
                    npoly=submit_conf_dict["npoly"],
                    nstop=submit_conf_dict["nstop"],
                    nk=submit_conf_dict["nk"],
                    nm=submit_conf_dict["nm"],
                    multifile=multifile,
                    output=epack_path,
                )
            )

        # Shift mass of eigenvalues
        for mass_label in tasks.mass:
            if mass_label == "zero":
                continue
            mass = str(submit_config.mass[mass_label])
            modules.append(
                templates.epack_modify(
                    name=f"evecs_mass_{mass_label}", eigen_pack="epack", mass=mass
                )
            )

        if tasks.epack.save_evals:
            eval_path = outfile_config_list.eval.filestem.format(**submit_conf_dict)
            modules.append(
                templates.eval_save(
                    name="eval_save", eigen_pack="epack", output=eval_path
                )
            )

    if tasks.high_modes:

        def m1_eq_m2(x):
            return x[-2] == x[-1]

        def m1_ge_m2(x):
            return x[-2] >= x[-1]

        modules.append(templates.sink(name="sink", mom="0 0 0"))

        for tsource in run_tsources:
            modules.append(
                templates.noise_rw(
                    name=f"noise_t{tsource}",
                    nsrc=submit_conf_dict["noise"],
                    t0=tsource,
                    tstep=submit_conf_dict["time"],
                )
            )

        solver_labels = []
        if tasks.epack:
            solver_labels.append("ranLL")

        for res in tasks.resid:
            solver_labels.append(f"ama_{res}")

        high_path = outfile_config_list.high_modes.filestem

        for op in tasks.high_modes.operations:
            glabel = op.gamma.name.lower()
            mlabel = op.mass

            for tsource in run_tsources:
                guess = ""
                for slabel in solver_labels:
                    quark = f"quark_{slabel}_{glabel}_mass_{mlabel}_t{tsource}"
                    source = f"noise_t{tsource}"
                    solver = f"stag_{slabel}_mass_{mlabel}"

                    modules.append(
                        templates.quark_prop(
                            name=quark,
                            source=source,
                            solver=solver,
                            guess=guess,
                            gammas=op.gamma.gamma_string,
                            apply_g5="true",
                            gauge="" if op.gamma.local else "gauge",
                        )
                    )
                    guess = quark

            for tsource, slabel, m1label, m2label in filter(m1_ge_m2, quark_iter):
                quark1 = f"quark_{slabel}_{glabel}_mass_{m1label}_t{tsource}"
                quark2 = f"quark_{slabel}_pion_local_mass_{m2label}_t{tsource}"

                if m1label == m2label:
                    mass_label = f"mass_{m1label}"
                    mass_output = f"{submit_config.mass_out_label[m1label]}"
                else:
                    mass_label = f"mass_{m1label}_mass_{m2label}"
                    mass_output = (
                        f"{submit_config.mass_out_label[m1label]}"
                        f"_m{submit_config.mass_out_label[m2label]}"
                    )

                output = high_path.format(
                    mass=mass_output,
                    dset=slabel,
                    gamma_label=glabel,
                    tsource=tsource,
                    **submit_conf_dict,
                )

                modules.append(
                    templates.prop_contract(
                        name=f"corr_{slabel}_{glabel}_{mass_label}_t{tsource}",
                        source=quark1,
                        sink=quark2,
                        sink_func="sink",
                        source_shift=f"noise_t{tsource}_shift",
                        source_gammas=op.gamma.gamma_string,
                        sink_gammas=op.gamma.gamma_string,
                        apply_g5="true",
                        gauge="" if op.gamma.local else "gauge",
                        output=output,
                    )
                )

        for mass_label in tasks.high_modes.mass:
            for res in tasks.resid:
                if tasks.high_modes.solver == "rb":
                    modules.append(
                        templates.rb_cg(
                            name=f"stag_ama_mass_{mass_label}",
                            action=f"stag_mass_{mass_label}",
                            residual=str(res),
                        )
                    )
                else:
                    modules.append(
                        templates.mixed_precision_cg(
                            name=f"stag_ama_mass_{mass_label}",
                            outer_action=f"stag_mass_{mass_label}",
                            inner_action=f"istag_mass_{mass_label}",
                            residual=str(res),
                        )
                    )

            if "ranLL" in solver_labels:
                modules.append(
                    templates.lma_solver(
                        name=f"stag_ranLL_mass_{mass_label}",
                        action=f"stag_mass_{mass_label}",
                        low_modes=f"evecs_mass_{mass_label}",
                    )
                )

    schedule = []

    return modules, schedule


def catalog_files(
    task_config: TestResidTask,
    submit_config: SubmitHadronsConfig,
    outfile_config_list: OutfileList,
) -> pd.DataFrame:
    def generate_outfile_formatting():
        if task_config.epack:
            if task_config.epack.save_eigs:
                if task_config.epack.multifile:
                    yield (
                        {"eig_index": list(range(int(submit_config.eigs)))},
                        outfile_config_list.eigdir,
                    )
                else:
                    yield {}, outfile_config_list.eig
            if task_config.epack.save_eigs:
                yield {}, outfile_config_list.eval

        if task_config.high_modes:
            res = {"tsource": list(map(str, submit_config.tsource_range)), "dset": []}
            if task_config.epack:
                res["dset"].append("ranLL")
            if not task_config.high_modes.skip_cg:
                res["dset"].append("ama")

            for op in task_config.high_modes.operations:
                res["gamma_label"] = op.gamma.name.lower()
                res["mass"] = [submit_config.mass_out_label[m] for m in op.mass]
                yield res, outfile_config_list.high_modes

    def build_row(filepath: str, repls: t.Dict[str, str]) -> t.Dict[str, str]:
        repls["filepath"] = filepath
        return repls

    outfile_generator = generate_outfile_formatting()
    replacements = submit_config.string_dict()

    df = []
    for task_replacements, outfile_config in outfile_generator:
        outfile = outfile_config.filestem + outfile_config.ext
        filekeys = utils.formatkeys(outfile)
        replacements.update(task_replacements)
        files = utils.process_files(
            outfile,
            processor=build_row,
            replacements={k: v for k, v in replacements.items() if k in filekeys},
        )
        dict_of_rows = {
            k: [file[k] for file in files] for k in files[0] if len(files) > 0
        }

        new_df = pd.DataFrame(dict_of_rows)
        new_df["good_size"] = outfile_config.good_size
        new_df["exists"] = new_df["filepath"].apply(os.path.exists)
        new_df["file_size"] = None
        new_df.loc[new_df["exists"], "file_size"] = new_df[new_df["exists"]][
            "filepath"
        ].apply(os.path.getsize)
        df.append(new_df)

    df = pd.concat(df, ignore_index=True)

    return df


def bad_files(
    task_config: TestResidTask,
    submit_config: SubmitHadronsConfig,
    outfile_config_list: OutfileList,
) -> t.List[str]:
    df = catalog_files(task_config, submit_config, outfile_config_list)
    return list(df[(df["file_size"] >= df["good_size"]) != True]["filepath"])


def processing_params(
    task_config: TestResidTask,
    submit_config: SubmitHadronsConfig,
    outfile_config_list: OutfileList,
) -> t.Dict:
    proc_params = {"run": []}

    infile_stem = outfile_config_list.high_modes.filename
    outfile = outfile_config_list.high_modes.filestem
    filekeys = utils.formatkeys(infile_stem)
    outfile = outfile.replace("correlators", "dataframes")
    outfile = outfile.replace("_{series}", "")
    outfile = outfile.replace("_t{tsource}", "")
    outfile += ".h5"
    replacements = {
        k: v for k, v in submit_config.string_dict().items() if k in filekeys
    }
    replacements["tsource"] = list(map(str, submit_config.tsource_range))

    solver_labels = []
    if task_config.epack:
        solver_labels.append("ranLL")
    if task_config.high_modes:
        if not task_config.high_modes.skip_cg:
            solver_labels.append("ama")

        for op in task_config.high_modes.operations:
            gamma_label = op.gamma.name.lower()
            replacements["gamma_label"] = gamma_label
            for m, dset in itertools.product(op.mass, solver_labels):
                mass_label = submit_config.mass_out_label[m]
                file_label = f"{gamma_label}_{mass_label}_{dset}"
                proc_params["run"].append(file_label)
                replacements["mass"] = mass_label
                replacements["dset"] = dset

                h5_datasets = {
                    g: f"/meson/meson_{i}/corr"
                    for i, g in enumerate(op.gamma.gamma_list)
                }
                array_params = {
                    g: {"order": ["t"], "labels": {"t": f"0..{submit_config.time - 1}"}}
                    for g in h5_datasets.keys()
                }
                proc_params[file_label] = {
                    "logging_level": getattr(submit_config, "logging_level", "INFO"),
                    "load_files": {
                        "filestem": infile_stem,
                        "regex": {"series": "[a-z]", "cfg": "[0-9]+"},
                        "replacements": utils.deep_copy_dict(replacements),
                        "h5_params": {"name": "gamma", "datasets": h5_datasets},
                        "array_params": array_params,
                    },
                    "out_files": {"filestem": outfile, "type": "dataframe"},
                }

    return proc_params


def get_task_factory():
    return TestResidTask.from_dict
