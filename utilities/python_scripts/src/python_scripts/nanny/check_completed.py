#! /usr/bin/env python3

# Python 3 version

import sys
import os
import re
import subprocess
import todo_utils
import python_scripts.utils as utils
import typing as t
import python_scripts.config as config
import fileoutput

######################################################################
def job_still_queued(param, job_id):
    """Get the status of the queued job"""
    # This code is locale dependent

    scheduler = param['submit']['scheduler']

    user = os.environ['USER']
    if scheduler == 'LSF':
        cmd = " ".join(["bjobs", "-u", user, "|", "grep -w", job_id])
    elif scheduler == 'PBS':
        cmd = " ".join(["qstat", "-u", user, "|", "grep -w", job_id])
    elif scheduler == 'SLURM':
        cmd = " ".join(["squeue", "-u", user, "|", "grep -w", job_id])
    elif scheduler == 'INTERACTIVE':
        cmd = " ".join(["squeue", "-u", user, "|", "grep -w", job_id])
    elif scheduler == 'Cobalt':
        cmd = " ".join(["qstat", "-fu", user, "|", "grep -w", job_id])
    else:
        print("Don't recognize scheduler", scheduler)
        print("Quitting")
        sys.exit(1)

    # print(cmd)
    reply = ""
    try:
        reply = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        status = e.returncode
        # If status is other than 0 or 1, we have an squeue/bjobs problem
        # Treat job as unfinished
        if status != 1:
            print("ERROR", status, "Can't get the job status.  Skipping.")
            return True

    if len(reply) > 0:
        a = reply.decode().split()
        if scheduler == 'LSF':
            # The start time
            if a[2] == 'PEND':
                time = 'TBA'
            else:
                time = a[5] + " " + a[6] + " " + a[7]
            field = "start"
            jobstat = a[2]
        elif scheduler == 'PBS':
            time = a[8]
            field = "queue"
            jobstat = a[9]
        elif scheduler == 'SLURM':
            time = a[5]
            field = "run"
            jobstat = a[4]
        elif scheduler == 'INTERACTIVE':
            time = a[5]
            field = "run"
            jobstat = a[4]
        elif scheduler == 'Cobalt':
            time = a[5]
            field = "run"
            jobstat = a[8]
        else:
            print("Don't recognize scheduler", scheduler)
            print("Quitting")
            sys.exit(1)

        print("Job status", jobstat, field, "time", time)
        # If job is being canceled, jobstat = C (PBS).  Treat as finished.
        if jobstat == "C":
            return False
        else:
            return True

    return False


######################################################################
def mark_completed_todo_entry(series_cfg, prec_tsrc, todo_list):
    """Update the todo_list, change status to X"""

    key = series_cfg + "-" + prec_tsrc
    todo_list[key] = [series_cfg, prec_tsrc, "X"]
    print("Marked cfg", series_cfg, prec_tsrc, "completed")


######################################################################
def mark_checking_todo_entry(series_cfg, prec_tsrc, todo_list):
    """Update the todo_list, change status to X"""

    key = series_cfg + "-" + prec_tsrc
    todo_list[key] = [series_cfg, prec_tsrc, "C"]


######################################################################
def tar_input_path(stream, s06Cfg, prec_tsrc):
    """Where the data and logs are found"""
    return os.path.join(stream, s06Cfg, prec_tsrc)


######################################################################
def check_path(param, job_key, cfgno, complain, file_key=None):
    """Complete the file path and check that it exists and has the correct size
    """

    # Substute variables coded in file path
    if file_key is not None:
        filepaths = []
    else:
        filepaths = [
            os.path.join(param['files']['home'], fv)
            for fk, fv in param['files'][job_key].items()
            if fk != 'good_size'
        ]

    good = True

    for filepath in filepaths:
        for v in param['lmi_param'].keys():
            filepath = re.sub(v, param['lmi_param'][v], filepath)

        series, cfg = cfgno.split('.')
        filepath = re.sub('SERIES', series, filepath)
        filepath = re.sub('CFG', cfg, filepath)

        try:
            file_size = os.path.getsize(filepath)

            if file_size < param['files'][job_key]['good_size']:
                good = False
                if complain:
                    print("File", filepath, "not found or not of correct size")
        except OSError:
            good = False
            if complain:
                print("File", filepath, "not found or not of correct size")

    return good


######################################################################
def good_links(param, cfgno):
    """Check that the ILDG links look OK"""

    return check_path(param, 'fnlinks', cfgno, True)


######################################################################
def good_eigs(param, cfgno):
    """Check that the eigenvector file looks OK"""

    good = check_path(param, 'eigs', cfgno, False, 'eig')

    if not good:
        # Check file in subdir
        good = check_path(param, 'eigsdir', cfgno, True, 'eigdir')

    return good


######################################################################
def good_lma(param, cfgno):
    """Check that the LMA output looks OK"""

    lma = param['files']['lma']
    good = check_path(param, 'lma', cfgno, True, 'ama')
    good = good and check_path(param, 'lma', cfgno, True, 'ranLL')

    if not good and 'ama_alt' in lma.keys():
        good = check_path(param, 'lma', cfgno, True, 'ama_alt')
        good = good and check_path(param, 'lma', cfgno, True, 'ranLL_alt')

    return good


######################################################################
def good_a2a_local(param, cfgno):
    """Check that the A2A output looks OK"""

    return check_path(param, 'a2a_local', cfgno, True)


######################################################################
def good_meson(tasks, cfgno, param) -> bool:
    """Check that the A2A output looks OK"""

    check_path(param, 'a2a_onelink', cfgno, True)
    return

######################################################################


def good_contract_py(param, cfgno):
    """Check that the contrraction output looks OK"""

    return check_path(param, 'contract_py', cfgno, True)


######################################################################
# def get_high_modes_output(tasks: t.Dict, param: t.Dict) -> t.List[OutputCheck]:

#     epack_tasks = tasks['epack']

#     generate_eigs: bool
#     multifile: bool
#     save_evals: bool

#     generate_eigs = not epack_tasks['load']
#     multifile = epack_tasks.get('multifile', False)
#     save_evals = generate_eigs or epack_tasks.get('save_evals', False)

#     ret = []
#     if generate_eigs:
#         ret.append(files['eigdir' if multifile else 'eig'])
#     if save_evals:
#         ret.append(files['eval'])

#     return ret


# def get_meson_output(tasks: t.Dict, param: t.Dict) -> t.List[OutputCheck]:

#     meson_tasks = tasks['meson']

#     ret = [
#         OutputCheck(
#             filestem=files[meson],
#             good_size=files['good_size']
#         )
#         for meson in meson_tasks
#     ]

#     return ret


def get_task_outputs(task_config: t.List[config.ConfigBase], param: t.Dict) \
        -> t.List[str]:
    """Call `get_{task}_output` function for each task in list"""

    outfiles = []
    for cfg in task_config:
        if hasattr(cfg, 'output'):
            if callable(cfg.output):
                outfiles.append(cfg.output(**param))
            elif isinstance(cfg.output, t.List):
                for output in cfg.output:
                    outfiles.append(output.format(**param))
            else:
                raise ValueError(("Expecting output of all tasks to provide"
                                  "lists of strings to be formatted"))

    return sum(outfiles, [])

######################################################################


def next_finished(param, todo_list, entry_list):
    """Find the next well-formed entry marked "Q" whose job is no longer
    in the queue
    """
    a = ()
    nskip = 0
    while len(entry_list) > 0:
        cfgno = entry_list.pop(0)
        a = todo_list[cfgno]
        index, cfgno, step = todo_utils.find_next_queued_task(a)
        if index == 0:
            continue

        if step == "":
            nskip = 5

        # Skip entries to Avoid collisions with other check-completed processes
        if nskip > 0:
            nskip -= 1  # Count down from nonzero
            a = ()
            continue

        print("--------------------------------------------------------------")
        print("Checking cfg", todo_list[cfgno])
        print("--------------------------------------------------------------")

        # Is job still queued?
        job_id = a[index + 1]
        if job_still_queued(param, job_id):
            index = 0  # To signal no checking
            continue
        break

    return index, cfgno, step


######################################################################
def check_pending_jobs(YAML):
    """Process all entries marked Q in the todolist"""

    # Read primary parameter file
    param = utils.load_param(YAML)

    # Read the todo file
    todo_file = param['nanny']['todo_file']
    lock_file = todo_utils.lock_file_name(todo_file)

    # First, just get a list of entries
    todo_utils.wait_set_todo_lock(lock_file)
    todo_list = todo_utils.read_todo(todo_file)
    todo_utils.remove_todo_lock(lock_file)
    entry_list = sorted(todo_list, key=todo_utils.key_todo_entries)

    # Run through the entries. The entry_list is static, but the
    # todo file could be changing due to other proceses
    while len(entry_list) > 0:
        # Reread the todo file (it might have changed)
        todo_utils.wait_set_todo_lock(lock_file)
        todo_list = todo_utils.read_todo(todo_file)

        index, cfgno, step = next_finished(param, todo_list, entry_list)
        if index == 0:
            todo_utils.remove_todo_lock(lock_file)
            continue

        step = step[:-1]
        # Mark that we are checking this item and rewrite the todo list
        todo_list[cfgno][index] = step + "C"
        todo_utils.write_todo(todo_file, todo_list)
        todo_utils.remove_todo_lock(lock_file)

        if step not in param["job_setup"].keys():
            print("ERROR: unrecognized step key", step)
            sys.exit(1)

        # Check that the job completed successfully
        status = True
        if step == "S":
            status = status and good_links(param, cfgno)
        else:
            tasks: t.Dict[str, t.Any]
            try:
                tasks = param["job_setup"][step]['tasks']
            except KeyError:
                raise KeyError(("`tasks` not found in `job_setup` parameters"
                                f"for step `{step}`."))

            task_config: t.Dict[str, config.ConfigBase]
            task_config = {
                key: config.get_config(key)(task)
                for key, task in tasks.items()
            }
            run_config = config.get_config('lmi_param')(param['lmi_param'])
            task_outfiles = {
                key: output.get_outfile_iter(key)(task[key], param['files'], run_config)
            }

        # status = all([
            # call(f"good_{key}", tasks[key], cfgno, param)
            # for key in tasks.keys()
        # ])
        sys.stdout.flush()

        # Update the entry in the todo file
        todo_utils.wait_set_todo_lock(lock_file)
        todo_list = todo_utils.read_todo(todo_file)
        if status:
            todo_list[cfgno][index] = step + "X"
            print("Job step", step, "is COMPLETE")
        else:
            todo_list[cfgno][index] = step + "XXfix"
            print("Marking todo entry XXfix.  Fix before rerunning.")
        todo_utils.write_todo(todo_file, todo_list)
        todo_utils.remove_todo_lock(lock_file)

        # Take a cat nap (avoids hammering the login node)
        subprocess.check_call(["sleep", "1"])


############################################################
def call(func_name, *args, **kwargs):
    func = globals().get(func_name, None)
    if callable(func):
        return func(*args, **kwargs)
    else:
        raise AttributeError(
            f"Function '{func_name}' not found or is not callable.")


############################################################
def main():

    # Parameter file

    YAML = "params.yaml"

    check_pending_jobs(YAML)


############################################################
if __name__ == '__main__':
    main()
