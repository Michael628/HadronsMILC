#! /usr/bin/env python3

import logging
import sys
import os
import re
import subprocess

from python_scripts import utils
import typing as t
from python_scripts.nanny import (
    config,
    fileio,
    todo_utils
)

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
        for v in param['run_params'].keys():
            filepath = re.sub(v, param['run_params'][v], filepath)

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
def good_output(step: str, cfgno: str, param: t.Dict) -> bool:

    try:
        job_config = config.get_job_config(param, step)
    except KeyError:
        raise NotImplementedError("Todo: return fallback branching code for legacy output checker.")

    run_config = config.get_run_config(param, job_config)
    run_config.series = cfgno.split('.')[0]
    run_config.cfg = cfgno.split('.')[1]
    task_config = job_config.tasks
    outfile_config_list = config.get_outfile_config(param)

    bad_files = fileio.find_bad_files(task_config, outfile_config_list, run_config)
    if bad_files:
        logging.warning(f"File `{bad_files[0]}` not found or not of correct file size.")
        return False
    else:
        return True


######################################################################
def check_jobs(YAML):
    """Process all entries marked Q in the todolist"""

    # Read primary parameter file
    param = utils.load_param(YAML)

    # Read the to-do file
    todo_file = param['nanny']['todo_file']
    lock_file = todo_utils.lock_file_name(todo_file)

    # First, just get a list of entries
    todo_utils.wait_set_todo_lock(lock_file)
    todo_list = todo_utils.read_todo(todo_file)
    todo_utils.remove_todo_lock(lock_file)
    entry_list = sorted(todo_list, key=todo_utils.key_todo_entries)

    # Run through the entries. The entry_list is static, but the
    # to-do file could be changing due to other proceses
    while len(entry_list) > 0:
        # Reread the to-do file (it might have changed)
        todo_utils.wait_set_todo_lock(lock_file)
        todo_list = todo_utils.read_todo(todo_file)

        index, cfgno, step = next_finished(param, todo_list, entry_list)
        if index == 0:
            todo_utils.remove_todo_lock(lock_file)
            continue

        step = step[:-1]
        # Mark that we are checking this item and rewrite the to-do list
        todo_list[cfgno][index] = step + "C"
        todo_utils.write_todo(todo_file, todo_list)
        todo_utils.remove_todo_lock(lock_file)

        if step not in param["job_setup"].keys():
            print("ERROR: unrecognized step key", step)
            sys.exit(1)

        # Check that the job completed successfully
        status = good_output(step, cfgno, param)
        sys.stdout.flush()

        # Update the entry in the to-do file
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
def main():

    # Parameter file

    YAML = "params.yaml"

    check_jobs(YAML)


############################################################
if __name__ == '__main__':
    main()
