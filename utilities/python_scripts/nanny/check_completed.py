#! /usr/bin/env python3

# Python 3 version

import sys, os, yaml, re, subprocess, copy
from todo_utils import *

# Check job completion.  For any completed jobs, mark the todo list
# C. DeTar

# Usage

# From the ensemble directory containing the todo file
# ../../scripts/check_completed.py

# Requires a todo file with a list of configurations to be processed
# In addition to the modules imported above, requires the following YAML files:
# ../scripts/params-allHISQ-plus5.yaml
# ../scripts/params-launch.yaml
# params-machine.yaml
# params-ens.yaml

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
        reply = subprocess.check_output(cmd, shell = True)
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
                time = a[5] + " " +  a[6] + " " + a[7]
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
    todo_list[key] = [ series_cfg, prec_tsrc, "X" ]
    print("Marked cfg", series_cfg, prec_tsrc, "completed")


######################################################################
def mark_checking_todo_entry(series_cfg, prec_tsrc, todo_list):
    """Update the todo_list, change status to X"""

    key = series_cfg + "-" + prec_tsrc
    todo_list[key] = [ series_cfg, prec_tsrc, "C" ]

#######################################################################
def decode_series_cfg(series_cfg):
    """Decode series, cfg, as it appeaers in the todo file"""
    return series_cfg.split(".")

#######################################################################
def decode_prec_tsrc(series_cfg):
    """Decode prec, tsrc, as it appeaers in the todo file
       Takes P.nn -> [P, nnn]"""
    return series_cfg.split(".")

######################################################################
def purge_props(param,series_cfg):
    """Purge propagators for the specified configuration"""

    print("Purging props for", series_cfg)
    series, cfg = decode_series_cfg(series_cfg)
    config_id = code_cfg(series, cfg)
    prop = param['files']['prop']
    subdirs = prop['subdirs'] + [ config_id ]
    remote_path = os.path.join(*subdirs)
    cmd = ' '.join([ "nohup", "/bin/rm -r", remote_path, "> /dev/null 2> /dev/null &"])
    print(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: can't remove props.  Error code", e.returncode, ".")

######################################################################
def purge_rands(param,series_cfg):
    """Purge random sources for the specified configuration"""

    print("Purging rands for", series_cfg)
    series, cfg = decode_series_cfg(series_cfg)
    config_iD = code_cfg(series, cfg)
    rand = param['files']['rand']
    subdirs = rand['subdirs'] + [ config_iD ]
    remote_path = os.path.join(*subdirs)
    cmd = ' '.join([ "nohup", "/bin/rm -r", remote_path, "> /dev/null 2> /dev/null &"])
    print(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: can't remove rands.  Error code", e.returncode, ".")

######################################################################
def tar_input_path(stream, s06Cfg, prec_tsrc):
    """Where the data and logs are found"""
    return os.path.join(stream, s06Cfg, prec_tsrc)

######################################################################
def purge_sym_links(param, job_case):
    """Purge symlinks for the specified job_id"""

    (stream, series, cfg, prec, tsrc, s06Cfg, tsrc_id, job_id, job_seq_no)  = job_case

    print("Purging symlinks for job", job_id)

    io = param['files']['out']
    logs_path = os.path.join(tar_input_path(stream, s06Cfg, tsrc_id), io['subdir'])
    cmd = ' '.join([ "find -P", logs_path, "-lname '?*Job'"+ job_id + "'*' -exec /bin/rm '{}' \;"])
    print(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: rmdir exited with code", e.returncode, ".")

######################################################################
def good_logs(param, job_case):
    """Check that the log files are complete"""

    (stream, series, cfg, prec, tsrc, s06Cfg, tsrc_id, job_id, job_seq_no)  = job_case
    prec_tsrc_config_id = [ prec, tsrc, series, cfg ]

    for step in range(param['job']['steprange']['high']):
        expect_file = out_file_name(stream, prec_tsrc_config_id, job_seq_no, '', "step" + str(step))
        log_path = os.path.join(stream, s06Cfg, tsrc_id, "logs", expect_file)
        try:
            stat = os.stat(log_path)
        except OSError:
            print("ERROR: Can't find expected output file", path)
            return False

        # Check for "RUNNING COMPLETED"
        entries = count_phrase(log_path, 'RUNNING COMPLETED')
        if entries < 1:
            print("ERROR: did not find 'RUNNING COMPLETED' in", log_path)
            return False

        # Check for nonconvergence, signaled by lines with "NOT"
        entries = count_phrase(log_path, "NOT")
        if entries > 0:
            print("WARNING: ", entries, "lines with 'NOT' suggesting nonconvergence")
#            return False
    
    # Passed these tests
    print("Output files OK")
    print("COMPLETE")

    return True

######################################################################
def check_path(param, job_key, cfgno, complain, file_key=None):
    """Complete the file path and check that it exists and has the correct size"""

    # Substute variables coded in file path
    if file_key != None:
        filepaths = []
    else:
        filepaths = [os.path.join( param['files']['home'], fv) for fk,fv in param['files'][job_key].items() if fk != 'good_size']

    good = True

    for filepath in filepaths:
        for v in param['lmi_param'].keys():
            filepath = re.sub(v, param['lmi_param'][v], filepath)

        series, cfg = cfgno.split('.')
        filepath = re.sub('SERIES', series, filepath)
        filepath = re.sub('CFG', cfg, filepath)

        try:
            file_size =  os.path.getsize(filepath)

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
    good =          check_path(param, 'lma', cfgno, True, 'ama')
    good = good and check_path(param, 'lma', cfgno, True, 'ranLL')

    if not good and 'ama_alt' in lma.keys():
        good =          check_path(param, 'lma', cfgno, True, 'ama_alt')
        good = good and check_path(param, 'lma', cfgno, True, 'ranLL_alt')
    
    return good

######################################################################
def good_a2a_local(param, cfgno):
    """Check that the A2A output looks OK"""

    return check_path(param, 'a2a_local', cfgno, True)

######################################################################
def good_a2a_onelink(param, cfgno):
    """Check that the A2A output looks OK"""

    return check_path(param, 'a2a_onelink', cfgno, True)

######################################################################
def good_contract_local(param, cfgno):
    """Check that the contrraction output looks OK"""

    return check_path(param, 'contract_local', cfgno, True)

######################################################################
def good_contract_onelink(param, cfgno):
    """Check that the contrraction output looks OK"""
                     
    return check_path(param, 'contract_onelink', cfgno, True)

######################################################################
def good_contract_py(param, cfgno):
    """Check that the contrraction output looks OK"""

    return check_path(param, 'contract_py', cfgno, True)

######################################################################
def move_failed_outputs(job_case):
    """Move failed output to temporary failure archive"""
    
    (stream, series, cfg, prec, tsrc, s06Cfg, tsrc_id, job_id, job_seq_no)  = job_case

    bad_output_path = tar_input_path(stream, s06Cfg, tsrc_id)
    fail_path = os.path.join(stream, s06Cfg, "fail", job_id)

    # Move the failed output
    cmd = " ".join(["mkdir -p ", fail_path, "; mv", bad_output_path, fail_path])
    print(cmd)
    try:
        subprocess.check_output(cmd, shell = True).decode("ASCII")
    except subprocess.CalledProcessError as e:
        status = e.returncode

######################################################################
def next_finished(param, todo_list, entry_list):
    """Find the next well-formed entry marked "Q" whose job is no longer in the queue"""
    a = ()
    nskip = 0
    while len(entry_list) > 0:
        cfgno = entry_list.pop(0)
        a = todo_list[cfgno]
        index, cfgno, step = find_next_queued_task(a)
        if index == 0:
            continue
        
        if step == "":
            nskip = 5  # So we stay this many Q entries away from another check_completed process
            
        # Skip entries to Avoid collisions with other check-completed processes
        if nskip > 0:
            nskip -= 1 # Count down from nonzero
            a = ()
            continue
    
        print("------------------------------------------------------------------------------------")
        print("Checking cfg", todo_list[cfgno]                     )
        print("------------------------------------------------------------------------------------")
        
        # Is job still queued?
        job_id = a[index+1]
        if job_still_queued(param, job_id):
            index = 0  # To signal no checking
            continue
        break

    return index, cfgno, step

######################################################################
def check_pending_jobs(YAML):
    """Process all entries marked Q in the todolist"""

    # Read primary parameter file
    param = load_param(YAML)

    # Read the todo file
    todo_file = param['nanny']['todo_file']
    lock_file = lock_file_name(todo_file)

    # First, just get a list of entries
    wait_set_todo_lock(lock_file)
    todo_list = read_todo(todo_file)
    remove_todo_lock(lock_file)
    entry_list = sorted(todo_list,key=key_todo_entries)

    # Run through the entries. The entry_list is static, but the
    # todo file could be changing due to other proceses
    while len(entry_list) > 0:
        # Reread the todo file (it might have changed)
        wait_set_todo_lock(lock_file)
        todo_list = read_todo(todo_file)

        index, cfgno, step = next_finished(param, todo_list, entry_list)
        if index == 0:
            remove_todo_lock(lock_file)
            continue

        step = step[:-1]
        # Mark that we are checking this item and rewrite the todo list
        todo_list[cfgno][index] = step + "C"
        write_todo(todo_file, todo_list)
        remove_todo_lock(lock_file)

        if step not in param["job_setup"].keys():
            print("ERROR: unrecognized step key", step)
            sys.exit(1)

        # Check that the job completed successfully
        sfx = ""
        status = True
        if step == "S":
            status = status and good_links(param, cfgno)
        if step == "E":
            status = status and good_eigs(param, cfgno)
        if step in ["L","A","B","N","D"]:
            status = status and good_lma(param, cfgno)
        if step in ["L","A","M","M1","D"]:
            status = status and good_a2a_onelink(param, cfgno)
        if step in ["L","A","M", "M2"]:
            status = status and good_a2a_local(param, cfgno)
        if step == "H":
            status = status and good_contract_local(param, cfgno)
        if step == "I":
            status = status and good_contract_onelink(param, cfgno)
        if step == "G":
            status = status and good_contract_py(param, cfgno)

        sys.stdout.flush()

        # Update the entry in the todo file
        wait_set_todo_lock(lock_file)
        todo_list = read_todo(todo_file)
        if status:
            todo_list[cfgno][index] = step+"X"
            print("Job step", step, "is COMPLETE")
        else:
            todo_list[cfgno][index] = step+"XXfix"
            print("Marking todo entry XXfix.  Fix before rerunning.")
        write_todo(todo_file, todo_list)
        remove_todo_lock(lock_file)

        # Take a cat nap (avoids hammering the login node)
        subprocess.check_call(["sleep", "1"])

############################################################
def main():

    # Parameter file

    YAML = "params.yaml"

    check_pending_jobs(YAML)


############################################################
main()
