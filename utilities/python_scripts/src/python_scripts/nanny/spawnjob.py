#! /usr/bin/env python3

# Python 3 version

import sys, os, yaml, re, subprocess
from todo_utils import *
from functools import reduce

from dict2xml import dict2xml as dxml
import importlib
import hadrons_templates

sys.path.append("../templates")

# Nanny script for managing job queues
# C. DeTar 7/11/2022

# Usage

# From the ensemble directory containing the todo file
# ../scripts/spawnjob.py

# Requires a todo file with a list of configurations to be processed

# The todo file contains the list of jobs to be done.
# The line format of the todo file is
# <cfgno> <task code and flag> <task jobid> <task code and flag> <task jobid> etc.
# Example: a.1170 SX 0 EX 2147965 LQ 2150955 A 0 H 0
# Where cfgno is tne configuration number in the format x.nnn where x is the series letter and nnn
# is the configuration number in the series
# The second letter in the task code and flag is the flag
# The task codes are
# S links smearing job
# E eigenvector generation job
# L LMA job
# M meson job
# A A2A + meson job
# H contraction local job
# I contraction one-link job
# They must be run in this sequence (but L and A can be concurrent)
# If flag is "X", the job has been finished
# If it is "Q", the job was queued with the given <jobid>
# If it is "C", the job is finished and is undergoing checking and tarring
# If the flag letter is empty, the job needs to be run

# Requires Todo_utils.py and params-launch.yaml with definitions of variables needed here

######################################################################
def count_queue( scheduler,  myjob_name_pfx ):
    """Count my jobs in the queue"""

    user = os.environ['USER']

    if scheduler == 'LSF':
        cmd = ' '.join(["bjobs -u", user, "| grep", user, "| grep ", myjob_name_pfx, "| wc -l"])
    elif scheduler == 'PBS':
        cmd = ' '.join(["qstat -u", user, "| grep", user, "| grep ", myjob_name_pfx, "| wc -l"])
    elif scheduler == 'SLURM':
        cmd = ' '.join(["squeue -u", user, "| grep", user, "| grep ", myjob_name_pfx, "| wc -l"])
    elif scheduler == 'INTERACTIVE':
        cmd = ' '.join(["squeue -u", user, "| grep", user, "| grep ", myjob_name_pfx, "| wc -l"])
    elif scheduler == 'Cobalt':
        cmd = ' '.join(["qstat -fu", user, "| grep", user, "| grep ", myjob_name_pfx, "| wc -l"])
    else:
        print("Don't recognize scheduler", scheduler)
        print("Quitting")
        sys.exit(1)

    nqueued = int(subprocess.check_output(cmd,shell=True))

    return nqueued

######################################################################
def next_cfgno_steps( max_cases, todo_list ):
    """Get next sets of cfgnos / job steps from the todo file"""

    # Return a list of cfgnos and indices to be submitted in the next job
    # All subjobs in a single job must do the same step

    step = "none"
    cfgno_steps = []
    for line in sorted(todo_list,key=key_todo_entries):
        a = todo_list[line]
        if len(a) < 2:
            print("ERROR: bad todo line format");
            print(a)
            sys.exit(1)
    
        index, cfgno, new_step = find_next_unfinished_task(a)
        if index > 0:
            if step == "none":
                step = new_step
            elif step != new_step:
                # Ensure only one step per job
                break
            cfgno_steps.append([cfgno, index])
            # We don't bundle the S (links) or H (contraction) steps
            if step in ['S','H','I']:
                break
        # Stop when we have enough for a bundle
        if len(cfgno_steps) >= max_cases:
            break

    ncases = len(cfgno_steps)
    
    if ncases > 0:
        print("Found", ncases, "cases...", cfgno_steps)
        sys.stdout.flush()

    return step, cfgno_steps

######################################################################
def sub_iOTemplate(param, io_stem_template):
    """Replace keys with values in io_stem
       e.g. converts lma-eEIGS-nNOISE to lma-e2000-n1"""

    s = io_stem_template

    for k in param["lmi_param"].keys():
        s = re.sub(k,  param['lmi_param'][k],  s)
    
    return s

######################################################################
def set_env(param, series, cfgno):
    """Set environment variables"""

    # Set environment parameters for the job script
    lmi_param = param['lmi_param']
    for key in lmi_param.keys():
        os.environ[key] = str(lmi_param[key])
    
    # These will be ignored in the bundled job script
    os.environ['SERIES'] = series
    os.environ['CFG']    = str(cfgno)

    # Compute starting time for loose and fine
    dt     = int(lmi_param['DT'])   # Spacing of source times

    cfg0   = int(param['precess']['loose']['cfg0'])  # Base configuration
    dcfg   = int(param['precess']['loose']['dcfg'])  # Interval between cfgnos
    tstep  = int(param['precess']['loose']['tstep']) # Precession interval
    t0loose = ( (int(cfgno) - cfg0)//dcfg * tstep ) % dt  # Starting time

    cfg0   = int(param['precess']['fine']['cfg0'])
    dcfg   = int(param['precess']['fine']['dcfg'])
    tstep  = int(param['precess']['fine']['tstep'])
    nt     = int(param['precess']['fine']['nt'])
    t0fine = ( t0loose + (int(cfgno) - cfg0)//dcfg * tstep * dt ) % nt

    os.environ['T0LOOSE']  = str(t0loose)
    os.environ['T0FINE']  = str(t0fine)


######################################################################
def make_inputs(param, step, cfgno_steps):
    """Create input XML files for this job"""

    ncases = len(cfgno_steps)
    INPUTXMLLIST = ""

    if step == 'S':

        # Special treatment for link smearing. No bundling. Self-generated input.
        if ncases > 0:
            print("WARNING: No bundling of smearing jobs")
            print("Will submit only one case")
        ( cfgno_series, _ ) = cfgno_steps[0]
        ( series, cfgno ) = cfgno_series.split(".")
        set_env(param, series, cfgno)

    else:

        for i in range(ncases):
            ( cfgno_series, _ ) = cfgno_steps[i]

            # Extract series and cfgno  a.1428 -> a 1428
            ( series, cfgno ) = cfgno_series.split(".")
            
            # Define common environment variables
            set_env(param, series, cfgno)
            
            # Name of the input XML file
            io_stem_template = param['job_setup'][step]['io']
            # Replace variables in io_stem if need be
            io_stem = sub_iOTemplate(param, io_stem_template)
            input_xml = io_stem + "-" + cfgno_series + ".xml"
        
            if 'param_file' in param['job_setup'][step]:
                # String naming one of the files in ../templates without the .py
                param_module = f"{param['job_setup'][step]['param_file']}_params"

                # import paramModule as pm
                pm = importlib.import_module(param_module)

                # Template dictionary fir the input XML
                templates = hadrons_templates.generate_templates()

                # Generate the input XML file
                fp = open("in/" + input_xml, "w")
                print(dxml(pm.build_params(**templates)), file=fp)
                fp.close()

            INPUTXMLLIST = INPUTXMLLIST + " " + input_xml

    os.environ['INPUTXMLLIST'] = INPUTXMLLIST

######################################################################
def submit_job(param, step, cfgno_steps, max_cases):
    """Submit the job"""

    ncases = len(cfgno_steps)

    job_script = param['job_setup'][step]['run']
    wall_time = param['job_setup'][step]['wall_time']

    layout = param['submit']['layout']
    basenodes = layout[step]['nodes']
    ppj = reduce((lambda x,y: x*y),layout[step]['geom'])
    ppn = layout['ppn'] if "ppn" not in layout[step].keys() else layout[step]['ppn']
    jpn = int(ppn/ppj)
    basetasks = basenodes * ppn if basenodes > 1 or jpn <= 1 else ppj
    nodes = basenodes*ncases if jpn <= 1 else int((basenodes*ncases + jpn -1)/jpn)
    NP = str(nodes * ppn)
    geom = ".".join([str(i) for i in layout[step]['geom']])

    # Append the number of cases to the step tag, as in A -> A3
    job_name   = param['submit']['job_name_pfx'] + "-" + step + str(ncases)
    os.environ['NP'] = NP
    os.environ['PPN'] = str(ppn)
    os.environ['PPJ'] = str(ppj)
    os.environ['BASETASKS'] = str(basetasks)
    os.environ['BASENODES'] = str(basenodes)
    os.environ['LAYOUT']  = geom

    # Check that the job script exists
    try:
        stat = os.stat(job_script)
    except OSError:
        print("Can't find the job script:", job_script)
        print("Quitting")
        sys.exit(1)

    # Job submission command depends on locale
    scheduler = param['submit']['scheduler']
    if scheduler == 'LSF':
        cmd = [ "bsub", "-nnodes", str(nodes), "-J", job_name, job_script ]
    elif scheduler == 'PBS':
        cmd = [ "qsub", "-l", ",".join(["nodes="+str(nodes)]), "-N", job_name, job_script ]
    elif scheduler == 'SLURM':
        # NEEDS UPDATING
        cmd = [ "sbatch", "-N", str(nodes), "-n", NP, "-J", job_name, "-t", wall_time, job_script ]
    elif scheduler == 'INTERACTIVE':
        cmd = [ "./"+job_script ]
    elif scheduler == 'Cobalt':
        # NEEDS UPDATING IF WE STILL USE Cobalt
        cmd = [ "qsub", "-n", str(nodes), "--jobname", job_name, archflags, "--mode script", "--env LATS="+LATS+":NCASES="+NCASES+":NP="+NP, job_script ]
    else:
        print("Don't recognize scheduler", scheduler)
        print("Quitting")
        sys.exit(1)

    # Run the job submission command
    cmd = " ".join(cmd)
    print(cmd)
    reply = ""
    try:
        reply = subprocess.check_output(cmd, shell=True).decode().splitlines()
    except subprocess.CalledProcessError as e:
        print('\n'.join(reply))
        print("Job submission error.  Return code", e.returncode)
        print("Quitting");
        sys.exit(1)

    print('\n'.join(reply))

    # Get job ID
    if scheduler == 'LSF':
        # a.2100 Q Job <99173> is submitted to default queue <batch>
        jobid = reply[0].split()[1].split("<")[1].split(">")[0]
        if type(jobid) is bytes:
            jobid = jobid.decode('ASCII')
    elif scheduler == 'PBS':
        # 3314170.kaon2.fnal.gov submitted
        jobid = reply[0].split(".")[0]
    elif scheduler == 'SLURM':
        # Submitted batch job 10059729
        jobid = reply[len(reply)-1].split()[3]
    elif scheduler == 'INTERACTIVE':
        jobid = os.environ['SLURM_JOBID']
    elif scheduler == 'Cobalt':
        # ** Project 'semileptonic'; job rerouted to queue 'prod-short'
        # ['1607897']
        jobid = reply[-1]
    if type(jobid) is bytes:
        jobid = jobid.decode('ASCII')

    cfgnos = ""
    for cfgno, index in cfgno_steps:
        cfgnos = cfgnos + cfgno
    date = subprocess.check_output("date",shell=True).rstrip().decode()
    print(date, "Submitted job", jobid, "for", cfgnos, "step", step)

    return (0, jobid)

######################################################################
def mark_queued_todo_entries(step, cfgno_steps, jobid, todo_list):
    """Update the todo_file, change status to "Q" and mark the job number"""

    for k in range(len(cfgno_steps)):
        c, i = cfgno_steps[k]

        todo_list[c][i]   = step+"Q"
        todo_list[c][i+1] = jobid

######################################################################                                                                     
def check_complete():
    """Check completion of queued jobs and purge scratch files"""
    cmd = "../scripts/check_completed.py"
    print(cmd)
    sys.stdout.flush()
    reply = ""
    try:
        reply = subprocess.check_output(cmd, shell=True).decode().splitlines()
    except subprocess.CalledProcessError as e:
        print("Error checking job completion.  Return code", e.returncode)
        for line in reply:
            print(line)

    for line in reply:
        print(line)

    return

######################################################################
def nanny_loop(YAML):
    """Check job periodically and submit to the queue"""
    
    date = subprocess.check_output("date",shell=True).rstrip().decode()
    hostname = subprocess.check_output("hostname",shell=True).rstrip().decode()
    print(date, "Spawn job process", os.getpid(), "started on", hostname)
    sys.stdout.flush()

    param = load_param(YAML)

    # Keep going until
    #   we see a file called "STOP" OR
    #   we have exhausted the list OR
    #   there are job submission or queue checking errors

    check_count = int(param['nanny']['check_interval'])
    while True:
        if os.access("STOP", os.R_OK):
            print("Spawn job process stopped because STOP file is present")
            break

        todo_file   = param['nanny']['todo_file']
        max_cases   = param['nanny']['max_cases']
        job_name_pfx = param['submit']['job_name_pfx']
        scheduler  = param['submit']['scheduler']

        lock_file = lock_file_name(todo_file)

        # Count queued jobs with our job name
        nqueued = count_queue( scheduler, job_name_pfx )
  
        # Submit until we have the desired number of jobs in the queue
        if nqueued < param['nanny']['max_queue']:
            wait_set_todo_lock(lock_file)
            todo_list = read_todo(todo_file)
            remove_todo_lock(lock_file)

            # List a set of cfgnos
            step, cfgno_steps = next_cfgno_steps(max_cases, todo_list)
            ncases = len(cfgno_steps)
        
            # Check completion and purge scratch files for complete jobs
            if check_count == 0:
                check_complete()
                check_count = int(param['nanny']['check_interval'])
            
            if ncases > 0:

                # Make input
                make_inputs(param, step, cfgno_steps)

                # Submit the job

                status, jobid = submit_job(param, step, cfgno_steps, max_cases)
            
                # Job submissions succeeded
                # Edit the todo_file, marking the lattice queued and
                # indicating the jobid
                if status == 0:
                    wait_set_todo_lock(lock_file)
                    todo_list = read_todo(todo_file)
                    mark_queued_todo_entries(step, cfgno_steps, jobid, todo_list)
                    write_todo(todo_file, todo_list)
                    remove_todo_lock(lock_file)
                else:
                    # Job submission failed
                    if status == 1:
                        # Fatal error
                        print("Quitting");
                        sys.exit(1)
                    else:
                        print("Will retry submitting", cfgno_steps, "later")

        sys.stdout.flush()
            
        subprocess.call(["sleep", str( param['nanny']['wait'] ) ])
        check_count -= 1

        # Reload parameters in case of hot changes
        param = load_param(YAML)

############################################################
def main():

    # Set permissions
    os.system("umask 022")
    print(sys.argv)

    YAML = "params.yaml"
        
    nanny_loop(YAML)

############################################################
main()

