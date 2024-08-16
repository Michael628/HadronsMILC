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

# Requires TodoUtils.py and params-launch.yaml with definitions of variables needed here

######################################################################
def countQueue( scheduler,  myjobNamePfx ):
    """Count my jobs in the queue"""

    user = os.environ['USER']

    if scheduler == 'LSF':
        cmd = ' '.join(["bjobs -u", user, "| grep", user, "| grep ", myjobNamePfx, "| wc -l"])
    elif scheduler == 'PBS':
        cmd = ' '.join(["qstat -u", user, "| grep", user, "| grep ", myjobNamePfx, "| wc -l"])
    elif scheduler == 'SLURM':
        cmd = ' '.join(["squeue -u", user, "| grep", user, "| grep ", myjobNamePfx, "| wc -l"])
    elif scheduler == 'INTERACTIVE':
        cmd = ' '.join(["squeue -u", user, "| grep", user, "| grep ", myjobNamePfx, "| wc -l"])
    elif scheduler == 'Cobalt':
        cmd = ' '.join(["qstat -fu", user, "| grep", user, "| grep ", myjobNamePfx, "| wc -l"])
    else:
        print("Don't recognize scheduler", scheduler)
        print("Quitting")
        sys.exit(1)

    nqueued = int(subprocess.check_output(cmd,shell=True))

    return nqueued

######################################################################
def nextCfgnoSteps( maxCases, todoList ):
    """Get next sets of cfgnos / job steps from the todo file"""

    # Return a list of cfgnos and indices to be submitted in the next job
    # All subjobs in a single job must do the same step

    step = "none"
    cfgnoSteps = []
    for line in sorted(todoList,key=keyToDoEntries):
        a = todoList[line]
        if len(a) < 2:
            print("ERROR: bad todo line format");
            print(a)
            sys.exit(1)
    
        index, cfgno, newStep = findNextUnfinishedTask(a)
        if index > 0:
            if step == "none":
                step = newStep
            elif step != newStep:
                # Ensure only one step per job
                break
            cfgnoSteps.append([cfgno, index])
            # We don't bundle the S (links) or H (contraction) steps
            if step in ['S','H','I']:
                break
        # Stop when we have enough for a bundle
        if len(cfgnoSteps) >= maxCases:
            break

    ncases = len(cfgnoSteps)
    
    if ncases > 0:
        print("Found", ncases, "cases...", cfgnoSteps)
        sys.stdout.flush()

    return step, cfgnoSteps

######################################################################
def subIOTemplate(param, ioStemTemplate):
    """Replace keys with values in ioStem
       e.g. converts lma-eEIGS-nNOISE to lma-e2000-n1"""

    s = ioStemTemplate

    for k in param["LMIparam"].keys():
        s = re.sub(k,  param['LMIparam'][k],  s)
    
    return s

######################################################################
def setEnv(param, series, cfgno):
    """Set environment variables"""

    # Set environment parameters for the job script
    LMIparam = param['LMIparam']
    for key in LMIparam.keys():
        os.environ[key] = str(LMIparam[key])
    
    # These will be ignored in the bundled job script
    os.environ['SERIES'] = series
    os.environ['CFG']    = str(cfgno)

    # Compute starting time for loose and fine
    dt     = int(LMIparam['DT'])   # Spacing of source times

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
def makeInputs(param, step, cfgnoSteps):
    """Create input XML files for this job"""

    ncases = len(cfgnoSteps)
    INPUTXMLLIST = ""

    if step == 'S':

        # Special treatment for link smearing. No bundling. Self-generated input.
        if ncases > 0:
            print("WARNING: No bundling of smearing jobs")
            print("Will submit only one case")
        ( cfgnoSeries, _ ) = cfgnoSteps[0]
        ( series, cfgno ) = cfgnoSeries.split(".")
        setEnv(param, series, cfgno)

    elif step == 'G':
        INPUTXMLLIST = " ".join([a[0] for a in cfgnoSteps])
        
    else:

        for i in range(ncases):
            ( cfgnoSeries, _ ) = cfgnoSteps[i]

            # Extract series and cfgno  a.1428 -> a 1428
            ( series, cfgno ) = cfgnoSeries.split(".")
            
            # Define common environment variables
            setEnv(param, series, cfgno)
            
            # Name of the input XML file
            ioStemTemplate = param['jobSetup'][step]['io']
            # Replace variables in ioStem if need be
            ioStem = subIOTemplate(param, ioStemTemplate)
            inputXML = ioStem + "-" + cfgnoSeries + ".xml"
        
            # String naming one of the files in ../templates without the .py
            paramModule = f"{param['jobSetup'][step]['params']}_params"

            # import paramModule as pm
            pm = importlib.import_module(paramModule)

            # Template dictionary fir the input XML
            templates = hadronsTemplates.generateTemplates()

            # Generate the input XML file
            fp = open("in/" + inputXML, "w")
            print(dxml(pm.buildParams(**templates)), file=fp)
            fp.close()

            INPUTXMLLIST = INPUTXMLLIST + " " + inputXML

    os.environ['INPUTXMLLIST'] = INPUTXMLLIST

######################################################################
def submitJob(param, step, cfgnoSteps, maxCases):
    """Submit the job"""

    ncases = len(cfgnoSteps)

    jobScript = param['jobSetup'][step]['run']
    wallTime = param['jobSetup'][step]['wallTime']

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
    jobName   = param['submit']['jobNamePfx'] + "-" + step + str(ncases)
    os.environ['NP'] = NP
    os.environ['PPN'] = str(ppn)
    os.environ['PPJ'] = str(ppj)
    os.environ['BASETASKS'] = str(basetasks)
    os.environ['BASENODES'] = str(basenodes)
    os.environ['LAYOUT']  = geom

    # Check that the job script exists
    try:
        stat = os.stat(jobScript)
    except OSError:
        print("Can't find the job script:", jobScript)
        print("Quitting")
        sys.exit(1)

    # Job submission command depends on locale
    scheduler = param['submit']['scheduler']
    if scheduler == 'LSF':
        cmd = [ "bsub", "-nnodes", str(nodes), "-J", jobName, jobScript ]
    elif scheduler == 'PBS':
        cmd = [ "qsub", "-l", ",".join(["nodes="+str(nodes)]), "-N", jobName, jobScript ]
    elif scheduler == 'SLURM':
        # NEEDS UPDATING
        cmd = [ "sbatch", "-N", str(nodes), "-n", NP, "-J", jobName, "-t", wallTime, jobScript ]
    elif scheduler == 'INTERACTIVE':
        cmd = [ "./"+jobScript ]
    elif scheduler == 'Cobalt':
        # NEEDS UPDATING IF WE STILL USE Cobalt
        cmd = [ "qsub", "-n", str(nodes), "--jobname", jobName, archflags, "--mode script", "--env LATS="+LATS+":NCASES="+NCASES+":NP="+NP, jobScript ]
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
    for cfgno, index in cfgnoSteps:
        cfgnos = cfgnos + cfgno
    date = subprocess.check_output("date",shell=True).rstrip().decode()
    print(date, "Submitted job", jobid, "for", cfgnos, "step", step)

    return (0, jobid)

######################################################################
def markQueuedTodoEntries(step, cfgnoSteps, jobid, todoList):
    """Update the todoFile, change status to "Q" and mark the job number"""

    for k in range(len(cfgnoSteps)):
        c, i = cfgnoSteps[k]

        todoList[c][i]   = step+"Q"
        todoList[c][i+1] = jobid

######################################################################                                                                     
def checkComplete():
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
def nannyLoop(YAML):
    """Check job periodically and submit to the queue"""
    
    date = subprocess.check_output("date",shell=True).rstrip().decode()
    hostname = subprocess.check_output("hostname",shell=True).rstrip().decode()
    print(date, "Spawn job process", os.getpid(), "started on", hostname)
    sys.stdout.flush()

    param = loadParam(YAML)

    # Keep going until
    #   we see a file called "STOP" OR
    #   we have exhausted the list OR
    #   there are job submission or queue checking errors

    checkCount = int(param['nanny']['checkInterval'])
    while True:
        if os.access("STOP", os.R_OK):
            print("Spawn job process stopped because STOP file is present")
            break

        todoFile   = param['nanny']['todoFile']
        maxCases   = param['nanny']['maxCases']
        jobNamePfx = param['submit']['jobNamePfx']
        scheduler  = param['submit']['scheduler']

        lockFile = lockFileName(todoFile)

        # Count queued jobs with our job name
        nqueued = countQueue( scheduler, jobNamePfx )
  
        # Submit until we have the desired number of jobs in the queue
        if nqueued < param['nanny']['maxQueue']:
            waitSetTodoLock(lockFile)
            todoList = readTodo(todoFile)
            removeTodoLock(lockFile)

            # List a set of cfgnos
            step, cfgnoSteps = nextCfgnoSteps(maxCases, todoList)
            ncases = len(cfgnoSteps)
        
            # Check completion and purge scratch files for complete jobs
            if checkCount == 0:
                checkComplete()
                checkCount = int(param['nanny']['checkInterval'])
            
            if ncases > 0:

                # Make input
                makeInputs(param, step, cfgnoSteps)

                # Submit the job

                status, jobid = submitJob(param, step, cfgnoSteps, maxCases)
            
                # Job submissions succeeded
                # Edit the todoFile, marking the lattice queued and
                # indicating the jobid
                if status == 0:
                    waitSetTodoLock(lockFile)
                    todoList = readTodo(todoFile)
                    markQueuedTodoEntries(step, cfgnoSteps, jobid, todoList)
                    writeTodo(todoFile, todoList)
                    removeTodoLock(lockFile)
                else:
                    # Job submission failed
                    if status == 1:
                        # Fatal error
                        print("Quitting");
                        sys.exit(1)
                    else:
                        print("Will retry submitting", cfgnoSteps, "later")

        sys.stdout.flush()
            
        subprocess.call(["sleep", str( param['nanny']['wait'] ) ])
        checkCount -= 1

        # Reload parameters in case of hot changes
        param = loadParam(YAML)

############################################################
def main():

    # Set permissions
    os.system("umask 022")
    print(sys.argv)

    YAML = "params.yaml"
        
    nannyLoop(YAML)

############################################################
main()

