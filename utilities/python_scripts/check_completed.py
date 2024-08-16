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
def jobStillQueued(param, jobID):
    """Get the status of the queued job"""
    # This code is locale dependent

    scheduler = param['submit']['scheduler']
    
    user = os.environ['USER']
    if scheduler == 'LSF':
        cmd = " ".join(["bjobs", "-u", user, "|", "grep -w", jobID])
    elif scheduler == 'PBS':
        cmd = " ".join(["qstat", "-u", user, "|", "grep -w", jobID])
    elif scheduler == 'SLURM':
        cmd = " ".join(["squeue", "-u", user, "|", "grep -w", jobID])
    elif scheduler == 'Cobalt':
        cmd = " ".join(["qstat", "-fu", user, "|", "grep -w", jobID])
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
def markCompletedTodoEntry(seriesCfg, precTsrc, todoList):
    """Update the todoList, change status to X"""

    key = seriesCfg + "-" + precTsrc
    todoList[key] = [ seriesCfg, precTsrc, "X" ]
    print("Marked cfg", seriesCfg, precTsrc, "completed")


######################################################################
def markCheckingTodoEntry(seriesCfg, precTsrc, todoList):
    """Update the todoList, change status to X"""

    key = seriesCfg + "-" + precTsrc
    todoList[key] = [ seriesCfg, precTsrc, "C" ]

#######################################################################
def decodeSeriesCfg(seriesCfg):
    """Decode series, cfg, as it appeaers in the todo file"""
    return seriesCfg.split(".")

#######################################################################
def decodePrecTsrc(seriesCfg):
    """Decode prec, tsrc, as it appeaers in the todo file
       Takes P.nn -> [P, nnn]"""
    return seriesCfg.split(".")

######################################################################
def purgeProps(param,seriesCfg):
    """Purge propagators for the specified configuration"""

    print("Purging props for", seriesCfg)
    series, cfg = decodeSeriesCfg(seriesCfg)
    configID = codeCfg(series, cfg)
    prop = param['files']['prop']
    subdirs = prop['subdirs'] + [ configID ]
    remotePath = os.path.join(*subdirs)
    cmd = ' '.join([ "nohup", "/bin/rm -r", remotePath, "> /dev/null 2> /dev/null &"])
    print(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: can't remove props.  Error code", e.returncode, ".")

######################################################################
def purgeRands(param,seriesCfg):
    """Purge random sources for the specified configuration"""

    print("Purging rands for", seriesCfg)
    series, cfg = decodeSeriesCfg(seriesCfg)
    configID = codeCfg(series, cfg)
    rand = param['files']['rand']
    subdirs = rand['subdirs'] + [ configID ]
    remotePath = os.path.join(*subdirs)
    cmd = ' '.join([ "nohup", "/bin/rm -r", remotePath, "> /dev/null 2> /dev/null &"])
    print(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: can't remove rands.  Error code", e.returncode, ".")

######################################################################
def tarInputPath(stream, s06Cfg, precTsrc):
    """Where the data and logs are found"""
    return os.path.join(stream, s06Cfg, precTsrc)

######################################################################
def purgeSymLinks(param, jobCase):
    """Purge symlinks for the specified jobID"""

    (stream, series, cfg, prec, tsrc, s06Cfg, tsrcID, jobID, jobSeqNo)  = jobCase

    print("Purging symlinks for job", jobID)

    io = param['files']['out']
    logsPath = os.path.join(tarInputPath(stream, s06Cfg, tsrcID), io['subdir'])
    cmd = ' '.join([ "find -P", logsPath, "-lname '?*Job'"+ jobID + "'*' -exec /bin/rm '{}' \;"])
    print(cmd)
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: rmdir exited with code", e.returncode, ".")

######################################################################
def goodLogs(param, jobCase):
    """Check that the log files are complete"""

    (stream, series, cfg, prec, tsrc, s06Cfg, tsrcID, jobID, jobSeqNo)  = jobCase
    precTsrcConfigId = [ prec, tsrc, series, cfg ]

    for step in range(param['job']['steprange']['high']):
        expectFile = outFileName(stream, precTsrcConfigId, jobSeqNo, '', "step" + str(step))
        logPath = os.path.join(stream, s06Cfg, tsrcID, "logs", expectFile)
        try:
            stat = os.stat(logPath)
        except OSError:
            print("ERROR: Can't find expected output file", path)
            return False

        # Check for "RUNNING COMPLETED"
        entries = countPhrase(logPath, 'RUNNING COMPLETED')
        if entries < 1:
            print("ERROR: did not find 'RUNNING COMPLETED' in", logPath)
            return False

        # Check for nonconvergence, signaled by lines with "NOT"
        entries = countPhrase(logPath, "NOT")
        if entries > 0:
            print("WARNING: ", entries, "lines with 'NOT' suggesting nonconvergence")
#            return False
    
    # Passed these tests
    print("Output files OK")
    print("COMPLETE")

    return True

######################################################################
def checkPath(param, jobKey, fileKey, cfgno, complain):
    """Complete the file path and check that it exists and has the correct size"""

    # Substute variables coded in file path
    filepath = os.path.join( param['files']['home'], param['files'][jobKey][fileKey] )
    for v in param['LMIparam'].keys():
        filepath = re.sub(v, param['LMIparam'][v], filepath)
    series, cfg = cfgno.split('.')
    filepath = re.sub('SERIES', series, filepath)
    filepath = re.sub('CFG', cfg, filepath)

    good = True
    try:
        fileSize =  os.path.getsize(filepath)
    except OSError:
        good = False

    if good and fileSize >= param['files'][jobKey]['goodSize']:
        return True
    
    if complain:
        print("File", filepath, "not found or not of correct size")

    return False
    
######################################################################
def goodLinks(param, cfgno):
    """Check that the ILDG links look OK"""

    good = checkPath(param, 'fnlinks', 'fat', cfgno, True)
    good = good and checkPath(param, 'fnlinks', 'lng', cfgno, True)

    return good

######################################################################
def goodEigs(param, cfgno):
    """Check that the eigenvector file looks OK"""

    good = checkPath(param, 'eigs', 'eig', cfgno, False)

    if not good:
        # Check file in subdir
        good = checkPath(param, 'eigsdir', 'eigdir', cfgno, True)

    return good


######################################################################
def goodLMA(param, cfgno):
    """Check that the LMA output looks OK"""

    lma = param['files']['lma']
    good =          checkPath(param, 'lma', 'ama', cfgno, True)
    good = good and checkPath(param, 'lma', 'ranLL', cfgno, True)

    if not good and 'ama_alt' in lma.keys():
        good =          checkPath(param, 'lma', 'ama_alt', cfgno, True)
        good = good and checkPath(param, 'lma', 'ranLL_alt', cfgno, True)
    
    return good

######################################################################
def goodA2ALocal(param, cfgno):
    """Check that the A2A output looks OK"""

    good =          checkPath(param, 'a2a_local', 'gamma5', cfgno, True)
    good = good and checkPath(param, 'a2a_local', 'gammaX', cfgno, True)
    good = good and checkPath(param, 'a2a_local', 'gammaY', cfgno, True)
    good = good and checkPath(param, 'a2a_local', 'gammaZ', cfgno, True)
                     
    return good

######################################################################
def goodA2AOnelink(param, cfgno):
    """Check that the A2A output looks OK"""

    good =          checkPath(param, 'a2a_onelink', 'gammaX', cfgno, True)
    good = good and checkPath(param, 'a2a_onelink', 'gammaY', cfgno, True)
    good = good and checkPath(param, 'a2a_onelink', 'gammaZ', cfgno, True)
                     
    return good

######################################################################
def goodContractLocal(param, cfgno):
    """Check that the contrraction output looks OK"""

    good =          checkPath(param, 'contract_local', 'pion', cfgno, True)
    good = good and checkPath(param, 'contract_local', 'vecX', cfgno, True)
    good = good and checkPath(param, 'contract_local', 'vecY', cfgno, True)
    good = good and checkPath(param, 'contract_local', 'vecZ', cfgno, True)
                     
    return good

######################################################################
def goodContractOnelink(param, cfgno):
    """Check that the contrraction output looks OK"""

    good =          checkPath(param, 'contract_onelink', 'vecX', cfgno, True)
    good = good and checkPath(param, 'contract_onelink', 'vecY', cfgno, True)
    good = good and checkPath(param, 'contract_onelink', 'vecZ', cfgno, True)
                     
    return good

######################################################################
def goodContractOnelinkPy(param, cfgno):
    """Check that the contrraction output looks OK"""

    good =          checkPath(param, 'contract_onelink_py', 'vec', cfgno, True)

    return good

######################################################################
def moveFailedOutputs(jobCase):
    """Move failed output to temporary failure archive"""
    
    (stream, series, cfg, prec, tsrc, s06Cfg, tsrcID, jobID, jobSeqNo)  = jobCase

    badOutputPath = tarInputPath(stream, s06Cfg, tsrcID)
    failPath = os.path.join(stream, s06Cfg, "fail", jobID)

    # Move the failed output
    cmd = " ".join(["mkdir -p ", failPath, "; mv", badOutputPath, failPath])
    print(cmd)
    try:
        subprocess.check_output(cmd, shell = True).decode("ASCII")
    except subprocess.CalledProcessError as e:
        status = e.returncode

######################################################################
def nextFinished(param, todoList, entryList):
    """Find the next well-formed entry marked "Q" whose job is no longer in the queue"""
    a = ()
    nskip = 0
    while len(entryList) > 0:
        cfgno = entryList.pop(0)
        a = todoList[cfgno]
        index, cfgno, step = findNextQueuedTask(a)
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
        print("Checking cfg", todoList[cfgno]                     )
        print("------------------------------------------------------------------------------------")
        
        # Is job still queued?
        jobID = a[index+1]
        if jobStillQueued(param, jobID):
            index = 0  # To signal no checking
            continue
        break

    return index, cfgno, step

######################################################################
def checkPendingJobs(YAML):
    """Process all entries marked Q in the todolist"""

    # Read primary parameter file
    param = loadParam(YAML)

    # Read the todo file
    todoFile = param['nanny']['todoFile']
    lockFile = lockFileName(todoFile)

    # First, just get a list of entries
    waitSetTodoLock(lockFile)
    todoList = readTodo(todoFile)
    removeTodoLock(lockFile)
    entryList = sorted(todoList,key=keyToDoEntries)

    # Run through the entries. The entryList is static, but the
    # todo file could be changing due to other proceses
    while len(entryList) > 0:
        # Reread the todo file (it might have changed)
        waitSetTodoLock(lockFile)
        todoList = readTodo(todoFile)

        index, cfgno, step = nextFinished(param, todoList, entryList)
        if index == 0:
            removeTodoLock(lockFile)
            continue

        step = step[:-1]
        # Mark that we are checking this item and rewrite the todo list
        todoList[cfgno][index] = step + "C"
        writeTodo(todoFile, todoList)
        removeTodoLock(lockFile)

        if step not in param["jobSetup"].keys():
            print("ERROR: unrecognized step key", step)
            sys.exit(1)

        # Check that the job completed successfully
        sfx = ""
        status = True
        if step == "S":
            status = status and goodLinks(param, cfgno)
        if step == "E":
            status = status and goodEigs(param, cfgno)
        if step in ["L","A","B","N","D"]:
            status = status and goodLMA(param, cfgno)
        if step in ["L","A","M","M1","D"]:
            status = status and goodA2AOnelink(param, cfgno)
        if step in ["L","A","M"]:
            status = status and goodA2ALocal(param, cfgno)
        if step == "H":
            status = status and goodContractLocal(param, cfgno)
        if step == "I":
            status = status and goodContractOnelink(param, cfgno)
        if step == "G":
            status = status and goodContractOnelinkPy(param, cfgno)

        sys.stdout.flush()

        # Update the entry in the todo file
        waitSetTodoLock(lockFile)
        todoList = readTodo(todoFile)
        if status:
            todoList[cfgno][index] = step+"X"
            print("Job step", step, "is COMPLETE")
        else:
            todoList[cfgno][index] = step+"XXfix"
            print("Marking todo entry XXfix.  Fix before rerunning.")
        writeTodo(todoFile, todoList)
        removeTodoLock(lockFile)

        # Take a cat nap (avoids hammering the login node)
        subprocess.check_call(["sleep", "1"])

############################################################
def main():

    # Parameter file

    YAML = "params.yaml"

    checkPendingJobs(YAML)


############################################################
main()
