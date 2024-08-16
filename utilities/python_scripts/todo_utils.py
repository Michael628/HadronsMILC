# Scripts supporting job queue management
# spawnjob.py and check_completed.py

# For Python 3 version

import sys, os, yaml, subprocess, time

######################################################################
def lockFileName(todoFile):
    """Directory entry"""
    return todoFile + ".lock"

######################################################################
def waitSetTodoLock(lockFile):
    """Set lock file"""

    while os.access(lockFile, os.R_OK):
        print("Lock file present. Sleeping.")
        sys.stdout.flush()
        time.sleep(600)

    subprocess.call(["touch", lockFile])
    
######################################################################
def removeTodoLock(lockFile):
    """Remove lock file"""
    subprocess.call(["rm", lockFile])


############################################################
def updateParam(param, paramUpdate):
    """Update the param dictionary according to terms in paramUpdate"""

    # Updating is recursive in the tree so we can update selected branches
    # leaving the remainder untouched
    for b in paramUpdate.keys():
        try:
            k = paramUpdate[b].keys()
            n = len(k)
        except AttributeError:
            n = 0

        if b in param.keys() and n > 0:
            # Keep descending until we run out of branches
            updateParam(param[b], paramUpdate[b])
        else:
            # Then stop, replacing just the last branch or creating a new one
            param[b] = paramUpdate[b]

    return param

######################################################################
def loadParam(file):
    """Read the YAML parameter file"""

    try:
        param = yaml.safe_load(open(file,'r'))
    except subprocess.CalledProcessError as e:
        print("WARNING: loadParam failed for", e.cmd)
        print("return code", e.returncode)
        sys.exit(1)

    return param

############################################################
def loadParamsJoin(YAMLEns, YAMLAll):
    """Concatenate two YAML parameter files and load
    We need this because YAMLEns defines a reference needed
    by YAMLAll"""
    
    # Initial parameter file
    try:
        ens = open(YAMLEns,'r').readlines()
        all = open(YAMLAll,'r').readlines()
        param = yaml.safe_load("".join(ens+all))
    except:
        print("ERROR: Error loading the parameter files", YAMLEns, YAMLAll)
        sys.exit(1)

    return param

######################################################################
def readTodo(todoFile):
    """Read the todo file"""
    
    todoList = dict()
    try:
        with open(todoFile) as todo:
            todoLines = todo.readlines()
    except IOError:
        print("Can't open", todoFile)
        sys.exit(1)

    for line in todoLines:
        if len(line) == 1:
            continue
        a = line.split()
        for i in range(len(a)):
            if type(a[i]) is bytes:
                a[i] = a[i].decode('ASCII')
        key = a[0]
        todoList[key] = a

    todo.close()
    return todoList

######################################################################
def keyToDoEntries(td):
    """Sort key for todo entries with format x.nnnn"""

    (stream, cfg) = td.split(".")
    return "{0:s}{1:010d}".format(stream, int(cfg))

######################################################################
def cmpToDoEntries(td1, td2):
    """Compare todo entries with format x.nnnn"""
    # Python 2.7 only

    (stream1, cfg1) = td1.split(".")
    (stream2, cfg2) = td2.split(".")

    # Sort first on stream, then on cfg
    order = cmp(stream1, stream2)
    if order == 0:
        order = cmp(int(cfg1), int(cfg2))

    return order

######################################################################
def writeTodo(todoFile, todoList):
    """Write the todo file"""

    # Back up the files
    subprocess.call(["mv", todoFile, todoFile + ".bak"])

    try:
        todo = open(todoFile, "w")

    except IOError:
        print("Can't open", todoFile, "for writing")
        sys.exit(1)
            
    for line in sorted(todoList, key=keyToDoEntries):
        print(" ".join(todoList[line]), file=todo)

    todo.close()

        
######################################################################
def findNextUnfinishedTask(a):
    """Examine todo line "a" to see if more needs to be done"""

    # Format
    # a.1170 SX 0 EX 2147965 LQ 2150955 A 0 H 0

    index = 0
    cfgno = a[0]
    step = ""
    
    for i in range(1,len(a),2):
        if ( 'Q' in a[i] ) or ( 'C' in a[i] ) or ( 'fix' in a[i] ):
            # If any entry for this cfg has a Q, we don't try to run a subsequent step
            # because of dependencies.
            # If it is being checked, (marked 'C'),  we also skip it.
            # If it is marked 'fix', we also skip it.
            break
        if not ( 'X' in a[i] ) and not ( 'C' in a[i] ) :
            # Found an unfinised task
            index = i
            cfgno = a[0]
            step = a[i]
            break

    return index, cfgno, step

######################################################################
def findNextQueuedTask(a):
    """Examine todo line "a" to see if more needs to be done"""

    # Format
    # a.1170 SX 0 EX 2147965 LQ 2150955 A 0 H 0

    index = 0
    cfgno = a[0]
    step = ""
    
    for i in range(1,len(a),2):
        if 'Q' in a[i]:
            index = i
            step = a[i]
            break
        elif 'C' in a[i]:
            # If we find a 'C', this entry is being checked by another process
            # Empty "step" flags this
            index = i
            break

    return index, cfgno, step

            
        
