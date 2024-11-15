# Scripts supporting job queue management
# spawnjob.py and check_completed.py

# For Python 3 version

import sys
import os
import subprocess
import time


######################################################################
def lock_file_name(todo_file):
    """Directory entry"""
    return todo_file + ".lock"


######################################################################
def wait_set_todo_lock(lock_file):
    """Set lock file"""

    while os.access(lock_file, os.R_OK):
        print("Lock file present. Sleeping.")
        sys.stdout.flush()
        time.sleep(600)

    subprocess.call(["touch", lock_file])


######################################################################
def remove_todo_lock(lock_file):
    """Remove lock file"""
    subprocess.call(["rm", lock_file])


######################################################################
def read_todo(todo_file):
    """Read the todo file"""

    todo_list = dict()
    try:
        with open(todo_file) as todo:
            todo_lines = todo.readlines()
    except IOError:
        print("Can't open", todo_file)
        sys.exit(1)

    for line in todo_lines:
        if len(line) == 1:
            continue
        a = line.split()
        for i in range(len(a)):
            if isinstance(a[i], bytes):
                a[i] = a[i].decode('ASCII')
        key = a[0]
        todo_list[key] = a

    todo.close()
    return todo_list


######################################################################
def key_todo_entries(td):
    """Sort key for todo entries with format x.nnnn"""

    (stream, cfg) = td.split(".")
    return "{0:s}{1:010d}".format(stream, int(cfg))


######################################################################
def cmp_todo_entries(td1, td2):
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
def write_todo(todo_file, todo_list):
    """Write the todo file"""

    # Back up the files
    subprocess.call(["mv", todo_file, todo_file + ".bak"])

    try:
        todo = open(todo_file, "w")

    except IOError:
        print("Can't open", todo_file, "for writing")
        sys.exit(1)

    for line in sorted(todo_list, key=key_todo_entries):
        print(" ".join(todo_list[line]), file=todo)

    todo.close()


######################################################################
def find_next_unfinished_task(a):
    """Examine todo line "a" to see if more needs to be done"""

    # Format
    # a.1170 SX 0 EX 2147965 LQ 2150955 A 0 H 0

    index = 0
    cfgno = a[0]
    step = ""

    for i in range(1, len(a), 2):
        if ('Q' in a[i]) or ('C' in a[i]) or ('fix' in a[i]):
            # If any entry for this cfg has a Q, we don't try to run a subsequent step
            # because of dependencies.
            # If it is being checked, (marked 'C'),  we also skip it.
            # If it is marked 'fix', we also skip it.
            break
        if not ('X' in a[i]) and not ('C' in a[i]):
            # Found an unfinised task
            index = i
            cfgno = a[0]
            step = a[i]
            break

    return index, cfgno, step


######################################################################
def find_next_queued_task(a):
    """Examine todo line "a" to see if more needs to be done"""

    # Format
    # a.1170 SX 0 EX 2147965 LQ 2150955 A 0 H 0

    index = 0
    cfgno = a[0]
    step = ""

    for i in range(1, len(a), 2):
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
