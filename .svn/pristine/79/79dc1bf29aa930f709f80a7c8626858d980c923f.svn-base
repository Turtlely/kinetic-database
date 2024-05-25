#!/usr/bin/env python3

# NOTE: to test MySQL/anything in remote, you must first
#       install MySQL-server (atleast until I setup the DB on Theory server)
#       then run remote_initialize.py

import os
import sys
import subprocess
import math
import numpy
import re

try:
    import sqlite3
except:
    sys.exit("SQLite could not be imported into Python.")

from kdb.aselite import Atoms, read_any
from kdb.common import Kdb


def get_args():
    args = sys.argv
    if len(args) <= 1:
        sys.exit("At least 1 argument is required. Please pass 'local', 'remote', and/or \
                    'clean' as parameters.")
    return args


def logger(line):
    with open("test.log", "a") as f:
        f.write(line)


def print_header(test_mode):
    print( "===========================")
    print(f"     {test_mode} KDB Test")
    print( "===========================")


def remove_files():
    # remove any old queries and/or local databases
    try:
        os.system("rm -rf kdbmatches")
        os.system("rm -f kdb.db")
        os.system("rm -f test.log")
    except:
        pass


def grade_test(output, *keywords):
    for keyword in keywords:
        if keyword not in output:
            return "FAIL"
    return "PASS"


def equivalent(a, b):  # used in comparisons to assess equivalence
    ta = type(a)
    tb = type(b)
    assert ta == tb
    if (ta == float) or (ta == numpy.float64) or (ta == numpy.float32):
        return math.isclose(a, b, abs_tol=1.0e-13)
    elif ta == numpy.ndarray:
        return numpy.all(abs(a - b) < 1e-13)
    elif (ta != Atoms):
        return a == b
    else:  # Atoms objects
        return \
                numpy.all(a.numbers == b.numbers) and \
                numpy.all(abs(a.positions - b.positions) < 1e-13) and \
                numpy.all(abs(a.cell - b.cell) < 1e-13) and \
                numpy.all(a.constraints[0].index == b.constraints[0].index) and \
                numpy.all(a.pbc == b.pbc)


def localdb_cmp(db1, db2): # returns None if the 2 db's are the same
    # open databases to compare
    conn1 = sqlite3.connect(db1)
    conn2 = sqlite3.connect(db2)

    # see if list of tables are the same
    cmd = """SELECT name FROM sqlite_master WHERE type='table';"""  # get list of all tables
    cursor1 = conn1.cursor()
    tbls1 = cursor1.execute(cmd).fetchall()
    cursor2 = conn2.cursor()
    tbls2 = cursor2.execute(cmd).fetchall()
    if set(tbls1) != set(tbls2):
        conn1.close()
        conn2.close()
        return f"List of tables is different\n{tbls1}\nvs\n{tbls2}"

    # see if all entries are equal
    for temp in tbls1:  # tbls1 should be equivalent to tbls2
        tablename = temp[0]
        cmd = f"""SELECT * from {tablename}"""
        entries1 = conn1.execute(cmd).fetchall()
        entries2 = conn2.execute(cmd).fetchall()
        N1 = len(entries1)
        N2 = len(entries2)
        if N1 != N2:
            conn1.close()
            conn2.close()
            return f"Different number of entries in table {tablename}: {N1} vs {N2}"
        for i in range(N1):  # rows
            row1 = entries1[i]
            row2 = entries2[i]
            nr1 = len(row1)
            nr2 = len(row2)
            if nr1 != nr2:
                conn1.close()
                conn2.close()
                return f"Different number of columns in row {i} of table {tablename}: {nr1} \
                            vs {nr2}"
            for j in range(nr1):  # columns
                entry1 = row1[j]
                entry2 = row2[j]
                if not equivalent(entry1, entry2):
                    conn1.close()
                    conn2.close()
                    return f"Different entry in position ({i}, {j}) of {tablename}: {entry1} \
                                vs {entry2}"
    
    # all tests passed
    conn1.close()
    conn2.close()
    return None


def localquery_cmp(target_saddles, new_saddles):  # returns None if 100% successful
    target_leftovers = target_saddles.copy()  # will hold unmatched target saddles
    new_leftovers = new_saddles.copy()  # will hole unmatched new saddles
    for target_saddle_file in target_saddles:
        target_saddle = read_any(target_saddle_file)  # to match
        target_product_file = target_saddle_file.replace("SADDLE", "PRODUCT")
        target_product = read_any(target_product_file)
        target_mode_file = target_saddle_file.replace("SADDLE", "MODE")
        target_mode = Kdb().load_mode(target_mode_file)

        for new_saddle_file in new_leftovers:
            new_saddle = read_any(new_saddle_file)
            if equivalent(target_saddle, new_saddle):
                # compare product, and mode
                new_product_file = new_saddle_file.replace("SADDLE", "PRODUCT")
                new_product = read_any(new_product_file)
                new_mode_file = new_saddle_file.replace("SADDLE", "MODE")
                new_mode = Kdb().load_mode(new_mode_file)
                if equivalent(target_product, new_product) and equivalent(target_mode, new_mode):
                    # remove from leftovers list
                    target_leftovers.remove(target_saddle_file)
                    new_leftovers.remove(new_saddle_file)
                    break

    if target_leftovers or new_leftovers:  # if there are leftovers that are unmatched
        return f"Some saddles (and their associated products and modes) could not be matched: \
                    \n\tReferences: {str(target_leftovers)} \
                    \n\tNew queries: {str(new_leftovers)}"
    return None



def local_test():
    print_header("Local")
    remove_files()

    # test insert with local sqlite3 database
    sys.stdout.write("Testing insert... ")
    sys.stdout.flush()
    insert_out = subprocess.check_output(["kdb_local_client.py","insert",
                    "test_vars/reactant.con","test_vars/saddle.con","test_vars/product.con",
                    "test_vars/mode.dat"])
    insert_out = insert_out.decode('ascii')
    print(grade_test(insert_out, "success"))
    logger("INSERT OUTPUT\n" + insert_out + "\n")


    # test check for duplicates with local sqlite3 database
    sys.stdout.write("Testing insert duplicate check... ")
    sys.stdout.flush()
    dupcheck_out = subprocess.check_output(["kdb_local_client.py","insert",
                    "test_vars/reactant.con","test_vars/saddle.con","test_vars/product.con",
                    "test_vars/mode.dat"])
    dupcheck_out = dupcheck_out.decode('ascii')
    print(grade_test(dupcheck_out, "duplicate", "FAILED"))
    logger("INSERT DUPLICATE OUTPUT\n" + dupcheck_out + "\n")

    # test query with local sqlite3 database
    sys.stdout.write("Testing query... ")
    sys.stdout.flush()
    query_out = subprocess.check_output(["kdb_local_client.py","query","test_vars/reactant.con"])
    query_out = query_out.decode('ascii')
    print(grade_test(query_out, "Finished", "written"))
    logger("QUERY OUTPUT\n" + query_out + "\n")

    
    # comparing new kdb.db with the expected kdb.db in test_vars/
    sys.stdout.write("Testing insert math... ")
    sys.stdout.flush()
    cmp_err = localdb_cmp('kdb.db', 'test_vars/kdb.db')
    if cmp_err is None:
        print("PASS")
    else:
        print("FAIL")
    logger("INSERT MATH ERRORS\n" + str(cmp_err) + "\n")

    # comparing new query results with the expected results in test_vars/kdbmatches/
    sys.stdout.write("Testing query math... ")
    sys.stdout.flush()
    target_saddles = subprocess.check_output("ls -1 test_vars/kdbmatches/SADDLE*", shell=True) \
                    .decode('ascii').strip().split("\n")  # target files in a list
    new_saddles = subprocess.check_output("ls -1 kdbmatches/SADDLE*", shell=True) \
                    .decode('ascii').strip().split("\n")
    cmp_err = localquery_cmp(target_saddles, new_saddles)
    if cmp_err is None:
        print("PASS")
    else:
        print("FAIL")
    logger("QUERY MATH ERRORS\n" + str(cmp_err) + "\n")

    print("See test.log to view KDB outputs and error messages\n")


def remote_test():
    print_header("Remote")
    remove_files()

    print("testing insert")
    out1 = subprocess.check_call(["kdb_remote_client.py","insert","test_vars/reactant.con",
                                    "test_vars/saddle.con","test_vars/product.con",
                                    "--mode test_vars/mode.dat"])
    print(out1.decode('ascii'))

    print("testing query")
    out2 = subprocess.check_call(["kdb_remote_client.py","query","test_vars/reactant.con"])
    print(out2.decode('ascii'))


if __name__ == "__main__":
    args = get_args()
    if "local" in args:
        local_test()
    if "remote" in args:
        remote_test()
    if "clean" in args:
        remove_files()
