#!/usr/bin/env python3

import sys

from kdb.common import Kdb
from kdb.local_insert import LocalInsert
from kdb.local_query import  LocalQuery
from kdb.local_db import LocalDB
from kdb.aselite import read_any
from kdb.config import *


def run(args):
    if not Kdb().check_version():
        sys.exit("Python 3.6 or greater is required.")

    if len(args) < 1:
        sys.exit("\nFirst parameter should be either: insert, query\n")
    if args[0] == 'insert':
        if len(args) < 4:
            sys.exit("\nParameters for insert should include reactant, saddle, and product files.\n")

        # Read files and args: reactant, saddle, product, [mode], [b_f], [b_r]
        try:
            reactant = read_any(args[1])
            saddle   = read_any(args[2])
            product  = read_any(args[3])
        except IOError:
            sys.exit("\nOne or more files could not be read. Exiting...\n")
        try:
            mode = Kdb().load_mode(args[4])
        except:
            mode = None

        if mode is None:  # no mode was given
            try:
               b_f=float(args[4])
            except:
               b_f=0.0
            try:
               b_r=float(args[5])
            except:
               b_r=0.0
        else:  # mode.dat was given
            try:
               b_f=float(args[5])
            except:
               b_f=0.0
            try:
               b_r=float(args[6])
            except:
               b_r=0.0
        db = LocalDB(KDB_NAME)
        params = db.get_params()
        #insert
        exitCode = LocalInsert().insert(reactant, saddle, product, mode=mode,forward_bar=b_f,reverse_bar=b_r, dc=params['dc'], nf=params['nf'], mac=params['mac'], kdbname=KDB_NAME)

        if exitCode:
            print("*** KDB INSERT FAILED ***")
        else:
            print("KDB insert success!")
    
    elif args[0] == 'query':
        if len(args) < 2:
            print("\nParameters for query should include a reactant file.\n")
            sys.exit()
        if len(args) > 2 and args[2] == '-c':
            config = parse_cfg()
        else:
            config = None
        #read file
        try:
            reactant = read_any(args[1])
        except IOError:
            print("\nReactant file could not be read. Exiting...\n")
            sys.exit()
        #grab params
        db = LocalDB(KDB_NAME)
        params = db.get_params()
        #query
        LocalQuery().query(reactant, OUTPUT_DIR, dc=params['dc'], nf=params['nf'], kdbname=KDB_NAME, custom_config=config, local_query=True)

if __name__ == "__main__":
    args = sys.argv[1:]
    run(args)
