#!/usr/bin/env python3

#from future.standard_library import install_aliases
#install_aliases()

import urllib.request, urllib.parse, urllib.error
import http.client
import pickle
from kdb.aselite import read_any, write_vasp
import sys, codecs
import os
import shutil
import getpass
import time
from kdb.config import parse_cfg
from kdb import remote_config
from os import listdir
from os.path import isfile, join
try:
    from kdb.common import Kdb
except ImportError as error:
    from common import Kdb

host = '146.6.145.114'
port = 8080

def server_create_account():
    #grab info from user
    info = remote_config.RemoteConfig().config()
    #create/populate dictionary
    params = {}
    params['first']    = info[0]
    params['last']     = info[1]
    params['email']    = info[2]
    params['password'] = info[3]
    #format for http post
    params  = urllib.parse.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn    = http.client.HTTPConnection(host=host, port=port)
    #send http POST request
    conn.request('POST', '/account_create', params, headers)
    #grab results
    response = conn.getresponse()
    #print(response.status, response.reason)
    data = response.read()
    #if account created store email/password for later use
    if data.decode('ascii') == "account added":
#        print(data)
#        print("info2",info[2],' ', "info3", info[3])
        set_account_info([info[2],info[3]])
        return 1
    else:
#        print(data)
        answer = input('Would you like to try again? [y/n] ')
        if 'y' in answer.lower()    :
            out = server_create_account()
        else:
            return 0
    return out


def get_account_info():  # get account info (either user input or .kdb file)
    try:
        with open('.kdb', 'rb') as infile:
            email_pass = pickle.load(infile)
    except (EOFError, IOError) as e:
        answer = input("No account information found. Do you have an account? [y/n] ")    
        if 'y' in answer.lower():
            email_pass = []
            email_pass.append(input('email: '))
            email_pass.append(getpass.getpass('password: '))
            set_account_info(email_pass)
            print("Attempting to insert process...")
        elif 'n' in answer.lower():
            print("No problem. Let's create one.")
            if server_create_account():
                with open('.kdb', 'rb') as infile:
                    email_pass = pickle.load(infile)
            else:
                sys.exit("Account could not be created.")
    return email_pass


def set_account_info(info):  # write account information to .kdb file
    with open('.kdb', 'wb') as outfile:
        pickle.dump([info[0], info[1]], outfile)

def server_insert(reactant, saddle, product, mode,b_f,b_r):
    email_pass = get_account_info()

    params = {}
    #print (" params['email'] and params['password']: ", email_pass[0],email_pass[1])
    params['email']    = email_pass[0]
    params['password'] = email_pass[1]
    params['reactant'] = codecs.encode(pickle.dumps(reactant), "base64").decode()
    params['saddle']   = codecs.encode(pickle.dumps(saddle),"base64").decode()
    params['product']  = codecs.encode(pickle.dumps(product),"base64").decode()
    if mode is not None:
        params['mode'] = codecs.encode(pickle.dumps(mode), "base64").decode()
    if b_f != 0.0:
        params['barrier_forward'] = b_f
    if b_r != 0.0:
        params['barrier_reverse'] = b_r
    params  = urllib.parse.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn    = http.client.HTTPConnection(host=host, port=port)
    try:
        conn.request('POST', '/insert', params, headers)
    except ConnectionRefusedError:
        sys.exit("Connection refused by server. Remote KDB might be down.")

    response = conn.getresponse()
    print(response.status, response.reason)
    data = response.read().decode()
    if data == "invalid account info":
        os.remove('.kdb')

def server_query(reactant, **kwargs):
    """params = {}
    params['reactant'] = reactant
    if kwargs['config']:
        params['config'] = parse_cfg()
    params = pickle.dumps(params)
    headers = {'Content-type': 'bin-data', 'Accept': 'text/plain'}
    conn = http.client.HTTPConnection(host=host, port=port)
    conn.request('POST', '/query', params, headers)
    print ("Test")
    response = conn.getresponse()
    thread_id = response.read()
    if thread_id == b'busy':
        print('Too many queries running..')
    else:
        print(f'Task running on PID {thread_id.decode()}')
        params = {}
        params['thread_number'] = thread_id
        params = pickle.dumps(params)
        conn.request('POST', '/request', params)
        response = conn.getresponse().read()
        try:
            response = pickle.loads(response)
        except pickle.UnpicklingError:
            response = response.decode()
        except:
            print('Uncaught error.')
            exit(1)
        if response == 'null':
            print('Invalid request.')
            exit(1)
        elif type(response) == dict:
            write_matches(response)
        else:
            print('unknown error...')
            exit(1)"""
    #print(f'Task running on PID {thread_id.decode()}')
    params = {}        #params['thread_number'] = thread_id
    #params['reactant'] = pickle.dumps(reactant)
    #params = pickle.dumps(params)
    params['reactant'] = codecs.encode(pickle.dumps(reactant), "base64").decode()
    params = urllib.parse.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn = http.client.HTTPConnection(host=host, port=port)
    try:
        conn.request('POST', '/query', params, headers)
    except ConnectionRefusedError:
        sys.exit("Connection refused by server. Remote KDB might be down.")
    #response = conn.getresponse().read()
    response = conn.getresponse().read()
    data = pickle.loads(response)
    if data is None:
        print("No matches found")
    else:
        write_matches(data) 

def write_matches(suggestion_dict):
    # TODO: the output directory should come from the user's parameters, which should replace all the "kdbmatches" below
    try:
        os.mkdir('kdbmatches')
    except:
        shutil.rmtree('kdbmatches')
        os.mkdir('kdbmatches')
    matches_counter = 0
    for key in suggestion_dict:
        if key == 'scores':
            continue
        matches_counter += 1
        write_vasp('kdbmatches' + "/SADDLE_%d_%d" % (key[0],key[1]), suggestion_dict[key][0])
        write_vasp('kdbmatches' + "/PRODUCT_%d_%d" % (key[0],key[1]), suggestion_dict[key][1])
        try:
            Kdb().save_mode('kdbmatches' + "/MODE_%d_%d" % (key[0],key[1]), suggestion_dict[key][2])
        except:
            pass
    '''
    onlyfiles = [f for f in listdir("kdbmatches") if isfile(join("kdbmatches", f))]
    SaddleFiles = []
    ProductFiles = []
    for file in onlyfiles:
        if "SADDLE" in file:
            SaddleFiles.append(file)
        if "PRODUCT" in file:
            ProductFiles.append(file)
    saddleCounter = 0
    productCounter = 0
    for sadFile in SaddleFiles:
        for prodFile in ProductFiles:
            if (sadFile.replace('SADDLE', '') == prodFile.replace('PRODUCT', '')):
                os.rename((os.getcwd()+"/kdbmatches/" + sadFile), (os.getcwd()+"/kdbmatches/"+"SADDLE_%s" % saddleCounter))
                os.rename((os.getcwd()+"/kdbmatches/" + prodFile), (os.getcwd()+"/kdbmatches/"+"PRODUCT_%s" % saddleCounter))
        saddleCounter = saddleCounter + 1
    print("Number of kdb matches", saddleCounter)
    '''
    print("Number of kdb matches", matches_counter)

    print("done, output now in kdbmatches/")

def run(args):
    if not Kdb().check_version():
        sys.exit()
    if len(args) < 1:
        print("first parameter should be either: insert or query")
        sys.exit()
    if args[0] == 'insert':
        print ("all args : ", args)
        if len(args) < 4:
            print("parameters for insert should include reactant, saddle, and product files.")
            sys.exit()
        try:
            reactant = read_any(args[1])
            saddle   = read_any(args[2])
            product  = read_any(args[3])
        except IOError:
            print("One or more files could not be read.")
            sys.exit()
        try: #mode should be 4, BUG?
            mode = Kdb().load_mode(args[4])
        except:
            mode = None

        if mode is None:
            try:
               b_f=float(args[4])
            except:
               b_f=0.0
            try:
               b_r=float(args[5])
            except:
               b_r=0.0
        else: # if no mode given, shift by 1 arg
            try:
               b_f=float(args[5])
            except:
               b_f=0.0
            try:
               b_r=float(args[6])
            except:
               b_r=0.0
        server_insert(reactant, saddle, product, mode, b_f, b_r)
    elif args[0] == 'query':
        if len(args) < 2:
            print("parameters for query should include a reactant file.")
            sys.exit()
        try:
            reactant = read_any(args[1])
        except IOError:
            print("could not read reactant file.")
            sys.exit()
        if len(args) > 2 and args[2] == '-c':
            config = True
        else:
            config = False
        server_query(reactant, config=config)


if __name__ == "__main__":
    args = sys.argv[1:]
    run(args)
