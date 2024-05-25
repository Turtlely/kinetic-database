import urllib.parse
import http.client as httplib
import pickle
from aselite import read_any, write_vasp, write_xyz, write_con
import sys
import os
import shutil
import getpass
import remote_config
import codecs
try:
    from kdb.common import Kdb
except (ImportError, ModuleNotFoundError) as error:
    from common import Kdb
from zipfile import ZipFile
from os import listdir
from os.path import isfile, join
#import cgitb

#cgitb.enable()
host = '146.6.145.114'
#host = theory.cm.utexas.edu
port = 8080
#print host
def get_all_file_paths(directory): 
  
    # initializing empty file paths list 
    file_paths = []
    # crawling through directory and subdirectories 
    for path, directories, files in os.walk(directory): 
        for filename in files: 
            # join the two strings in order to form the full filepath. 
            filepath = os.path.join(path, filename) 
            file_paths.append(filepath) 
    return file_paths
def create_zip():
     # path to folder which needs to be zipped 
    directory = '/home/www/KDB/kdb_user_info/kdbmatches'
    #print directory
    # calling function to get all file paths in the directory 
    file_paths = get_all_file_paths(directory)
    #print file_paths 
  
    # printing the list of all files to be zipped 
    print('Following files will be zipped:') 
    for file_name in file_paths: 
        print(file_name)
    # writing files to a zipfile
    with ZipFile('/kdb_user_info/kdbmatches.zip','w') as zip:
        # writing each file one by one
        for file in file_paths:
            zip.write(file)
    print('All files zipped successfully!')
def server_query(args, **kwargs):
    reactant = args # [0]
    email = kwargs['user']
    params = {}
    directory = '/tmp/' + email # "/home/www/KDB/kdb_user_info/kdbmatches/" + email
    params = pickle.dumps({'reactant':reactant, 'config':{'USE_GRAPH':False}})
    #params['reactant'] = codecs.encode(pickle.dumps(reactant), 'base64').decode()
    #params  = urllib.parse.urlencode(params)
    headers = {'Content-type': 'bin-data', 'Accept': 'text/plain'}
    conn = httplib.HTTPConnection(host=host, port=port)
    try:
        conn.request('POST', '/query', params, headers)
    except ConnectionRefusedError:
        sys.stderr.write("Bottle is not running.")
        exit(1)
    response = conn.getresponse()
    data = response.read().decode()
    print(data)
    #suggestion_dict = pickle.loads(data)
    #try:
    #    os.mkdir(directory)
	#print "directory made"
    #except:
    #    shutil.rmtree(directory)
    #    os.mkdir(directory)
    #for key in suggestion_dict:
    #    if key == 'scores':
    #        continue
	#print "NO WAYY"
        #write('/kdb_user_info/kdbmatches' + f'vasp_{suggestion_dict[key][0]}', 'vasp')
    #    write_vasp(directory + "/VASP_SADDLE_%d" % key, suggestion_dict[key][0])
    #    write_vasp(directory + "/VASP_PRODUCT_%d" % key, suggestion_dict[key][1])
    #    write_con(directory + "/CON_SADDLE_%d" % key, suggestion_dict[key][0])
    #    write_con(directory + "/CON_PRODUCT_%d" % key, suggestion_dict[key][1])
    #    write_xyz(directory + "/SADDLE_%d" % key, suggestion_dict[key][0])
    #    write_xyz(directory + "/PRODUCT_%d" % key, suggestion_dict[key][1])
    #    try:
    #        Kdb().save_mode(directory + "/MODE_%d" % key, suggestion_dict[key][2])
    #    except:
    #        pass
    # Return the ordered array for the website
    #try:
    #    print(str(suggestion_dict['scores']), end='')
    #except:
    #    pass
    #create_zip()
    #print "done, output now in kdbmatches/"
def run(args):
    if not Kdb().check_version():
        sys.exit()
    if len(args) < 1:
        print("parameters for query should include a reactant file.")
        sys.exit()
    try:
        reactant = read_any(args[0])
        email = args[1]
    except IOError:
        print("could not read reactant file.")
        sys.exit(1)
    except IndexError:
        sys.stderr.write("User queried with invalid email.")
        sys.exit(1)
    if type(reactant) is list: # Movie file submission?
        reactant = reactant[0]
    server_query(reactant, user=email)

if __name__ == "__main__":
    args = sys.argv[1:]
    run(args)

