import urllib
import http.client as httplib
import pickle
import codecs
from aselite import read_any, write_vasp
import sys
import os
import shutil
import getpass
import remote_config
try :
    from common import Kdb
except ImportError as error:
    from kdb.common import Kdb
import cgitb
import cgi
# get the info from the html form
form = cgi.FieldStorage()
host = '146.6.145.114'
port = 8080

def get_account_info():
    with open('/kdb/kdb/.kdb', 'rb') as infile:
        output = pickle.load(infile)
        print(output)
    return output

def set_account_info(info):
    with open('.kdb', 'wb') as outfile:    
	#info =  [(sys.argv[0]),(sys.argv[1])]
        pickle.dump([info[0], info[1]], outfile)
def server_insert(reactant, saddle, product, mode):
    email_pass = []
    email_pass.append(args[0])
    email_pass.append(args[1])
    #print email_pass
    set_account_info(email_pass)
    params = {}
    params['email']    = email_pass[0]
    params['password'] = email_pass[1]
    params['reactant'] = codecs.encode(pickle.dumps(reactant), 'base64').decode()
    params['saddle']   = codecs.encode(pickle.dumps(saddle), 'base64').decode()
    params['product']  = codecs.encode(pickle.dumps(product), 'base64').decode()
    if mode is not None:
        params['mode'] = pickle.dumps(mode)
    params  = urllib.parse.urlencode(params)
    headers = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain'}
    conn    = httplib.HTTPConnection(host=host, port=port)
    conn.request('POST', '/insert', params, headers)
    response = conn.getresponse()
    #print response.status, response.reason
    data = response.read()
    print(data)
    if data == "invalid account info":
        os.remove('.kdb')

def run(args):
    #if not Kdb().check_version():
    #	sys.exit()
    #print len(args)
    if len(args) < 4:
        print("parameters for insert should include reactant, saddle, and product files.")
        sys.exit()
    try:
	#print "hello"
    	reactant = read_any(args[2])
    	saddle   = read_any(args[3])
    	product  = read_any(args[4])
    except IOError:
	 #print "hello"
         print("One or more files could not be read.")
         sys.exit()
    try:
        mode = Kdb().load_mode(args[5])
#        mode = read_any(args[5])
    except:
        mode = None    
    server_insert(reactant, saddle, product, mode)

if __name__ == "__main__":
   args = sys.argv[1:]
   run(args)

