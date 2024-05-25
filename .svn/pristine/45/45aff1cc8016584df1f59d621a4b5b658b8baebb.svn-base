import http.client as httplib
import sys
import os
import pickle
import shutil
from aselite import read_any, write_vasp, write_xyz, write_con

host = '146.6.145.114'
port = 8080

def request(tid, directory):
    # Sends a request to the bottle server, checks if thread is finished
    data = {'thread_number':int(tid)}
    data = pickle.dumps(data)
    directory = '/tmp/' + directory
    headers = {'Content-type' : 'bin-data', 'Accept' : 'text/plain'}
    conn = httplib.HTTPConnection(host=host, port=port)
    try:
        conn.request('POST', '/request', data, headers)
    except ConnectionRefusedError:
        exit(1)
    response = conn.getresponse().read()
    if response != b'null':
        suggestion_dict = pickle.loads(response)
        # todo: make directory, place files etc
    else:
        print('inc')
        return
    try:
        os.mkdir(directory)
    except:
        shutil.rmtree(directory)
        os.mkdir(directory)
    for key in suggestion_dict:
        if key == 'scores':
            continue
        write_vasp(directory + "/VASP_SADDLE_%d" % key, suggestion_dict[key][0])
        write_vasp(directory + "/VASP_PRODUCT_%d" % key, suggestion_dict[key][1])
        write_con(directory + "/CON_SADDLE_%d" % key, suggestion_dict[key][0])
        write_con(directory + "/CON_PRODUCT_%d" % key, suggestion_dict[key][1])
        write_xyz(directory + "/SADDLE_%d" % key, suggestion_dict[key][0])
        write_xyz(directory + "/PRODUCT_%d" % key, suggestion_dict[key][1])
        try:
            Kdb().save_mode(directory + "/MODE_%d" % key, suggestion_dict[key][2])
        except:
            pass
    # Return the ordered array for the website
    try:
        print(str(suggestion_dict['scores']), end='')
    except:
        pass


def run(args):
    if len(args) != 2:
        print('Usage: python html_request.py thread_id directory')
        exit(1)
    request(*args)

if __name__ == '__main__':
    run(sys.argv[1:])
