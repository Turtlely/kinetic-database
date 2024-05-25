from kdb.bottle import route, post, run, request
from kdb.remote_db import RemoteDB
from kdb.remote_insert import RemoteInsert
from kdb.remote_query import RemoteQuery
#from bottle import route, post, run, request
#from remote_db import RemoteDB
#from remote_insert import RemoteInsert
#from remote_query import RemoteQuery

from random import randrange
import pickle, codecs
import io
import multiprocessing

host = '146.6.145.114'
port = 8080
max_threads = 5
current_threads = []

class ThreadObject():
    def __init__(self, thread, query):
        self.thread = thread
        self.query = query

@post('/account_create')
def create_account():
    email    = request.forms.get('email')
    password = request.forms.get('password')
    first    = request.forms.get('first')
    last     = request.forms.get('last')

    #print('account created')

    db = RemoteDB()
    output = db.add_user(first, last, email, password)
    return output

@post('/login')
def login():
    email = request.forms.get('email')
    password = request.forms.get('password')

    db = RemoteDB()
    #User exists?
    return "route working"

@post('/insert')
def insert():
    db = RemoteDB()
    if not db.is_user(request.forms.get('email'),request.forms.get('password')):
        return "invalid account info"
    reactant = pickle.loads(codecs.decode(request.forms.get('reactant').encode(),"base64")) #NK bytes required but got string
    saddle = pickle.loads(codecs.decode(request.forms.get('saddle').encode(),"base64"))
    product = pickle.loads(codecs.decode(request.forms.get('product').encode(),"base64"))
    try:
        mode = pickle.loads(codecs.decode(request.forms.get('mode').encode(),"base64")) 
    except:
        mode = None
    try:
        forward_bar=request.forms.get('barrier_forward')
    except:
        forward_bar=0.0
    try:
        reverse_bar=request.forms.get('barrier_reverse')
    except:
        reverse_bar=0.0
    insert_class = RemoteInsert()
    insert_class.email = request.forms.get('email')
    insert_class.password = request.forms.get('password')
    #output = insert_class.insert(reactant, saddle, product, mode) #should be kdbinsert's insert function
    output = insert_class.insert(reactant, saddle, product, mode,forward_bar,reverse_bar)
    #output = insert_class.insert_into_db(reactant, saddle, product, mode) #NK commented this line out
    return str(output)

@post('/query')
def query():
    """
if request.get_header('Content-type') == 'bin-data': # We sent some config options
        recv = pickle.loads(request.body.read())
        reactant = recv['reactant']
        try:
            config = recv['config']
        except KeyError:
            config = None
    else:
        reactant = pickle.loads(codecs.decode(request.forms.get('reactant').encode(), 'base64'))
        config = None
    query_class = RemoteQuery()
    parent, child = multiprocessing.Pipe()
    t = multiprocessing.Process(target=query_class.query, args=(reactant,"./kdbmatches", 0.2, 0.3, False, "kdb.db", config, child))
    #query_class.query(reactant, custom_config=config)
    if len(current_threads) >= max_threads:
        # todo: server busy
        return 'busy'
    else:
        #t.id = len(current_threads)
        #to = ThreadObject(t, query_class)
        to = ThreadObject(t, parent)
        current_threads.append(to)
        t.start()
        return str(t.pid)
    """
    #reactant = None
    #reactant = pickle.loads(request.forms.get('reactant'))
    #config = None
    reactant = pickle.loads(codecs.decode(request.forms.get('reactant').encode(), 'base64'))
    #if reactant is None:
    #reactant = pickle.loads(request.forms.get('reactant'))
    query_class = RemoteQuery()
    output = query_class.query(reactant)
    #output = pickle.dumps(query_class.return_dict)
    return pickle.dumps(output)


@post('/request')
def thread_request():
    # Request to see if a thread is finished, if so, return the return dict
    data = pickle.loads(request.body.read())
    #thread = data['thread_number']

    """if len(current_threads) < 1:
        # No threads running
        return 'null'
    else:
        t = None
        for to in current_threads:
            if to.thread.pid == int(thread):
                t = to
                break
        else:
            # Thread id not found
            return 'null'"""
    return_dict = t.query.recv()
    #t.thread.join()
    #current_threads.remove(t)
    return pickle.dumps(return_dict)

    #if t.thread.is_alive():
    #    # is thread still running?
    #    print('null2')
    #    return 'null'
    #else:
    #    current_threads.remove(t)
    #    print('DUMPED')
    #    print(t.query.get())
    #    return pickle.dumps(t.query.get())

def start():
    #run(host=host, port=port, server='paste')
    run(host=host, port=port)

if __name__ == "__main__":
    start()
