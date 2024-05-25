from kdb.kdbquery import KdbQuery
from optparse import OptionParser
import sys
from kdb.aselite import read_any
from kdb.remote_db import RemoteDB
from kdb.server_config import *

class RemoteQuery(KdbQuery):
    def __init__(self):
        #print ("Dictionary Reset")
        self.return_dict = {}
        self.scores = list()
    # overloads KdbQuery.query_db()
    def query_db(self, **args):
        db = RemoteDB()
        name = db.get_name(args['reactant'].get_chemical_symbols())
        entries = db.get_process(name)
        return entries, name

    def output_query(self, outputdir, numMatches, suggestion, sugproduct, entry, modeTemp=None, score=None):
        #print (numMatches)
        keyTup = (entry['id'], numMatches)
        stringSugID = str(entry['id'])
        self.return_dict[keyTup] = [suggestion, sugproduct]
        if modeTemp is not None:
            self.return_dict[keyTup].append(modeTemp)
        if score is not None:
            self.scores.append(score)

    def update_scores(self):
        # Packs the scores (best matches) in order for website
        def rank(x):
            return self.scores[x]

        ordered = [n for n in range(len(self.scores))]
        ordered = sorted(ordered, key=rank)
        self.return_dict['scores'] = ordered + self.scores
        #print('actual scores: {}\nsorted scores: {}'.format(self.scores, ordered))

if __name__ == "__main__":

    # Parse command line options.
    parser = OptionParser(usage = "%prog [options] reactant.con")
    parser.add_option("--nodupes", dest = "nodupes", action="store_true",
                      help = "detect and remove duplicate suggestions (can be expensive)")
    options, args = parser.parse_args()

    # Make sure we get the reactant file name.
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    # Load the reactant con file.
    reactant = read_any(args[0])
    db = RemoteDB()
    params = db.get_params()
    query_sub_class = RemoteQuery()
    query_sub_class.query(reactant, "./kdbmatches", dc=params['dc'], nf=params['nf'], nodupes=options.nodupes, kdbname=db_name)
