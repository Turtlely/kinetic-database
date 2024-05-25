# This file holds variables that will never be changed during execution.

# The below parameters are only used when a database has not been created.
MOBILE_ATOM_CUTOFF = 0.7
NEIGHBOR_FUDGE = 0.2
DISTANCE_CUTOFF = 0.3
PBC_MAPPING_CHECK = True
REBOX_SUGGESTIONS = True
KEEP_BARRIERS = False
ENERGY_CEILING = 0.0
ENERGY_FLOOR = 0.0
USE_GROUP_SIMILARITY = False
KDB_NAME = 'kdb.db'
OUTPUT_DIR = "./kdbmatches"
USE_GRAPH = True

# Can easily be made into a non-python file
# Sent over bottle and compared to server-side values
# Overridden on server side
def parse_cfg():
    try:
        file = open('config.ini')
    except:
        print('Config file not found.')
        exit()
    ret = {}
    lines = file.readlines()
    file.close()

    for line in lines:
        try:
            if line.strip()[0] == '#': # comment line
                continue
        except IndexError: # Empty line
            continue
        splitted = line.split('=')
        for j in range(len(splitted)):
            splitted[j] = splitted[j].strip()
        if '#' in line:
            ret[splitted[0]] = splitted[1][:splitted[1].index('#')]
        else:
            ret[splitted[0]] = splitted[1]

    return ret
