#!/usr/bin/env python3

# authors/licensing {{{

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,argparse,yaml,os,glob,shutil,re,time,csv,zipfile
from collections import defaultdict
from lib import *
#from db import *
from clades import *

# }}}
# setup notes {{{

# beyond the modules specified in the lib areas, two env variables are
# required (ie: can be set in ~/.bashrc):

# REDUX_PATH - where these bin and conf files are
# REDUX_ENV - where the data will be

# }}}

# trace

start_time = time.clock()
trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

# env

try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.environ['REDUX_PATH']+'/config.yaml'
except:
    trace(0,"Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
try:
    REDUX_ENV = os.environ['REDUX_ENV']
except:
    trace(0,"Missing environment variable REDUX_ENV. Aborting.")
    sys.exit()

config = yaml.load(open(REDUX_CONF))
#print ('\n'+sys.argv[0]+' version: '+config['VERSION']+'\n')

# arg parser (we can replace this -if useful- with getOpts or whatever)

parser = argparse.ArgumentParser()

#note: probably needs to be rethought
parser.add_argument('-A', '--all', help='perform all possible steps (prob best not to use for now)', action='store_true')

#note: redux2.bash refactoring area (part 1)
parser.add_argument('-c', '--backup', help='do a "backup" (redux.bash stuff)', action='store_true')

#note: redux2.bash refactoring area (part 2)
parser.add_argument('-p', '--prep', help='prep file structure (redux.bash stuff)', action='store_true')

#note: redux2.py refactoring area
parser.add_argument('-r', '--redux2', help='redux2.py stuff (v1 schema) ', action='store_true')

#note: sort prototype stuff
parser.add_argument('-o', '--sort', help='sort data prototype (s_ schema currently)', action='store_true')

#note: Jef's new v2 schema
parser.add_argument('-n', '--new', help='new v2 schema', action='store_true')

#note: clades.py args
parser.add_argument('-v', '--verbose', action='count')
parser.add_argument('action', nargs='*')
parser.add_argument('-s', '--snp', nargs=1)
parser.add_argument('-i', '--implications', action='store_true')
parser.add_argument('-t', '--tree', action='store_true')
parser.add_argument('-b', '--badlist', action='store_true')
parser.add_argument('-k', '--kits', action='store_true')

#TODO: what I had
args = parser.parse_args()
#TODO: clades way of doing it
namespace = parser.parse_args(sys.argv[1:])
verbose = vars(namespace)['verbose']

#TODO: new clades code
if not verbose:
    verbose = config['DEBUG']
if namespace.snp:
    querysnp = vars(namespace)['snp'][0]
else:
    querysnp = False

#TODO: new clades code
# create: new database from all of the .bed and vcf files
# fastest if you remove the old .db file first
create = stats1 = stats2 = docalls = listfiles = listbed = mergeup = updatesnps = False
for a in namespace.action:
    if a == 'create':
        create = True
    elif a == 'stats1':
        stats1 = True
    elif a == 'stats2':
        stats2 = True
    elif a == 'docalls':
        docalls = True
    elif a == 'listfiles':
        listfiles = True
    elif a == 'listbed':
        listbed = True
    elif a == 'updatesnps':
        updatesnps = True
    elif a == 'mergeup':
        mergeup = True
    else:
        print('unknown:', action, 'exiting')

#TODO: clades line
t0 = time.time()

#TODO: clades stuff
if create:
    try:
        os.unlink(config['DB_FILE'])
    except:
        pass

#TODO: clades line
#dbconn = sqlite3.connect(config['DB_FILE'])
#dbcurs = dbconn.cursor()

#args

if args.all:
    go_backup()
    go_prep()
    go_db()
else:
    # redux.bash stuff
    if args.backup:
        go_backup()
    # redux.bash stuff
    if args.prep:
        go_prep()
    # redux2.py stuff (v1 schema)
    if args.redux2:
        go_v1_db()
    # sort prototype
    if args.sort:
        go_sort_db()
    # (v2 schema)
    if args.new:
        go_db()

if create:
    c_create()
if docalls:
    c_docalls()
if stats1:
    c_stats1()
if stats2:
    c_ctats2()
if listfiles:
    c_listfiles()
if listbed:
    c_listbed()
if updatesnps:
    c_updatesnps()
if querysnp:
    c_querysnp(namespace)
if mergeup: #incomplete
    c_mergeup()

trace(0, "** script complete.\n")
#TODO: clades line
#dbconn.commit()
#dbcurs.close()
#TODO: clades line
trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))
