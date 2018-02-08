#!/usr/bin/env python3

# license {{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,argparse,yaml,os,glob,shutil,re,time,csv,zipfile
from collections import defaultdict
from controller import *
#from lib import *
import lib,controller
#from clades import *

# }}}

#debugging {{{

def trace (level, msg):
    print(msg)
    #if level <= config['verbosity']:
    #    print(msg)
    #TODO: below line in clades.py
    #sys.stderr(flush)
    
def debug_chk(TYPE,msg):
    if config[TYPE]:
        print(msg)

#}}}
# conf {{{

try:
    config = yaml.load(open(os.environ['REDUX_CONF']))
    #conf = Dict2Obj(config)
except:
    trace(0,"Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

#}}}

# trace (commenting out){{{

#start_time = time.clock()
#trace (1, "Beginning run [%s]" % time.strftime("%H:%M:%S"))

# }}}

# arg parser {{{

# redux2 {{{

parser = argparse.ArgumentParser()
# note: probably needs to be rethought
parser.add_argument('-A', '--all', help='perform all possible steps (prob best not to use for now)', action='store_true')
# note: redux2.bash refactoring area (part 1)
parser.add_argument('-c', '--backup', help='do a "backup" (redux.bash stuff)', action='store_true')
# note: redux2.bash refactoring area (part 2)
parser.add_argument('-p', '--prep', help='prep file structure (redux.bash stuff)', action='store_true')
# note: redux2.py refactoring area
parser.add_argument('-r', '--redux2', help='redux2.py stuff (v1 schema) ', action='store_true')

#}}}
# clades {{{

parser.add_argument('-v', '--verbose', action='count')
parser.add_argument('action', nargs='*')
parser.add_argument('-s', '--snp', nargs=1)
parser.add_argument('-i', '--implications', action='store_true')
parser.add_argument('-t', '--tree', action='store_true')
parser.add_argument('-b', '--badlist', action='store_true')
parser.add_argument('-k', '--kits', action='store_true')

#}}}
# v2 {{{

parser.add_argument('-n', '--new', help='new v2 schema', action='store_true')

# }}}
# sort {{{

parser.add_argument('-o', '--sort', help='sort matrix data prototype (s_ schema currently)', action='store_true')
parser.add_argument('-om', '--sortmatrix', help='sort matrixdata prototype (s_ schema currently)', action='store_true')
parser.add_argument('-ot', '--sorttree', help='sort tree data prototype (s_ schema currently)', action='store_true')

#}}}
# variant {{{

#parser.add_argument('-gv', '--variant', help='variant', type=int, action='store_true')
parser.add_argument('-gv', '--variant', help='variant', type=int)
parser.add_argument('-uv', '--update_variant', help='update variant', action='store_true')

#}}}

#}}}
# arg/namespace {{{

# TODO: what I had
args = parser.parse_args()
# TODO: clades way of doing it
namespace = parser.parse_args(sys.argv[1:])
verbose = vars(namespace)['verbose']

# }}}
# TODO: new clades code {{{
if not verbose:
    verbose = config['DEBUG']
if namespace.snp or len(namespace.action):
    cladesO = Clades();
    cladesO.dbo = DB()
    cladesO.dbo.db = cladesO.dbo.db_init()
    cladesO.namespace = namespace
if namespace.snp:
    cladesO.querysnp = vars(namespace)['snp'][0]

#}}}
# TODO: new clades code {{{

# create: new database from all of the .bed and vcf files
# fastest if you remove the old .db file first
for a in namespace.action:
    if a == 'create':
        cladesO.create = True
    elif a == 'stats1':
        cladesO.stats1 = True
    elif a == 'stats2':
        cladesO.stats2 = True
    elif a == 'docalls':
        cladesO.docalls = True
    elif a == 'listfiles':
        cladesO.listfiles = True
    elif a == 'listbed':
        cladesO.listbed = True
    elif a == 'updatesnps':
        cladesO.updatesnps = True
    elif a == 'mergeup':
        cladesO.mergeup = True
    else:
        print('unknown:', action, 'exiting')

# }}}

# TODO: clades line
#t0 = time.time()

# TODO: clades stuff {{{

if namespace.snp or len(namespace.action):
    if cladesO.create:
        try:
            os.unlink(config['DB_FILE'])
        except:
            pass

# }}}

if args.all:

    # all {{{

    c_r2_backup()
    c_r2_prep()
    c_r2_db()

    # }}}

else: #this area calls controllers

    # redux2 {{{

    # redux.bash stuff
    if args.backup:
        c_r2_backup()

    # redux.bash stuff
    if args.prep:
        c_r2_prep()

    # redux2.py stuff (v1 schema)
    if args.redux2:
        c_r2_db()

    # }}}
    # v2 {{{

    if args.new:
        c_v2_db()

    # }}}
    # sort {{{

    # sort prototype - matrix
    if args.sort or args.sortmatrix:
        c_sort_sample_db()
        c_sort_db_matrix()

    # sort prototype - tree
    if args.sorttree:
        #c_sort_sample_db()
        c_sort_db_tree()

    # }}}
    # variant {{{

    if args.variant:
        #print(args)
        #c_sort_sample_db()
        c_variant(args.variant)

    if args.update_variant:
        c_upd_variant()

    # }}}

# clades (handled in a diff way)  {{{

if namespace.snp or len(namespace.action):
    if cladesO.create:
        cladesO.c_create()
    if cladesO.docalls:
        cladesO.c_docalls()
    if cladesO.stats1:
        cladesO.c_stats1()
    if cladesO.stats2:
        cladesO.c_ctats2()
    if cladesO.listfiles:
        cladesO.c_listfiles()
    if cladesO.listbed:
        cladesO.c_listbed()
    if cladesO.updatesnps:
        cladesO.c_updatesnps()
    if cladesO.querysnp:
        cladesO.c_querysnp()
    if cladesO.mergeup: #incomplete
        cladesO.c_mergeup()

# }}}

# trace (commenting out) {{{

# trace(0, "** script complete.\n")
# TODO: clades line
# trace(1, 'done at {:.2f} seconds'.format(time.time() - t0))

# }}}
