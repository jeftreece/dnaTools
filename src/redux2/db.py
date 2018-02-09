# license {{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,os,sqlite3,yaml,time,csv,json,numpy as np
from beautifultable import BeautifulTable
from anytree import Node, RenderTree
#import string
import copy #only used for STASHprint (debugging)
from collections import OrderedDict

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
except:
    trace(0,"Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

#}}}

#TODO: fix this
#start_time = 0 

class DB(object):
    
    def __init__(self):
        #proper db class attributes
        self.db = None
        self.dc = None
        
    def db_init(self):
        #trace (1, "Initialising database...")
        return sqlite3.connect('variant.db')
        
    def cursor(self):
        return self.db.cursor()
        
    def sql_exec_file(self,FILE):
        fh = open(config['REDUX_SQL']+'/'+FILE,'r');
        try:
            #print(fh.read())
            self.dc.executescript(fh.read())
        finally:
            fh.close()
        
    def sql_exec(self,sql):
        if config['DBG_SQL']:
            print("[SQL] "+sql)
        self.dc.execute(sql)
        if config['COMMITS_ON']:
            self.commit()
        
    def sql_exec_many(self,sql,itms):
        if config['DBG_SQL']:
            print("[SQL] "+sql)
        self.dc.executemany(sql,itms)
        if config['COMMITS_ON']:
            self.commit()
        
    def fetchall(self):
        return self.dc.fetchall()
        
    def fetchone(self):
        return self.dc.fetchone()
        
    def commit(self):
        self.db.commit()
