# license {{{

# Purpose: Y-DNA NGS analytics
# Github repo: https://github.com/jazdrv/dnaTools/tree/master
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
# yaml {{{

try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.environ['REDUX_PATH']+'/config.yaml'
except:
    trace(0,"Missing environment variable REDUX_PATH. Aborting.")
    sys.exit()

config = yaml.load(open(REDUX_CONF))

#}}}

#TODO: fix this
start_time = 0 

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
        
    def run_sql_file(self,FILE):
        fh = open(config['REDUX_SQL']+'/'+FILE,'r');
        try:
            #print(fh.read())
            self.dc.executescript(fh.read())
        finally:
            fh.close()
        
    def exec_sql_exec(self,sql):
        self.dc.executemany(sql)
        self.commit()
        
    def exec_sql_exec_many(self,sql):
        self.dc.executemany(sql)
        self.commit()
        
    def commit(self):
        self.db.commit()
