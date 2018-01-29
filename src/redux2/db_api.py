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
# conf {{{

try:
    config = yaml.load(open(os.environ['REDUX_CONF']))
except:
    trace(0,"Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

#}}}

#TODO: needed?
#start_time

#redux2

def r2_schema(self):
    self.run_sql_file('redux2-schema.sql')
        
def r2_ins_variants(self,variant_array):
    self.sql_exec_many('''INSERT INTO v1_variants(id,ref,alt) VALUES (?,?,?)''', variant_array)
        
def r2_ins_calls(self):
    self.sql_exec('''INSERT INTO v1_calls(variant,person) SELECT id, person FROM v1_variants CROSS JOIN v1_people''')
        
def r2_ins_hg19(self,snp_reference):
    self.sql_exec_many("INSERT INTO v1_hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)",
        ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
          
def r2_ins_hg38(self,snp_reference):
    self.sql_exec_many("INSERT INTO v1_hg38(grch38,grch38end,name,anc,der) VALUES (?,?,?,?,?)",
        ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))

#clades

def clades_schema(self):
    self.run_sql_file('clades-schema.sql')

#v2

def v2_schema(self):
    self.run_sql_file('schema-v2.sql')

