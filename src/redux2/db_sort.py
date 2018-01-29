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
# yaml {{{

try:
    sys.path.append(os.environ['REDUX_PATH'])
    REDUX_CONF = os.environ['REDUX_PATH']+'/config.yaml'
except:
    trace(0,"Missing environment variable REDUX_PATH. Aborting.")
    sys.exit()

config = yaml.load(open(REDUX_CONF))

#}}}

#TODO: needed?
#start_time = 0

def sort_schema(self):
    self.run_sql_file('sort-schema.sql')
        
def sort_insert_sample_data(self):

    #hide-me {{{

    #sample data: 3019783,M343,1,1,Null,1,1,1,1,1,1,1
    #{'v': '3019783', 'n': 'M343', 'k1': '1', 'k2': '1', 'k3': 'Null', 'k4': '1', 'k5': '1', 'k6': '1', 'k7': '1', 'k8': '1', 'k9': '1', 'k10': '1', 'k11': None, 'k12': None}

    #create table s_kits(
    # kit_id  int,  -- later this can be person_id
    # sort_order int
    #);

    #create table s_variants (
    # -- variant_id int, -- not needed for prototype
    # variant_loc int,  -- PK
    # name varchar(20)
    # -- old_reference varchar(2), -- commenting out right now cuz not part of ian's doc
    #);

    #create table s_calls(
    # kit_id int,
    # variant_loc int,
    # assigned boolean
    #);

    # }}}

    cols=10
    #for k in range(1,cols+1):
    kits = "A B C D E F G H I J"
    for k in kits.split():
    self.dc.execute("insert into s_kits (kit_id) values ('"+str(k)+"');")

    with open(REDUX_DATA+'/sample-sort-data.csv','r') as FILE:
        #for row in csv.DictReader(FILE,'v n k1 k2 k3 k4 k5 k6 k7 k8 k9 k10'.split()):
        for row in csv.DictReader(FILE,'vi v n A B C D E F G H I J'.split()):
            row = json.loads(json.dumps(row).replace('\\ufeff','')) #hack: remove byte order mark
            self.dc.execute("insert into s_variants (variant_id,variant_loc,name) values ("+row['vi']+",'"+row['v']+"','"+row['n']+"');")
            #print(' - inserting sample variant data: '+str(row['v']))
            #for k in range(1,cols+1):
            for k in kits.split():
                #kv = str(row['k'+str(k)])
                kv = str(row[str(k)])
                #'null' if kv == "None" else kv
                vv = str(row['v'])
                #print (kv)
                self.dc.execute("insert into s_calls (kit_id,variant_loc,assigned) values ('k"+str(k)+"','"+vv+"',"+kv+");")

    #hide-me {{{
                #self.commit()
                #print (kv+":"+vv)
            #break;
            #sys.exit()
            #print(row)
            #print(row.encode('utf-8-sig'))
        #for (l,n,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12) in dr]
        #    print (n)
        #self.dc.executemany("insert into s_kits (col1, col2) VALUES (?, ?);", to_db)
        #(variant_loc,name,) = t_db
        #con.commit()
        #con.close()
    #}}}

    self.commit()
        

