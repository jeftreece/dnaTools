# license {{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,os,sqlite3,yaml,time,csv,json,numpy as np
from beautifultable import BeautifulTable
import itertools
from anytree import Node, RenderTree
import copy #only used for STASHprint (debugging)
from collections import OrderedDict
from lib import *
#import string

# }}}

#debugging {{{

def trace (level, msg):
    print(msg)
    #if level <= config['verbosity']:
    #    print(msg)
    #TODO: below line in clades.py
    #sys.stderr(flush)
    
def debug_chk(var,msg,lev=0):
    if config[var] > lev:
        print(msg)

def sql_on():
    config['DEBUG_SQL'] = True
    
def sql_off():
    config['DEBUG_SQL'] = False

#}}}
# conf {{{

global config
try:
    config = yaml.load(open(os.environ['REDUX_CONF']))
    #print(conf)
except:
    trace(0,"Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

#}}}

class Sort(object):
    
    def __init__(self):

        #db attributes
        self.dbo = None #db object
        #self.db = None #sqlite db instance
        #self.dc = None #sqlite db cursor

        #tree attributes    
        self.TREE = {}
        self.REF = None
        self.MATRIX_MODE = config['MATRIX_MODE']
        self.TREE_MODE = config['TREE_MODE'] #(sort tree presentation) 1=letters, 2=names, 3=letters+names 
        self.TDATA = None

        #matrix attributes
        self.KITS = None
        self.VARIANTS = None
        self.DATA = None
        self.CNTS = {}
        self.NP = None
        self.NONES = []
        self.MDATA = None

    # schema / sample data

    def sort_schema(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.dbo.sql_exec_file('sort-schema.sql')
        
    def sort_ins_sample_data(self):

        cols=10
        kits = "A B C D E F G H I J"
        for k in kits.split():
            #s_kits
            self.dbo.sql_exec("insert into s_kits (kit_id) values ('"+str(k)+"');")

        with open(config['REDUX_DATA']+'/sample-sort-data.csv','r') as FILE:
            for row in csv.DictReader(FILE,'vi v n A B C D E F G H I J'.split()):
                row = json.loads(json.dumps(row).replace('\\ufeff','')) #hack: remove byte order mark
                #s_variants
                self.dbo.sql_exec("insert into s_variants (variant_id,variant_loc,name) values ("+row['vi']+",'"+row['v']+"','"+row['n']+"');")
                for k in kits.split():
                    kv = str(row[str(k)])
                    vv = str(row['v'])
                    #s_calls
                    self.dbo.sql_exec("insert into s_calls (kit_id,variant_loc,assigned) values ('k"+str(k)+"','"+vv+"',"+kv+");")

    # matrix

    def sort_matrix(self):

        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #get data
        self.get_matrix_data()

        #step 0
        debug_chk('DEBUG_MATRIX',"data - step 0 (default)",1)
        self.stdout_tbl_matrix()

        #step 1
        debug_chk('DEBUG_MATRIX',"data - step 1",1)
        self.sort_step1()
        self.stdout_tbl_matrix()

        #step 2
        debug_chk('DEBUG_MATRIX',"data - step 2",1)
        self.sort_step2()
        self.stdout_tbl_matrix()

        #step 3
        debug_chk('DEBUG_MATRIX',"data - step 3",1)
        self.sort_step3()
        self.stdout_tbl_matrix()

        sys.exit()

    def sort_step1(self):

        #print("...")
        #print('kB')
        #print(self.get_matrix_col_data(name='kB'))
        #print("...")
        #print('Z381')
        #print(self.get_matrix_row_data(name='Z381'))
        #print("...")
        #print(self.get_matrix_row_indices_by_val(1,name='Z381'))
        #print("...")
        #print(self.get_lowest_superset_variant(name='Z381'))
        #sys.exit()

        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        for K,V in self.get_axis('variants'):
            #if 0 not in self.get_numpy_matrix_row_as_list(V[1]):
            if -1 not in self.get_numpy_matrix_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_numpy_matrix_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('vp'):
            #if 0 in self.get_numpy_matrix_row_as_list(V[1]) and 'None' not in self.get_numpy_matrix_row_as_list(V[1]):
            if -1 in self.get_numpy_matrix_row_as_list(V[1]) and 0 not in self.get_numpy_matrix_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_numpy_matrix_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('variants'):
            #if 0 in self.get_numpy_matrix_row_as_list(V[1]) and 'None' in self.get_numpy_matrix_row_as_list(V[1]):
            if -1 in self.get_numpy_matrix_row_as_list(V[1]) and 0 in self.get_numpy_matrix_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_numpy_matrix_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for NO in new_orders:
            self.set_new_order(NO[0],NO[1],variantType=True)
        self.NP = np.matrix(list(DATA.values()))
        
    def sort_step2(self):
        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        self.NP = np.transpose(self.NP)
        for K,V in self.get_axis('kp'):
            new_orders.append([K,cnt])
            DATA[K] = self.get_numpy_matrix_row_as_list(V[1],noneToStr=False)
            cnt = cnt + 1
        for NO in new_orders:
            self.set_new_order(NO[0],NO[1],kitType=True)
        self.NP = np.matrix(list(DATA.values()))
        self.NP = np.transpose(self.NP)
        
    def sort_step3(self):
        print("")
        #print(self.get_matrix_row(1))
        #sys.exit()
        print("Processing Nones:")
        print("----------------")
        #variant list that have kits with negative (zero) values
        zlist = np.unique(np.argwhere(self.NP == 0)[:,0]).tolist()
        #iterate all None situations
        for non in ((np.argwhere(self.NP == -1)).tolist()):
            if non[0] in zlist:
                print(self.stdout_coord(non[1],non[0]) + "("+str(self.get_lowest_superset_variant(order=non[0]))+")")
                #for itm in list(self.VARIANTS.items()):
                #    if itm[1][1] == non[0]:
                #        variant = itm[0]
                #        break
                #for itm in list(self.KITS.items()):
                #    if itm[1][1] == non[1]:
                #        kit = itm[0]
                #        break
                #print("kit:"+str(kit)+",variant:"+str(variant))

        print("")
        print(self.MIXA)
        print("")
        sys.exit()

    def stdout_tbl_matrix(self):
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',"big_matrix view{{"+"{",1)
        table = BeautifulTable()
        table.column_headers = ['top']+self.get_cur_kit_list()
        for K,V in self.get_axis('variants'):
            table.append_row([K]+self.get_numpy_matrix_row_as_list(V[1]))
        #table = table.replace(" 0 ","\033[1;31m 0 \033[1;37m")
        debug_chk('DEBUG_MATRIX',table,1)
        debug_chk('DEBUG_MATRIX',"}}"+"}",1)
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',self.get_axis('kits',keysOnly=True),1)
        debug_chk('DEBUG_MATRIX',self.get_axis('variants',keysOnly=True),1)
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',self.NP,1)
        debug_chk('DEBUG_MATRIX',"",1)
        
    def stdout_coord(self,X,Y,moreInfo=False):
        buf = ""
        if moreInfo:
            buf = "coord: "+str(X)+","+str(Y)
        kit = self.get_kit_name_by_order(X)
        variant = self.get_variant_name_by_order(Y)
        buf = buf +  "k|v: "+str(kit)+"|"+(variant)
        if moreInfo:
            buf = buf + "value:"+str(self.NP[Y,X])
        return buf

    def get_cur_kit_list(self):
        return self.get_axis('kits',keysOnly=True)
        
    def get_numpy_matrix_row_as_list(self,rownum,noneToStr=True):
        if noneToStr:
            return ['None' if v is None else v for v in self.NP[rownum,:].tolist()[0]]
        else:
            return self.NP[rownum,:].tolist()[0]
        
    def get_cur_variant_list(self):
        return self.get_axis('variants',keysOnly=True)
        
    def get_axis(self,orderByType=None,keysOnly=False): # gets the variant/col or kit/row names (and optionally order info too)
        if orderByType in ['variants','kits']:
            if orderByType == 'variants' : SCH = self.VARIANTS
            if orderByType == 'kits' : SCH = self.KITS
            if keysOnly:
                return [i[0] for i in sorted(SCH.items(), key=lambda e: e[1][1])]
            else:
                return sorted(SCH.items(), key=lambda e: e[1][1])
        if orderByType in ['kp','kn','kx','vp','vn','vx']:
            if keysOnly:
                return list(OrderedDict(sorted(self.CNTS[orderByType].items(), key=lambda item: item[1],reverse=True)).keys())
            else:
                listByCount = list(OrderedDict(sorted(self.CNTS[orderByType].items(), key=lambda item: item[1],reverse=True)).keys())
                if orderByType in ['vp','vn','vx']:
                    return [(key, self.VARIANTS[key]) for key in listByCount]
                if orderByType in ['kp','kn','kx']:
                    return [(key, self.KITS[key]) for key in listByCount]
        
    def set_new_order(self,val,cnt,kitType=False,variantType=False):
        if kitType:
            self.KITS[val][1] = cnt
        if variantType:
            self.VARIANTS[val][1] = cnt
        
    def get_variant_name_by_order(self,Y): #get variant name
        variant = None
        for itm in list(self.VARIANTS.items()):
            if itm[1][1] == Y:
                variant = itm[0]
                break
        return variant
        
    def get_kit_name_by_order(self,X): #get kit name
        kit = None
        for itm in list(self.KITS.items()):
            if itm[1][1] == X:
                kit = itm[0]
                break
        return kit
        
    def get_variant_order_by_name(self,y_val): #get variant order from its name
        return self.VARIANTS[y_val][1]
        
    def get_kit_order_by_name(self,x_val): #get kit order from its name
        return self.KITS[x_val][1]
        
    def get_matrix_row_data(self,row_order=None,name=None): #get same type variant data for each kit
        if name is not None:
            row_order = self.get_variant_order_by_name(name)
        if row_order is not None:
            return self.NP[row_order,]
        
    def get_matrix_col_data(self,col_order=None,name=None): #get all type variant data for one kit
        if name is not None:
            col_order = self.get_kit_order_by_name(name)
        if col_order is not None:
            return self.NP[:,col_order].T
        
    def get_matrix_row_indices_by_val(self,val,row_order=None,name=None): #like get_matrix_row_data but retrieves index info for given val
        if name is not None:
            row_order = self.get_variant_order_by_name(name)
        if row_order is not None:
            return np.argwhere(self.NP[row_order,] == val).T[1,]
        
    def get_matrix_col_indices_by_val(self,val,col_order=None,name=None): #like get_matrix_col_data but retrieves index info for given val
        if name is not None:
            col_order = self.get_kit_order_by_name(name)
        if col_order is not None:
            return np.argwhere(self.NP[:,col_order] == val)

    def get_lowest_superset_variant(self,order=None,name=None): #order is variant's order in matrix, name is variant name
        if name is not None:
            order = self.get_variant_order_by_name(name)
        pos_conditions = self.get_matrix_row_indices_by_val(1,row_order=order)
        VAR1 = np.argwhere(self.NP[:,pos_conditions]==1)[:,0]
        unique_elements, counts_elements = np.unique(VAR1, return_counts=True)
        VAR2 = np.asarray((unique_elements, counts_elements)).T
        VAR3 = VAR2[VAR2[:,1]==len(pos_conditions)][:,0]
        idx = np.argwhere(VAR3==order) # idx - make sure we exclude the incoming variant
        min_superset_pos_cnt = 0 # #default
        min_superset_variant = None #default
        for super_v in np.delete(VAR3, idx): # here we use idx
            tmp_name = self.get_variant_name_by_order(super_v)
            #print(self.CNTS['vp'][name])
            if self.CNTS['vp'][tmp_name] > len(pos_conditions): #has to be bigger than the default
                if min_superset_pos_cnt == 0 or self.CNTS['vp'][tmp_name] < min_superset_pos_cnt:
                    min_superset_pos_cnt = self.CNTS['vp'][tmp_name]
                    min_superset_variant = tmp_name
        return min_superset_variant

    def get_matrix_data(self):

        #sql 
        sql2 = '''
            SELECT C.kit_id, V.name, C.assigned, V.variant_id
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc
            ORDER by 4
            '''

        #get data
        self.dbo.sql_exec(sql2)
        F = self.dbo.fetchall()

        #retrieve data from sqlite like so: [V][K] [x,x,x,x,x,x,...]
        DATA = OrderedDict()
        for row in F:
            if row[1] not in DATA:
                DATA[row[1]] = []
            DATA[row[1]].append(row[2])
        
        #create numpy version of data
        for key,value in DATA.items():
            self.NP = np.matrix(list(DATA.values()))

        #all unique kits
        self.KITS = {}
        cnt=0
        for k in sorted(list(set([itm[0] for itm in F]))):
            self.KITS[k]=[cnt,cnt] #default sort,custom sort 
            cnt = cnt + 1

        if 1 == 1: # hack - to get variants default sort just like Iain's PDF
            sql1 = '''
                SELECT distinct V.name, V.variant_id
                FROM s_calls C, s_variants V
                WHERE C.variant_loc = V.variant_loc
                ORDER by 2
                '''
            self.dbo.sql_exec(sql1)
            F = self.dbo.fetchall()

        #all unique variants
        self.VARIANTS = {}
        cnt=0
        for itm in F:
            self.VARIANTS[itm[0]]=[cnt,cnt] #default sort,custom sort 
            cnt = cnt + 1

        #get count data
        self.get_matrix_count_data()

        #get relations data
        self.get_matrix_relations_data()
        
    def get_matrix_count_data(self):

        #vars
        self.CNTS = {}
        sqlc = {}

        #sql - cnt variants
        sqlc['vp'] = '''
            SELECT count(V.name), V.name
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc and C.assigned = 1 
            group by 2;
            '''
        sqlc['vn'] = '''
            SELECT count(V.name), V.name
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc and C.assigned = -1
            group by 2;
            '''
        sqlc['vx'] = '''
            SELECT count(V.name), V.name
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc and C.assigned = 0
            group by 2;
            '''

        #sql - cnt kits
        sqlc['kp'] = '''
            SELECT count(C.kit_id), C.kit_id
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc and C.assigned = 1 
            group by 2;
            '''
        sqlc['kn'] = '''
            SELECT count(C.kit_id), C.kit_id
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc and C.assigned = -1
            group by 2;
            '''
        sqlc['kx'] = '''
            SELECT count(C.kit_id), C.kit_id
            FROM s_calls C, s_variants V
            WHERE C.variant_loc = V.variant_loc and C.assigned = 0
            group by 2;
            '''
        #get all cnts
        for key, sql in sqlc.items():
            self.CNTS[key] = {}
            self.dbo.sql_exec(sql)
            F = self.dbo.fetchall()
            for itm in F:
                self.CNTS[key][itm[1]] = itm[0]
        
    def get_matrix_relations_data(self):

        #sql - get negatives (without kits){{{

        sql = '''
            SELECT distinct V1.name, V2.name
            FROM s_calls C1, s_calls C2,
            s_variants V1, s_variants V2
            WHERE
            C1.kit_id = C2.kit_id AND
            C1.assigned = 1 AND
            C2.assigned = -1 AND
            C1.variant_loc = V1.variant_loc AND
            C2.variant_loc = V2.variant_loc
            ORDER by 1,2;
            '''

        self.dbo.sql_exec(sql)
        self.NEGA = self.dbo.fetchall()

        #}}}
        #sql - get positives (without kits) {{{

        sql = '''
            SELECT distinct V1.name, V2.name
            FROM s_calls C1, s_calls C2,
            s_variants V1, s_variants V2
            WHERE
            C1.kit_id = C2.kit_id AND
            C1.assigned = 1 AND
            C2.assigned = 1 AND
            C1.variant_loc = V1.variant_loc AND
            C2.variant_loc = V2.variant_loc AND
            V1.variant_loc != V2.variant_loc AND
            V1.name||'|'||V2.name NOT IN (
                SELECT distinct QV1.name||'|'||QV2.name as name
                FROM s_calls QC1, s_calls QC2,
                s_variants QV1, s_variants QV2
                WHERE
                QC1.kit_id = QC2.kit_id AND
                QC1.assigned = 1 AND
                QC2.assigned = -1 AND
                QC1.variant_loc = QV1.variant_loc AND
                QC2.variant_loc = QV2.variant_loc)
            ORDER by 1,2;
            '''
        self.dbo.sql_exec(sql)
        self.POSA = self.dbo.fetchall()

        #}}}
        #sql - get mixes (without kits){{{

        sql = '''
            SELECT distinct V1.name, V2.name
            FROM s_calls C1, s_calls C2,
            s_variants V1, s_variants V2
            WHERE
            C1.kit_id = C2.kit_id AND
            C1.assigned = 1 AND
            C2.assigned = 1 AND
            C1.variant_loc = V1.variant_loc AND
            C2.variant_loc = V2.variant_loc AND
            V1.variant_loc != V2.variant_loc AND
            V1.name||'|'||V2.name IN (
                SELECT distinct QV1.name||'|'||QV2.name as name
                FROM s_calls QC1, s_calls QC2,
                s_variants QV1, s_variants QV2
                WHERE
                QC1.kit_id = QC2.kit_id AND
                QC1.assigned = 1 AND
                QC2.assigned = -1 AND
                QC1.variant_loc = QV1.variant_loc AND
                QC2.variant_loc = QV2.variant_loc)
            ORDER by 1,2;
            '''

        self.dbo.sql_exec(sql)
        self.MIXA = self.dbo.fetchall()

        #}}}

    def _bak_sort_step3(self):
        print("")
        print("Processing Nones:")
        print("----------------")
        #variant list that have kits with negative (zero) values
        zlist = np.unique(np.argwhere(self.NP == 0)[:,0]).tolist()
        #iterate all None situations
        for non in ((np.argwhere(self.NP == -1)).tolist()):
            if non[0] in zlist:
                for itm in list(self.VARIANTS.items()):
                    if itm[1][1] == non[0]:
                        variant = itm[0]
                        break
                for itm in list(self.KITS.items()):
                    if itm[1][1] == non[1]:
                        kit = itm[0]
                        break
                print("kit:"+str(kit)+",variant:"+str(variant))

        print("")

        #print(np.nonzero(self.NP == 1)[:1])
        #print(np.nonzero(self.NP == 1))
        #print(((np.argwhere(self.NP == 1)))
        #sys.exit()
        #print(np.where(self.NP == 0))
        #np.multiply(a,b)
        #>>> a = np.array([[1,2],[3,4]])
        #>>> b = np.array([[5],[7]])
        #>>> np.multiply(a,b)

        #class myarray(np.ndarray):
        #    def __new__(cls, *args, **kwargs):
        #        return np.array(*args, **kwargs).view(myarray)
        #    def index(self, value):
        #        return np.where(self==value)

        #prep
        NP = np.copy(self.NP.T)
        N10 = np.arange(1,16)
        NX = N10*NP
        #NXt = np.copy(NX.T)
        print(list(itertools.permutations(NX, 2)))

        sys.exit()
        print(NX)
        sys.exit()
        print("")

        #(beg)is this really getting me combinations?
        m,n = NX.shape
        #x = np.array([x for i in range(m) for x in itertools.product(NX[i, 0 : 1], NX[i, 1 : n])])
        x = np.array([x for i in range(m) for x in itertools.product(NX[i, 0 : n], NX[i, 1 : n])])
        #(end)is this really getting me combinations?

        print(x)
        sys.exit()

        NA = x[np.all(x != 0, axis=1)] #remove zeros

        #print(np.unique(NA), axis=0) #dupes
        print(np.unique(NA)) #dupes

        #(beg)what is this doing?
        #b = np.ascontiguousarray(NA).view(np.dtype((np.void, NA.dtype.itemsize * NA.shape[1])))
        #_, idx = np.unique(b, return_index=True)
        #unique_a = NA[idx]
        #print(unique_a)
        #(end)what is this doing?

        #y = np.ascontiguousarray(x).view(np.dtype((np.void, x.dtype.itemsize * x.shape[1])))
        #_, idx = np.unique(y, return_index=True)
        #print(y)

        #https://exceptionshub.com/permutations-with-unique-values.html
        #https://stackoverflow.com/questions/19744542/itertools-product-eliminating-repeated-elements
        #https://stackoverflow.com/questions/20764926/combinations-without-using-itertools-combinations
        #https://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        #https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-87.php
        #https://stackoverflow.com/questions/38187286/find-unique-pairs-in-list-of-pairs
        #https://stackoverflow.com/questions/29585379/efficient-way-of-making-a-list-of-pairs-from-an-array-in-numpy    

        #print(NXt)
        #print(NP.shape)
        #print(NP.shape[1])

        #m * c[:, np.newaxis]
        #https://stackoverflow.com/questions/18522216/multiplying-across-in-a-numpy-array
        #NP
        #NP[NP > 0] = 5
        #print(NP)

        #self.get_matrix_relations_data()
        #Note: these return in two arrays Array-Y(Variant-Axis),Array-X(kitAxis)
        #print(np.nonzero(self.NP == 1))
        #print(np.where(self.NP == 0))
        #print(self.NP)
        sys.exit()
        
    def _bak_get_matrix_relations_data(self):
        self.MDATA = {}
        
        '''
        for K,V in self.get_axis('variants'):

            posY = list(list(np.nonzero(self.NP == 1))[0])

            (1) grab all ones (x,y) coords
            (2) what are the unique y's in that (these are the known positive variants)
            (3) per row ... isolate those situations where we see a
            (4) what are the unique combos of 1,1's (per row)
            (5) what are the unique combos of 1,0's (per row)

            3 is k1
            4 is k2

            6 is A col
            7 is B col

            0 is false
            1 is true

            k = A+, k = B-
            k = A+, k = B-

            convert matrix to -1,1,None 
            then multiple 1/-1 by variant order num (+1)
            may not need -> P = get uniq k+
            P = get uniq k1+,k2+ combos (as names or order ids) --> at least one pos num (no Nones)
            rP = get reverse P combos (as names or order ids)
            uP = uniq(P + rP)
            N = get uniq k1+,k2- combos (as names or order ids)
            intersection(uP + N) = mix
            outside(uP + N - favor uP) = pos
            outside(uP + N - favor N) = neg
            
            ----
            v,k,b
            ----
            get uniq x|6|0
             A:-
            get uniq x|6|1
             A:3|6|1 > what 3|(not 6)|0 exist?
             A:4|6|1 > what 4|(not 6)|0 exist? (first case, ... becomes mix for 6 -- and 6 locked down)
             A:5|6|1 > what 5|(not 6)|0 exist?
            get uniq x|7|0
             A:3|7|0
            get uniq x|7|1
             A:4|7|1
            3,6,1  4,6,1    3,7,1    4,7,0
            A+     A+       B+       B-
            ----
            A mix{B}
            filter on all 6's ... is there a 0|1?

            ----
            v,k,b
            ----
            3,6,1  4,6,1    3,7,1    4,7,1
            A+     A+       B+       B+
            ----
            A pos{B}

            https://stackoverflow.com/questions/38187286/find-unique-pairs-in-list-of-pairs
            https://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array/16973510#16973510
            https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays

            find unique(pos-pos in same kits)
            find unique(pos-neg in same kits)
        '''

        #mix ->A+ on k5,k6,k7 .... B- k5, B+ k7
        #for each variant:
        #    get me the set of kit rows where we see positive vals
        #        for these combo of rows:
        #            what are the unique other pos variants
        #              - is there at least one pos, no neg: pos
        #              - is there at least one pos, one neg: mix
        #            what are the unique neg variants that we haven't seen yet?
        #              - these are all negs

        #{{{
        #Note: these return in two arrays Array-Y(Variant-Axis),Array-X(kitAxis)
        #like this: (array([ 0,x,x,x...]),array([x,x,,...]))
        #posX = list(list(np.nonzero(self.NP == 1))[1])
        #posY = list(list(np.nonzero(self.NP == 1))[0])
        #negX = list(list(np.where(self.NP == 0))[1])
        #negY = list(list(np.where(self.NP == 0))[0])
        #loop pos values 
        #for cnt1 in range(len(posX)):
        #    V1 = self.get_variant_name_by_order(posY[cnt1])
        #    X1 = posX[cnt]
        #    for cnt2 in range(len(negX)):
        #        V2 = self.get_variant_name_by_order(negY[cnt2])
        #        X2 = negX[cnt2]
        #        if X1==X2 and V1 != V2:
        #            self.MDATA[V1]['mix'].append(V2)
        #}}}

    # tree

    def sort_tree(self):

        #prep data
        self.sort_tree_prep_data()
                    
        #end collapse vim marker
        debug_chk('DEBUG_TREE',"}"+"}}",2)

        #build unsorted tree with all nodes under top
        self.TREE = {}
        self.TREE['_'] = Node("_")
        self.TREE['top'] = Node('top', parent=self.TREE['_'])
        self.TREE['dupes'] = Node('dupes', parent=self.TREE['_'])
        for key, value in self.TDATA.items():
            self.TREE[key] = Node(key, parent=self.TREE['top'])

        #sort it
        self.sort_variant_tree(run_mix=True,run=1)
        self.sort_variant_tree(run_mix=True,run=2)
        self.sort_variant_tree(run_pos=True,run=3)

        #self.sort_variant_tree(run_neg=True,run=3)
        #self.sort_variant_tree(run_all=True,run=4)
        #self.sort_variant_tree(run_neg=True,run=5)
        #self.sort_variant_tree(run_all=True,run=6)
        #self.sort_variant_tree(run_all=True,run=7)
        #cnt = 2 
        #while cnt < config['TREE_SORT_RUN_CNT']:
        #    self.sort_variant_tree(run_all=True,run=cnt)
        #    cnt = cnt + 1

        sys.exit()

        #json/stdout a variable (debugging)
        self.stdout_dump_var(newV2)
        sys.exit()
        
    def sort_variant_tree (self,run_mix=False,run_pos=False,run_neg=False,run_all=None,run=1):

        #init sort logging and var prep
        debug_chk('DEBUG_TREE',"===",1)
        debug_chk('DEBUG_TREE',"RUN:"+str(run),1)

        #show pre-proc default data 
        if config['DEBUG_TREE']>3:
            self.stdout_variant_relations_data(self.TDATA,'STASH pre-proc',run)

        #prep ref
        if self.REF is None:
            self.REF = {}
            for key, value in self.TDATA.items():
                self.REF[key] = {'mix':[],'pos':[],'neg':[],'dup':[]}

        self.mix_rule_chks(run_mix or run_all,run)
        self.pos_rule_chks(run_pos or run_all,run)
        self.neg_rule_chks(run_neg or run_all,run)

        #show the final tree diagram after run completion
        debug_chk('DEBUG_TREE',"---",1)
        debug_chk('DEBUG_TREE',"RUN:"+str(run)+" DONE",1)
        for pre, fill, node in RenderTree(self.TREE['top']):
            #print("%s%s" % (pre, node.name))
            debug_chk('DEBUG_TREE',"%s%s" % (pre, node.name),1)
        if len(self.TREE['dupes'].children):
            for pre, fill, node in RenderTree(self.TREE['dupes']):
                #print("%s%s" % (pre, node.name))
                debug_chk('DEBUG_TREE',"%s%s" % (pre, node.name),1)

        #show post-proc remaining data that didn't get complete (now in STASH)
        if config['DEBUG_TREE']>3:
            self.stdout_variant_relations_data(self.TDATA,'STASH post-proc',run)
            self.stdout_variant_relations_data(self.REF,'REF post-proc',run)
        #debug_chk('DEBUG_TREE',"---",3)
        

    def mix_rule_chks(self,run_flg,run,hardFlg=False):
        debug_chk('DEBUG_TREE',"---",2)
        debug_chk('DEBUG_TREE',"mix-checks {"+"{{",2) #beg collapse vim marker
        debug_chk('DEBUG_TREE',"---",2)
        if run_flg:
            for key, value in self.TDATA.items():
                for Vz in value['mix']:
                    #chk1 - if the two nodes are on the same level, then: create a default parental relation + don't STASH
                    if self.TREE[Vz].parent == self.TREE[key].parent and self.TREE[key].parent != self.TREE['dupes'] :
                        debug_chk('DEBUG_TREE',"MIX-CHK1 - "+key+"|"+Vz+" - parents are same level - so put "+Vz+" under "+key,4)
                        ch1 = list(self.TREE[Vz].parent.children).remove(self.TREE[Vz])
                        ch2 = self.TREE[Vz].children
                        if ch1 is None:
                            self.TREE[Vz].parent = None
                        else:
                            self.TREE[Vz].parent.children = tuple(ch1)
                        self.TREE[Vz] = Node(Vz, parent=self.TREE[key])
                        self.TREE[Vz].children = ch2
                        self.TDATA[key]['mix'].remove(Vz)
                        if key not in self.REF[Vz]['mix']:
                            self.REF[Vz]['mix'].append(key)
                    #chk2 - if there is already a good direct line established, then: don't STASH
                    elif self.TREE[Vz] in self.TREE[key].descendants:
                        debug_chk('DEBUG_TREE',"MIX-CHK2 - "+key+"|"+Vz+" - condition satisfied - "+Vz+" already under "+key,4)
                        self.TDATA[key]['mix'].remove(Vz)
                        if key not in self.REF[Vz]['mix']:
                            self.REF[Vz]['mix'].append(key)
                    #chk3 - if there's one dupe and there is already a good direct line established, then: don't STASH
                    elif len(self.REF[key]['dup']) > 0 and self.TREE[Vz] in self.TREE[self.REF[key]['dup'][0]].descendants:
                        debug_chk('DEBUG_TREE',"MIX-CHK3 - "+key+"|"+Vz+" - (dupe) condition satisfied - "+Vz+" already under "+key,4)
                        self.TDATA[key]['mix'].remove(Vz)
                        if key not in self.REF[Vz]['mix']:
                            self.REF[Vz]['mix'].append(key)
                    #chk4 - if there's one dupe and there is already a good direct line established, then: don't STASH
                    elif len(self.REF[Vz]['dup']) > 0 and self.TREE[self.REF[Vz]['dup'][0]] in self.TREE[key].descendants:
                        debug_chk('DEBUG_TREE',"MIX-CHK4 - "+key+"|"+Vz+" - (dupe) condition satisfied - "+Vz+" already under "+key,4)
                        self.TDATA[key]['mix'].remove(Vz)
                        if key not in self.REF[Vz]['mix']:
                            self.REF[Vz]['mix'].append(key)
                    #chk5 - if both are dupes and there is already a good direct line established, then: don't STASH
                    elif len(self.REF[Vz]['dup']) > 0 and len(self.REF[key]['dup']) > 0 and self.TREE[self.REF[Vz]['dup'][0]] in self.TREE[self.REF[key]['dup'][0]].descendants:
                        debug_chk('DEBUG_TREE',"MIX-CHK5 - "+key+"|"+Vz+" - (dupe) condition satisfied - "+Vz+" already under "+key,4)
                        self.TDATA[key]['mix'].remove(Vz)
                        if key not in self.REF[Vz]['mix']:
                            self.REF[Vz]['mix'].append(key)
                    #chk6 - if anything else, then: STASH
                    else:
                        debug_chk('DEBUG_TREE',"MIX-CHK6 - "+key+"|"+Vz+" - parents not same level and not direct lineage - put in STASH",2)

        #self.stdout_variant_relations_data(self.TDATA,'DEBUG','5')
        #sys.exit()
        debug_chk('DEBUG_TREE',"}"+"}}",2)  #end collapse vim marker
        
    def pos_rule_chks(self,run_flg,run):
        debug_chk('DEBUG_TREE',"pos-checks {"+"{{",2)  #beg collapse vim marker
        debug_chk('DEBUG_TREE',"---",2)
        if run_flg:
            for key, value in self.TDATA.items():
                for Vz in value['pos']:
                    #chk1 - if the two nodes have direct line relation, then: don't STASH
                    if self.TREE[Vz] in self.TREE[key].descendants or self.TREE[key] in self.TREE[Vz].descendants:
                        debug_chk('DEBUG_TREE',"POS-CHK1 - "+key+"|"+Vz+" - direct lineage relation found",4)
                        self.TDATA[key]['pos'].remove(Vz)
                        if key not in self.REF[Vz]['pos']:
                            self.REF[Vz]['pos'].append(key)
                    #chk2a - if one of the two nodes have direct line relation via dupe, then: don't STASH
                    elif len(self.REF[key]['dup']) > 0 and self.TREE[Vz] in self.TREE[self.REF[key]['dup'][0]].descendants:
                        debug_chk('DEBUG_TREE',"POS-CHK2a - "+key+"|"+Vz+" - (via dupe) direct lineage relation found",4)
                        self.TDATA[key]['pos'].remove(Vz)
                        if key not in self.REF[Vz]['pos']:
                            self.REF[Vz]['pos'].append(key)
                    #chk2b - if one of the two nodes have direct line relation via dupe, then: don't STASH
                    elif len(self.REF[key]['dup']) > 0 and self.TREE[self.REF[key]['dup'][0]] in self.TREE[Vz].descendants:
                        debug_chk('DEBUG_TREE',"POS-CHK2b - "+key+"|"+Vz+" - (via dupe) direct lineage relation found",4)
                        self.TDATA[key]['pos'].remove(Vz)
                        if key not in self.REF[Vz]['pos']:
                            self.REF[Vz]['pos'].append(key)
                    #chk3a - if one of the two nodes have direct line relation via dupe, then: don't STASH
                    elif len(self.REF[Vz]['dup']) > 0 and self.TREE[key] in self.TREE[self.REF[Vz]['dup'][0]].descendants:
                        debug_chk('DEBUG_TREE',"POS-CHK3a - "+key+"|"+Vz+" - (via dupe) direct lineage relation found",4)
                        self.TDATA[key]['pos'].remove(Vz)
                        if key not in self.REF[Vz]['pos']:
                            self.REF[Vz]['pos'].append(key)
                    #chk3b - if one of the two nodes have direct line relation via dupe, then: don't STASH
                    elif len(self.REF[Vz]['dup']) > 0 and self.TREE[self.REF[Vz]['dup'][0]] in self.TREE[key].descendants:
                        debug_chk('DEBUG_TREE',"POS-CHK3b - "+key+"|"+Vz+" - (via dupe) direct lineage relation found",4)
                        self.TDATA[key]['pos'].remove(Vz)
                        if key not in self.REF[Vz]['pos']:
                            self.REF[Vz]['pos'].append(key)
                    #chk4 - if both nodes have direct line relation via dupe, then: don't STASH
                    elif len(self.REF[key]['dup']) > 0 and len(self.REF[Vz]['dup']) > 0:
                        if self.TREE[self.REF[key]['dup'][0]] in self.TREE[self.REF[Vz]['dup'][0]].descendants:
                            debug_chk('DEBUG_TREE',"POS-CHK4a - "+key+"|"+Vz+" - (via dupe) direct lineage relation found",4)
                            self.TDATA[key]['pos'].remove(Vz)
                            if key not in self.REF[Vz]['pos']:
                                self.REF[Vz]['pos'].append(key)
                        if self.TREE[self.REF[Vz]['dup'][0]] in self.TREE[self.REF[key]['dup'][0]].descendants:
                            debug_chk('DEBUG_TREE',"POS-CHK4b - "+key+"|"+Vz+" - (via dupe) direct lineage relation found",4)
                            self.TDATA[key]['pos'].remove(Vz)
                            if key not in self.REF[Vz]['pos']:
                                self.REF[Vz]['pos'].append(key)
                    #chk5 - if the two nodes don't have direct line relation ... are they perhaps dupes? if so, don't STASH 
                    elif run>1 and self.dupe_variant_check(key,Vz):
                        #self.TREE[Vz] not in self.TREE[key].descendants and not in self.TREE[key] in self.TREE[Vz].descendants:
                        debug_chk('DEBUG_TREE',"POS-CHK5 - "+key+"|"+Vz+" - dupe relation found",4)
                        self.TDATA[key]['pos'].remove(Vz)
                        if key not in self.REF[Vz]['pos']:
                            self.REF[Vz]['pos'].append(key)
                    #chk6 - other situations 
                    else:
                        debug_chk('DEBUG_TREE',"POS-CHK6 - "+key+"|"+Vz+" - other situations - put in STASH",2)

        debug_chk('DEBUG_TREE',"}"+"}}",2)  #end collapse vim marker
        
    def neg_rule_chks(self,run_flg,run):
        debug_chk('DEBUG_TREE',"neg-checks {"+"{{",2)  #beg collapse vim marker
        debug_chk('DEBUG_TREE',"---",2)
        if run_flg:
            for key, value in self.TDATA.items():
                for Vz in value['neg']:
                    #chk1 - if the two nodes don't have direct line relation, then: don't STASH
                    if self.TREE[Vz] not in self.TREE[key].descendants and self.TREE[key] not in self.TREE[Vz].descendants:
                        debug_chk('DEBUG_TREE',"NEG-CHK1 - "+key+"|"+Vz+" - no direct lineage relation found",4)
                        self.TDATA[key]['neg'].remove(Vz)
                        if key not in self.REF[Vz]['neg']:
                            self.REF[Vz]['neg'].append(key)
                    #chk2 - if the two nodes don't have direct line relation, then: don't STASH
                    elif self.TREE[Vz] in self.TREE[key].descendants:
                        debug_chk('DEBUG_TREE',"NEG-CHK2 - "+key+"|"+Vz+" - anc to dec relation found",4)
                        self.TDATA[key]['neg'].remove(Vz)
                        if key not in self.REF[Vz]['neg']:
                            self.REF[Vz]['neg'].append(key)
                    #chk3 - if anything else, then: STASH
                    else:
                        debug_chk('DEBUG_TREE',"NEG-CHK3 - "+key+"|"+Vz+" - direct lineage found - put in STASH",2)
            else:
                #since we're skipping "neg" relations, they need to go to STASH
                self.TDATA[key]['neg'] = value['neg']
        debug_chk('DEBUG_TREE',"}"+"}}",2) #end collapse vim marker
        #return STASH
        
    def dupe_variant_check(self,variant1,variant2): #requires run > 1
        if self.REF[variant2]['dup'] == [variant1]:
            return True #don't stash
        if self.REF[variant1]['dup'] == [variant2]:
            return True #don't stash
        #when siblings
        if self.TREE[variant1].parent == self.TREE[variant2].parent:
            mixList1 = []
            mixList2 = []
            chk = 0
            #check processed variants
            for k,v in self.REF.items():
                #nothing processed previously had variant1 as a child
                if k not in [variant1,variant2] and variant1 in v['mix']:
                    mixList1.append(k)
                #nothing processed previously had variant2 as a child
                if k not in [variant1,variant2] and variant2 in v['mix']:
                    mixList2.append(k)
            if len(mixList1) == 0 and len(mixList2)>0: #variant1 references
                for k in mixList2:
                    chk = 1
                    if self.TREE[variant2] not in self.TREE[k].ancestors:
                        chk = 0
                #chk = 1 : v2 is part of all k's ancestors (so v1 is a dupe of v2)
            if len(mixList2) == 0 and len(mixList1)>0: #variant2 references 
                for k in mixList1:
                    chk = 2
                    #print(str(variant1)+":"+str(k))
                    if self.TREE[variant1] not in self.TREE[k].ancestors:
                        chk = 0
                #chk == 2: v1 is part of all k's ancestors (so v2 is a dupe of v1)
            if chk == 1:
                #self.REF[variant2]['dup'] = list(set(self.REF[variant2]['DUP']+[variant1]))
                self.REF[variant1]['dup'] = [variant2]
                self.TREE[variant1].parent = self.TREE['dupes']
            if chk == 2:
                #self.REF[variant1]['dup'] = list(set(self.REF[variant1]['DUP']+[variant2]))
                self.REF[variant2]['dup'] = [variant1]
                self.TREE[variant2].parent = self.TREE['dupes']
        if self.REF[variant2]['dup'] == [variant1]:
            return True #don't stash
        elif self.REF[variant1]['dup'] == [variant2]:
            return True #don't stash
        else:
            return False #stash
        #HIDE-ME {{{
        #print(chk)
        #sys.exit()
        #print(mixList1) #['B', 'C', 'D', 'E', 'F', 'G', 'I', 'J', 'K', 'L', 'N', 'O']
        #for k in mixList1:
        #    print(k+":"+str(self.TREE[k].ancestors))
        #print(mixList2) #[]
        #for k in mixList2:
        #    print(k+":"+str(self.TREE[k].ancestors))
        #...
        #print(variant1)
        #print("stash:"+str(STASH[variant1]))
        #print("ref:"+str(self.REF[variant1]))
        #...
        #print(variant2)
        #print("stash:"+str(STASH[variant2]))
        #print("ref:"+str(self.REF[variant2]))
        #}}}
        #...
        sys.exit()
        #self.TREE[Vz] not in self.TREE[key].descendants and not in self.TREE[key] in self.TREE[Vz].descendants:

    def stdout_variant_relations_data(self,DATA,dataStr,run=1):

        mixlen = 0
        poslen = 0
        neglen = 0
        duplen = 0

        #print the counts
        print("---")
        print(dataStr+"{{"+"{") #beg vim marker
        for key, value in DATA.items():
            mixlen = mixlen+len(value['mix'])
            poslen = poslen+len(value['pos'])
            neglen = neglen+len(value['neg'])
            duplen = neglen+len(value['dup'])
            #DUPlen = neglen+len(value['DUP'])
        print("RUN:"+str(run)+"("+dataStr+") - mix cnt:"+str(mixlen))
        print("RUN:"+str(run)+"("+dataStr+") - pos cnt:"+str(poslen))
        print("RUN:"+str(run)+"("+dataStr+") - neg cnt:"+str(neglen))
        print("RUN:"+str(run)+"("+dataStr+") - dup cnt:"+str(duplen))
        #print("RUN:"+str(run)+"("+dataStr+") - DUP cnt:"+str(DUPlen))
        print("")

        #print the data
        print("RUN:"+str(run)+"("+dataStr+") - data")
        for key, value in DATA.items():
            print(key+'|'+str(value).replace("'","").replace(" ",""))
        print("}}"+"}") #end vim marker
        

    def sort_tree_prep_data(self):

        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #beg collapse vim marker
        debug_chk('DEBUG_TREE',"PREP {"+"{{",2)

        #all kits, variant, assignment mixes 

        #Letters 
        if self.TREE_MODE == 1:
            sql = "select C.kit_id,C.variant_loc,C.assigned from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 1,2,3"
        #Names
        if self.TREE_MODE == 2:
            sql = "select C.kit_id,V.name,C.assigned from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 1,2,3"
        #Combo - Names+Letters
        if self.TREE_MODE == 3:
            sql = "select C.kit_id,'('||C.variant_loc||') '||V.name,C.assigned from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 1,2,3"

        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        #print(F)
        #sys.exit()

        #all unique kits
        #print("---")
        #print("all unique kits")
        #KITS = sorted(list(set([itm[0] for itm in F])))
        #print(KITS)

        #all unique variants
        VARIANTS = sorted(list(set([itm[1] for itm in F])))
        #VARIANTSp = ['+'+itm[1] for itm in F]
        #VARIANTSn = ['-'+itm[1] for itm in F]
        #VARIANTSa = sorted(list(set(VARIANTSp+VARIANTSn)))
        debug_chk('DEBUG_TREE',"---",5)
        debug_chk('DEBUG_TREE',"all unique variants",5)
        debug_chk('DEBUG_TREE',VARIANTS,5)
        #print("---")
        #print("all unique variants - pos+neg")
        #print(VARIANTSa)
        #sys.exit()
        
        #kits with positive assignments
        Fp = sorted(list(set([i[0] for i in list(filter(lambda x: x[2]==1, F))])))
        debug_chk('DEBUG_TREE',"---",5)
        debug_chk('DEBUG_TREE',"kits with positive assignment variant calls",5)
        debug_chk('DEBUG_TREE',Fp,5)

        #kits with negative assignments
        Fn = sorted(list(set([i[0] for i in list(filter(lambda x: x[2]==0, F))])))
        debug_chk('DEBUG_TREE',"---",5)
        debug_chk('DEBUG_TREE',"kits with positive negative variant calls",5)
        debug_chk('DEBUG_TREE',Fn,5)

        #per all the kits with positive variants (build new dict)
        debug_chk('DEBUG_TREE',"---",5)
        debug_chk('DEBUG_TREE',"dict of kits with their positive assignment variant calls",5)
        KA={}
        for k in Fp:
            Kp = sorted(list(set(['+'+i[1] for i in list(filter(lambda x: x[0]==k and x[2]==1, F))])))
            #['A+', 'D+', 'F+', 'H+', 'M+']
            debug_chk('DEBUG_TREE',k+" "+str(Kp),5)
            KA[k] = {'len':len(Kp),'plen':len(Kp),'sort':0,'variants':Kp}

        #per all the kits with negative variants (build new dict)
        debug_chk('DEBUG_TREE',"---",5)
        debug_chk('DEBUG_TREE',"dict of kits with their negative assignment variant calls",5)
        for k in Fn:
            Kn = sorted(list(set(['-'+i[1] for i in list(filter(lambda x: x[0]==k and x[2]==0, F))])))
            debug_chk('DEBUG_TREE',k+" "+str(Kn),5)
            if k in KA.keys():
                KA[k]['len'] = len(KA[k]['variants'])+len(Kn)
                KA[k]['variants'] = sorted(KA[k]['variants']+Kn)
            else:
                KA[k] = {'len':len(Kn),'plen':0,'sort':0,'variants':Kn}

        #loop dict to create list version of the data
        newV1 = []
        for key, value in KA.items():
            newV1.append({'kit':key,'variants':value['variants'],'sort':value['sort'],'len':value['len'],'plen':value['plen']})

        #sort this new list
        cnt = 0
        for d in sorted(newV1, key=lambda k: (k['plen'],k['len']), reverse=True):
            d.update((k, cnt) for k, v in d.items() if k == "sort")
            cnt = cnt + 1

        #create a var for the sorted version (not necessary)
        newV2 = sorted(newV1, key=lambda k: (k['sort']))

        #print to stdout so I can see what I'm doing
        debug_chk('DEBUG_TREE',"---",5)
        debug_chk('DEBUG_TREE',"combined dict of kits with pos+neg variant calls - sorted",5)
        #newV3 = {}
        for d in newV2:
            #newV3[d['kit']] = d['variants']
            STR = d['kit']+':'+str(d['variants'])
            debug_chk('DEBUG_TREE',STR.replace("'",""),5)

        #build variant relationship data that we need for sorting
        debug_chk('DEBUG_TREE',"---",2)
        self.TDATA = {}
        for VX in VARIANTS:
            #DATA[VX] = {'mix':[],'pos':[],'neg':[],'dup':[],'DUP':[]}
            self.TDATA[VX] = {'mix':[],'pos':[],'neg':[],'dup':[]}
            VXP = '+'+VX
            for VY in VARIANTS:
                VYP = '+'+VY
                if VXP != VYP:
                    VYN = '-'+VY
                    chk1 = False
                    chk2 = False
                    for d in newV2:
                        if VXP in d['variants']:
                            if chk1 is False and VYP in d['variants']:
                                chk1 = True
                            if chk2 is False and VYN in d['variants']:
                                chk2 = True
                        if chk1 is True and chk2 is True:
                            self.TDATA[VX]['mix'].append(VY)
                            break
                    if chk1 is True and chk2 is False:
                        self.TDATA[VX]['pos'].append(VY)
                    if chk2 is True and chk1 is False:
                        self.TDATA[VX]['neg'].append(VY)
          

    # TODO: set up a random approach to pushing data into the sort. troubleshoot results

    #NOTES{{{
        # HIDE-ME: rules {{{
        # -----------------------------------
        # mix: (1|2) means 1 is above 2
        # pos: (1|2) means 1 is a direct ancestor or direct descendant or dupe, not a "cousin", "uncle", or 
        # neg: (1|2) means 1 is a "cousin" or "uncle" or "sibling" or direct ancestor of 2
        # -----------------------------------
        #note: next -- attempt to automate the "rules" to sort the "blocks" data
        #sample data:A|{mix: [B,C,D,E,F,G,I,J,K,L,N,O], pos: [H,M], neg: []}
        # -----------------------------------
        #TODO: need to track +/- in the tree nodes???
        # ----------------------------------- }}}
        # HIDE-ME: results when applying these rules manually:{{{

        # A|{mix: [B,C,D,E,F,G,I,J,K,L,N,O], pos: [H,M], neg: []}
        # =M|{mix: [B,C,D,E,F,G,I,J,K,L,N,O], pos: [A,H], neg: []}
        # =H|{mix: [B,C,D,F,G,I,J,K,L,N,O], pos: [A,E,M], neg: []}
        #     1-J|{mix: [], pos: [A,H,M], neg: [B,C,D,F,G,I,K,L,N,O]}
        #     2-D|{mix: [B,C,E,F,G,I,K,L,N,O], pos: [A,H,M], neg: [J]}
        #         1-F|{mix: [], pos: [A,D,H,M],neg: [B,C,E,G,I,J,K,L,N,O]}
        #         2-C|{mix: [B,G,I,L,N,O], pos: [A,D,E,H,M], neg: [F,J,K]}
        #             1-E|{mix: [], pos: [A,C,D,H,M,O], neg: [B,F,G,I,J,K,L,N]}
        #                 1-L|{mix: [B,G,I,K,?O], pos: [A,C,D,H,M], neg: [F,J,N]}
        #                     1-G|{mix: [N], pos: [A,C,D,H,L,M], neg: [B,F,I,J,K,O]}
        #                         1-N|{mix: [], pos: [A,C,D,G,H,M], neg: [B,F,I,J,K,O]}
        #                     2-O|{mix: [B,I,K,?L], pos: [A,C,D,E,H,M], neg: [F,G,J,N]}
        #                         1-I|{mix: [K], pos: [A,B,C,D,H,L,M,O], neg: [F,G,J,N]}
        #                             1-B|{mix: [K], pos: [A,C,D,H,I,L,M,O], neg: [F,G,J,N]}
        #                                 1-K|{mix: [], pos: [A,B,D,H,I,L,M,O], neg: [F,G,J,N]}

        #}}}
        # HIDE-ME: disputes: {{{
        # (1) L:mix-O vs. O:mix-L 
        #
        # resolutions: 
        # (1) winning rule is L:mix-O
        #
        # reasons for (1) resolution:
        #   (1) L-mix-I
        #   (2) L-mix-B
        #   (3) L-mix-K
        #   (4) B-pos-L
        #   (5) K-pos-L
        #   (6) I-pos-L

        #}}}
        # HIDE-ME: sample output (at the moment) {{{

        # top
        # ├── (A) M343(=)
        # │   ├── (D) U106 (ok)
        # │   │   ├── (C) Z381 (ok)
        # │   │   │   └── (L) A297 (!!!) <--- recurrent rule (see Iain notes)
        # │   │   │       ├── (G) Z156 (ok)
        # │   │   │       │   └── (N) Z306 (ok)
        # │   │   │       └── (O) L48 (ok)
        # │   │   │           ├── (B) Z9 (ok)
        # │   │   │           │   └── (K) Z8 (ok)
        # │   │   │           └── (I) Z28 (!!!) <-- Iain says this is equiv to Z9
        # │   │   ├── (E) Z301 (!!!) <-- should be under (C) Z381 (known Problem#1)
        # │   │   └── (F) Z18 (ok)
        # │   └── (J) P312 (ok)
        # ├── (H) L11 (=)
        # └── (M) M269 (=)

        # TODO: Z301+ exists, but you don't have any Z381+ Z301- tests in the example, 
        # so you have to treat Z381 and Z301 as equivalent]

        # TODO: A297 is a recurrent SNP that doesn't fit into the phylogeny, so
        # could either be dumped or listed as recurrent.

        # TODO: Z28 is equivalent to Z9

        #}}}
        # HIDE-ME: PROBLEM 1: {{{

        # POS C|E
        # POS O|E

        # C,O have positive E values (and not seeing E in direct line when under D directly)
        # the C could actually be ok, if it's a dupe. but the O -- no.

        # how to come up with other ideas?

        # (1) what is the last mix node for E 
        #     (answer: D)
        # (2) what are other desc of D that have direct lines (or are dupes # with) with C and O
        #     (answer: C,L,O,B,I,K)
        #     Arch Need: keep ref of what's been processed like so:
        #     E {'ref-mix': [A,D], 'ref-pos': [], 'ref-neg': [E]}
        # (3) are there any neg for E in those results? (if so -- it can't be
        #     one of them ... or any of their direct lines)
        #     (answer: No)
        # (4  what is the first possibility? ... use dupe options last
        #     (answer: E is C's parent)
        # (5) (ISSUE) what about E's children if it has any -- and their possible conflicts 
        #     if move E?

        # }}}
        # HIDE-ME: the ones I missed: {{{

        # ------------------------
        # OK - F should be under D
        # #1 - E should be under C like this D>C>E
        # (fixed by #1) L should be under E like this D>C>E>L ...
        # B should be under I (B and I should be dupes -- my manual work was wrong)
        # ------------------------

        #}}}
    #}}}

    # misc

    def stdout_dump_var(self,var):
        #TODO: put this somewhere else
        print(json.dumps(var, indent=4, sort_keys=True))

