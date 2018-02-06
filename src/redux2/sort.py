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
#colors {{{

RED = '\033[31m' #Z=1,N=1
GREEN = '\033[32m'
WHITE = '\033[37m'

# }}}

#TODO: order of priority
'''
1. A297 junk rule
2. redo the horiz/vert sort
3. get_matrix_data should exlude all neg variants too
4. hide top from matrix stdout
5. do I still need self.perfect_variants, etc?
6. write good notes and commenting for the rules
'''

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

        self.perfect_variants = None
        self.imperfect_variants = None
        self.imperfect_known_variants = None
        self.imperfect_unknown_variants = None
        self.resolved_variants = []

    # schema / sample data

    def sort_schema(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.dbo.sql_exec_file('sort-schema.sql')
        
    def sort_ins_sample_data(self):

        #kits
        kits = "A B C D E F G H I J"
        for k in kits.split():
            sql = "insert into s_kits (kit_id) values ('%s');" % k
            self.dbo.sql_exec(sql)

        #artificial top
        sql = "INSERT into s_variants (variant_id,variant_loc,name) VALUES (%s,'%s','%s');" % (-999,'top','top')
        self.dbo.sql_exec(sql)
        for k in kits.split():
            sql = "INSERT into s_calls (kit_id,variant_loc,assigned) VALUES ('k%s','%s',%s);" % (k,'top',1)
            self.dbo.sql_exec(sql)

        #variants + calls
        with open(config['REDUX_DATA']+'/sample-sort-data.csv','r') as FILE:
            for row in csv.DictReader(FILE,'vi v n A B C D E F G H I J'.split()):
                row = json.loads(json.dumps(row).replace('\\ufeff','')) #hack: remove byte order mark
                #s_variants
                sql = "INSERT into s_variants (variant_id,variant_loc,name) VALUES (%s,'%s','%s');" % (row['vi'],row['v'],row['n'])
                self.dbo.sql_exec(sql)
                for k in kits.split(): #kit_id
                    kv = str(row[str(k)]) #assigned
                    vv = str(row['v']) #variant_loc
                    #s_calls
                    sql1 = "INSERT into s_calls (kit_id,variant_loc,assigned) VALUES ('k%s','%s',%s);" % (k,vv,kv)
                    self.dbo.sql_exec(sql1)

    # matrix

    def sort_matrix(self):

        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #get data
        self.get_matrix_data()

        #stdout relations data
        if config['DEBUG_RELATIONS']:
            self.stdout_matrix_relations_data()

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

        #step 4
        #debug_chk('DEBUG_MATRIX',"data - step 4",1)
        #self.sort_step4()

        sys.exit()

    def sort_step1(self):
        self.matrix_vertical_sort_new()
        
    def sort_step2(self):
        self.matrix_horizontal_sort()
        
    def sort_step3(self):
        #(beg)stash these 
        self.perfect_variants = self.get_perfect_variants_idx()
        self.imperfect_variants = self.get_imperfect_variants_idx()
        self.imperfect_known_variants = self.get_imperfect_known_variants_idx()
        self.imperfect_unknown_variants = self.get_imperfect_unknown_variants_idx()
        #(end)stash these
        #variant list that have kits with negative (zero) values
        zlist = np.unique(np.argwhere(self.NP == -1)[:,0]).tolist()
        #iterate all None situations
        self.unk_variants = ((np.argwhere(self.NP == 0)).tolist())
        unk_variants = self.unk_variants[:]
        #unresolved list
        print("")
        print("unresolved")
        print(self.get_coord_name_by_order(self.unk_variants))
        print("...")
        print("")
        print("Processing Nones:")
        #loop unk variants
        for unk in unk_variants:
            if unk[0] in zlist: #unk[0] = variant_order, unk[1] = kit_order
                print("----------------")
                coord = self.get_coord(unk[1],unk[0])
                print("%s:[%s,%s] " % (coord,unk[0],unk[1]))
                print("{{"+"{") #beg vim marker
                print("")
                supsetsP = self.get_supset_variants(override_val=1,variant_order=unk[0],kit_order=unk[1],convertToNames=False)
                supsets = self.get_supset_variants(variant_order=unk[0],kit_order=unk[1],convertToNames=False)
                subsetsP = self.get_subset_variants(override_val=1,variant_order=unk[0],kit_order=unk[1],convertToNames=False)
                subsets = self.get_subset_variants(variant_order=unk[0],kit_order=unk[1],convertToNames=False)
                print("%s"%self.test_rule1_subsets(unk_variant=unk,subsets=subsets,supsets=supsets))
                print("")
                print("}}"+"}") #end vim marker

        #unresolved list
        print("----------------")
        print("")
        print("unresolved")
        print(self.get_coord_name_by_order(self.unk_variants))
        #resolved list
        print("resolved")
        print(self.resolved_variants)
        #stdout relations data
        if config['DEBUG_RELATIONS']:
            self.stdout_matrix_relations_data()
        print("")
        #sys.exit()

        self.matrix_vertical_sort()
        self.matrix_horizontal_sort()

    def matrix_vertical_sort_new(self):
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
        
    def matrix_horizontal_sort_new(self):
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
        
    def matrix_vertical_sort(self):
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
        
    def matrix_horizontal_sort(self):
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
        

    def stdout_tbl_matrix(self):
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',"big_matrix view{{"+"{",1)
        debug_chk('DEBUG_MATRIX',"",1)
        table = BeautifulTable()
        table.column_headers = ['top']+self.get_cur_kit_list()
        for K,V in self.get_axis('variants'):
            table.append_row([K]+self.get_numpy_matrix_row_as_list(V[1]))
        debug_chk('DEBUG_MATRIX',table,1)
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',"}}"+"}",1)
        debug_chk('DEBUG_MATRIX',"small_matrix view{{"+"{",1)
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX','kits: '+str(self.get_axis('kits',keysOnly=True)),1)
        debug_chk('DEBUG_MATRIX','variants: '+str(self.get_axis('variants',keysOnly=True)),1)
        debug_chk('DEBUG_MATRIX',"",1)
        if config['MATRIX_COLORS']:
            debug_chk('DEBUG_MATRIX',str(self.NP).replace("-1"," -").replace("0",'%s0%s'%(RED,WHITE)),1)
        else:
            debug_chk('DEBUG_MATRIX',str(self.NP).replace("-1"," -"),1)
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',"}}"+"}",1)
        debug_chk('DEBUG_MATRIX',"",1)
        
    def stdout_matrix_relations_data(self,dataStr='',run=0):
        if dataStr != '' and run != 0:
            print("---")
            print(str(dataStr)+"{{"+"{") #beg vim marker
        print("")
        #print counts
        print("relation counts:")
        print("-------------------------")
        print("mix:"+str(len(self.MIXA)))
        print("pos:"+str(len(self.POSA)))
        print("neg:"+str(len(self.NEGA)))
        print("unk:"+str(len(self.UNKA)))
        print("")
        #print data
        print("relation data:")
        print("-------------------------")
        for K in sorted(list(set([itm1[0] for itm1 in self.MIXA]+[itm1[0] for itm1 in self.POSA]+[itm1[0] for itm1 in self.NEGA]))):
            M = ",".join(sorted([itm2[1] for itm2 in self.MIXA if itm2[0] == K]))
            P = ",".join(sorted([itm2[1] for itm2 in self.POSA if itm2[0] == K]))
            N = ",".join(sorted([itm2[1] for itm2 in self.NEGA if itm2[0] == K]))
            U = ",".join(sorted([itm2[1] for itm2 in self.UNKA if itm2[0] == K]))
            print (str(K.lower())+"| mix:["+str(M.lower())+"], pos:["+str(P.lower())+"], neg:["+str(N.lower())+"], unk:["+str(U.lower())+"]")
        print("")
        if dataStr != '' and run != 0:
            print("}}"+"}") #end vim marker

    def set_new_order(self,val,cnt,kitType=False,variantType=False):
        if kitType:
            self.KITS[val][1] = cnt
        if variantType:
            self.VARIANTS[val][1] = cnt
        
    def set_new_axis(self,vals,cnts,kitType=False,variantType=False):
        self.VARIANTS = {}
        for x in range(len(vals)):
            self.VARIANTS[vals[x]] = [0,cnts[x]]

    def get_coord(self,kit_order,variant_order,moreInfo=False):
        buf = ""
        if moreInfo:
            buf = "coord: "+str(kit_order)+","+str(variant_order)
        kit = self.get_kit_name_by_order(kit_order)
        variant = self.get_variant_name_by_order(variant_order)
        buf = buf +  "k|v: "+str(kit)+"|"+(variant)
        if moreInfo:
            buf = buf + "value:"+str(self.NP[variant_order,kit_order])
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
                    #print(self.VARIANTS)
                    return [(key, self.VARIANTS[key]) for key in listByCount]
                if orderByType in ['kp','kn','kx']:
                    return [(key, self.KITS[key]) for key in listByCount]
        
    def get_kit_name_by_order(self,kit_order,listFlg=False):
        #listFlg: force listFlg as return data type
        intFlg = True
        try:
            value = int(kit_order)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            kit = None
            for itm in list(self.KITS.items()):
                if itm[1][1] == kit_order:
                    #print (".1..")
                    #print (itm[0])
                    #print (self.imperfect_known_variants)
                    #print (".2..")
                    if listFlg:
                            kit = itm[0]
                    else:
                            kit = itm[0]
                    break
            if listFlg:
                if kit == None:
                    return []
                else:
                    return [kit]
            else:
                return kit
        else: #assume it's a list/set/numpy array (whatever) > that I can cast to a list if need be
            kitList = []
            for ko in list(kit_order):
                for itm in list(self.KITS.items()):
                    #print (".1..")
                    #print (itm[0])
                    #print (self.imperfect_known_variants)
                    #print (".2..")
                    if itm[1][1] == ko:
                        kitList.append(itm[0])
                        break
            return(kitList)
        
    def get_coord_name_by_order(self,coord_order):
        coord_name = []
        for C in coord_order:
            variant_name = self.get_variant_name_by_order(variant_order=C[0])
            kit_name = self.get_kit_name_by_order(kit_order=C[1])
            coord_name.append((variant_name,kit_name))
        return coord_name
        
    def get_variant_name_by_order(self,variant_order,listFlg=False,impKnownFlg=False): #get variant name (can also take a list)
        #listFlg: force listFlg as return data type
        #impKnownFlg: force it to be a name that's part of the imperfect known list of variants
        intFlg = True
        try:
            value = int(variant_order)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            variant = None
            #hack to have an artifial top
            #if variant_order == -999: #top
            #    return 'top'
            #normal variants in the matrix
            for itm in list(self.VARIANTS.items()):
                if itm[1][1] == variant_order:
                    if listFlg:
                        if impKnownFlg is False or self.imperfect_known_variants is None:
                            variant = itm[0]
                        elif itm[0] in self.imperfect_known_variants:
                            variant = itm[0]
                    else:
                        if impKnownFlg is False or self.imperfect_known_variants is None:
                            variant = itm[0]
                        elif itm[1][1] in self.imperfect_known_variants:
                            variant = itm[0]
                    break
            if listFlg:
                if variant == None:
                    return []
                else:
                    return [variant]
            else:
                return variant
        else: #assume it's a list/set/numpy array (whatever) > that I can cast to a list if need be
            variantList = []
            for vo in list(variant_order):
                #hack to have an artificial top
                #if vo == -999: #top
                #    variantList.append('top')
                #normal variants in the matrix
                if 1 == 2:
                    foo = 1 #TODO: can fix this later
                else:
                    for itm in list(self.VARIANTS.items()):
                        if itm[1][1] == vo:
                            if impKnownFlg is False or self.imperfect_known_variants is None:
                                variantList.append(itm[0])
                            elif itm[1][1] in self.imperfect_known_variants:
                                variantList.append(itm[0])
                            break
            return(variantList)
        
        
    def get_variant_order_by_name(self,variant_name): #get variant order from its name
        return self.VARIANTS[variant_name][1]
        
    def get_kit_order_by_name(self,kit_name): #get kit order from its name
        return self.KITS[kit_name][1]
        
    def get_matrix_row_data(self,variant_order=None,variant_name=None): #get same type variant data for each kit
        if variant_name is not None:
            variant_order = self.get_variant_order_by_name(variant_name)
        if variant_order is not None:
            return self.NP[variant_order,]
        
    def get_matrix_col_data(self,kit_order=None,kit_name=None): #get all type variant data for one kit
        if kit_name is not None:
            kit_order = self.get_kit_order_by_name(kit_name)
        if kit_order is not None:
            return self.NP[:,kit_order].T
        
    def get_matrix_row_indices_by_val(self,val,variant_order=None,variant_name=None,overrideData=None): #like get_matrix_row_data but retrieves index info for given val
        if variant_name is not None:
            variant_order = self.get_variant_order_by_name(variant_name)
        #if variant_order == -999: #top hack
        #    return np.asarray((list(range(len(self.KITS))))) #TODO: working testing!!! (zak)
        if variant_order is not None and overrideData is not None: # we're sending in a custom evaluation
            return np.argwhere(overrideData[0,] == val).T[1,] #with override data, there's only one line evaluated - 1d datset
        if variant_order is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[variant_order,] == val).T[1,] #default data, it's the entire matrix - 2d dataset 
        
    def get_matrix_col_indices_by_val(self,val,kit_order=None,kit_name=None,overrideData=None): #like get_matrix_col_data but retrieves index info for given val
        if kit_name is not None:
            kit_order = self.get_kit_order_by_name(kit_name)
        if kit_order is not None and overrideData is not None:
            return np.argwhere(overrideData[:,0] == val).T[0,] #with override data, there's only one line evaluated - 1d dataset
        if kit_order is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[:,kit_order] == val).T[0,] #default data, it's the entire matrix - 2d dataset

    def use_imperfect_known_variants_only(self,variants):
        if self.imperfect_known_variants == None:
            return variants
        newList = []
        for v in variants:
            if v in self.imperfect_known_variants:
                newList.append(v)
        return newList

    def get_perfect_variants(self):
        neg_idx = np.argwhere(self.NP==-1) #get index to negative data
        neg_idx_r = np.unique(neg_idx[:,0]) #get unique rows of those indices
        print(neg_idx_r)
        print(len(self.NP))
        sys.exit()
        #Note: https://stackoverflow.com/questions/25330959/how-to-select-inverse-of-indexes-of-a-numpy-array
        #(beg) technique to delete things other than the idx values
        mask = np.ones(len(self.NP), np.bool)
        mask[neg_idx_r] = 0
        variant_data = self.NP[mask]
        #(end) 
        print(variant_data)
        sys.exit()
        variant_names = self.get_variant_name_by_order(variant_idx)
        self.set_new_axis(variant_names,variant_idx_r,variantType=True) #reset variant axis
        self.NP = variant_data #reset data
        
    def get_imperfect_variants(self):
        variant_idx = np.argwhere(self.NP==-1) #get index to negative data
        variant_idx_r = np.unique(variant_idx[:,0]) #get unique rows of those indices
        variant_data = self.NP[variant_idx_r]
        #print(variant_data)
        #sys.exit()
        variant_names = self.get_variant_name_by_order(variant_idx_r)
        self.set_new_axis(variant_names,variant_idx_r,variantType=True) #reset variant axis
        #sys.exit()
        #self.NP = None
        #print(variant_data)
        self.NP = variant_data #reset data
        #print(self.NP)
        #sys.exit()

    def get_perfect_variants_idx(self):
        allrows_idx = list(range(len(self.VARIANTS)))
        neg_idx_r = list(self.get_imperfect_variants_idx())
        return np.unique(np.delete(allrows_idx, neg_idx_r, axis=0))
        
    def get_imperfect_variants_idx(self):
        neg_idx = np.argwhere(self.NP==-1)
        return np.unique(neg_idx[:,0])
        
    def get_imperfect_known_variants_idx(self):
        imp_idx = list(self.get_imperfect_variants_idx())
        unk_idx = list(np.unique(np.argwhere(self.NP==0)[:,0]))
        for x in unk_idx:
            if x in imp_idx:
                imp_idx.remove(x)
        return imp_idx
        
    def get_imperfect_unknown_variants_idx(self):
        unk_idx = list(np.unique(np.argwhere(self.NP==0)[:,0]))
        prf_idx = list(self.get_perfect_variants_idx())
        for x in prf_idx:
            if x in unk_idx:
                unk_idx.remove(x)
        return unk_idx

    def get_subset_variants(self,override_val=None,variant_order=None,variant_name=None,kit_order=None,kit_name=None, convertToNames=True,impKnownFlg=False):
        #variant_order: is variant's order in matrix, name is variant name
        #override_val: is the override val (ie: check what conditions are after setting a coord to be 1 and not 0, for example)

        def get_subsets(pc,vo):
            VAR1p = np.argwhere(self.NP[:,pc]==1)[:,0] #looking for variants w/pos assignments like the incoming variant's pos conditions
            idxP = np.argwhere(VAR1p==vo) #idx make sure we exclude the incoming variant
            VAR2p = np.delete(VAR1p, idxP)
            unique_elements_p, counts_elements_p = np.unique(VAR2p, return_counts=True)
            VAR3p = np.asarray((unique_elements_p, counts_elements_p)).T
            #...
            if 1 == 2: # this code allows for unks to be considered
                VAR1u = np.argwhere(self.NP[:,pc]==0)[:,0] #looking for variants w/unk assignments like the incoming variant's pos conditions
                idxU = np.argwhere(VAR1u==vo) #idx make sure we exclude the incoming variant
                VAR2u = np.delete(VAR1u, idxU)
                unique_elements_u, counts_elements_u = np.unique(VAR2u, return_counts=True)
                VAR3u = np.asarray((unique_elements_u, counts_elements_u)).T
                VAR3x = np.concatenate((VAR3p,VAR3u), axis=0)
            #...
            else: #this is just positives (I think more what we're looking for)
                VAR3x = VAR3p
            #...
            #Note: for the following "adding technique" -- we need to exclude unk situations for the comparison (special handling)
            #beg - adding technique - got this idea here: https://stackoverflow.com/questions/30041286/sum-rows-where-value-equal-in-column
            unq, unq_inv = np.unique(VAR3x[:,0], return_inverse=True)
            out = np.zeros((len(unq), VAR3x.shape[1]), dtype=VAR3x.dtype) #create empty array to put the added values
            out[:, 0] = unq #fill the first column
            np.add.at(out[:, 1:], unq_inv, VAR3x[:, 1:])
            #end - adding technique
            #Note: sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
            out1 = out[out[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
            VAR4 = np.argwhere(out1[:,1]<len(pc)) #[:,0]
            #end - sorting technique
            out2a = out1[:,0]
            out2b = out1[:,1]
            VAR5a = out2a[list(VAR4.T[0])] #these are the superset variant orders ids (in order, max first)
            VAR5b = out2b[list(VAR4.T[0])]#these are the superset variant order ids (in order, max first)
            VAR6 = np.asarray((VAR5a,VAR5b)).T #merged for return
            return VAR6
        
        if variant_name is not None:
            variant_order = self.get_variant_order_by_name(variant_name)
        if kit_name is not None:
            kit_order = self.get_kit_order_by_name(kit_name)
        if override_val is not None and kit_order is not None:
            overrideData = self.get_row_when_value_override_coord(override_val,kit_order=kit_order,variant_order=variant_order)
            pc = self.get_matrix_row_indices_by_val(1,variant_order=variant_order,overrideData=overrideData) #pos conditions when override coord with a value
        else:
            pc = self.get_matrix_row_indices_by_val(1,variant_order=variant_order) #default pos conditions
        if convertToNames is True:
            subs = self.get_variant_name_by_order(variant_order=get_subsets(pc,variant_order)[:,0],impKnownFlg=impKnownFlg,listFlg=1)
        else:
            subs = get_subsets(pc,variant_order)[:,0]
        if override_val is not None and kit_order is not None:
            if override_val==1:
                #print("- subsetsP: "+",".join(subs)+" pc:"+str(pc))
                if config['DBG_SUBS_SUPS'] == True:
                    print("- subsetsP: "+",".join([str(i) for i in subs]) +" pc:"+str(pc))
            else:
                if config['DBG_SUBS_SUPS'] == True:
                    print("- subsetsN: "+",".join([str(i) for i in subs]) +" pc:"+str(pc))
        else:
            if config['DBG_SUBS_SUPS'] == True:
                print("- subsets: "+",".join([str(i) for i in subs]) +" pc:"+str(pc))
        return subs
        
    def get_supset_variants(self,override_val=None,variant_order=None,variant_name=None,kit_order=None,kit_name=None, convertToNames=True,impKnownFlg=False):
        #variant order: is variant's order in matrix, name is variant name
        #override_val: is the override val (ie: check what conditions are after setting a coord to be 1 and not 0, for example)

        def get_supsets(pc,vo):
            VAR1p = np.argwhere(self.NP[:,pc]==1)[:,0] #looking for variants w/pos assignments like the incoming variant condition
            unqP, cntP = np.unique(VAR1p, return_counts=True)
            VAR2p = np.asarray((unqP, cntP)).T
            #...
            if 1 == 2: # this code allows for unks to be considered (not sure this is good code)
                VAR1u = np.argwhere(self.NP[:,pc]==0)[:,0] #looking for variants w/pos assignments like the incoming variant condition
                unqU, cntU = np.unique(VAR1u, return_counts=True)
                VAR2u = np.asarray((unqU, cntU)).T
                VAR2x = np.concatenate((VAR2p,VAR2u), axis=0)
            else:
                VAR2x = VAR2p
                VAR2y = VAR2x[VAR2x[:,1]==len(pc)] #has to have at least what the incoming variant had in count
                idxU = np.argwhere(VAR2y[:,0]==vo) #idx make sure we exclude the incoming variant
                VAR3 = np.delete(VAR2y[:,0], idxU) #idx again/delete
            #...
            if len(VAR3) == 0: return [] #there are no supsets according to the filter
            #(beg) master list of all positives
            allPos = np.argwhere(self.NP==1)[:,0]

            unqA, cntA = np.unique(allPos, return_counts=True)
            AP = np.asarray((unqA, cntA))[1,]
            #(end) master list of all positives
            VAR5 = AP[list(VAR3),] #extrapolate the right subset mix to master list of all positives
            VAR6 = np.asarray((VAR3,VAR5)).T
            #Note: sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
            VAR7 = VAR6[VAR6[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
            #VAR7 = VAR6[VAR6[:,1].argsort()] #normal sort (2nd col) - so we can get the minimum ones first
            #(end) sorting
            return VAR7

        if variant_name is not None:
            variant_order = self.get_variant_order_by_name(variant_name)
        if kit_name is not None:
            kit_order = self.get_kit_order_by_name(kit_name)
        if override_val is not None and kit_order is not None:
            overrideData = self.get_row_when_value_override_coord(override_val,kit_order=kit_order,variant_order=variant_order)
            pc = self.get_matrix_row_indices_by_val(1,variant_order=variant_order,overrideData=overrideData) #pos conditions when override coord with a value
        else:
            pc = self.get_matrix_row_indices_by_val(1,variant_order=variant_order) #default pos conditions
        supsX = get_supsets(pc,variant_order)
        if len(supsX) == 0:
            sups = supsX
        elif convertToNames is True:
            sups = self.get_variant_name_by_order(variant_order=supsX[:,0],impKnownFlg=impKnownFlg,listFlg=1)
        else:
            sups = supsX[:,0]
                
        if override_val is not None and kit_order is not None:
            if override_val==1:
                if config['DBG_SUBS_SUPS'] == True:
                    print("- supsetsP: "+",".join([str(i) for i in sups]) +" pc:"+str(pc))
            else:
                if config['DBG_SUBS_SUPS'] == True:
                    print("- supsetsN: "+",".join([str(i) for i in sups]) +" pc:"+str(pc))
        else:
            if config['DBG_SUBS_SUPS'] == True:
                print("- supsets: "+",".join([str(i) for i in sups]) +" pc:"+str(pc))
        return sups
        
    def get_min_superset_variant(self,variant_order=None,variant_name=None): #order is variant's order in matrix, name is variant name
        if variant_name is not None:
            variant_order = self.get_variant_order_by_name(variant_name)
        sups = self.get_supset_variants(variant_order=variant_order)
        if len(sups) is 0:
            return None
        else:
            return sups[0]

    def get_coord_value(self,kit_order=None,variant_order=None,kit_name=None,variant_name=None):
        if kit_order is not None and variant_order is not None:
            return self.NP[variant_order][kit_order]
        if kit_name is not None and variant_name is not None:
            get_kit_order_by_name(kit_name)
            get_variant_order_by_name(variant_name)
            return self.NP[variant_order][kit_order]
        
    def get_row_when_value_override_coord(self,override_val,kit_order=None,variant_order=None,kit_name=None,variant_name=None):
        #override_val -- is the override val (ie: check what conditions are after setting a coord to be 1 and not 0, for example)
        row = self.get_matrix_row_data(variant_order=variant_order)
        #print("*****beg:This is the val we're overriding -kit/variant")
        #print(kit_order)
        #print(variant_order)
        #print("*****end:This is the val we're overriding")
        #(beg) this technique found here: https://stackoverflow.com/questions/6431973/how-to-copy-data-from-a-numpy-array-to-another
        #note: this is necessary because otherwise, it seems to be working from
        #      the same memory pointer when I push the override test value in
        rowO = np.empty_like(row)
        rowO[:] = row
        #(end)
        #print(rowO)
        rowO[0,kit_order] = override_val
        #print(rowO)
        #sys.exit()
        return rowO

    def get_matrix_data(self):

        #sql - exclude perfect variants
        sql0 = '''
            SELECT C.kit_id, V.name, C.assigned, V.variant_id
            FROM s_calls C, s_variants V,
            (SELECT DISTINCT V.variant_loc
             FROM s_calls C, s_variants V
             WHERE (C.assigned = -1 OR V.name = 'top') AND
             V.variant_loc = C.variant_loc
            ) VX
            WHERE C.variant_loc = V.variant_loc AND
            C.variant_loc = VX.variant_loc
            ORDER by 4;
            '''

        #get data
        self.dbo.sql_exec(sql0)
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
                FROM s_calls C, s_variants V ,
                (SELECT distinct C.variant_loc
                FROM s_calls C, s_variants V
                WHERE (C.assigned = -1 OR V.name = 'top') AND
                V.variant_loc = C.variant_loc) VX
                WHERE C.variant_loc = V.variant_loc AND
                C.variant_loc = VX.variant_loc
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

        #chk matrix (debugging)
        #self.stdout_tbl_matrix()
        #sys.exit()
        
    def get_matrix_count_data(self):

        #NOTE: this could be done with numpy - which is better? does it matter?

        #vars
        self.CNTS = {}
        sqlc = {}

        sql = '''
            FROM s_calls C, s_variants V,
            (SELECT DISTINCT C.variant_loc
            FROM s_calls C, s_variants V
            WHERE (C.assigned = -1 OR V.name = 'top') AND
            V.variant_loc = C.variant_loc
            )VX
            WHERE C.variant_loc = VX.variant_loc AND
            C.variant_loc = V.variant_loc AND
            '''

        #sql - cnt variants
        sqlc['vp'] = "SELECT count(V.name), V.name %s C.assigned = 1 GROUP BY 2;" % sql
        sqlc['vn'] = "SELECT count(V.name), V.name %s C.assigned = -1 GROUP BY 2;" % sql
        sqlc['vx'] = "SELECT count(V.name), V.name %s C.assigned = 0 GROUP BY 2;" % sql

        #sql - cnt kits
        sqlc['kp'] = "SELECT count(C.kit_id), C.kit_id %s C.assigned = 1 GROUP BY 2;" % sql
        sqlc['kn'] = "SELECT count(C.kit_id), C.kit_id %s C.assigned = -1 GROUP BY 2;" % sql
        sqlc['kx'] = "SELECT count(C.kit_id), C.kit_id %s C.assigned = 0 GROUP BY 2;" % sql

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
        #sql - get unknowns (without kits) {{{

        sql = '''
            SELECT distinct V1.name, V2.name
            FROM s_calls C1, s_calls C2,
            s_variants V1, s_variants V2
            WHERE
            C1.kit_id = C2.kit_id AND
            C1.assigned = 1 AND
            C2.assigned = 0 AND
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
                (QC2.assigned = -1 OR QC2.assigned == 1) AND
                QC1.variant_loc = QV1.variant_loc AND
                QC2.variant_loc = QV2.variant_loc)
            ORDER by 1,2;
            '''
        self.dbo.sql_exec(sql)
        self.UNKA = self.dbo.fetchall()

        #}}}

    def test_rule1_subsets(self,unk_variant,subsets,supsets):
        kit_order = unk_variant[1]
        pc = self.get_matrix_col_indices_by_val(1,kit_order=kit_order)
        if config['DEBUG_RULE1']:
            print("!!!pc: positive conditions")
            print(self.get_variant_name_by_order(variant_order=pc))

        #standard subset test
        for sub in subsets:
            if sub in pc:
                if config['DEBUG_RULE1']:
                    print("!!!sub: ")
                    print(self.get_variant_name_by_order(variant_order=sub))
                self.unk_variants.remove(unk_variant)
                self.NP[unk_variant[0],unk_variant[1]]=1
                self.resolved_variants.append((unk_variant,True))
                if config['DEBUG_RULE1']:
                    print("")
                return 'RULE1: True: sub of %s' % self.get_variant_name_by_order(variant_order=sub)
        #Note: if get here -- go to rule 2
        return self.test_rule2_supsets(unk_variant=unk_variant,subsets=subsets,supsets=supsets)
        
    def test_rule2_supsets(self,unk_variant,subsets,supsets):
        #Note: I think the rule goes like this:  if there is an
        #existing perfect variant (now that I understand the definition) has the
        #same (or superset of the) known +/-'s to then imperfect variant when
        #the unresolved ? is turned to a +, then it's a +.

        variant_order = unk_variant[0]
        kit_order = unk_variant[1]
        overrideData = self.get_row_when_value_override_coord(1,kit_order=kit_order,variant_order=variant_order)

        kpc = self.get_matrix_row_indices_by_val(1,variant_order=variant_order,overrideData=overrideData)
        knc = self.get_matrix_row_indices_by_val(-1,variant_order=variant_order)

        if config['DEBUG_RULE2']:
            print("!!!kpc: kit positive conditions")
            print(kpc)
            print(self.get_kit_name_by_order(kit_order=kpc))
            print("!!!knc: kit negative conditions")
            print(knc)
            print(self.get_kit_name_by_order(kit_order=knc))

        for sup in supsets:
            if config['DEBUG_RULE2']:
                print("!!!sup: ")
                print(self.get_variant_name_by_order(variant_order=sup))
            kpc4sup = self.get_matrix_row_indices_by_val(1,variant_order=sup)
            if config['DEBUG_RULE2']:
                print("!!!kpc4sup: ")
                print(self.get_kit_name_by_order(kit_order=kpc4sup))
            knc4sup = self.get_matrix_row_indices_by_val(-1,variant_order=sup)
            if config['DEBUG_RULE2']:
                print("!!!knc4sup: ")
                print(self.get_kit_name_by_order(kit_order=knc4sup))
            lenPos = len(list(set(kpc).intersection(set(kpc4sup))))
            lenNeg = len(list(set(knc).intersection(set(knc4sup))))
            if config['DEBUG_RULE2']:
                print("lenPos: %s" % lenPos)
                print("lenNeg: %s" % lenNeg)
            if lenPos == len(kpc) and lenNeg == len(knc):
                self.unk_variants.remove(unk_variant)
                self.NP[unk_variant[0],unk_variant[1]]=1
                self.resolved_variants.append((unk_variant,True))
                if config['DEBUG_RULE2']:
                    print("")
                return 'RULE2: True - ambiguous: equivalent/subset of %s'%(self.get_variant_name_by_order(variant_order=sup))

        #if we get here ... go to rule 3
        return self.test_rule3_diff_branches(unk_variant,subsets,supsets)
        
    def test_rule3_diff_branches(self,unk_variant,subsets,supsets):
        variant_order = unk_variant[0]
        kit_order = unk_variant[1]
        #supsets = self.get_supset_variants(variant_order=unk_variant[0],convertToNames=False) #superset of given coord
        #supsets.append(-999) #top
        #subsets = self.use_imperfect_known_variants_only(self.get_subset_variants(variant_order=unk_variant[0],convertToNames=False)) #superset of given coord
        #vi = self.use_imperfect_known_variants_only(self.get_matrix_col_indices_by_val(1,kit_order=kit_order))
        vi = self.get_matrix_col_indices_by_val(1,kit_order=kit_order)
        #vi.append(-999) #top
        ki = self.get_matrix_row_indices_by_val(1,variant_order=variant_order).tolist()
        if ki is None:
            ki = []
        rule_p1_list = []

        if config['DEBUG_RULE3']:
            print("[1]!!! variant_order: %s" %self.get_variant_name_by_order(variant_order))
            print("[2]!!! supsets: %s" %self.get_variant_name_by_order(supsets))
            print("[3]!!! vi (imperfect knowns when set coord to pos): %s" %self.get_variant_name_by_order(vi))
            print("[4]!!! ki (related kits when set coord to pos): %s" %self.get_kit_name_by_order(ki))
            print("")

        #Note: first deal with the variant that might might share the unknown variant we're wondering about
        #Q: Is it not a direct relation?
        #Q: does it have a superset?
        for V in vi:

            #sups4V = self.use_imperfect_known_variants_only(self.get_supset_variants(variant_order=V,convertToNames=False)) #superset of related coord
            sups4V = self.get_supset_variants(variant_order=V,convertToNames=False) #superset of related coord
            #sups4V.append(-999) #top

            if config['DEBUG_RULE3']:
                print("[P1.5]!!! V: %s" %self.get_variant_name_by_order(V))
                print("[P1.6]!!! sups4V(%s): %s" %(self.get_variant_name_by_order(V),self.get_variant_name_by_order(sups4V)))

            #Is it a direct relation? 
            if V not in supsets and V not in subsets: # and variant_order not in sups4V:
                if config['DEBUG_RULE3']:
                    print("")
                directRelation_chk = False
                if config['DEBUG_RULE3']:
                    print("directRelation Chk is: %s (continue with V=%s)" %(directRelation_chk,self.get_variant_name_by_order(V)))
                    print("")
            else:
                directRelation_chk = True
                if config['DEBUG_RULE3']:
                    print("directRelation Chk is: %s (don't continue with V=%s)" %(directRelation_chk,self.get_variant_name_by_order(V)))
                    print("")
                continue # try the next V

            if directRelation_chk == False:

                if config['DEBUG_RULE3']:
                    print("[P1.7]!!! V in: %s" %self.get_variant_name_by_order(V))
                k4V = self.get_matrix_row_indices_by_val(1,variant_order=V).tolist() #other kits per V
                k4V.remove(kit_order)
                if config['DEBUG_RULE3']:
                    print("[P1.8]!!! k4V(%s): %s" %(self.get_variant_name_by_order(V),self.get_kit_name_by_order(k4V)))

                #Does it have any supersets?
                if len(sups4V) == 0:
                    #sups4V.append(-999) #top
                    if config['DEBUG_RULE3']:
                        print("sups Chk: sups not in sups4V (don't continue with V=%s)"%self.get_variant_name_by_order(V))
                        print("")
                    continue # try the next V

                for sup4V in sups4V:
                    if config['DEBUG_RULE3']:
                        print("")
                        print("sups Chk: sups in sups4V (continue with V=%s)"%self.get_variant_name_by_order(V))
                        print("")
                        print("[P1.9]!!! sup4V(%s): %s" %(self.get_variant_name_by_order(V),self.get_variant_name_by_order(sup4V)))
                    k4sup4V = self.get_matrix_row_indices_by_val(1,variant_order=sup4V).tolist() #other kits per V
                    if config['DEBUG_RULE3']:
                        print("[P1.10]!!! k4sup4V(%s)(%s)(bef.rem) in: %s" %(self.get_variant_name_by_order(V),self.get_variant_name_by_order(sup4V),self.get_kit_name_by_order(k4sup4V)))
                    if kit_order in k4sup4V: #(1)can't be the given coord's kit
                        k4sup4V.remove(kit_order) #(2)can't be the given coord's kit
                    if config['DEBUG_RULE3']:
                        print("[P1.11]!!! k4sup4V(%s)(%s)(aft.rem) in: %s" %(self.get_variant_name_by_order(V),self.get_variant_name_by_order(sup4V),self.get_kit_name_by_order(k4sup4V)))
                    k_in_k4V_and_k4sup4V = list(set(k4V).intersection(set(k4sup4V)))

                    #Is there additional overlap btw this other variant and its superset?
                    if len(k_in_k4V_and_k4sup4V) == 0:
                        if config['DEBUG_RULE3']:
                            print("P1.k4V + k4sup4V intersection chk: no additional overlap (don't continue with V=%s)"%self.get_variant_name_by_order(V))
                            print("")
                        continue

                    #Everything seems ok for part one of this rule
                    if config['DEBUG_RULE3']:
                        print("k4V + k4sup4V intersection chk: additional overlap (continue with V=%s)"%self.get_variant_name_by_order(V))
                        print("")
                        print("[P1.12]!!! k_in_k4V_and_k4sup4V: %s" %self.get_kit_name_by_order(k_in_k4V_and_k4sup4V))
                        print("")
                    rule_p1_list.append((sup4V,k_in_k4V_and_k4sup4V))

        if len(rule_p1_list):
            rule_p2_list = []

            #Loop the known supersets of the given variant
            for sup4vo in supsets:

                #What positive variant relations do those supersets have? 
                k4sup4vo = self.get_matrix_row_indices_by_val(1,variant_order=sup4vo).tolist() #other kits per V
                #k4sup4vo.append(-999) #top
                
                if config['DEBUG_RULE3']:
                    print("[P2.13]!!! k4sup4vo(%s): %s"%(self.get_variant_name_by_order(sup4vo),self.get_kit_name_by_order(k4sup4vo)))
                    print("[P2.14]!!! ki: %s"%self.get_kit_name_by_order(ki))
                    print("")

                #Is there any overlap with the given coord?
                k_in_k4sup4vo_and_ki = set(k4sup4vo).intersection(set(ki))
                if len(k_in_k4sup4vo_and_ki) == 0:
                    print("k4sup4vo + ki intersection chk: no additional overlap (don't continue with sup4vo=%s)"%self.get_variant_name_by_order(sup4vo))
                    print("")
                    continue

                #If so, everything seems ok for part two of this rule
                if config['DEBUG_RULE3']:
                    print("k4sup4vo + ki intersection chk: additional overlap (continue with sup4vo=%s)"%self.get_variant_name_by_order(sup4vo))
                    print("")
                    print("[P2.15]!!! k_in_k4sup4vo_and_ki in: %s" %self.get_kit_name_by_order(k_in_k4sup4vo_and_ki))
                rule_p2_list.append((sup4vo,k_in_k4sup4vo_and_ki))

            #Check that both rules are satisfied and not shared btw both lists
            if len(rule_p2_list):
                for itm1 in rule_p1_list:
                    if itm1 not in rule_p2_list:
                        for itm2 in rule_p2_list:
                            if itm2 not in rule_p1_list:
                                #If here, we have a winner!
                                msg = "%s:%s" % (itm1,itm2)
                                #self.resolved_variants.append((unk_variant,False))
                                self.unk_variants.remove(unk_variant)
                                self.NP[unk_variant[0],unk_variant[1]]=-1
                                self.resolved_variants.append((unk_variant,False))
                                #if config['DEBUG_RULE3']:
                                if config['DEBUG_RULE3']:
                                    print("")
                                #print("}}"+"}")
                                msg1a = self.get_variant_name_by_order(itm1[0])
                                msg1b = self.get_variant_name_by_order(itm1[1])
                                msg2a = self.get_variant_name_by_order(itm2[0])
                                msg2b = self.get_variant_name_by_order(itm2[1])
                                return "RULE3: False - 2 diff branches (%s,%s) + (%s,%s)"% (msg1a,msg1b,msg2a,msg2b)

        #print("}}"+"}")
        if config['DEBUG_RULE3']:
            print("")
        return "RULEX: Unk"

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


