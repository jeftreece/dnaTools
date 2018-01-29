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

class Sort(object):
    
    def __init__(self):
        self.TREE = {}
        self.REF = None
        self.DBG = 1 #>1 = means show don't stash msgs
        self.MODE = 1 #(sort tree presentation) 1=letters, 2=names, 3=letters+names 
        self.KITS = None
        self.VARIANTS = None
        self.DATA = None
        self.KDATA = None
        self.VDATA = None
        self.CNTS = {}
        self.NP = None
        self.db = None #db init
        self.dc = None #db cursor

    #TOGGLE BTW SORT TOOLS

    def sort_data(self):
        #self.sort_data_matrix()
        self.sort_data_tree()

    #SORT - MATRIX FORMAT

    def sort_data_matrix(self):

        #get counts
        self.sort_cnts()

        #sql
        sql = "select C.kit_id,V.name,C.assigned,V.variant_id from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 4"
        sql1 = "select distinct V.name,V.variant_id from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 2"

        #all unique variants
        self.dc.execute(sql1)
        F = self.dc.fetchall()
        self.VARIANTS = {}
        cnt=0
        for itm in F:
            self.VARIANTS[itm[0]]=[cnt,cnt] #default sort,custom sort 
            cnt = cnt + 1

        #get data
        self.dc.execute(sql)
        F = self.dc.fetchall()
        
        #all unique kits
        self.KITS = {}
        cnt=0
        for k in sorted(list(set([itm[0] for itm in F]))):
            self.KITS[k]=[cnt,cnt] #default sort,custom sort 
            cnt = cnt + 1

        #retrieve data from sqlite like so: [V][K] [x,x,x,x,x,x,...]
        DATA = OrderedDict()
        for row in F:
            if row[1] not in DATA:
                DATA[row[1]] = []
            DATA[row[1]].append(row[2])

        #numpy version of data
        for key,value in DATA.items():
            self.NP = np.matrix(list(DATA.values()))

        print("")
        print("data - default")

        #1st tbl out
        self.stdout_tbl_matrix()

        #2nd tbl prep 
        self.sort_step1()

        print("data - step 1")

        #2nd tbl out
        self.stdout_tbl_matrix()

        print("data - step 2")

        #2nd tbl prep
        self.sort_step2()

        #2nd tbl out
        self.stdout_tbl_matrix()

        sys.exit()
        
    def sort_cnts(self):
        #vars
        self.MODE = 2
        self.CNTS = {}
        sqlc = {}
        #sql - cnt variants
        sqlc['vp'] = "select count(V.name),V.name from s_calls C,s_variants V where C.variant_loc=V.variant_loc and C.assigned=1 group by 2"
        sqlc['vn'] = "select count(V.name),V.name from s_calls C,s_variants V where C.variant_loc=V.variant_loc and C.assigned=0 group by 2"
        sqlc['vx'] = "select count(V.name),V.name from s_calls C,s_variants V where C.variant_loc=V.variant_loc and C.assigned is null group by 2"
        #sql - cnt kits
        sqlc['kp'] = "select count(C.kit_id),C.kit_id from s_calls C,s_variants V where C.variant_loc=V.variant_loc and C.assigned=1 group by 2"
        sqlc['kn'] = "select count(C.kit_id),C.kit_id from s_calls C,s_variants V where C.variant_loc=V.variant_loc and C.assigned=0 group by 2"
        sqlc['kx'] = "select count(C.kit_id),C.kit_id from s_calls C,s_variants V where C.variant_loc=V.variant_loc and C.assigned is null group by 2"
        #get all cnts
        for key, sql in sqlc.items():
            self.CNTS[key] = {}
            self.dc.execute(sql)
            F = self.dc.fetchall()
            for itm in F:
                self.CNTS[key][itm[1]] = itm[0]
        
    def sort_step1(self):

        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        for K,V in self.get_axis('variants'):
            if 0 not in self.get_numpy_matrix_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_numpy_matrix_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('vp'):
            if 0 in self.get_numpy_matrix_row_as_list(V[1]) and 'None' not in self.get_numpy_matrix_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_numpy_matrix_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('variants'):
            if 0 in self.get_numpy_matrix_row_as_list(V[1]) and 'None' in self.get_numpy_matrix_row_as_list(V[1]):
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
        
    def get_cur_kit_list(self):
        return self.get_axis('kits',keysOnly=True)
        
    def get_cur_variant_list(self):
        return self.get_axis('variants',keysOnly=True)
        
    def get_numpy_matrix_row_as_list(self,rownum,noneToStr=True):
        if noneToStr:
            return ['None' if v is None else v for v in self.NP[rownum,:].tolist()[0]]
        else:
            return self.NP[rownum,:].tolist()[0]
        
    def get_axis(self,orderByType=None,keysOnly=False):
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
        
    def stdout_tbl_matrix(self):
        print("")
        table = BeautifulTable()
        table.column_headers = ['top']+self.get_cur_kit_list()
        for K,V in self.get_axis('variants'):
            table.append_row([K]+self.get_numpy_matrix_row_as_list(V[1]))
        print(table)
        print("")
        

    #SORT - TREE FORMAT

    #TODO: set up a random approach to pushing data into the sort. troubleshoot results

    def sort_data_tree(self):

        #beg collapse vim marker
        print("PREP {"+"{{")

        self.MODE = 3

        #all kits, variant, assignment mixes 
        #self.commit()

        #Letters 
        if self.MODE == 1:
            sql = "select C.kit_id,C.variant_loc,C.assigned from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 1,2,3"
        #Names
        if self.MODE == 2:
            sql = "select C.kit_id,V.name,C.assigned from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 1,2,3"
        #Combo - Names+Letters
        if self.MODE == 3:
            sql = "select C.kit_id,'('||C.variant_loc||') '||V.name,C.assigned from s_calls C,s_variants V where C.variant_loc=V.variant_loc order by 1,2,3"

        self.dc.execute(sql)
        F = self.dc.fetchall()
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
        print("---")
        print("all unique variants")
        print(VARIANTS)
        #sys.exit()
        #print("---")
        #print("all unique variants - pos+neg")
        #print(VARIANTSa)
        #sys.exit()
        
        #kits with positive assignments
        Fp = sorted(list(set([i[0] for i in list(filter(lambda x: x[2]==1, F))])))
        print("---")
        print("kits with positive assignment variant calls")
        print(Fp)

        #kits with negative assignments
        Fn = sorted(list(set([i[0] for i in list(filter(lambda x: x[2]==0, F))])))
        print("---")
        print("kits with negative assignment variant calls")
        print(Fn)
        #sys.exit()

        #per all the kits with positive variants (build new dict)
        print("---")
        print("dict of kits with their positive assignment variant calls")
        KA={}
        for k in Fp:
            Kp = sorted(list(set(['+'+i[1] for i in list(filter(lambda x: x[0]==k and x[2]==1, F))])))
            #['A+', 'D+', 'F+', 'H+', 'M+']
            print(k+" "+str(Kp))
            #sys.exit()
            KA[k] = {'len':len(Kp),'plen':len(Kp),'sort':0,'variants':Kp}

        #per all the kits with negative variants (build new dict)
        print("---")
        print("dict of kits with their negative assignment variant calls")
        for k in Fn:
            Kn = sorted(list(set(['-'+i[1] for i in list(filter(lambda x: x[0]==k and x[2]==0, F))])))
            print(k+" "+str(Kn))
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
        print("---")
        print("combined dict of kits with pos+neg variant calls - sorted")
        #newV3 = {}
        for d in newV2:
            #newV3[d['kit']] = d['variants']
            STR = d['kit']+':'+str(d['variants'])
            print(STR.replace("'",""))

        #build variant relationship data that we need for sorting
        print("---")
        DATA = {}
        for VX in VARIANTS:
            DATA[VX] = {'mix':[],'pos':[],'neg':[]}
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
                            DATA[VX]['mix'].append(VY)
                            break
                    if chk1 is True and chk2 is False:
                        DATA[VX]['pos'].append(VY)
                    if chk2 is True and chk1 is False:
                        DATA[VX]['neg'].append(VY)
                    
        #print to stdout so I can see what I'm doing
        print("variant relationship data needed for sorting")
        for key, value in DATA.items():
            print(key+'|'+str(value))

        #end collapse vim marker
        print("}"+"}}")

        #build unsorted tree with all nodes under top
        self.TREE = {}
        self.TREE['top'] = Node("top")
        #STASH = {}
        for key, value in DATA.items():
            self.TREE[key] = Node(key, parent=self.TREE['top'])

        #sort it
        #STASH = self.sort_variant_tree(DATA,run_mix=True,run_pos=True,run_neg=True,run=1)
        STASH = self.sort_variant_tree(DATA,run_all=True,run=1)
        STASH = self.sort_variant_tree(STASH,run_all=True,run=2)

        sys.exit()

        #json/stdout a variable (debugging)
        self.stdout_dump_var(newV2)
        sys.exit()
        
    def sort_variant_tree (self,DATA,run_mix=False,run_pos=False,run_neg=False,run_all=None,run=1):

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
        # HIDE-ME: sample data{{{

        # variant,name,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
        # A,M343,1,1,Null,1,1,1,1,1,1,1
        # B,Z9,0,0,0,1,0,0,0,0,1,0
        # C,Z381,0,1,1,Null,0,0,1,1,1,0
        # D,U106,1,1,1,1,0,1,1,1,1,0
        # E,Z301,Null,Null,1,Null,Null,0,Null,1,Null,Null
        # F,Z18,1,0,0,0,0,1,0,0,0,0
        # G,Z156,0,1,0,0,0,0,1,0,0,0
        # H,L11,1,1,1,1,1,Null,1,1,1,1
        # I,Z28,0,0,0,1,0,Null,0,0,1,0
        # J,P312,0,0,0,0,1,0,0,0,0,1
        # K,Z8,0,0,0,1,0,0,0,0,0,0
        # L,A297,0,Null,0,1,0,0,1,0,Null,0
        # M,M269,1,1,1,1,1,1,1,1,1,1
        # N,Z306,0,1,0,0,0,0,0,0,0,0
        # O,L48,0,0,1,1,0,0,0,1,1,0

        #}}}
        # HIDE-ME: processed raw results: {{{

        # A|{mix: [B,C,D,E,F,G,I,J,K,L,N,O], pos: [H,M], neg: []}
        # B|{mix: [K], pos: [A,C,D,H,I,L,M,O], neg: [F,G,J,N]}
        # C|{mix: [B,G,I,L,N,O], pos: [A,D,E,H,M], neg: [F,J,K]}
        # D|{mix: [B,C,E,F,G,I,K,L,N,O], pos: [A,H,M], neg: [J]}
        # E|{mix: [], pos: [A,C,D,H,M,O], neg: [B,F,G,I,J,K,L,N]}
        # F|{mix: [], pos: [A,D,H,M], neg: [B,C,E,G,I,J,K,L,N,O]}
        # G|{mix: [N], pos: [A,C,D,H,L,M], neg: [B,F,I,J,K,O]}
        # H|{mix: [B,C,D,F,G,I,J,K,L,N,O], pos: [A,E,M], neg: []}
        # I|{mix: [K], pos: [A,B,C,D,H,L,M,O], neg: [F,G,J,N]}
        # J|{mix: [], pos: [A,H,M], neg: [B,C,D,F,G,I,K,L,N,O]}
        # K|{mix: [], pos: [A,B,D,H,I,L,M,O], neg: [F,G,J,N]}
        # L|{mix: [B,G,I,K,O], pos: [A,C,D,H,M], neg: [F,J,N]}
        # M|{mix: [B,C,D,E,F,G,I,J,K,L,N,O], pos: [A,H], neg: []}
        # N|{mix: [], pos: [A,C,D,G,H,M], neg: [B,F,I,J,K,O]}
        # O|{mix: [B,I,K,L], pos: [A,C,D,E,H,M], neg: [F,G,J,N]}

        #}}}
        # HIDE-ME: rules {{{

        # mix: (1|2) means 1 is above 2
        # pos: (1|2) means 1 is a direct ancestor or direct descendant or dupe, not a "cousin", "uncle", or 
        #      "sibling" - to 2. (differ to dupe in cases where there's no other clues)
        # neg: (1|2) means 1 is a "cousin" or "uncle" or "sibling" or direct ancestor of 2

        # }}}
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

        #init sort logging and var prep
        STASH = {}
        print("===")
        print("variant tree sort start")
        print("---")

        #show pre-proc default data 
        self.stdout_dump_variant_relations_data(DATA,'pre-proc',run)
        #prep stash
        for key, value in DATA.items():
            STASH[key] = {'mix':[],'pos':[],'neg':[]}

        #prep ref
        if self.REF is None:
            self.REF = {}
            for key, value in DATA.items():
                self.REF[key] = {'r_mix':[],'r_pos':[],'r_neg':[]}

        #MIX RULE CHKS{{{

        print("mix-checks {"+"{{") #beg collapse vim marker
        for key, value in DATA.items():
            if run_all or run_mix:
                for Vz in value['mix']:

                    #chk1 - if the two nodes are on the same level, then: create a default parental relation + don't STASH
                    if self.TREE[Vz].parent == self.TREE[key].parent:
                        if self.DBG > 1:
                            print("---")
                            print("MIX-CHK1 - "+key+"|"+Vz+" - parents are same level - so put "+Vz+" under "+key)
                        ch1 = list(self.TREE[Vz].parent.children).remove(self.TREE[Vz])
                        ch2 = self.TREE[Vz].children
                        if ch1 is None:
                            self.TREE[Vz].parent = None
                        else:
                            self.TREE[Vz].parent.children = tuple(ch1)
                        self.TREE[Vz] = Node(Vz, parent=self.TREE[key])
                        self.TREE[Vz].children = ch2
                        if key not in self.REF[Vz]['r_mix']:
                            self.REF[Vz]['r_mix'].append(key)
                        if self.DBG > 1:
                            for pre, fill, node in RenderTree(self.TREE['top']):
                                print("%s%s" % (pre, node.name))

                    #chk2 - if there is already a good direct line established, then: don't STASH
                    elif self.TREE[Vz] in self.TREE[key].descendants:
                        if self.DBG > 1:
                            print("MIX-CHK2 - "+key+"|"+Vz+" - condition satisfied - "+Vz+" already under "+key)
                        if key not in self.REF[Vz]['r_mix']:
                            self.REF[Vz]['r_mix'].append(key)

                    #chk3 - if anything else, then: STASH
                    else:
                        if self.DBG > 1:
                            print("---")
                        print("MIX-CHK3 - "+key+"|"+Vz+" - parents not same level and not direct lineage - put in STASH")
                        STASH[key]['mix'].append(Vz)

            else:
                #since we're skipping "mix" relations, they need to go to STASH
                STASH[key]['mix'] = value['mix']
        print("}"+"}}") #end collapse vim marker

        #}}}
        #POS RULE CHKS{{{

        print("pos-checks {"+"{{") #beg collapse vim marker
        for key, value in DATA.items():
            if run_all or run_pos:
                for Vz in value['pos']:
                    #chk1 - if the two nodes have direct line relation, then: don't STASH
                    if self.TREE[Vz] in self.TREE[key].descendants or self.TREE[key] in self.TREE[Vz].descendants:
                        if self.DBG > 1:
                            print("---")
                            print("POS-CHK1 - "+key+"|"+Vz+" - direct lineage relation found - put in STASH")
                        if key not in self.REF[Vz]['r_pos']:
                            self.REF[Vz]['r_pos'].append(key)
                    #chk2 - if anything else, then: STASH
                    else:
                        if self.DBG > 1:
                            print("---")
                        print("POS-CHK2 - "+key+"|"+Vz+" - no direct lineage found - put in STASH")
                        STASH[key]['pos'].append(Vz)
                        #print('key-desc: '+str(self.TREE[key].descendants))
                        #print('Vz-desz:'+str(self.TREE[Vz].descendants))
                        #print(STASH)
                        #sys.exit()
            else:
                #since we're skipping "pos" relations, they need to go to STASH
                STASH[key]['pos'] = value['pos']
        print("}"+"}}") #end collapse vim marker

        #}}}
        #NEG RULE CHKS{{{

        print("neg-checks {"+"{{") #beg collapse vim marker
        for key, value in DATA.items():
            if run_all or run_neg:
                for Vz in value['neg']:
                    #chk1 - if the two nodes don't have direct line relation, then: don't STASH
                    if self.TREE[Vz] not in self.TREE[key].descendants and self.TREE[key] not in self.TREE[Vz].descendants:
                        if self.DBG > 1:
                            print("---")
                            print("NEG-CHK1 - "+key+"|"+Vz+" - no direct lineage relation found - put in STASH")
                        if key not in self.REF[Vz]['r_neg']:
                            self.REF[Vz]['r_neg'].append(key)
                    #chk2 - if the two nodes don't have direct line relation, then: don't STASH
                    elif self.TREE[Vz] in self.TREE[key].descendants:
                        if self.DBG > 1:
                            print("---")
                            print("NEG-CHK2 - "+key+"|"+Vz+" - anc to dec relation found - put in STASH")
                        if key not in self.REF[Vz]['r_neg']:
                            self.REF[Vz]['r_neg'].append(key)
                    #chk3 - if anything else, then: STASH
                    else:
                        if self.DBG > 1:
                            print("---")
                        print("NEG-CHK3 - "+key+"|"+Vz+" - direct lineage found - put in STASH")
                        STASH[key]['neg'].append(Vz)
                        #print('key-desc: '+str(self.TREE[key].descendants))
                        #print('Vz-desz:'+str(self.TREE[Vz].descendants))
                        #print(STASH)
                        #sys.exit()
            else:
                #since we're skipping "neg" relations, they need to go to STASH
                STASH[key]['neg'] = value['neg']
        print("}"+"}}") #end collapse vim marker

        #}}}

        #show the final tree diagram after run completion
        print("---")
        print("RUN:"+str(run)+" DONE")
        for pre, fill, node in RenderTree(self.TREE['top']):
            print("%s%s" % (pre, node.name))
        
        #show post-proc remaining data that didn't get complete (now in STASH)
        print("---")
        self.stdout_dump_variant_relations_data(STASH,'post-proc',run)
        print("---")

        return STASH
        
    def stdout_dump_variant_relations_data(self,DATA,dataStr,run=1):
        dataPrint = copy.deepcopy(DATA)
        refPrint = copy.deepcopy(self.REF)
        mixlen = 0
        poslen = 0
        neglen = 0
        chkRef1 = False
        if self.REF is not None:
            r_mixlen = 0
            r_poslen = 0
            r_neglen = 0
        for key, value in DATA.items():
            dataPrint[key]['mix'] = ','.join(map(str, value['mix']))
            dataPrint[key]['pos'] = ','.join(map(str, value['pos']))
            dataPrint[key]['neg'] = ','.join(map(str, value['neg']))
            mixlen = mixlen+len(value['mix'])
            poslen = poslen+len(value['pos'])
            neglen = neglen+len(value['neg'])
            if self.REF is not None:
                refPrint[key]['r_mix'] = ','.join(map(str, self.REF[key]['r_mix']))
                refPrint[key]['r_pos'] = ','.join(map(str, self.REF[key]['r_pos']))
                refPrint[key]['r_neg'] = ','.join(map(str, self.REF[key]['r_neg']))
                r_mixlen = r_mixlen+len(self.REF[key]['r_mix'])
                r_poslen = r_poslen+len(self.REF[key]['r_pos'])
                r_neglen = r_neglen+len(self.REF[key]['r_neg'])
        print("RUN:"+str(run)+"("+dataStr+") mix cnt:"+str(mixlen))
        print("RUN:"+str(run)+"("+dataStr+") pos cnt:"+str(poslen))
        print("RUN:"+str(run)+"("+dataStr+") neg cnt:"+str(neglen))
        if self.REF is not None:
            print("RUN:"+str(run)+"("+dataStr+") r_mix cnt:"+str(r_mixlen))
            print("RUN:"+str(run)+"("+dataStr+") r_pos cnt:"+str(r_poslen))
            print("RUN:"+str(run)+"("+dataStr+") r_neg cnt:"+str(r_neglen))
        print("---")
        #beg collapse vim marker
        print("RUN:"+str(run)+"("+dataStr+") DATA/STASH dump {"+"{{")
        self.stdout_dump_var(dataPrint)
        #end collapse vim marker
        #beg collapse vim marker
        if self.REF is not None:
            print("}"+"}}")
            print("RUN:"+str(run)+"("+dataStr+") REF dump {"+"{{")
            self.stdout_dump_var(refPrint)
        #end collapse vim marker
        print("}"+"}}")
        

    #MISC DEF

    def stdout_dump_var(self,var):
        #TODO: put this somewhere else
        print(json.dumps(var, indent=4, sort_keys=True))

