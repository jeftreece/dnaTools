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
1. split checker ... anything with a split needs to be flagged. doesn't get an ambiguous bump up
1. A297 junk rule
2. redo the horiz/vert sort
3. get_mx_data should exlude all neg variants too
4. hide top from matrix stdout
5. do I still need self.perfect_variants, etc?
6. write good notes and commenting for the rules
7. put debugging code into subset/supset area
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
        #self.perfect_variants = None
        #self.imperfect_variants = None
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
                    if kv != 0:
                        sql1 = "INSERT into s_calls (kit_id,variant_loc,assigned) VALUES ('k%s','%s',%s);" % (k,vv,kv)
                        self.dbo.sql_exec(sql1)

    # matrix

    def sort_matrix(self):

        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #get data
        self.get_mx_data()

        #stdout relations data
        if config['DEBUG_RELATIONS']:
            self.stdout_mx_relations_data()

        #step 0
        debug_chk('DEBUG_MATRIX',"data - step 0 (default)",1)
        self.stdout_tbl_mx()

        #step 1
        debug_chk('DEBUG_MATRIX',"data - step 1",1)
        self.sort_step1()
        self.stdout_tbl_mx()

        #step 2
        debug_chk('DEBUG_MATRIX',"data - step 2",1)
        self.sort_step2()
        self.stdout_tbl_mx()

        #step 3
        debug_chk('DEBUG_MATRIX',"data - step 3",1)
        self.sort_step3()
        self.stdout_tbl_mx()

        #step 4
        #debug_chk('DEBUG_MATRIX',"data - step 4",1)
        #self.sort_step4()

        sys.exit()

    def sort_step1(self):
        self.mx_vertical_sort_new()
        
    def sort_step2(self):
        self.mx_horizontal_sort()
        
    def sort_step3(self):
        print("Processing Imperfect Variants:")
        for impVix in self.get_imperfect_variants_idx():
            print("----------------")
            print("impVix: %s"%impVix)
            print("impVix(names): %s"%self.get_vname_by_vix(impVix))
            print("{{"+"{") #beg vim marker
            kzc = self.get_kixs_by_val(val=0,vix=impVix)
            print("kzc: %s"%kzc)
            print("kzc(names): %s"%self.get_kname_by_kix(kzc))
            results = []
            #sys.exit()
            print("")
            #supsetsP = self.get_supset_variants(override_val=1,vix=impVix,kix=unk[1],convertToNames=False,perfectFlg=True)
            #supsets = self.get_supset_variants(vix=impVix,kix=unk[1],convertToNames=False,perfectFlg=True)
            #subsetsP = self.get_subset_variants(override_val=1,vix=impVix,kix=unk[1],convertToNames=False,perfectFlg=True)
            #subsets = self.get_subset_variants(vix=impVix,kix=unk[1],convertToNames=False,perfectFlg=True)
            #rule0 - A297 - needs supsets{{{

            #results.append(['R0',self.test_rule0_consistency(vix=impVix,kzc=kzc,subsets=subsets,supsets=supsets)])

            #}}}
            #rule1 - M301 - needs supsetsP (ambiguous promotion){{{

            #good - k|v: kD|Z301:[10,0]
            #good - RULE2: True - ambiguous: equivalent/subset of L48
            #good - k|v: kI|Z301:[10,1]
            #good - RULE2: True - ambiguous: equivalent/subset of L48

            #bad - missing k|v: kB|Z301
            #bad - missing k|v: kG|Z301

            #Problem: ???

            results.append(['R1',self.test_rule1_supsets(vix=impVix,kzc=kzc)])

            #}}} 
            #rule1a - Z381 - needs subsets{{{

            #good - k|v: kD|Z381:[9,0]
            #good - RULE1: True: sup of L48

            #bad - k|v: kI|A297:[12,1]
            #bad - RULE1: True: sup of Z9
            #bad - k|v: kB|A297:[12,2]
            #bad - RULE1: True: sup of Z156

            #Problem: it doesn't check sup consistency

            #results.append(['R1a',self.test_rule1a_subsets(unk_variant=unk,subsets=subsets,supsets=supsets,supsetsP=supsetsP)])

            #}}}
            #rule3 - Z28 - needs supsets (2 diff branches - False){{{

            #bad - it's doing false for everything
            #Problem: ???
            #results.append(['R3',self.test_rule3_diff_branches(unk_variant=unk,subsets=subsets,supsets=supsets)])

            #}}}
            #print("")
            print("}}"+"}") #end vim marker
            for R in results:
                print(R[1])

        #unresolved list
        print("----------------")
        print("")
        #print("unresolved")
        #print(self.get_coord_name_by_order(self.unk_variants))
        #resolved list
        #print("resolved")
        #print(self.resolved_variants)
        #stdout relations data
        #if config['DEBUG_RELATIONS']:
        #    self.stdout_mx_relations_data()
        #print("")
        #sys.exit()

        self.get_mx_count_data()
        self.mx_vertical_sort()
        self.mx_horizontal_sort()

    def test_rule1_supsets(self,vix,kzc):

        #M301 rule (part of it){{{

        #Note: if there is an existing perfect variant that has the same 
        #(or superset of the) known +/-'s to the imperfect variant when 
        #the unresolved ? is turned to a +, then it's a +.
        #Note: this rule requires supsetsP!

        #find the the kits that perfectly match positives for filling any
        #combination of unk values

        #what is the kpc of this variant}}}

        if config['DEBUG_RULE1']:
            print("[R1.0] START RULE 1")
            print("")

        sups4vix = self.get_supset_variants(vix=vix,convertToNames=False,perfectFlg=True)

        if config['DEBUG_RULE1']:
            print("[R1.3] sups4vix: %s" %sups4vix)
            print("[R1.4] sups4vix(names): %s" % self.get_vname_by_vix(vix=sups4vix))
            print("")

        subs4vix = self.get_subset_variants(vix=vix,convertToNames=False,perfectFlg=True)

        if config['DEBUG_RULE1']:
            print("[R1.1] subs4vix: %s" %subs4vix)
            print("[R1.2] subs4vix(names): %s" % self.get_vname_by_vix(vix=subs4vix))
            print("")

        overrideData = self.get_row_when_override_kixs(1,vix=vix,kixs=kzc)
        kpuc = np.argwhere(overrideData[0,] == 1).T[1,]
        eqv1 = np.argwhere(self.NP[:,kpuc]==1)[:,0]
        eqv2 = np.unique(eqv1)

        if config['DEBUG_RULE1']:
            #print ("[R1.5] kpuc equiv (not quite): %s:"%eqv1)
            print("[R1.3] kpuc: %s:"%kpuc)
            print("[R1.4] kpuc(names): %s" % self.get_kname_by_kix(kix=kpuc))
            print("[R1.5] kpuc equiv/uniq (not quite): %s:"%eqv2)
            print("[R1.6] kpuc equiv/uniq/names: %s" % self.get_vname_by_vix(vix=eqv2))
            print("")

        #{{{
        #what could it be with all unks to positive?
        #self.get_kixs_by_val(1,
        #what have the min kpc of the bottom variant?
        #excluding to just the kpc of the max variant?

        #vix = unk_variant[0]
        #kix = unk_variant[1]
        #supsets = self.get_supset_variants(vix=vix,kix=kix,convertToNames=False,perfectFlg=True)

        #overrideData = self.get_row_when_override_coord(1,kix=kix,vix=vix)

        #kpc = self.get_kixs_by_val(1,vix=vix,overrideData=overrideData)
        #knc = self.get_kixs_by_val(-1,vix=vix)
        #}}}

        subs4kpuc = self.get_subset_variants(vix=vix,convertToNames=False,perfectFlg=True,kpc=kpuc)
        sups4kpuc = self.get_supset_variants(vix=vix,convertToNames=False,perfectFlg=True,kpc=subs4kpuc)

        if config['DEBUG_RULE1']:
            print("[R1.5] subs4kpuc: %s" %subs4kpuc)
            print("[R1.6] subs4kpuc(names): %s" % self.get_vname_by_vix(vix=subs4kpuc))
            print("")
            print("[R1.7] sups4kpuc: %s" %sups4kpuc)
            print("[R1.8] sups4kpuc(names): %s" % self.get_vname_by_vix(vix=sups4kpuc))
            print("")

        print("")

        #Need to redo how I'm approaching this 
        sys.exit()

        if config['DEBUG_RULE1']:
            print("")
            print("[R1.1]!!!kpc: kit positive conditions: %s" %kpc)
            print("[R1.2]!!!kpc(names): %s"%self.get_kname_by_kix(kix=kpc))
            print("[R1.3]!!!knc: kit negative conditions: %s" %knc)
            print("[R1.4]!!!knc(names): %s"% self.get_kname_by_kix(kix=knc))
            print("[R1.5]!!!supsetsP: %s" %supsetsP)
            print("[R1.6]!!!supsetsP(names): %s" % self.get_vname_by_vix(vix=supsetsP))

        for sup in reversed(supsetsP):

            if config['DEBUG_RULE1']:
                print("")
                print("[R1.5]!!!sup: %s" % sup)
                print("[R1.6]!!!sup(name): %s" % self.get_vname_by_vix(vix=sup))

            kpc4sup = self.get_kixs_by_val(1,vix=sup)

            if config['DEBUG_RULE1']:
                print("[R1.7]!!!kpc4sup: %s"%kpc4sup)
                print("[R1.8]!!!kpc4sup(name): %s"%self.get_kname_by_kix(kix=kpc4sup))

            knc4sup = self.get_kixs_by_val(-1,vix=sup)

            if config['DEBUG_RULE1']:
                print("[R1.9]!!!knc4sup: %s"%knc4sup)
                print("[R1.10]!!!knc4sup(name): %s"%self.get_kname_by_kix(kix=knc4sup))

            lenPos = len(list(set(kpc).intersection(set(kpc4sup))))
            lenNeg = len(list(set(knc).intersection(set(knc4sup))))

            if config['DEBUG_RULE1']:
                print("[R1.11]!!!lenPos: %s" % lenPos)
                print("[R1.12]!!!lenNeg: %s" % lenNeg)

            if lenPos == len(kpc) and lenNeg == len(knc):
                self.unk_variants.remove(unk_variant)
                self.NP[unk_variant[0],unk_variant[1]]=1
                self.resolved_variants.append((unk_variant,True))

                print("")
                return 'RULE1: True - ambiguous: equivalent/subset of %s'%(self.get_vname_by_vix(vix=sup))

        #if we get here ... go to rule 3
        print("")
        return 'RULE1: Unk'
        return self.test_rule3_diff_branches(unk_variant,subsets,supsets)
        
    def test_rule1a_subsets(self,unk_variant,subsets,supsets,supsetsP):
        #How this works:
        #1. check to see if this coord's related vpc variants are subsets of this coord's variant.

        #Z381 rule - needs subsets

        kix = unk_variant[1]
        vpc = self.get_vixs_by_val(1,kix=kix)

        if config['DEBUG_RULE1a']:
            print("[R1a.1] vpc: positive conditions: %s" %vpc)
            print("[R1a.2] vpc(names): %s"%self.get_vname_by_vix(vix=vpc))
            print("[R1a.3] subsets: %s"%subsets)
            print("[R1a.4] subsets(names): %s"%self.get_vname_by_vix(vix=subsets))

        #standard subset test
        for sub in subsets:
            if sub in vpc:

                if config['DEBUG_RULE1a']:
                    print("")
                    print("[R1a.5] sub: %s"%sub)
                    print("[R1a.6] sub(name):"%self.get_vname_by_vix(vix=sub))

                self.unk_variants.remove(unk_variant)
                self.NP[unk_variant[0],unk_variant[1]]=1
                self.resolved_variants.append((unk_variant,True))

                if config['DEBUG_RULE1a']:
                    print("")

                return 'RULE1a: True: sup of %s' % self.get_vname_by_vix(vix=sub)

        #Note: if get here -- go to rule 2
        return "RULE1a: Unk"
        return self.test_rule2_supsets(unk_variant=unk_variant,subsets=subsets,supsets=supsets,supsetsP=supsetsP)
        
    def test_rule3_diff_branches(self,unk_variant,subsets,supsets):

        #Z28 rule

        vix = unk_variant[0]
        kix = unk_variant[1]
        #supsets = self.get_supset_variants(vix=unk_variant[0],convertToNames=False) #superset of given coord
        #supsets.append(-999) #top
        #subsets = self.use_imperfect_variants_only(self.get_subset_variants(vix=unk_variant[0],convertToNames=False)) #superset of given coord
        #vi = self.use_imperfect_variants_only(self.get_vixs_by_val(1,kix=kix))
        vi = self.get_vixs_by_val(1,kix=kix)
        #vi.append(-999) #top
        ki = self.get_kixs_by_val(1,vix=vix).tolist()
        if ki is None:
            ki = []
        rule_p1_list = []

        if config['DEBUG_RULE3']:
            print("[R3a.1] vix: %s" %self.get_vname_by_vix(vix))
            print("[R3a.2] supsets: %s" %self.get_vname_by_vix(supsets))
            print("[R3a.3] vi (imperfect knowns when set coord to pos): %s" %self.get_vname_by_vix(vi))
            print("[R3a.4] ki (related kits when set coord to pos): %s" %self.get_kname_by_kix(ki))
            print("")

        #Note: first deal with the variant that might might share the unknown variant we're wondering about
        #Q: Is it not a direct relation?
        #Q: does it have a superset?
        for V in vi:

            #sups4V = self.use_imperfect_variants_only(self.get_supset_variants(vix=V,convertToNames=False)) #superset of related coord
            sups4V = self.get_supset_variants(vix=V,convertToNames=False) #superset of related coord
            #sups4V.append(-999) #top

            if config['DEBUG_RULE3']:
                print("[R3a.5] V: %s" %self.get_vname_by_vix(V))
                print("[R3a.6] sups4V(%s): %s" %(self.get_vname_by_vix(V),self.get_vname_by_vix(sups4V)))
                print("[R3a.6] supsets: %s" %supsets)
                print("[R3a.7] subsets: %s" %subsets)
                print("[R3a.8] supsets(names): %s" %self.get_vname_by_vix(supsets))
                print("[R3a.9] subsets(names): %s" %self.get_vname_by_vix(subsets))

            #Is it a direct relation? 
            if V not in supsets and V not in subsets: # and vix not in sups4V:
                if config['DEBUG_RULE3']:
                    print("")
                directRelation_chk = False
                if config['DEBUG_RULE3']:
                    print("[R3a.10] directRelation Chk is: %s (continue with V=%s)" %(directRelation_chk,self.get_vname_by_vix(V)))
                    print("")
            else:
                directRelation_chk = True
                if config['DEBUG_RULE3']:
                    print("[R3a.11] directRelation Chk is: %s (don't continue with V=%s)" %(directRelation_chk,self.get_vname_by_vix(V)))
                    print("")
                continue # try the next V

            if directRelation_chk == False:

                if config['DEBUG_RULE3']:
                    print("[R3a.12] V in: %s" %self.get_vname_by_vix(V))
                k4V = self.get_kixs_by_val(1,vix=V).tolist() #other kits per V
                k4V.remove(kix)
                if config['DEBUG_RULE3']:
                    print("[R3a.13] k4V(%s): %s" %(self.get_vname_by_vix(V),self.get_kname_by_kix(k4V)))

                #Does it have any supersets?
                if len(sups4V) == 0:
                    #sups4V.append(-999) #top
                    if config['DEBUG_RULE3']:
                        print("[R3a.14] sups Chk: sups not in sups4V (don't continue with V=%s)"%self.get_vname_by_vix(V))
                        print("")
                    continue # try the next V

                for sup4V in sups4V:
                    if config['DEBUG_RULE3']:
                        print("")
                        print("sups Chk: sups in sups4V (continue with V=%s)"%self.get_vname_by_vix(V))
                        print("")
                        print("[R3a.15] sup4V(%s): %s" %(self.get_vname_by_vix(V),self.get_vname_by_vix(sup4V)))
                    k4sup4V = self.get_kixs_by_val(1,vix=sup4V).tolist() #other kits per V
                    if config['DEBUG_RULE3']:
                        print("[R3a.16] k4sup4V(%s)(%s)(bef.rem) in: %s" %(self.get_vname_by_vix(V),self.get_vname_by_vix(sup4V),self.get_kname_by_kix(k4sup4V)))
                    if kix in k4sup4V: #(1)can't be the given coord's kit
                        k4sup4V.remove(kix) #(2)can't be the given coord's kit
                    if config['DEBUG_RULE3']:
                        print("[R3a.17] k4sup4V(%s)(%s)(aft.rem) in: %s" %(self.get_vname_by_vix(V),self.get_vname_by_vix(sup4V),self.get_kname_by_kix(k4sup4V)))
                    k_in_k4V_and_k4sup4V = list(set(k4V).intersection(set(k4sup4V)))

                    #Is there additional overlap btw this other variant and its superset?
                    if len(k_in_k4V_and_k4sup4V) == 0:
                        if config['DEBUG_RULE3']:
                            print("[R3a.18] k4V + k4sup4V intersection chk: no additional overlap (don't continue with V=%s)"%self.get_vname_by_vix(V))
                            print("")
                        continue

                    #Everything seems ok for part one of this rule
                    if config['DEBUG_RULE3']:
                        print("[R3a.19] k4V + k4sup4V intersection chk: additional overlap (continue with V=%s)"%self.get_vname_by_vix(V))
                        print("")
                        print("[R3a.20] k_in_k4V_and_k4sup4V: %s" %self.get_kname_by_kix(k_in_k4V_and_k4sup4V))
                        print("")
                    rule_p1_list.append((sup4V,k_in_k4V_and_k4sup4V))

        if len(rule_p1_list):
            rule_p2_list = []

            #Loop the known supersets of the given variant
            for sup4vo in supsets:

                #What positive variant relations do those supersets have? 
                k4sup4vo = self.get_kixs_by_val(1,vix=sup4vo).tolist() #other kits per V
                #k4sup4vo.append(-999) #top
                
                if config['DEBUG_RULE3']:
                    print("[R3b.21] k4sup4vo(%s): %s"%(self.get_vname_by_vix(sup4vo),self.get_kname_by_kix(k4sup4vo)))
                    print("[R3b.22] ki: %s"%self.get_kname_by_kix(ki))
                    print("")

                #Is there any overlap with the given coord?
                k_in_k4sup4vo_and_ki = set(k4sup4vo).intersection(set(ki))
                if len(k_in_k4sup4vo_and_ki) == 0:
                    print("[R3b.23] k4sup4vo + ki intersection chk: no additional overlap (don't continue with sup4vo=%s)"%self.get_vname_by_vix(sup4vo))
                    print("")
                    continue

                #If so, everything seems ok for part two of this rule
                if config['DEBUG_RULE3']:
                    print("[R3b.24] k4sup4vo + ki intersection chk: additional overlap (continue with sup4vo=%s)"%self.get_vname_by_vix(sup4vo))
                    print("")
                    print("[R3b.25] k_in_k4sup4vo_and_ki in: %s" %self.get_kname_by_kix(k_in_k4sup4vo_and_ki))
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
                                #print("}}"+"}")
                                msg1a = self.get_vname_by_vix(itm1[0])
                                msg1b = self.get_vname_by_vix(itm1[1])
                                msg2a = self.get_vname_by_vix(itm2[0])
                                msg2b = self.get_vname_by_vix(itm2[1])
                                print("")
                                return "RULE3: False - 2 diff branches (%s,%s) + (%s,%s)"% (msg1a,msg1b,msg2a,msg2b)

        #print("}}"+"}")
        if config['DEBUG_RULE3']:
            print("")

        return "RULE3: Unk"
        return "RULEX: Unk"
    def test_rule0_consistency(self,vix,subsets,supsets):
        '''{{{

        [step 1 - rule1 - M301]{{{ 

        If we make the ? a +, are there any +/- combos (superset candidates) like this already

         [b] Yes - go to standard Rule 2 (convert "?" to "+" - and flag as ambiguuous)
        x[a] No - go to step 2

        }}}
        [step 2]{{{

        (step 2a) ok, let's try the + for the unk.
        (step 2b) ok, let's try it w/o the pos for unk. (ph.consistency)

        }}}
        [step 3]{{{
        
        What is its minimum superset [(a) with (b) w/o] considering the unk call?

        x[A] U106

        (step 3) What are the unique highest subset combo children of v1 (other than v5) that
        touch +'s (in equal or different ways that v5 does to the superset)?

        v2--> L48
        v3--> Z156

        Are they subsets? No

        Would changing unknowns to positive, make any all of them subsets? Yes
        If can change all of them ... then promote.

        If no...still can't ... step 4
        }}}
        [step 4]{{{
        
        With these subsets, do the following study:

        (a) One for they way they intersect the superset (v1).
        (b) One for they way they intersect the superset (v1) and the variant being studied (v5).
        (c) Add up the total unique connections btw both lists
        (d) Note the unique difference in kits btw both lists

        Intersection List1          Intersection List2
        ------------------          ------------------
        L48:  v1+v2 [D,I,C,H]       L48: v1+v2+v5 [C,H] D,I <- possible
        Z156: v1+v3 [B,G]           Z156: v1+v3+v5 []   B,G <- possible
        -------                     --------
        total = 6                   total = [2] 4<- possible (due to unks)

        (step 5)
        
        Make conclusions for split/no split

        Conclusions:

        (1) No diff btw two lists - ambiguous  promote all shared points btw the two lists that might be unk
        (2) Intersection List 1 is greater == split
        (3) Intersection List 2 is greater == flag for study

        }}}
        [step 5 - split process/loop]{{{

        (step 1) move Zxx to Zxx.tmp
        (step 2) of v1's highest children, which intersects v5 the best? v2 (or other v1 child)
        (step 3) splice Zxx.tmp to create a Zxx.# and if anything left unexplained, a new Zxx.tmp that
                 accomodates v2 (or other v1 child) as a parent
        (step 4) repeat steps 2-3 until no more Zxx.tmp

        }}}

        }}}'''

        #A297 rule

        vix = unk_variant[0]
        kix = unk_variant[1]
        kpc = self.get_kixs_by_val(val=1,vix=vix)

        if config['DEBUG_RULE0']:
            print("[R0.1] kix: %s" % kix)
            print("[R0.2] vix: %s" % vix)
            print("[R0.8] kpc: %s" % kpc)
            print("[R0.8] kpc(names): %s" % self.get_kname_by_kix(kpc))

        if config['DEBUG_RULE0']:
            print("[R0.3] supsets: %s" % supsets)
            print("[R0.4] supsets(names): %s" % self.get_vname_by_vix(vix=supsets))
            #print("[R0.5] !!! - subsets: %s" % subsets)
            #print("[R0.6] !!! - subsets(names): %s" % self.get_vname_by_vix(vix=subsets))

        if len(supsets):
            sup = supsets[len(supsets)-1]
            kpc4sup = self.get_kixs_by_val(val=1,vix=sup)

            if config['DEBUG_RULE0']:
                print("")
                print("[R0.7] sup: %s" % self.get_vname_by_vix(vix=sup))
                print("[R0.8] kpc4sup: %s" % kpc4sup)
                print("[R0.8] kpc4sup(names): %s" % self.get_kname_by_kix(kpc4sup))

            subs4sup = self.get_subset_variants(vix=sup,perfectFlg=True,convertToNames=False)

            if config['DEBUG_RULE0']:
                print("[R0.9] subs4sup: %s" % subs4sup)
                print("[R0.10] subs4sup(names): %s" % self.get_subset_variants(vix=sup,perfectFlg=True))
                print("")

        #for sup in supsets:
        #    if config['DEBUG_RULE0']:
        #        print("")
        #        print("[cons.7] !!! - sup: %s" % self.get_vname_by_vix(vix=sup))
        #    kpc = self.get_kixs_by_val(val=1,vix=sup)
        #    if config['DEBUG_RULE0']:
        #        print("[cons.8] !!! - kpc: %s" % kpc)
        #        print("[cons.8] !!! - kpc(names): %s" % self.get_kname_by_kix(kpc))

        return 'RULE0: Unk'
        

    def get_subset_variants(self,override_val=None,vix=None,kix=None,convertToNames=True,perfectFlg=False,kpc=None):

        #-------------------------------------------------------------------------------------------------
        # How subsets are determined:
        #-------------------------------------------------------------------------------------------------
        # get_subsets (private routine)
        # - vix - the row-idx of the reference variant
        # - kpc - the col kit-idxs where this vix is "actively" seen as having the given "val" value
        # - VAR1 - evaluate the incoming kpc for positives along the matrix - create an index with the results
        # - VAR2 - delete any mention of the vix in these index results (since we are looking for other variants)
        # - VAR3 - get a count of how many times these other variants hit the kpc we're looking at
        # - out/out1 - can't remember what the out logic does again, check that
        # - VAR4 - then filter out those variants that have less kpc hits than the vix
        # - VAR5 - the supersets
        # - VAR6 - the superset counts to vpc
        # - VAR7 - merge of VAR5 and VAR6
        #-------------------------------------------------------------------------------------------------

        def get_subsets(kpc,vix):
            VAR1 = np.argwhere(self.NP[:,kpc]==1)[:,0] #looking for variants w/pos assignments like the incoming variant's pos conditions
            if config['DBG_SUBS_IN']:
                print("[subin.1] VAR1: %s"%VAR1)
            idx = np.argwhere(VAR1==vix) #idx make sure we exclude the incoming variant
            VAR2 = np.delete(VAR1, idx)
            if config['DBG_SUBS_IN']:
                print("[subin.2] VAR2p: %s"%VAR2)
            VAR3 = self.get_uniq_counts(VAR2)
            if config['DBG_SUBS_IN']:
                print("[subin.3] VAR3: %s"%VAR3)
            #...
            #Note: for the following "adding technique" -- we need to exclude unk situations for the comparison (special handling)
            #beg - adding technique - got this idea here: https://stackoverflow.com/questions/30041286/sum-rows-where-value-equal-in-column
            unq, unq_inv = np.unique(VAR3[:,0], return_inverse=True)
            out = np.zeros((len(unq), VAR3.shape[1]), dtype=VAR3.dtype) #create empty array to put the added values
            out[:, 0] = unq #fill the first column
            np.add.at(out[:, 1:], unq_inv, VAR3[:, 1:])
            if config['DBG_SUBS_IN']:
                print("[subin.4] out: %s"%out)
            #end - adding technique
            #Note: sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
            out1 = out[out[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
            if config['DBG_SUBS_IN']:
                print("[subin.4a] out1: %s"%out1)
                print("[subin.5] kpc: %s"%kpc)
                print("[subin.6] kpc(names): %s"%self.get_kname_by_kix(kpc))
            VAR4 = np.argwhere(out1[:,1]<len(kpc)) #[:,0]
            if config['DBG_SUBS_IN']:
                print("[subin.7] VAR4: %s"%VAR4)
            #end - sorting technique
            out2a = out1[:,0]
            if config['DBG_SUBS_IN']:
                print("[subin.8] out2a: %s"%out2a)
            out2b = out1[:,1]
            if config['DBG_SUBS_IN']:
                print("[subin.9] out2b: %s"%out2b)
            VAR5 = out2a[list(VAR4.T[0])] #these are the subset vixs ids (in order, max first)
            if config['DBG_SUBS_IN']:
                print("[subin.10] VAR5: %s"%VAR5)
                print("[subin.11] VAR5(names): %s"%self.get_vname_by_vix(VAR5))
            VAR6 = out2b[list(VAR4.T[0])]#these are the subset variant counts
            if config['DBG_SUBS_IN']:
                print("[subin.12] VAR6: %s"%VAR6)
            VAR7 = np.asarray((VAR5,VAR6)).T #merged for return
            if config['DBG_SUBS_IN']:
                print("[subin.13] VAR7: %s"%VAR7)
            #if len(VAR6) == 0:
            #    return []
            return VAR7 #[:,0]

        #Note: override one specific coord
        if override_val is not None and kix is not None:
            overrideData = self.get_row_when_override_coord(override_val,kix=kix,vix=vix)
            kpc = self.get_kixs_by_val(1,vix=vix,overrideData=overrideData) #pos conditions when override coord with a value
        elif kpc is None:
            kpc = self.get_kixs_by_val(1,vix=vix) #default pos conditions
    
        if convertToNames is True:
            subs = self.get_vname_by_vix(vix=self.filter_perfect_variants(get_subsets(kpc,vix)[:,0].tolist()),listFlg=1)
        else:
            subs = self.filter_perfect_variants(vix=get_subsets(kpc,vix)[:,0].tolist())

        #Debugging stdout msgs
        if config['DBG_SUBS']:
            print("[sub.1] subs: %s" % subs)
        if config['DBG_SUBS']:
            suffix=''
            if override_val is not None and kix is not None:
                suffix = 'P' if override_val == 1 else 'N'
            print("[sub.2] subsets%s: %s kpc: %s" % (suffix,",".join([str(i) for i in subs]),kpc))

        return subs
        
    def get_supset_variants(self,override_val=None,vix=None,kix=None,convertToNames=True,perfectFlg=False,kpc=None):

        #-------------------------------------------------------------------------------------------------
        # How supsets are determined:
        #-------------------------------------------------------------------------------------------------
        # get_subsets (private routine)
        # - vix - the row-idx of the reference variant
        # - kpc - the col kit-idxs where this vix is "actively" seen as having the given "val" value
        # - VAR1 - evaluate the incoming kpc for positives along the matrix - create an index with the results
        # - VAR2 - get a count of how many times these other variants hit the kpc we're looking at
        # - VAR3 - has to have at least what the incoming variant had in count
        # - VAR4 - delete any mention of the vix in these index results (since we are looking for other variants)
        # -      - give an opportunity to return [] if no results
        # - allPos - all idx situations in the matrix where there are positives (we just want the col-variant idxs)
        # - AP   - a unique count on how many positives there were per variant of allPos
        # - VAR5 - use the VAR4 and AP idxs together to the hit counts on those variants that can be considered supersets
        # - VAR6 - merge VAR4 and VAR5 together to combine superset idxs with their hit counts
        # - VAR7 - reverse sort VAR6 - so the bigger ones are first
        #-------------------------------------------------------------------------------------------------

        def get_supsets(kpc,vix):
            if config['DBG_SUPS_IN']:
                print("[supin.0] kpc: %s"%kpc)
                print("[supin.0] kpc(names): %s"%self.get_kname_by_kix(kpc))
            VAR1 = np.argwhere(self.NP[:,kpc]==1)[:,0] #looking for variants w/pos assignments like the incoming variant condition
            if config['DBG_SUPS_IN']:
                print("[supin.1] VAR1: %s"%VAR1)
            VAR2 = self.get_uniq_counts(VAR1)
            if config['DBG_SUPS_IN']:
                print("[supin.2] VAR2: %s"%VAR2)
            VAR3 = VAR2[VAR2[:,1]==len(kpc)] #has to have at least what the incoming variant had in count
            if config['DBG_SUPS_IN']:
                print("[supin.3] VAR3: %s"%VAR3)
            idx = np.argwhere(VAR3[:,0]==vix) #idx make sure we exclude the incoming variant
            VAR4 = np.delete(VAR3[:,0], idx) #idx again/delete
            if config['DBG_SUPS_IN']:
                print("[supin.4] VAR4: %s"%VAR4)
            #...
            if len(VAR4) == 0: return [] #there are no supsets according to the filter
            #(beg) master list of all positives
            allPos = np.argwhere(self.NP==1)[:,0]
            if config['DBG_SUPS_IN']:
                print("[supin.5] allPos: %s"%allPos)
            unqA, cntA = np.unique(allPos, return_counts=True)
            AP = np.asarray((unqA, cntA))[1,]
            if config['DBG_SUPS_IN']:
                print("[supin.6] AP: %s"%AP)
            #(end) master list of all positives
            VAR5 = AP[list(VAR4),] #extrapolate the right subset mix to master list of all positives
            if config['DBG_SUPS_IN']:
                print("[supin.7] VAR5: %s"%VAR5)
            VAR6 = np.asarray((VAR4,VAR5)).T
            if config['DBG_SUPS_IN']:
                print("[supin.8] VAR6: %s"%VAR6)
            #Note: sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
            VAR7 = VAR6[VAR6[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
            if config['DBG_SUPS_IN']:
                print("[supin.9] VAR7: %s"%VAR7)
            #(end) sorting
            return VAR7

        #Note: override one specific coord
        if override_val is not None and kix is not None:
            overrideData = self.get_row_when_override_coord(override_val,kix=kix,vix=vix)
            kpc = self.get_kixs_by_val(1,vix=vix,overrideData=overrideData) #pos conditions when override coord with a value
        elif kpc is None:
            kpc = self.get_kixs_by_val(1,vix=vix) #default pos conditions

        #Note: get supsets
        supsY = get_supsets(kpc,vix)

        #Debugging
        if config['DBG_SUPS']:
            print("")
            print("[sup.1] kpc: %s"%kpc)
            print("[sup.2] supsY: %s"%supsY)

        #Perfect variants
        if len(supsY) == 0:
            return []
        else:
            supsX = self.filter_perfect_variants(supsY.tolist()) #.tolist())

        #Debugging
        if config['DBG_SUPS']:
            print("[sup.3] supsX: %s"%supsX)

        #Convert Names or not
        if len(supsX) == 0:
            sups = supsX
        elif convertToNames is True:
            sups = self.get_vname_by_vix(vix=supsX[:,0],perfectFlg=perfectFlg,listFlg=1)
        else:
            sups = np.asarray(supsX)[:,0]

        #Debugging
        if config['DBG_SUPS']:
            print("[sup.4] sups: %s"%sups)
        if override_val is not None and kix is not None:
            if override_val==1:
                if config['DBG_SUPS']:
                    print("[sup.5] supsetsP: "+",".join([str(i) for i in sups]) +" kpc:"+str(kpc))
            else:
                if config['DBG_SUPS']:
                    print("[sup.6] supsetsN: "+",".join([str(i) for i in sups]) +" kpc:"+str(kpc))
        else:
            if config['DBG_SUPS']:
                print("[sup.7] supsets: "+",".join([str(i) for i in sups]) +" kpc:"+str(kpc))

        return sups
        
    def use_perfect_known_variants_only(self,variants):
        if self.perfect_known_variants == None:
            return variants
        newList = []
        for v in variants:
            if v in self.perfect_known_variants:
                newList.append(v)
        return newList

    def get_row_when_override_kixs(self,override_val,vix,kixs):
        row = self.get_mx_kdata(vix=vix)
        if len(kixs) == 0:
            return row
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        for kx in kixs.tolist():
            rowO[0,kx] = override_val
        return rowO

    def get_vname_by_vix(self,vix,listFlg=False,perfectFlg=False): #get vname (can also take a list)
        #listFlg: force listFlg as return data type
        #perfectFlg: force it to be a name that's part of the imperfect known list of variants
        intFlg = True
        try:
            value = int(vix)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            variant = None
            #normal variants in the matrix
            for itm in list(self.VARIANTS.items()):
                if itm[1][1] == vix:
                    if listFlg:
                        if perfectFlg is False or self.perfect_variants is None:
                            variant = itm[0]
                        elif itm[0] in self.perfect_variants:
                            variant = itm[0]
                    else:
                        if perfectFlg is False or self.perfect_variants is None:
                            variant = itm[0]
                        elif itm[1][1] in self.perfect_variants:
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
            #print("here1am2")
            for vo in list(vix):
                #normal variants in the matrix
                for itm in list(self.VARIANTS.items()):
                    if itm[1][1] == vo:
                        if perfectFlg is False or self.perfect_variants is None:
                            variantList.append(itm[0])
                        elif itm[1][1] in self.perfect_known_variants:
                            variantList.append(itm[0])
                        break
            return(variantList)
        
    def get_kname_by_kix(self,kix,listFlg=False):
        #listFlg: force listFlg as return data type
        intFlg = True
        try:
            value = int(kix)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            kit = None
            for itm in list(self.KITS.items()):
                if itm[1][1] == kix:
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
            for ko in list(kix):
                for itm in list(self.KITS.items()):
                    if itm[1][1] == ko:
                        kitList.append(itm[0])
                        break
            return(kitList)
        
    def get_vix_by_name(self,vname): #get vix from its name
        return self.VARIANTS[vname][1]
        
    def get_kix_by_name(self,kname): #get kix from its name
        return self.KITS[kname][1]
        

    def get_kixs_by_val(self,val,vix=None,vname=None,overrideData=None): #like get_mx_kdata but retrieves index info for given val
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None and overrideData is not None: # we're sending in a custom evaluation
            return np.argwhere(overrideData[0,] == val).T[1,] #with override data, there's only one line evaluated - 1d datset
        if vix is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[vix,] == val).T[1,] #default data, it's the entire matrix - 2d dataset 
        
    def get_vixs_by_val(self,val,kix=None,kname=None,overrideData=None): #like get_mx_variant_data but retrieves index info for given val
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None and overrideData is not None:
            return np.argwhere(overrideData[:,0] == val).T[0,] #with override data, there's only one line evaluated - 1d dataset
        if kix is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[:,kix] == val).T[0,] #default data, it's the entire matrix - 2d dataset

    #...

    def mx_vertical_sort_new(self):
        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        for K,V in self.get_axis('variants'):
            if -1 not in self.get_mx_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('vp'):
            if -1 in self.get_mx_row_as_list(V[1]) and 0 not in self.get_mx_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('variants'):
            if -1 in self.get_mx_row_as_list(V[1]) and 0 in self.get_mx_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for NO in new_orders:
            self.set_new_order(NO[0],NO[1],variantType=True)
        self.NP = np.matrix(list(DATA.values()))
        
    def mx_horizontal_sort_new(self):
        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        self.NP = np.transpose(self.NP)
        for K,V in self.get_axis('kp'):
            new_orders.append([K,cnt])
            DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
            cnt = cnt + 1
        for NO in new_orders:
            self.set_new_order(NO[0],NO[1],kitType=True)
        self.NP = np.matrix(list(DATA.values()))
        self.NP = np.transpose(self.NP)
        
    def mx_vertical_sort(self):
        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        for K,V in self.get_axis('variants'):
            if -1 not in self.get_mx_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('vp'):
            if -1 in self.get_mx_row_as_list(V[1]) and 0 not in self.get_mx_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for K,V in self.get_axis('variants'):
            if -1 in self.get_mx_row_as_list(V[1]) and 0 in self.get_mx_row_as_list(V[1]):
                new_orders.append([K,cnt])
                DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
                cnt = cnt + 1
        for NO in new_orders:
            self.set_new_order(NO[0],NO[1],variantType=True)
        self.NP = np.matrix(list(DATA.values()))
        
    def mx_horizontal_sort(self):
        DATA = OrderedDict()
        cnt = 0 
        new_orders = []
        self.NP = np.transpose(self.NP)
        for K,V in self.get_axis('kp'):
            new_orders.append([K,cnt])
            DATA[K] = self.get_mx_row_as_list(V[1],noneToStr=False)
            cnt = cnt + 1
        for NO in new_orders:
            self.set_new_order(NO[0],NO[1],kitType=True)
        self.NP = np.matrix(list(DATA.values()))
        self.NP = np.transpose(self.NP)

    def stdout_variant(self):
        foo = 1

    def stdout_tbl_mx(self):
        debug_chk('DEBUG_MATRIX',"",1)
        debug_chk('DEBUG_MATRIX',"big_matrix view{{"+"{",1)
        debug_chk('DEBUG_MATRIX',"",1)
        table = BeautifulTable()
        table.column_headers = ['']+self.get_cur_kit_list()
        for K,V in self.get_axis('variants'):
            table.append_row([K]+self.get_mx_row_as_list(V[1]))
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
        
    def stdout_mx_relations_data(self,dataStr='',run=0):
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

    def get_coord(self,kix,vix,moreInfo=False):
        buf = ""
        if moreInfo:
            buf = "coord: "+str(kix)+","+str(vix)
        kit = self.get_kname_by_kix(kix)
        variant = self.get_vname_by_vix(vix)
        buf = buf +  "k|v: "+str(kit)+"|"+(variant)
        if moreInfo:
            buf = buf + "value:"+str(self.NP[vix,kix])
        return buf
        
    def get_cur_kit_list(self):
        return self.get_axis('kits',keysOnly=True)
        
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
        
    def get_coord_name_by_order(self,coord_order):
        coord_name = []
        for C in coord_order:
            vname = self.get_vname_by_vix(vix=C[0])
            kname = self.get_kname_by_kix(kix=C[1])
            coord_name.append((vname,kname))
        return coord_name

    def get_mx_row_as_list(self,rownum,noneToStr=True):
        if noneToStr:
            return ['None' if v is None else v for v in self.NP[rownum,:].tolist()[0]]
        else:
            return self.NP[rownum,:].tolist()[0]
        
    def get_mx_kdata(self,vix=None,vname=None): #get same type variant data for each kit
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None:
            return self.NP[vix,]
        
    def get_mx_variant_data(self,kix=None,kname=None): #get all type variant data for one kit
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None:
            return self.NP[:,kix].T
        

    def get_imperfect_variants(self):
        #TODO: need to fix this. I had the def wrong.
        variant_idx = np.argwhere(self.NP==-1) #get index to negative data
        variant_idx_r = np.unique(variant_idx[:,0]) #get unique rows of those indices
        variant_data = self.NP[variant_idx_r]
        #print(variant_data)
        #sys.exit()
        vnames = self.get_vname_by_vix(variant_idx_r)
        self.set_new_axis(vnames,variant_idx_r,variantType=True) #reset variant axis
        #sys.exit()
        #self.NP = None
        #print(variant_data)
        self.NP = variant_data #reset data
        #print(self.NP)
        #sys.exit()
        
    def get_perfect_variants_idx(self):
        idx = list(range(len(self.VARIANTS)))
        prf_idx = idx[:] #copy idx so I have this one for deleting
        unk_idx = list(np.unique(np.argwhere(self.NP==0)[:,0]))
        for x in idx:
            if x in unk_idx:
                prf_idx.remove(x)
        return idx
        
    def get_imperfect_variants_idx(self):
        #print(list(np.unique(np.argwhere(self.NP==0)[:,0]))[:]) #make it a copy
        #sys.exit()
        return list(np.unique(np.argwhere(self.NP==0)[:,0]))[:] #make it a copy
        
    def filter_perfect_variants(self,vix):
        #print("filter(bef) -  vix: %s" % vix)
        if any(isinstance(el, list) for el in vix):
            for itm in reversed([[n,v] for (n,(v,c)) in enumerate(vix)]):
                #print(vix)
                if itm[1] in self.get_imperfect_variants_idx():
                    #print(itm)
                    #print(itm[0])
                    vix.remove(vix[itm[0]])
            #print("filter(aft.1) -  vix: %s" % vix)
            return vix
        else:
            for itm in reversed([[n,v] for n,v in enumerate(vix)]):
                #print(itm)
                if itm[1] in self.get_imperfect_variants_idx():
                    #print(itm)
                    #print(itm[0])
                    vix.remove(vix[itm[0]])
            #print("filter(aft.2) -  vix: %s" % vix)
            return vix

    def get_coord_value(self,kix=None,vix=None,kname=None,vname=None):
        if kix is not None and vix is not None:
            return self.NP[vix][kix]
        if kname is not None and vname is not None:
            get_kix_by_name(kname)
            get_vix_by_name(vname)
            return self.NP[vix][kix]
        

    def get_row_when_override_coord(self,override_val,kix=None,vix=None,kname=None,vname=None):
        #override_val -- is the override val (ie: check what conditions are after setting a coord to be 1 and not 0, for example)
        row = self.get_mx_kdata(vix=vix)
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        rowO[0,kix] = override_val
        return rowO

    def get_mx_data(self):

        #sql - exclude perfect variants
        sql = "select * from perfect_variants_with_kits_assignments_and_unk;"

        #get data
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        #retrieve data from sqlite like so: [V][K] [x,x,x,x,x,x,...]
        DATA = OrderedDict()
        self.KITS = {}
        self.VARIANTS = {}
        cntV = 0
        cntK = 0
        for row in F:
            if row[1] not in DATA:
                DATA[row[1]] = []
            DATA[row[1]].append(row[2])
            #TODO: numpy way to do this? still needed? get rid of dupe info?
            if row[0] not in self.KITS:
                self.KITS[row[0]] = [cntK,cntK]
                cntK = cntK + 1
            #TODO: numpy way to do this? still needed? get rid of dupe info?
            if row[1] not in self.VARIANTS:
                self.VARIANTS[row[1]] = [cntV,cntV]
                cntV = cntV + 1
      
       #create numpy version of data
        for key,value in DATA.items():
            self.NP = np.matrix(list(DATA.values()))

        #chk matrix (debugging)
        #self.stdout_tbl_mx()
        #sys.exit()

        #get count data
        self.get_mx_count_data()

        #get relations data
        self.get_mx_relations_data()
    def get_variant_data_by_vname(self,vname):

        #sql - exclude perfect variants
        sql = '''
            SELECT C.kit_id, V.name, C.assigned, V.variant_id
            FROM s_calls C, s_variants V,
            WHERE C.variant_loc = V.variant_loc AND
            V.vname = %s
            ORDER by 4;
            ''' %vname

        #get data
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()

        #retrieve data from sqlite like so: [V][K] [x,x,x,x,x,x,...]
        DATA = OrderedDict()
        self.KITS = {}
        self.VARIANTS = {}
        cntV = 0
        cntK = 0
        for row in F:
            if row[1] not in DATA:
                DATA[row[1]] = []
            DATA[row[1]].append(row[2])
            #TODO: numpy way to do this? still needed? get rid of dupe info?
            if row[0] not in self.KITS:
                self.KITS[row[0]] = [cntK,cntK]
                cntK = cntK + 1
            #TODO: numpy way to do this? still needed? get rid of dupe info?
            if row[1] not in self.VARIANTS:
                self.VARIANTS[row[1]] = [cntV,cntV]
                cntV = cntV + 1
      
       #create numpy version of data
        for key,value in DATA.items():
            self.NP = np.matrix(list(DATA.values()))

        #chk matrix (debugging)
        #self.stdout_tbl_mx()
        #sys.exit()

        #get count data
        self.get_mx_count_data()

        #get relations data
        self.get_mx_relations_data()

    def get_mx_count_data(self):

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
        
    def get_mx_relations_data(self):
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

    def get_uniq_counts(self,arr):
        unique_elements, counts_elements = np.unique(arr, return_counts=True)
        return np.asarray((unique_elements, counts_elements)).T #get uniques

    # not used

    def get_min_superset_variant(self,vix=None,vname=None): #order is variant's order in matrix, name is vname
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        sups = self.get_supset_variants(vix=vix,perfectFlg=perfectFlg)
        if len(sups) is 0:
            return None
        else:
            return sups[0]

