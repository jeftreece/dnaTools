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
    
def il2s(lst):
    return ",".join(str(x) for x in lst)
    
def l2s(lst):
    return ",".join(lst)

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

#TODO: order of priority
'''
1. (x) split checker ... anything with a split needs to be flagged. doesn't get an ambiguous bump up
2. redo the horiz/vert sort
'''

class Variant(object):

    def __init__(self):
        self.dbo = None

    def proc(self,vname):
        self.sort.get_mx_data(recreateFlg = False)
        if vname.isdigit():
            vix = int(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        self.vix = vix
        self.set_info(lev=2)
        self.chk_var(allowImperfect=False)
        #self.sort.test_rule1_supsets(vix=vix)
        
    def info(self,vname):
        allowImperfect = config['allowImperfectWithVInfo']
        self.sort.get_mx_data(recreateFlg = False)
        if vname.isdigit():
            vix = int(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        self.vix = vix
        self.set_info(lev=2,allowImperfect=allowImperfect)
        self.stdout_info()
        
    def set_info(self,vix=None,lev=1,allowImperfect=False):

        lev = lev-1
        if vix is not None:
            self.vix = vix
        self.name = self.sort.get_vname_by_vix(self.vix)
        self.vixn = "%s - [%s]"%(self.name,self.vix)
        self.kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        self.knc = self.sort.get_kixs_by_val(val=-1,vix=self.vix)
        self.kuc = self.sort.get_kixs_by_val(val=0,vix=self.vix)
        self.kpcn = "kpc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kpc)),il2s(self.kpc))
        self.kncn = "knc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.knc)),il2s(self.knc))
        self.kucn = "kuc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kuc)),il2s(self.kuc))
        self.kpnc = sorted(self.kpc + self.knc)

        self.sups = self.get_rel(relType=1,allowImperfect=allowImperfect)
        self.sups.remove(0) #remove 'top'
        self.subs = self.get_rel(relType=-1,allowImperfect=allowImperfect)
        self.eqv = self.get_rel(relType=0,allowImperfect=allowImperfect)
        self.subsn = "subs: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.subs)),il2s(self.subs))
        self.supsn = "sups: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.sups)),il2s(self.sups))
        self.eqvn = "eqv: %s [%s]"%(l2s(self.sort.get_vname_by_vix(vix=self.eqv)),il2s(self.eqv))

        if lev > 0:
            self.supOs = []
            self.subOs = []
            for sup in self.sups:
                supO = Variant()
                supO.sort = self.sort
                supO.dbo = self.dbo
                supO.set_info(sup,lev,allowImperfect)
                self.supOs.append(supO)
            for sub in self.subs:
                subO = Variant()
                subO.sort = self.sort
                subO.dbo = self.dbo
                subO.set_info(sub,lev,allowImperfect)
                self.subOs.append(subO)
        
    def stdout_info(self,tp=None):

        print("")
        if tp is None:
            print("---------------------------------------------------------------------")

        if tp==1:
            print("[+] %s" %self.vixn)
            sp = "    "
        elif tp==-1:
            print("[-] %s" %self.vixn)
            sp = "    "
        else:
            print("vix: %s" %self.vixn)
            sp = ""
        print("%s%s" %(sp,self.kpcn))
        print("%s%s" %(sp,self.kncn))
        print("%s%s" %(sp,self.kucn))

        print("%s%s" %(sp,self.supsn))
        print("%s%s" %(sp,self.subsn))
        print("%s%s" %(sp,self.eqvn))

        noSupOs = 0
        noSubOs = 0
        try:
            self.supOs
        except:
            noSupOs = 1
        try:
            self.subOs
        except:
            noSubOs = 1

        if noSupOs == 0:
            for sup in self.supOs:
                if sup.name != 'top':
                    sup.stdout_info(tp=1)
        if noSubOs == 0:
            for sub in self.subOs:
                sub.stdout_info(tp=-1)

        if tp is None:
            print("---------------------------------------------------------------------")
            print("")
        
    def chk_var(self,allowImperfect):

        '''{{{
        [9] Z381:{{{
        ---------------------------------------------------------------------------------------------------------------

        +'s: B,C,G,H,I
        U: D

        sups are:

        U106:     pP:D              mN:E,J                   subBef: Yes  rel: sub (cuz of A-F-)

        So: with the sup issue alone, it's an ambiguous promotion.

        any other sup U106 subs (that overlap Z381) - that Z381 doesn't see as eq/sub?:
        Yes:

        L48:   mP(Z381):C,H,I   mP(U106):C,D,H,I       xP:D   commonP: C,H,I
        Z9:    mP(Z381):I       mp(U106):D,I           xP:D   commonP: I

        These xP's ... resolved with any unk's? Yes -- so then the D is promotion confirmed.

        }}}
        [10] Z301:{{{
        ---------------------------------------------------------------------------------------------------------------

        +'s: C,H
        -'s: F
        U: A,B,D,E,G,I,J

        sups are:

        L48:  pP:D       pN:A,B,E,J mU:-  mP:C,H  mN: F      subBef: Yes  rel: Subset/Equal+
        U106: pP:A,B,D   pN:E,J     mU:-  mP:C,H  mN: -      subBef: Yes  rel: Subset (cuz of F-)

        So: sub2U106 means: whatever u106 is negative for ... Z301 has to be neg too: E,J
        So: sub/eq+L48 means -- outside of E and J (unk>neg) and C,H (already +), the others could be - or + (in this case, we put ambiguous positve for all)

        any other sup L48/U106 subs (that overlap Z301) - that Z301 doesn't see as eq/sub?:
        No

        }}}
        [11] Z28:{{{
        ---------------------------------------------------------------------------------------------------------------
         
        +'s: D,I (Z9 eq)
        U: F

        sups are:


        L48:               pN:F     mN:A,B,E,G,J             subBef: Yes  rel: sub (cuz of C-H-)
        U106:     pP:F     pN:      mN:E,J                   subBef: Yes  rel: eq

        So: sub2L48 means: whatever L48 is negative for ... Z28 has to be neg too: F

        any other sup L48/U106 subs (that overlap Z28) - that Z28 doesn't see as eq/sub?:
        No

        }}}
        [12] A297:{{{
        ---------------------------------------------------------------------------------------------------------------

        +'s: D,G
        U: B,I

        sups are:

        U106:   pP:B,I              mN:E,J                 subBef: Yes  rel: sub (cuz of A-C-F-H-)

        So: with the sup issue alone: it's an ambiguous promotion for B,I

        any other sup U106 subs (that overlap A297) - that Z381 doesn't see as eq/sub?:
       
        L48:   mP(A297):D     mP(U106):C,D,H,I     xP:C,H,I       commonP: D (split.1)
        Z9:    mP(A297):D     mP(U106):I           xP:I           commonP: D

        Can the unk's resolve the XP's? No
        What commonP can be split? D (split.1)
        What remains? G (split.2)

        Everything else becomes neg

        }}}
        }}}'''

        print("")
        print("---------------------------------------------------------------------")
        print("")

        #starting unk list
        unkL = self.kuc
        unkD = {}
        for k in unkL:
            unkD[k] = [0,0]
        othK = {}

        #1st arg of set : sup test1 result
        #0 - nothing observed yet
        #1 = ambiguous promotion to positive
        #2 = hard promotion to positive
        #3 = hard negative
        
        #2nd arg of set : sup test2 result
        #0 - nothing observed yet
        #1 - means (so far it's seen as filling a need)

        #look for sups
        if len(self.sups):

            print("vix: %s [%s]"%(self.sort.get_vname_by_vix(self.vix),self.vix))
            print("unks: %s [%s]" % (l2s(self.sort.get_kname_by_kix(self.kuc)),il2s(self.kuc)))
            print("")
            for sup in reversed(self.sups):

                splitL = []
                print("sup: %s [%s]" % (self.sort.get_vname_by_vix(sup),sup))

                #(beg)sup check
                knc4sup = self.sort.get_kixs_by_val(val=-1,vix=sup)
                diffKnc = list(set(self.knc)-set(knc4sup))
                if len(diffKnc) > 0:
                    pN = list(set(self.kuc).intersection(set(knc4sup)))
                    if len(pN) > 0:
                        print(" - [1] vix unk %s [%s] is/are neg" % (il2s(self.sort.get_kname_by_kix(pN)),il2s(pN)))
                        for k in pN:
                            if k in unkL:
                                unkL.remove(k)
                            unkD[k][0] = 3 #hard failure
                        print(" - [1] remaining unk: %s [%s]" %(l2s(self.sort.get_kname_by_kix(unkL)),il2s(unkL)))
                    else:
                        print(" - [1] sup is truly sup to target variant, but all unks open to promotion")
                else:
                    print(" - [1] sup is sub/eq/sup to target variant")
                #for whatever hasn't been noted as a failure ... promote 
                for k in unkL:
                    if unkD[k] != 3:
                        unkD[k][0] = 1 #ambiguous promotions
                #(end)sup check

                #(beg)consistency check 
                v1 = Variant()
                v1.dbo = self.dbo
                v1.sort = self.sort
                v1.vix = sup
                v1.set_info()

                #remove any target subs or target equivalent variants from consideration
                v1subs = list(set(v1.subs)-set(self.subs)-set(self.eqv))

                #do consistency checks on remaining variants
                if len(v1subs):

                    VAR1a = np.argwhere(self.sort.NP[:,self.kpc] == 1)
                    mask = np.isin(VAR1a[:,0], v1subs)
                    idx = list(list(np.where(mask))[0])
                    VAR2a = np.unique(VAR1a[:,0][idx]) #these are the overlapping v1subs with target v
                    for v in VAR2a:
                        akpc = self.sort.get_kixs_by_val(val=1,vix=v)
                        #(mtP) common kpc btw v1sub and target variant
                        at_common_kpc = list(set(akpc).intersection(set(self.kpc)))
                        #(msP) common kpc btw v1sub and common sup
                        as_common_kpc = list(set(akpc).intersection(set(v1.kpc)))
                        #(xP) diff_common_kpc = [msP-mtP]
                        diff_common_kpc = list(set(as_common_kpc)-set(at_common_kpc))
                        #(cP) common_kp c = intersection btw msP and mtP
                        common_kpc = list(set(as_common_kpc).intersection(set(at_common_kpc)))
                        print(" - [2] v1sub: %s [%s]"%(self.sort.get_vname_by_vix(v),v))
                        print("       (mtP) shared btw vix + v1sub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(at_common_kpc)),il2s(at_common_kpc)))
                        print("       (msP) shared btw sup + v1sub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(as_common_kpc)),il2s(as_common_kpc)))
                        print("       (xP) msP-mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(diff_common_kpc)),il2s(diff_common_kpc)))
                        print("       (cP) common btw msP+mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(common_kpc)),il2s(common_kpc)))
                        for k in diff_common_kpc:
                            if k in unkL:
                                unkD[k][1] = unkD[k][1] + 1
                        if len(common_kpc) > 0:
                            splitL.append(common_kpc)
                        
                unique_splits = [list(x) for x in set(tuple(x) for x in splitL)]
                if len(unique_splits) > 0:
                    split_intersects = list(set.intersection(*map(set, unique_splits)))
                    if (len(unique_splits) > 1 and len(split_intersects) == 0):
                        print(" - [2] split required: btw %s" % unique_splits)
                #(end)consistency check 

                print("")

        print(unkD)

        '''{{{
        [9] Z381{{{

        any other sup U106 subs (that overlap Z381) - that Z381 doesn't see as eq/sub?:
        Yes:

               *                                       *
        L48:   mtP(Z381):C,H,I   msP(U106):C,D,H,I       xP:D   commonP: C,H,I
        Z9:    mtP(Z381):I       msP(U106):D,I           xP:D   commonP: I

        These xP's ... resolved with any unk's? Yes -- so then the D is promotion confirmed.

        if whatever is in XP is the only think missing and it's und ... that's
        the one.
        }}}
        [12] A297{{{

        any other sup U106 subs (that overlap A297) - that Z381 doesn't see as eq/sub?:
       
               *                                       *
        L48:   mtP(A297):D       msP(U106):C,D,H,I       xP:C,H,I       commonP: D (split.1)
        Z9:    mtP(A297):D       msP(U106):D,I           xP:I           commonP: D
        Z156   mtP(A297):G       msP(U106):B,G           xP:B           commonP: G

        Can the unk's resolve the XP's? No
        What commonP can be split? D (split.1)
        What remains? G (split.2)

        Everything else becomes neg
        }}}
        }}}'''

        sys.exit()

        '''{{{
        #target vix's combined kpc + knc idxs
        kpnc = sorted(self.kpc+self.knc)
        kpnc = self.kpc

        #NP data based on target vix's kpnc
        kpncNP = self.sort.NP[:,kpnc]

        ##NP data based on  target vix's kuc
        ##kucNP = self.sort.NP[:,self.kuc]
        ##print("")
        ##print(kucNP)
        ##print("")

        #target vix's data for its kpnc
        VkpncNP = kpncNP[self.vix]

        #other vixes that match target vix's kpnc based data
        VAR1 = np.argwhere(np.all(self.sort.NP[:,kpnc]==VkpncNP,axis=1)==True)[:,0]

        #other vixes (w/cnts) that have 1's somewhere in target vix's kuc 
        VAR2 = np.argwhere(self.sort.NP[:,self.kuc] == 1)

        #uniq VAR2 (just the vixes)
        VAR3 = np.unique(VAR2[:,0])

        #vix's are in both VAR1 and VAR2 lists
        VAR4 = list(set(list(VAR1)).intersection(list(VAR3)))
        #print(self.sort.get_vname_by_vix(VAR4))


        vList = []
        vnList = []
        if len(VAR4):
            for v in reversed(VAR4):
                #pos variants that 2nd variant has, that target vix doesn't 
                x = list(set(self.sort.get_kixs_by_val(val=1,vix=v))-set(self.kpc))
                #neg variants that 2nd variant has, that target vix doesn't 
                y = list(set(self.sort.get_kixs_by_val(val=-1,vix=v))-set(self.knc))
                if v != 0:
                    vList.append((v,x,y))
                    vnList.append((self.sort.get_vname_by_vix(v), l2s(self.sort.get_kname_by_kix(x)), l2s(self.sort.get_kname_by_kix(y))))
                #vkpc = self.sort.get_kix_by
                #self.sort.NP[:,kpnc]P
                #self.sort.NP[:,kpnc]

        print("")
        print("vix: %s]" % self.sort.get_vname_by_vix(self.vix))
        print("kpc: %s [%s]" % (il2s(self.kpc),l2s(self.sort.get_kname_by_kix(self.kpc))))
        print("knc: %s [%s]" % (il2s(self.knc),l2s(self.sort.get_kname_by_kix(self.knc))))
        print("kuc: %s [%s]" % (il2s(self.kuc),l2s(self.sort.get_kname_by_kix(self.kuc))))
        print("")
        print("vList(#'s)")
        print("")
        print(vList)
        print("")
        print("vList(names) - vixes that have all the +'s tvix has and then some, extra +'s to unks, shared -'s")
        print("")
        print(vnList)
        print("")

        
        #[(9, [1, 2, 4]), (2, [0, 1])]
        #...
        #what needs to happen:  is there any overlap btw results here? 
        #...
        #scenarios:
        #...
        #(1) completely separate (and target variant is parent to them both)
        #(2) they completely match - so they're dupes
        #(3) one is a superset of another
        #(4) they overlap but also have distinct things too. in this case ...  ???
        #...
        #how target variant fits in: 
        #...
        #(1) assume target variant is a parent - so matches all  positive vals
        #(2) assume target variant is a dupe as well ... matches all positive vals
        #(3) assume target variant is a dupe with the bigger variant - matches all +'s
        #...
        #is there something I need to do about checking the commonality with a
        #superset variant too? (like with subsets?)
        #what to do with other unks? (that don't fit?)

        #print(VAR4)
        #print(VAR2)

        #sorter = np.argsort(VAR2[:,0)
        #        idx = sorter[np.searchsorted(VAR2[:,0],,sorter=sorter)]
        #...
        #kpncNP = NP[:,kpnc]
        #VAR3 = np.argwhere(VAR2[:,self.kuc] == 1)

        #print(VAR1)
        #print(self.sort.get_vname_by_vix(VAR1))
        #(x) I want to just filter the kpc + knc columns
        #(x) vals = [-1,0]
        #(x) sorter = np.argsort(self.sort.NP)
        #(x) idx = sorter[np.searchsorted(self.sort.NP,vals,sorter=sorter)]
        #(x) [1] kpc+knc of the variant when not optimizing unk to pos  
        #(x) [1a] what other variants have an exact match at this kpc+knc?

        #[2] then further filter out those variants that have additional +'s for when the
        #    target has unks.
        #    if so, the variants could be at least this.(or a subset of these)
        #[2] do they have any common supsets where these "additional variants"
        #    are also shared?
        #    if so, there's an opportunity for even more +'s?
        #[3] are there any irregularities that prevent a promotion?
        }}}'''

    def get_rel(self,relType,override_val=None,kix=None,kpc=None,allowImperfect=False):

        #-------------------------------------------------------------------------------------------------{{{
        # 3 Ways to use this routine
        #-------------------------------------------------------------------------------------------------
        # 1. can use it with an override_val and a kix -- in this situation, it will override a given kit 
        #    (ie: one with an unkown val, to a positive val) and then determine the subsets of this 
        #    variant -- again with a routine generated kpc
        # 2. can use it with a vix (it gets the default subs that way), the routine creates a kpc
        # 3. can use it with a given kpc to override the kpc generated by this routine if you use methods 
        #    1 or 2. this is useful if you want to check subsets when overriding multiple kit unknown 
        #    values for a given variant
        #-------------------------------------------------------------------------------------------------
        # 3 Types of relations
        #-------------------------------------------------------------------------------------------------
        # 1. relType ==  1 supset
        # 2. relType == -1 subset
        # 3. relType ==  0 equivalent
        #-------------------------------------------------------------------------------------------------}}}

        #method 1
        if override_val is not None and kix is not None:
            if config['DBG_RELS']:
                print("")
                print("[rels.0] MODE (get_vix_rel) - ONE (override val w/vix)")
                print("")
            overrideData = self.sort.get_row_when_override_coord(override_val,kix=kix,vix=self.vix)
            #pos conditions when override coord with a value
            kpc = self.sort.get_kixs_by_val(1,vix=self.vix,overrideData=overrideData)

        #method 2
        elif kpc is None:
            if config['DBG_RELS']:
                print("")
                print("[rels.0] MODE (get_vix_rel) - TWO (only use vix)")
                print("")
            #default pos conditions when just use default target vix criteria
            kpc = self.sort.get_kixs_by_val(1,vix=self.vix) #default pos conditions
    
        #method3
        #requires only sending in a kpc
        else:
            if config['DBG_RELS']:
                print("")
                print("[rels.0] MODE (get_vix_rel) - THREE (sending in a kpc)")
                print("")

        if relType == 1: #supsets
            rels_ = self.get_supsets_or_eqv(kpc)
        elif relType == -1: #subsets
            rels_ = self.get_subsets(kpc)
        else: # relType == 0: subsets (eq_override flag)
            #Note: this creates data that needs to be compared to supsets, to
            #get equal variants ... this is def a hack!
            rels_ = self.get_supsets_or_eqv(kpc=kpc,eq_override=True)
        
        #if there's nothing here, best to just return
        if len(rels_) == 0:
            return []
        #two ways to return the data. with or without counts.
        #rels = self.filter_perfect_variants(vix=rel_.tolist()).tolist()
        if allowImperfect is False:
            rels = self.sort.filter_perfect_variants(vix=rels_[:,0].tolist()) #.tolist()
        else:
            rels = rels_[:,0].tolist()

        #debugging
        if config['DBG_RELS']:
            print("[rels.1] rels: %s" % rels)
        if config['DBG_RELS']:
            suffix=''
            if override_val is not None and kix is not None:
                suffix = 'P' if override_val == 1 else 'N'
            print("[rels.2] rels%s: %s kpc: %s" % (suffix,",".join([str(i) for i in rels]),kpc))

        return list(rels)
        
    def get_subsets(self,kpc):

        #-------------------------------------------------------------------------------------------------{{{
        # How subsets are determined:
        #-------------------------------------------------------------------------------------------------
        # - VAR1 - evaluate the incoming kpc for positives along the matrix - create an index with the results
        # - VAR2 - delete any mention of the vix in these index results (since we are looking for other variants)
        # - VAR3 - get a count of how many times these other variants hit the kpc we're looking at
        # - out/out1 - can't remember what the out logic does again, check that
        # - VAR4 - then filter out those variants that have less kpc hits than the vix
        # - VAR5 - the supersets
        # - VAR6 - the superset counts to vpc
        # - VAR7 - merge of VAR5 and VAR6
        #-------------------------------------------------------------------------------------------------}}}

        #of the top ones ... are they supsets? subsets? (start with top 2 - top)
        #what additional positives might they have that my variant has an unk for?
        #perfect fit? if so .. done

        if config['DBG_SUBS_IN']:
            print("")
            print("---------------------------------------------------------------------")
            print("")
            print("[subin.0] vix: %s"%self.vix)
            print("[subin.0] vix(name): %s"%self.sort.get_vname_by_vix(self.vix))
            print("[subin.0] kpc: %s"%kpc)
            print("[subin.0] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #looking for variants w/pos assignments that overlap the target variant
        VAR1 = np.argwhere(self.sort.NP[:,kpc]==1)[:,0]
        if config['DBG_SUBS_IN']:
            print("[subin.1] VAR1: %s"%VAR1)

        #idx make sure we exclude the target variant
        idx = np.argwhere(VAR1==self.vix)
        VAR2 = np.delete(VAR1, idx)
        if config['DBG_SUBS_IN']:
            print("[subin.2] VAR2: %s"%VAR2)

        #print (set(x for l in VAR1.tolist() for x in l))
        #get uniques of VAR2 with counts
        unique_elements, counts_elements = np.unique(VAR2, return_counts=True)
        VAR3 = np.asarray((unique_elements, counts_elements)).T #get uniques
        if config['DBG_SUBS_IN']:
            print("")
            print("[subin.3] VAR3: %s"%VAR3)

        #for the following "adding technique" -- we need to exclude unk situations for the comparison (special handling)
        #adding technique - got this idea here: https://stackoverflow.com/questions/30041286/sum-rows-where-value-equal-in-column
        unq, unq_inv = np.unique(VAR3[:,0], return_inverse=True)
        out = np.zeros((len(unq), VAR3.shape[1]), dtype=VAR3.dtype) #create empty array to put the added values
        out[:, 0] = unq #fill the first column
        np.add.at(out[:, 1:], unq_inv, VAR3[:, 1:])
        if config['DBG_SUBS_IN']:
            print("[subin.4] out: %s"%out)

        #sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
        out1 = out[out[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
        if config['DBG_SUBS_IN']:
            print("[subin.4a] out1: %s"%out1)
            print("[subin.5] kpc: %s"%kpc)
            print("[subin.6] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #subsets need to have less kpc than the target variant (cuz subsets)
        VAR4 = np.argwhere(out1[:,1]<len(kpc)) #[:,0]
        if config['DBG_SUBS_IN']:
            print("[subin.7] VAR4: %s"%VAR4)

        #separate the cols - for debugging (this is the subset variants)
        out2a = out1[:,0]
        if config['DBG_SUBS_IN']:
            print("[subin.8] out2a: %s"%out2a)

        #and this is their kpc count data
        out2b = out1[:,1]
        if config['DBG_SUBS_IN']:
            print("[subin.9] out2b: %s"%out2b)

        #the subset variants
        VAR5 = out2a[list(VAR4.T[0])]
        if config['DBG_SUBS_IN']:
            print("[subin.10] VAR5: %s"%VAR5)
            print("[subin.11] VAR5(names): %s"%self.sort.get_vname_by_vix(VAR5))

        #the subset variant counts
        VAR6 = out2b[list(VAR4.T[0])]
        if config['DBG_SUBS_IN']:
            print("[subin.12] VAR6: %s"%VAR6)

        #subsets shouldn't have diff intersecting positives with common supersets{{{
        subsc = list(VAR5)
        invalid_subs = []
        for v in subsc:
            v1 = Variant()
            v1.dbo = self.dbo
            v1.sort = self.sort
            v1.vix = v
            v1sups = v1.get_rel(relType=1)
            valid_sub = 1
            if len(v1sups):
                common_sups = set(v1sups).intersection(set(self.sups))
                v1kpc = self.sort.get_kixs_by_val(val=1,vix=v)
                for sup in common_sups:
                    v2 = Variant()
                    v2.dbo = self.dbo
                    v2.sort = self.sort
                    v2.vix = sup
                    v2kpc = self.sort.get_kixs_by_val(val=1,vix=sup)
                    common_kpc1 = set(self.kpc).intersection(set(v2kpc))
                    common_kpc2 = set(v1kpc).intersection(set(v2kpc))
                    for c in common_kpc2:
                        if c not in common_kpc1:
                            valid_sub = 0
                            invalid_subs.append(v)
                            break 
                    if valid_sub == 0 : break
        if len(invalid_subs):
            for v in invalid_subs:
                subsc.remove(v)
        search = np.array(invalid_subs)
        #(beg) technique found here: https://stackoverflow.com/questions/32191029/getting-the-indices-of-several-elements-in-a-numpy-array-at-once/32191125#32191125
        sorter = np.argsort(VAR5)
        idx = sorter[np.searchsorted(VAR5,invalid_subs,sorter=sorter)]
        #(end)#}}}
        VAR5a = np.delete(VAR5, idx)
        VAR6a = np.delete(VAR6, idx)

        #merged for return
        VAR7 = np.asarray((VAR5a,VAR6a)).T

        if config['DBG_SUBS_IN']:
            print("[subin.13] VAR7: %s"%VAR7)

        #if len(VAR6) == 0:
        #    return []

        #question is does it share other things in common with common supersets
        #of the target variant. if so, not subet.

        return VAR7 #[:,0]
        
    def get_supsets_or_eqv(self,kpc,eq_override=False):

        #-------------------------------------------------------------------------------------------------{{{
        # How supsets are determined:
        #-------------------------------------------------------------------------------------------------
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
        #-------------------------------------------------------------------------------------------------}}}

        if config['DBG_SUPS_IN']:
            print("")
            print("---------------------------------------------------------------------")
            print("")
            print("[supin.0] vix: %s"%self.vix)
            print("[supin.0] vix(name): %s"%self.sort.get_vname_by_vix(self.vix))
            print("[supin.0] kpc: %s"%kpc)
            print("[supin.0] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #looking for variants w/pos assignments that overlap the target variant
        VAR1 = np.argwhere(np.all(self.sort.NP[:,kpc]==[1]*len(kpc),axis=1)==True)[:,0]
        if config['DBG_SUPS_IN']:
            print("[supin.1] VAR1: %s"%VAR1)

        #get uniques of VAR1 with counts
        unique_elements, counts_elements = np.unique(VAR1, return_counts=True)
        VAR2 = np.asarray((unique_elements, counts_elements)).T #get uniques
        if config['DBG_SUBS_IN']:
            print("")
            print("[supin.2] VAR3: %s"%VAR2)

        #idx make sure we exclude the incoming variant
        idx = np.argwhere(VAR2[:,0] == self.vix)
        VAR4 = np.delete(VAR2[:,0], idx)
        if config['DBG_SUPS_IN']:
            print("[supin.4] VAR4: %s"%VAR4)

        #if there are no supsets, end here
        if len(VAR4) == 0: return []

        #get master idx of all positives for all data (+deal with equivalent variants)
        allPos = np.argwhere(self.sort.NP==1)[:,0]
        if config['DBG_SUPS_IN']:
            print("[supin.5] allPos: %s"%allPos)
        unqA, cntA = np.unique(allPos, return_counts=True)
        if len(VAR4):
            for idx, val in enumerate(unqA):
                try:
                    kpc_ = self.kpc
                except:
                    kpc_ = self.sort.get_kixs_by_val(val=1,vix=self.vix)
                #this logic drops equivalents
                if eq_override is False and val in VAR4 and cntA[idx] == len(kpc_):
                    idx1 = np.argwhere(VAR4==val)
                    VAR4 = np.delete(VAR4, idx1)
                #this logic drops everything but equivalents
                if eq_override is True and val in VAR4 and cntA[idx] != len(kpc_):
                    idx1 = np.argwhere(VAR4==val)
                    VAR4 = np.delete(VAR4, idx1)
        AP = np.asarray((unqA, cntA))[1,]
        if config['DBG_SUPS_IN']:
            print("[supin.6] AP: %s"%AP)

        #extrapolate the right supset mix using the master list idx
        VAR5 = AP[list(VAR4),]
        if config['DBG_SUPS_IN']:
            print("[supin.7] VAR5: %s"%VAR5)

        VAR6 = np.asarray((VAR4,VAR5)).T
        if config['DBG_SUPS_IN']:
            print("[supin.8] VAR6: %s"%VAR6)

        #sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
        VAR7 = VAR6[VAR6[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
        if config['DBG_SUPS_IN']:
            print("[supin.9] VAR7: %s"%VAR7)

        return VAR7

class Sort(object):

    def __init__(self):

        #db attributes
        self.dbo = None #db object

        #tree attributes    
        #self.TREE = {}
        #self.REF = None
        #self.MATRIX_MODE = config['MATRIX_MODE']
        #self.TREE_MODE = config['TREE_MODE'] #(sort tree presentation) 1=letters, 2=names, 3=letters+names 
        #self.TDATA = None

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
            sql = "INSERT into s_calls (kit_id,variant_loc,assigned) VALUES ('%s','%s',%s);" % (k,'top',1)
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
                        sql1 = "INSERT into s_calls (kit_id,variant_loc,assigned) VALUES ('%s','%s',%s);" % (k,vv,kv)
                        self.dbo.sql_exec(sql1)

    # argparser special routines

    def stdout_matrix(self,vix=None,refreshDbFlg=False):
        if refreshDbFlg:
            self.dbo.db = self.dbo.db_init()
            self.dbo.dc = self.dbo.cursor()
            self.get_mx_data(recreateFlg = False)
            #self.stdout_matrix()

        debug_chk('DBG_MATRIX',"",1)
        debug_chk('DBG_MATRIX',"matrix{{"+"{",1)
        debug_chk('DBG_MATRIX',"",1)

        #(beg)matrix
        table = BeautifulTable()
        table.column_headers = ['c']+['v']+self.get_cur_kit_list()
        table.append_row(['']+['']+[str(x) for x in list(range(10))])
        table.append_row(['','','','','','','','','','','',''])
        cntV = 0
        for K,V in self.get_axis('variants',idx=vix):
            table.append_row([cntV]+[K]+[str(x).replace('-1','') for x in self.get_mx_row_as_list(V[1])])
            table.row_seperator_char = ''
            table.column_seperator_char = ''
            table.column_alignments['v'] = BeautifulTable.ALIGN_LEFT
            cntV = cntV + 1
        debug_chk('DBG_MATRIX',table,1)
        #(end)matrix

        debug_chk('DBG_MATRIX',"",1)
        debug_chk('DBG_MATRIX',"}}"+"}",1)
        debug_chk('DBG_MATRIX',"",1)
        
    def stdout_unknowns(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.get_mx_data(recreateFlg = False)
        print("")
        print("---------------------------------------------------------------------")
        print("")
        cnt = 1
        for vix in self.get_imperfect_variants_idx():
            vt = Variant()
            vt.sort = self
            vt.dbo = self.dbo
            vt.set_info(vix)
            print("[%s]  var: %s" %(vt.vix,vt.name))
            print("     %s" %(vt.kucn))
            print("")
            cnt = cnt + 1
        print("---------------------------------------------------------------------")
        print("")

    # matrix calls

    def sort_matrix(self):

        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #get data
        self.get_mx_data()

        #step 0
        debug_chk('DBG_MATRIX',"data - step 0 (default)",1)
        self.stdout_matrix()

        #step 1
        debug_chk('DBG_MATRIX',"data - step 1",1)
        self.sort_step1()
        self.stdout_matrix()

        #step 2
        debug_chk('DBG_MATRIX',"data - step 2",1)
        self.sort_step2()
        self.stdout_matrix()

        #step 3
        debug_chk('DBG_MATRIX',"data - step 3",1)
        self.sort_step3()
        self.stdout_matrix()

        sys.exit()

    def sort_step1(self):
        self.mx_vertical_sort_new()
        self.save_mx_to_db()
        
    def sort_step2(self):
        self.mx_horizontal_sort()
        self.save_mx_to_db()
        
    def sort_step3(self):
        print("---------------------------------------------------------------------")
        print("Processing Imperfect Variants:")
        print("disabled. exiting.")
        sys.exit()
        for impVix in self.get_imperfect_variants_idx():
            print("---------------------------------------------------------------------")
            print("{{"+"{") #beg vim marker
            results = []
            print("}}"+"}") #end vim marker
            for R in results:
                print(R[1])

        print("----------------")
        print("")
        self.get_mx_count_data()
        self.mx_vertical_sort()
        self.mx_horizontal_sort()
        self.save_mx_to_db()

    def get_row_when_override_kixs(self,override_val,vix,kixs):
        row = self.get_mx_kdata(vix=vix)
        if len(kixs) == 0:
            return row
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        for kx in kixs: #.tolist():
            rowO[0,kx] = override_val
        return rowO

    def get_vname_by_vix(self,vix,listFlg=False): #get vname (can also take a list)
        #listFlg: force listFlg as return data type
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
                        variant = itm[0]
                    else:
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
            return list(np.argwhere(overrideData[0,] == val).T[1,]) #with override data, there's only one line evaluated - 1d datset
        if vix is not None: #no override -- use self.NP (all data)
            return list(np.argwhere(self.NP[vix,] == val).T[1,]) #default data, it's the entire matrix - 2d dataset 
        
    def get_knames_by_val(self,val,vix=None,vname=None,overrideData=None,toStr=True): #like get_mx_kdata but retrieves index info for given val
        if toStr:
            return l2s(self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData)))
        else:
            return self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData))
        
    def get_vixs_by_val(self,val,kix=None,kname=None,overrideData=None): #like get_mx_variant_data but retrieves index info for given val
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None and overrideData is not None:
            return np.argwhere(overrideData[:,0] == val).T[0,] #with override data, there's only one line evaluated - 1d dataset
        if kix is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[:,kix] == val).T[0,] #default data, it's the entire matrix - 2d dataset
        
    def get_vnames_by_val(self,val,kix=None,kname=None,overrideData=None,toStr=True): #like get_mx_variant_data but retrieves index info for given val
        if toStr:
            return l2s(self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData)))
        else:
            return self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData))

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

    def set_new_order(self,val,cnt,kitType=False,variantType=False):
        if kitType:
            self.KITS[val][1] = cnt
        if variantType:
            self.VARIANTS[val][1] = cnt
        
    def set_new_axis(self,vals,cnts,kitType=False,variantType=False):
        self.VARIANTS = {}
        for x in range(len(vals)):
            self.VARIANTS[vals[x]] = [0,cnts[x]]
        
    def save_mx_to_db(self):
        sql = "delete from s_mx_idxs;"
        #sql = sql+"delete from s_mx_kits;"
        #sql = sql+"delete from s_mx_calls;"
        #sql = sql+"delete from s_mx_idxs;"
        self.dbo.sql_exec(sql)
        itms = [(k,c2) for (n,(k,(c1,c2))) in enumerate(self.get_axis('variants'))]
        sql = "insert into s_mx_idxs (type_id,axis_id,idx) values (0,?,?);"
        self.dbo.sql_exec_many(sql,itms)
        itms = [[k,c2] for (n,(k,(c1,c2))) in enumerate(self.get_axis('kits'))]
        sql = "insert into s_mx_idxs (type_id,axis_id,idx) values (1,?,?);"
        self.dbo.sql_exec_many(sql,itms)

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
        
    def get_axis(self,orderByType=None,keysOnly=False,idx=None): # gets the variant/col or kit/row names (and optionally order info too)
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
        vnames = self.get_vname_by_vix(variant_idx_r)
        self.set_new_axis(vnames,variant_idx_r,variantType=True) #reset variant axis
        self.NP = variant_data #reset data
        
    def get_perfect_variants_idx(self):
        idx = list(range(len(self.VARIANTS)))
        prf_idx = idx[:] #copy idx so I have this one for deleting
        unk_idx = list(np.unique(np.argwhere(self.NP==0)[:,0]))
        for x in idx:
            if x in unk_idx:
                prf_idx.remove(x)
        return idx
        
    def get_imperfect_variants_idx(self):
        return list(np.unique(np.argwhere(self.NP==0)[:,0]))[:] #make it a copy
        
    def filter_perfect_variants(self,vix):
        if any(isinstance(el, list) for el in vix):
            for itm in reversed([[n,v] for (n,(v,c)) in enumerate(vix)]):
                if itm[1] in self.get_imperfect_variants_idx():
                    vix.remove(vix[itm[0]])
            return vix
        else:
            for itm in reversed([[n,v] for n,v in enumerate(vix)]):
                if itm[1] in self.get_imperfect_variants_idx():
                    vix.remove(vix[itm[0]])
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

    def get_mx_data(self,recreateFlg=True):

        #recreating from scratch
        if recreateFlg:
            sql = "select * from perfect_assignments_with_unk;"
        #going with the saved data
        else:
            sql = "select * from saved_assignments_with_unk;"

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
            #calls
            DATA[row[1]].append(row[2])
            #kits
            if row[0] not in self.KITS:
                self.KITS[row[0]] = [cntK,cntK]
                if recreateFlg:
                    sql = "insert into s_mx_kits (kit_id) values('%s');" % row[0]
                    self.dbo.sql_exec(sql)
                    sql = "insert into s_mx_idxs (type_id, axis_id, idx) values (%s,'%s',%s);" % (1,row[0],cntK)
                cntK = cntK + 1
            #variants
            if row[1] not in self.VARIANTS:
                self.VARIANTS[row[1]] = [cntV,cntV]
                if recreateFlg:
                    sql = "insert into s_mx_variants (variant_id,variant_loc,name) values('%s','%s','%s');" % (row[4],row[3],row[1])
                    self.dbo.sql_exec(sql)
                    sql = "insert into s_mx_idxs (type_id, axis_id, idx) values (%s,'%s',%s);" % (0,row[1],cntV)
                    self.dbo.sql_exec(sql)
                cntV = cntV + 1

        #create axes (saved version)
        #elif 1==2:
        #    sql = "select * from s_mx_idxs;"
        #    self.dbo.sql_exec(sql)
        #    F = self.dbo.fetchall()
        #    for row in F:
        #        if row[0] == 0: #variant
        #            self.VARIANTS[row[1]] = [row[2],row[2]]
        #        else: #kit
        #            self.KITS[row[1]] = [row[2],row[2]]
      
        #create numpy version of data
        for key,value in DATA.items():
            self.NP = np.matrix(list(DATA.values()))

        #print(self.VARIANTS)
        #print(self.KITS)
        #sys.exit()

        #create numpy version of data
        for key,value in DATA.items():
            self.NP = np.matrix(list(DATA.values()))

        #cpush this new stuff into saved/matrix tbls (recreation)
        if recreateFlg:
            sql = "insert into s_mx_calls (kit_id,variant_loc,assigned) select kit_id,variant_loc,assigned from perfect_assignments;"
            self.dbo.sql_exec(sql)

        #chk matrix (debugging)
        #self.stdout_matrix()
        #sys.exit()

        #get count data
        self.get_mx_count_data()
        
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

        #get count data
        self.get_mx_count_data()

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
        
