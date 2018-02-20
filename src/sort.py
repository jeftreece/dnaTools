# license + libs{{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

import sys,os,yaml,csv,json,h5py,numpy as np
from beautifultable import BeautifulTable
from collections import OrderedDict
import pandas as pd

#}}}

try:
    config = yaml.load(open(os.environ['REDUX_CONF_ZAK']))
except:
    print("Missing environment variable REDUX_CONF_ZAK. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

def l2s(lst):
    return ",".join(str(x) for x in lst)

class Variant(object):

    def __init__(self):
        self.dbo = None

    def proc(self,vname):
        self.proc_vname(vname)
        self.set_info(lev=2)
        self.proc_chk(allowImperfect=False)
        
    def info(self,vname):
        allowImperfect = config['allowImperfectWithVInfo']
        self.proc_vname(vname)
        self.set_info(lev=2,allowImperfect=allowImperfect)
        self.stdout_info()
        
    def proc_vname(self,vname):
        self.sort.restore_mx_data()
        if vname.isdigit() and int(vname)<=len(self.sort.VARIANTS):
            print("\n** Assuming, you're providing a vix ID...")
            vix = int(vname)
        elif vname.isdigit():
            vix = self.sort.get_vix_by_name(vname)
        else:
            vix = self.sort.get_vix_by_name(vname.upper())
        self.vix = vix
        
    def set_info(self,vix=None,lev=1,allowImperfect=False):

        lev = lev-1
        if vix is not None:
            self.vix = vix

        self.imperfect_vix = self.sort.get_imperfect_variants_idx()
        self.name = self.sort.get_vname_by_vix(self.vix)
        self.vixn = "%s - [%s]"%(self.name,self.vix)
        self.kpc = self.sort.get_kixs_by_val(val=1,vix=self.vix)
        self.knc = self.sort.get_kixs_by_val(val=-1,vix=self.vix)
        self.kuc = self.sort.get_kixs_by_val(val=0,vix=self.vix)
        self.kpcn = "kpc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kpc)),l2s(self.kpc))
        self.kncn = "knc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.knc)),l2s(self.knc))
        self.kucn = "kuc: %s [%s]"%(l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc))
        self.kpnc = sorted(self.kpc + self.knc)

        self.sups = self.get_rel(relType=1,allowImperfect=allowImperfect)
        #self.sups.remove(0) #remove 'top'
        self.subs = self.get_rel(relType=-1,allowImperfect=allowImperfect)
        self.eqv = self.get_rel(relType=0,allowImperfect=allowImperfect)
        self.subsn = "subs: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.subs)),l2s(self.subs))
        self.supsn = "sups: %s [%s]" %(l2s(self.sort.get_vname_by_vix(self.sups)),l2s(self.sups))
        self.eqvn = "eqv: %s [%s]"%(l2s(self.sort.get_vname_by_vix(vix=self.eqv)),l2s(self.eqv))

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
            print("")

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
            print("")
            print("---------------------------------------------------------------------")
            print("")
        
    def proc_chk(self,allowImperfect):

        print("")
        print("---------------------------------------------------------------------")
        print("")

        #starting unk list
        unkL = self.kuc
        unkD = {}
        for k in unkL:
            unkD[k] = [0,0,0] #0's are defaults
        othK = {}
        splitReq = False

        #1st arg of set : sup test1 result
        #1 - ambiguous promotion to positive
        #2 - hard promotion to positive
        #3 - hard negative
        
        #2nd arg of set : consistency test2 result
        #2+ - if it gets to 2+ ... it means this unk should be positive to appease other variants

        #3rd arg of set : consistency test2 result (split override)
        #1 - means neg (due to split)

        #look for sups
        if len(self.sups) == 0:
            print("vix: %s [%s]"%(self.sort.get_vname_by_vix(self.vix),self.vix))
            print("\nThis variant has no supsets, can't do anything with it!")
        else:

            print("vix: %s [%s]"%(self.sort.get_vname_by_vix(self.vix),self.vix))
            print("unks: %s [%s]" % (l2s(self.sort.get_kname_by_kix(self.kuc)),l2s(self.kuc)))
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
                        print(" - [1] vix unk %s [%s] is/are neg" % (l2s(self.sort.get_kname_by_kix(pN)),l2s(pN)))
                        for k in pN:
                            if k in unkL:
                                unkL.remove(k)
                            unkD[k][0] = 3 #hard failure
                    else:
                        print(" - [1] sup is truly sup to target variant, all unks open to promotion")
                else:
                    print(" - [1] sup is sub/eq/sup to target variant")
                #ambiguous promotions
                for k in unkL:
                    if unkD[k] != 3:
                        unkD[k][0] = 1 
                print(" - [1] remaining unk - ambig pos: %s [%s]" %(l2s(self.sort.get_kname_by_kix(unkL)),l2s(unkL)))
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
                        print("       (mtP) shared btw vix + v1sub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(at_common_kpc)),l2s(at_common_kpc)))
                        print("       (msP) shared btw sup + v1sub: %s [%s]"%(l2s(self.sort.get_kname_by_kix(as_common_kpc)),l2s(as_common_kpc)))
                        print("       (xP) msP-mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(diff_common_kpc)),l2s(diff_common_kpc)))
                        print("       (cP) common btw msP+mtP: %s [%s]"%(l2s(self.sort.get_kname_by_kix(common_kpc)),l2s(common_kpc)))
                        for k in diff_common_kpc:
                            if k in unkL:
                                unkD[k][1] = unkD[k][1] + 1
                                if unkD[k][1] > 1:
                                    print(" - [2] %s [%s] is a req positive" % (self.sort.get_kname_by_kix(k),k))
                        if len(common_kpc) > 0:
                            splitL.append(common_kpc)
                        
                #splits
                uniq_splits = [list(x) for x in set(tuple(x) for x in splitL)]
                if len(uniq_splits) > 0:
                    split_intersects = list(set.intersection(*map(set, uniq_splits)))
                    if (len(uniq_splits) > 1 and len(split_intersects) == 0):
                        print(" - [2] split required: btw %s" % uniq_splits)
                        print(" - [2] all unk to negative")
                        splitReq = True

                print("")

        #resolution vars
        pos = []
        pos_a = []
        neg = []
        s_ = []
        lenS = 1

        #resolutions
        if splitReq:
            for s in uniq_splits:
                if len(s) > 1 : lenS = len(s)
                else: s_.append(s[0])
            if len(uniq_splits)>1:
                print("splits: %s [%s]" % (l2s(self.sort.get_kname_by_kix(s_)),l2s(s_)))
            else:
                print("splits: %s" % l2s(uniq_splits))
            for k in self.kuc:
                unkD[k][2] = 1

        for k,v in unkD.items():
            if v[0] == 3 or v[2] == 1:
                neg.append(k)
            elif v[0] == 2 or v[1] > 1:
                pos.append(k)
            elif v[0] == 1:
                pos_a.append(k)
        if len(pos):
            print("pos: %s [%s]" % (l2s(self.sort.get_kname_by_kix(pos)),l2s(pos)))
        if len(pos_a):
            print("pos_a: %s [%s]" % (l2s(self.sort.get_kname_by_kix(pos_a)),l2s(pos_a)))
        if len(neg):
            print("neg: %s [%s]" %  (l2s(self.sort.get_kname_by_kix(neg)),l2s(neg)))
        
        print("")
        print("---------------------------------------------------------------------")
        print("")

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
            if config['DBG_RELS']: print("\n[rels.0] MODE (get_vix_rel) - ONE (override val w/vix)\n")
            overrideData = self.sort.get_row_when_override_coord(override_val,kix=kix,vix=self.vix)
            #pos conditions when override coord with a value
            kpc = self.sort.get_kixs_by_val(1,vix=self.vix,overrideData=overrideData)

        #method 2
        elif kpc is None:
            if config['DBG_RELS']: print("\n[rels.0] MODE (get_vix_rel) - TWO (only use vix)\n")
            #default pos conditions when just use default target vix criteria
            kpc = self.sort.get_kixs_by_val(1,vix=self.vix) #default pos conditions
    
        #method3
        #requires only sending in a kpc
        else:
            if config['DBG_RELS']: print("\n[rels.0] MODE (get_vix_rel) - THREE (sending in a kpc)\n")

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
        if config['DBG_RELS']: print("[rels.1] rels: %s" % rels)
        if config['DBG_RELS']:
            suffix=''
            if override_val is not None and kix is not None:
                suffix = 'P' if override_val == 1 else 'N'
            print("[rels.2] rels%s: %s kpc: %s" % (suffix,",".join([str(i) for i in rels]),kpc))

        return list(rels)
        
    def get_subsets(self,kpc):

        if config['DBG_SUBS_IN']:
            print("\n---------------------------------------------------------------------\n")
            print("[subin.0] vix: %s"%self.vix)
            print("[subin.0] vix(name): %s"%self.sort.get_vname_by_vix(self.vix))
            print("[subin.0] kpc: %s"%kpc)
            print("[subin.0] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #looking for variants w/pos assignments that overlap the target variant
        VAR1 = np.argwhere(self.sort.NP[:,kpc]==1)[:,0]
        if config['DBG_SUBS_IN']: print("[subin.1] VAR1: %s"%VAR1)

        #idx make sure we exclude the target variant
        idx = np.argwhere(VAR1==self.vix)
        VAR2 = np.delete(VAR1, idx)
        if config['DBG_SUBS_IN']: print("[subin.2] VAR2: %s"%VAR2)

        #print (set(x for l in VAR1.tolist() for x in l))
        #get uniques of VAR2 with counts
        unique_elements, counts_elements = np.unique(VAR2, return_counts=True)
        VAR3 = np.asarray((unique_elements, counts_elements)).T #get uniques
        if config['DBG_SUBS_IN']: print("\n[subin.3] VAR3: %s"%VAR3)

        #for the following "adding technique" -- we need to exclude unk situations for the comparison (special handling)
        #adding technique - got this idea here: https://stackoverflow.com/questions/30041286/sum-rows-where-value-equal-in-column
        unq, unq_inv = np.unique(VAR3[:,0], return_inverse=True)
        out = np.zeros((len(unq), VAR3.shape[1]), dtype=VAR3.dtype) #create empty array to put the added values
        out[:, 0] = unq #fill the first column
        np.add.at(out[:, 1:], unq_inv, VAR3[:, 1:])
        if config['DBG_SUBS_IN']: print("[subin.4] out: %s"%out)

        #sorting without fields - https://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
        out1 = out[out[:,1].argsort()[::-1]] #reverse sort (2nd col) -- to get the max ones first
        if config['DBG_SUBS_IN']:
            print("[subin.4a] out1: %s"%out1)
            print("[subin.5] kpc: %s"%kpc)
            print("[subin.6] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #subsets need to have less kpc than the target variant (cuz subsets)
        VAR4 = np.argwhere(out1[:,1]<len(kpc)) #[:,0]
        if config['DBG_SUBS_IN']: print("[subin.7] VAR4: %s"%VAR4)

        #separate the cols - for debugging (this is the subset variants)
        out2a = out1[:,0]
        if config['DBG_SUBS_IN']: print("[subin.8] out2a: %s"%out2a)

        #and this is their kpc count data
        out2b = out1[:,1]
        if config['DBG_SUBS_IN']: print("[subin.9] out2b: %s"%out2b)

        #the subset variants
        VAR5 = out2a[list(VAR4.T[0])]
        if config['DBG_SUBS_IN']:
            print("[subin.10] VAR5: %s"%VAR5)
            print("[subin.11] VAR5(names): %s"%self.sort.get_vname_by_vix(VAR5))

        #the subset variant counts
        VAR6 = out2b[list(VAR4.T[0])]
        if config['DBG_SUBS_IN']: print("[subin.12] VAR6: %s"%VAR6)

        #subsets shouldn't have diff intersecting positives with common supersets
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
        #(end)
        VAR5a = np.delete(VAR5, idx)
        VAR6a = np.delete(VAR6, idx)

        #merged for return
        VAR7 = np.asarray((VAR5a,VAR6a)).T

        if config['DBG_SUBS_IN']: print("[subin.13] VAR7: %s"%VAR7)

        return VAR7
        
    def get_supsets_or_eqv(self,kpc,eq_override=False):

        if len(kpc) == 0: #if there's no positive kpc to compare
            return [] #then return nothing

        if config['DBG_SUPS_IN']:
            print("\n---------------------------------------------------------------------\n")
            print("[supin.0] vix: %s"%self.vix)
            print("[supin.0] vix(name): %s"%self.sort.get_vname_by_vix(self.vix))
            print("[supin.0] kpc: %s"%kpc)
            print("[supin.0] kpc(names): %s"%self.sort.get_kname_by_kix(kpc))

        #look for variants that intersect target variant along its kpc (exclude target
        #variant from results) -- these are the superset candidates
        vixIntersects = np.argwhere(np.all(self.sort.NP[:,kpc]==[1]*len(kpc),axis=1)==True)[:,0]
        if len(vixIntersects) == 0: return [] #if there are no kpc intersects, end here
        if config['DBG_SUPS_IN']: print("[supin.1] vixIntersects: %s"%vixIntersects)

        #get the kpc for each of these superset candidates
        allPosVix = np.argwhere(self.sort.NP[vixIntersects]==1)[:,0]
        unqA, cntA = np.unique(allPosVix, return_counts=True)
        IC = np.asarray((vixIntersects,cntA)).T #IC = intersect counts
        if config['DBG_SUPS_IN']: print("[supin.2] intersectCnts: %s"%IC)

        #now check to see if these superset candidate vix have more/equal kpc to target variant
        for itm in reversed([(idx,val) for idx,val in enumerate(vixIntersects.tolist())]):
            #this logic drops equivalents (eq_override arg)
            if eq_override is False and itm[1] == IC[itm[0]][0] and IC[itm[0]][1] == len(kpc):
                IC = np.delete(IC,(itm[0]),axis=0)
            #this logic drops everything but equivalents
            elif eq_override is True and itm[1] == IC[itm[0]][0] and IC[itm[0]][1] != len(kpc):
                IC = np.delete(IC,(itm[0]),axis=0)

        #return filtered results ... includes counts (in case they're needed)
        if config['DBG_SUPS_IN']: print("[supin.3] intersectCnts: %s"%IC)
        return IC

class Sort(object):

    def __init__(self):
        self.dbo = None #db object
        self.KITS = None
        self.VARIANTS = None
        self.NP = None
        self.imperfect_vix = None

    # argparser special routines

    def stdout_matrix(self,vix=None,refreshDbFlg=False):
        if refreshDbFlg:
            self.dbo.db = self.dbo.db_init()
            self.dbo.dc = self.dbo.cursor()
            self.restore_mx_data()

        print("")
        print("---------------------------------------------------------------------")
        print("")


        #(beg)matrix
        table = BeautifulTable(max_width=95)
        table.column_headers = ['c']+['v']+[str(x) for x in self.get_axis('kits',keysOnly=True)]
        #table.append_row(['']+['']+[str(x) for x in list(range(10))])
        #table.append_row(['','','','','','','','','','','',''])
        cntV = 0
        for K,V in self.get_axis('variants',idx=vix):
            table.append_row([cntV]+[str(K)]+[str(x).replace('-1','') for x in self.get_mx_row_as_list(V[1])])
            table.row_seperator_char = ''
            table.column_seperator_char = ''
            table.column_alignments['v'] = BeautifulTable.ALIGN_LEFT
            cntV = cntV + 1
        print(table)
        #(end)matrix

        print("")
        print("---------------------------------------------------------------------")
        print("")

        
    def stdout_unknowns(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.restore_mx_data()
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

    # matrix 

    def sort_matrix(self):

        #db
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()

        #get data
        self.create_mx_data()

        #proc
        #self.stdout_matrix()
        self.mx_vsort()
        self.mx_hsort()
        self.stdout_matrix()
        #sys.exit()
        self.save_mx()
        sys.exit()

        #step x
        self.sort_step() # this can be expanded upon
        self.stdout_matrix()
        self.save_mx()
        sys.exit()
        
    def sort_step(self):
        print("Processing Imperfect Variants:")
        print("disabled. exiting.")
        sys.exit()
        for vix in self.get_imperfect_variants_idx():
            results = []
            for R in results:
                print(R[1])

        print("")
        #self.get_mx_count_data()
        self.mx_vsort()
        self.mx_hsort()

    def get_vid_by_vix(self,vix,listFlg=False):
        intFlg = True
        try:
            value = int(vix)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            vid = None
            #normal variants in the matrix
            for itm in list(self.VARIANTS.items()):
                if itm[1][1] == vix:
                    vid = itm[1][0]
                    break
            if listFlg:
                if vid == None:
                    return []
                else:
                    return [vid]
            else:
                return vid
        else: #assume it's a list/set/numpy array (whatever) > that I can cast to a list if need be
            vidList = []
            for vo in list(vix):
                #normal variants in the matrix
                for itm in list(self.VARIANTS.items()):
                    if itm[1][1] == vo:
                        vidList.append(itm[1][0])
                        break
            return(vidList)
        
    def get_kid_by_kix(self,kix,listFlg=False):
        intFlg = True
        try:
            value = int(kix)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            kid = None
            for itm in list(self.KITS.items()):
                if itm[1][1] == kix:
                    kid = itm[1][0]
                    break
            if listFlg:
                if kid == None:
                    return []
                else:
                    return [kid]
            else:
                return kid
        else: #assume it's a list/set/numpy array (whatever) > that I can cast to a list if need be
            kidList = []
            for ko in list(kix):
                for itm in list(self.KITS.items()):
                    if itm[1][1] == ko:
                        kidList.append(itm[1][0])
                        break
            return(kidList)
        
    def get_vname_by_vix(self,vix,listFlg=False):
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
            for vo in list(vix):
                #normal variants in the matrix
                for itm in list(self.VARIANTS.items()):
                    if itm[1][1] == vo:
                        variantList.append(itm[0])
                        break
            return(variantList)
        
    def get_kname_by_kix(self,kix,listFlg=False):
        intFlg = True
        try:
            value = int(kix)
        except:
            intFlg = False
        if intFlg: #typically, it's just the order number it's placed in the matrix
            kit = None
            for itm in list(self.KITS.items()):
                if itm[1][1] == kix:
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
        
    def get_vix_by_name(self,vname):
        return self.VARIANTS[vname][1]
        
    def get_kix_by_name(self,kname):
        return self.KITS[kname][1]
        
    def get_kid_by_name(self,kname):
        return self.KITS[kname][0]
        
    def get_vid_by_name(self,vname):
        return self.VARIANTS[vname][0]
        

    def get_kixs_by_val(self,val,vix=None,vname=None,overrideData=None):
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None and overrideData is not None: # we're sending in a custom evaluation
            return list(np.argwhere(overrideData[0,] == val).T[1,]) #with override data, there's only one line evaluated - 1d datset
        if vix is not None: #no override -- use self.NP (all data)
            return list(np.argwhere(self.NP[vix,] == val).T[1,]) #default data, it's the entire matrix - 2d dataset 
        
    def get_vixs_by_val(self,val,kix=None,kname=None,overrideData=None):
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None and overrideData is not None:
            return np.argwhere(overrideData[:,0] == val).T[0,] #with override data, there's only one line evaluated - 1d dataset
        if kix is not None: #no override -- use self.NP (all data)
            return np.argwhere(self.NP[:,kix] == val).T[0,] #default data, it's the entire matrix - 2d dataset
        
    def get_knames_by_val(self,val,vix=None,vname=None,overrideData=None,toStr=True):
        if toStr:
            return l2s(self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData)))
        else:
            return self.get_kname_by_kix(self.get_kix_by_val(val,vix,vname,overrideData))
        
    def get_vnames_by_val(self,val,kix=None,kname=None,overrideData=None,toStr=True):
        if toStr:
            return l2s(self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData)))
        else:
            return self.get_vname_by_vix(self.get_vixs_by_val(val,kix,kname,overrideData))

    def old_mx_vsort(self):

        PI = self.get_perfect_variants_idx()
        print(l2s(PI))
        print(self.NP[PI,])
        PIX = np.argwhere(self.NP[PI]==1)
        unqP, cntP = np.unique(PIX[:,0], return_counts=True)
        CNT = np.asarray((unqP,cntP)).T
        pdMRG = pd.merge(pd.DataFrame(CNT,columns=['A','B']),pd.DataFrame(PIX,columns=['A','C']), on='A').sort_values(['C','B','A'], ascending=[False, False, True])
        print(pdMRG)
        newO = pdMRG.as_matrix(columns=pdMRG.columns[0:1]).T
        _, idx = np.unique(newO, return_index=True)

        #{{{
        print("---")
        print(newO[:,idx])
        print("---")
        print(newO)
        print("---")
        print(newO[:,np.sort(idx)].tolist()[0])
        print("---")

        #print(inv)
        #print("---")
        #def return_counts(idx, inv):
        #    count = np.zeros(len(idx), np.int)
        #    np.add.at(count, inv, 1)
        #    return count
        #cnts = return_counts(idx,inv)
        #for itm in [(i,j) for i,j in enumerate(idx) if cnts[i] > 1]:
        #    print(itm)
        #print(np.argwhere(self.NP[prfIdx]==1))
        #}}}
        #self.NP[newO[:,np.sort(idx)].tolist()[0]],)

        allPL = newO[:,np.sort(idx)].tolist()[0]

        #prfPos = np.argwhere(self.NP[prfIdx]==1)[:,0]
        #unqP, cntP = np.unique(prfPos, return_counts=True)
        #allP = np.asarray((prfIdx,cntP))
        #allP = allP[:,np.argsort(-allP[1])]
        #allPL = allP.T[:,0]

        #imperfect variants
        impIdx = self.get_imperfect_variants_idx()
        impPos = np.argwhere(self.NP[impIdx]==1)[:,0]
        unqI, cntI = np.unique(impPos, return_counts=True)
        allI = np.asarray((impIdx,cntI))
        allI = allI[:,np.argsort(-allI[1])]
        allIL = allI.T[:,0]
        #names
        namesL = self.get_vname_by_vix(np.concatenate((allPL,allIL)))
        #new matrix
        self.NP = self.NP[np.concatenate((allPL,allIL))]
        #variants + kits
        for ix,v in enumerate(namesL):
            self.VARIANTS[v][1] = ix
        
    def mx_vsort(self):

        PI = self.get_perfect_variants_idx()
        print(l2s(PI))
        print(self.NP[PI,])
        PIX = np.argwhere(self.NP[PI]==1)
        unqP, cntP = np.unique(PIX[:,0], return_counts=True)
        CNT = np.asarray((unqP,cntP)).T
        pdMRG = pd.merge(pd.DataFrame(CNT,columns=['A','B']),pd.DataFrame(PIX,columns=['A','C']), on='A').sort_values(['C','B','A'], ascending=[False, False, True])
        print(pdMRG)
        newO = pdMRG.as_matrix(columns=pdMRG.columns[0:1]).T
        _, idx = np.unique(newO, return_index=True)

        #{{{
        print("---")
        print(newO[:,idx])
        print("---")
        print(newO)
        print("---")
        print(newO[:,np.sort(idx)].tolist()[0])
        print("---")

        #print(inv)
        #print("---")
        #def return_counts(idx, inv):
        #    count = np.zeros(len(idx), np.int)
        #    np.add.at(count, inv, 1)
        #    return count
        #cnts = return_counts(idx,inv)
        #for itm in [(i,j) for i,j in enumerate(idx) if cnts[i] > 1]:
        #    print(itm)
        #print(np.argwhere(self.NP[prfIdx]==1))
        #}}}
        #self.NP[newO[:,np.sort(idx)].tolist()[0]],)

        allPL = newO[:,np.sort(idx)].tolist()[0]

        #prfPos = np.argwhere(self.NP[prfIdx]==1)[:,0]
        #unqP, cntP = np.unique(prfPos, return_counts=True)
        #allP = np.asarray((prfIdx,cntP))
        #allP = allP[:,np.argsort(-allP[1])]
        #allPL = allP.T[:,0]

        #imperfect variants
        impIdx = self.get_imperfect_variants_idx()
        impPos = np.argwhere(self.NP[impIdx]==1)[:,0]
        unqI, cntI = np.unique(impPos, return_counts=True)
        allI = np.asarray((impIdx,cntI))
        allI = allI[:,np.argsort(-allI[1])]
        allIL = allI.T[:,0]
        #names
        namesL = self.get_vname_by_vix(np.concatenate((allPL,allIL)))
        #new matrix
        self.NP = self.NP[np.concatenate((allPL,allIL))]
        #variants + kits
        for ix,v in enumerate(namesL):
            self.VARIANTS[v][1] = ix
        
    def mx_hsort(self):
        allPosKix = np.argwhere(self.NP==1)[:,1]
        unqA, cntA = np.unique(allPosKix, return_counts=True)
        allA = np.asarray((unqA,cntA))
        allA = allA[:,np.argsort(-allA[1])]
        allL = allA.T[:,0]
        namesL = self.get_kname_by_kix(allL)
        self.NP = self.NP[:,allL]
        for ix,v in enumerate(namesL):
            self.KITS[v][1] = ix

    def get_axis(self,orderByType=None,keysOnly=False,idx=None):
        #Note: this is a useful function for being able to enumerate the VARIANTS/KITS dictionaries
        if orderByType in ['variants','kits']:
            if orderByType == 'variants' : SCH = self.VARIANTS
            if orderByType == 'kits' : SCH = self.KITS
            if keysOnly:
                return [i[0] for i in sorted(SCH.items(), key=lambda e: e[1][1])]
            else:
                return sorted(SCH.items(), key=lambda e: e[1][1])

    def get_mx_row_as_list(self,rownum,noneToStr=True):
        if noneToStr:
            return ['None' if v is None else v for v in self.NP[rownum,:].tolist()[0]]
        else:
            return self.NP[rownum,:].tolist()[0]
        
    def get_mx_kit_data(self,vix=None,vname=None):
        if vname is not None:
            vix = self.get_vix_by_name(vname)
        if vix is not None:
            return self.NP[vix,]
        
    def get_mx_variant_data(self,kix=None,kname=None):
        if kname is not None:
            kix = self.get_kix_by_name(kname)
        if kix is not None:
            return self.NP[:,kix].T

    def get_perfect_variants_idx(self):
        idx = list(range(len(self.VARIANTS)))
        prf_idx = idx[:]
        unk_idx = list(np.unique(np.argwhere(self.NP==0)[:,0]))
        for x in idx:
            if x in unk_idx:
                prf_idx.remove(x)
        return prf_idx
        
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

    def get_row_when_override_kixs(self,override_val,vix,kixs):
        row = self.get_mx_kit_data(vix=vix)
        if len(kixs) == 0:
            return row
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        for kx in kixs:
            rowO[0,kx] = override_val
        return rowO
        
    def get_row_when_override_coord(self,override_val,kix=None,vix=None,kname=None,vname=None):
        row = self.get_mx_kit_data(vix=vix)
        rowO = np.empty_like(row)
        rowO[:] = row #make duplicate copy - important!
        rowO[0,kix] = override_val
        return rowO

    # data 

    def sort_schema(self):
        self.dbo.db = self.dbo.db_init()
        self.dbo.dc = self.dbo.cursor()
        self.dbo.sql_exec_file('sort-schema.sql')
        
    def save_mx(self):
        #push kit/variant/numpy data into saved/matrix tbls + h5py file
        #deletes
        sql = "delete from x_mx_variants;"
        self.dbo.sql_exec(sql)
        sql = "delete from x_mx_kits;"
        self.dbo.sql_exec(sql)
        sql = "delete from x_mx_idxs;"
        self.dbo.sql_exec(sql)
        #save matrix variants
        sql = "insert into x_mx_variants (ID,name) values (?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([vid,nm])) for (n,(nm,(vid,idx))) in enumerate(self.get_axis('variants'))])
        #save matrix variants (idx order)
        sql = "insert into x_mx_idxs (type_id,axis_id,idx) values (0,?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([vid,idx])) for (n,(nm,(vid,idx))) in enumerate(self.get_axis('variants'))])
        #save matrix kits
        sql = "insert into x_mx_kits (ID,kitId) values (?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([kid,nm])) for (n,(nm,(kid,idx))) in enumerate(self.get_axis('kits'))])
        #save matrix kits (idx order)
        sql = "insert into x_mx_idxs (type_id,axis_id,idx) values (1,?,?);"
        self.dbo.sql_exec_many(sql,[(tuple([kid,idx])) for (n,(nm,(kid,idx))) in enumerate(self.get_axis('kits'))])
        #save numpy data
        #https://stackoverflow.com/questions/20928136/input-and-output-numpy-arrays-to-h5py
        h5f = h5py.File('data.h5', 'w')
        h5f.create_dataset('dataset_1', data=self.NP)
        h5f.close()
        
    def restore_mx_data(self):
        self.VARIANTS = {}
        self.KITS = {}
        #variants
        sql = '''
            SELECT V.name, V.ID, IX.idx
            FROM x_mx_idxs IX, x_mx_variants V
            WHERE IX.type_id = 0 AND
            IX.axis_id = V.ID
            ORDER BY IX.idx;
            ''';
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        for row in F:
            self.VARIANTS[row[0]] = (row[1],row[2]) #self.VARIANTS[name] = [vID,idx]
        #kits
        sql = '''
            SELECT K.kitId, K.ID, IX.idx
            FROM x_mx_idxs IX, x_mx_kits K
            WHERE IX.type_id = 1 AND
            IX.axis_id = K.ID
            ORDER BY IX.idx;
            ''';
        self.dbo.sql_exec(sql)
        F = self.dbo.fetchall()
        for row in F:
            self.KITS[row[0]] = (row[1],row[2]) #self.VARIANTS[name] = [vID,idx]
        #numpy data
        h5f = h5py.File('data.h5','r')
        self.NP = np.asmatrix(h5f['dataset_1'][:])
        h5f.close()
        
    def create_mx_data(self):

        #sql to fetch data
        sql = "select * from x_perfect_assignments_with_unk;"
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
            if row[5] not in self.KITS.keys():
                self.KITS[row[5]] = [row[0],cntK] # self.KITS[name] = [pID,idx]
                cntK = cntK + 1
            #variants
            if row[1] not in self.VARIANTS.keys():
                self.VARIANTS[row[1]] = [row[4],cntV] #self.VARIANTS[name] = [vID,idx]
                cntV = cntV + 1

        self.NP = np.matrix(list(DATA.values()))

        #create variables that we'll need for finding duplicates among the variants
        b = np.ascontiguousarray(self.NP).view(np.dtype((np.void, self.NP.dtype.itemsize * self.NP.shape[1])))
        _, idx, inv = np.unique(b, return_index=True, return_inverse=True)

        def return_counts(idx, inv):
            count = np.zeros(len(idx), np.int)
            np.add.at(count, inv, 1)
            return count
        cnts = return_counts(idx,inv)

        print("\nThere are %s dupes: moving them" % len([(i,j) for i,j in enumerate(idx) if cnts[i] > 1]))

        #idx_uniq = [i for i, j in enumerate(idx) if cnts[i] == 1]
        #idx_dupe = [i for i, j in enumerate(idx) if cnts[i] > 1]
        #print("idx: %s"%idx)
        #print("idx_uniq: %s"%idx_uniq)
        #print("idx_dupe: %s"%idx_dupe)
        #print("inv: %s"%inv)

        #https://stackoverflow.com/questions/34123122/finding-indices-of-non-unique-elements-in-numpy-array
        #find and separate dupe variants
        for itm in [(i,j) for i,j in enumerate(idx) if cnts[i] > 1]:
            #this gets the duplicate variants (based on idx order) 
            itms = list(np.delete(np.argwhere(inv==itm[0]),0))
            #this puts these duplicate variants into the dupes tbl next to the one in their series we're keeping
            sql = "insert into x_mx_dupe_variants(vID,dupe_vID) values (%s,?);" % self.get_vid_by_vix(itm[1])
            self.dbo.sql_exec_many(sql,[tuple([l]) for l in self.get_vid_by_vix(itms)])
            #remove dupe variants from the self.VARIANTS var
            for k in self.get_vname_by_vix(itms):
                self.VARIANTS.pop(k,None)

        #reset the matrix index for the remaining non-dupe variants
        for ix,vn in enumerate(self.get_vname_by_vix(idx)):
            self.VARIANTS[vn][1] = ix
            
        #reset the matrix (now that the dupes are out)
        self.NP = self.NP[idx,]

        self.stdout_matrix()

        #save matrix data
        #self.save_mx()

