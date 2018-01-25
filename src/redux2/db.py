# authors/licensing{{{

# @author: Iain McDonald
# Contributors: Jef Treece, Harald Alvestrand
# Purpose: Reduction and comparison script for Y-chromosome NGS test data
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007)
# https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import sys,os,sqlite3,yaml,time,csv,json,numpy as np
from anytree import Node, RenderTree
#import string

# }}}

REDUX_CONF = os.environ['REDUX_CONF']
config = yaml.load(open(REDUX_CONF))
start_time = 0 # need to fix this

#TODO: need to put this info properly into yaml file (for now, I'm hacking bashrc)
REDUX_ENV = os.environ['REDUX_ENV']
REDUX_SQL = os.environ['REDUX_SQL']
REDUX_DATA = os.environ['REDUX_DATA']

class DB(object):
    
    def __init__(self):
        self.db = None
        self.dc = None
        
    def db_init(self):
        #trace (1, "Initialising database...")
        return sqlite3.connect('variant.db')
        
    def cursor(self):
        return self.db.cursor()
        
    def run_sql_file(self,FILE):
        fh = open(REDUX_SQL+'/'+FILE,'r');
        try:
            #print(fh.read())
            self.dc.executescript(fh.read())
        finally:
            fh.close()
    def commit(self):
        self.db.commit()

    #redux2 ddl+dml

    def redux2_schema(self):
        self.run_sql_file('redux2-schema.sql')
        
    def insert_v1_variants(self,variant_array):
        self.dc.executemany('''INSERT INTO v1_variants(id,ref,alt) VALUES (?,?,?)''', variant_array)
        self.commit()
        
    def insert_v1_calls(self):
        self.dc.execute('''INSERT INTO v1_calls(variant,person)
            SELECT id, person
            FROM v1_variants CROSS JOIN v1_people''')
        self.commit()
        
    def insert_v1_hg19(self,snp_reference):
        self.dc.executemany("INSERT INTO v1_hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)",
            ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
        self.commit()
            
    def insert_v1_hg38(self,snp_reference):
        self.dc.executemany("INSERT INTO v1_hg38(grch38,grch38end,name,anc,der) VALUES (?,?,?,?,?)",
            ((rec[3], rec[4], rec[8], rec[10], rec[11]) for rec in snp_reference))
        self.commit()

    #v2 db schema ddl+dml

    def v2_schema(self):
        self.run_sql_file('schema-v2.sql')

    #clades db schema ddl+dml

    def clades_schema(self):
        self.run_sql_file('clades-schema.sql')

    #tree sort prototype ddl+dml
    
    def sort_schema(self):
        self.run_sql_file('sort-schema.sql')
        
    def insert_sample_sort_data(self):

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

        cols=10
        for k in range(1,cols+1):
            self.dc.execute("insert into s_kits (kit_id) values ("+str(k)+");")

        with open(REDUX_DATA+'/sample-sort-data.csv','r') as FILE:
            for row in csv.DictReader(FILE,'v n k1 k2 k3 k4 k5 k6 k7 k8 k9 k10'.split()):
                row = json.loads(json.dumps(row).replace('\\ufeff','')) #hack: remove byte order mark
                self.dc.execute("insert into s_variants (variant_loc,name) values ('"+row['v']+"','"+row['n']+"');")
                #print(' - inserting sample variant data: '+str(row['v']))
                for k in range(1,cols+1):
                    kv = str(row['k'+str(k)])
                    #'null' if kv == "None" else kv
                    vv = str(row['v'])
                    #print (kv)
                    self.dc.execute("insert into s_calls (kit_id,variant_loc,assigned) values ('k"+str(k)+"','"+vv+"',"+kv+");")
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
        self.commit()
        
    def sort_data(self):

        print("===")
        print("FILTER, step: A ")
        print("===")
        sql = "select distinct variant_loc from s_calls;"
        self.dc.execute(sql)
        A = [itm[0] for itm in self.dc.fetchall()]
        print("A - distinct variants")
        print(A) #A - distinct variants
        print("---")
        print('numA - num distinct variants')
        print(len(A)) #numA - num distinct variants
        #sys.exit()

        print("===")
        print("FILTER - step: B")
        print("===")
        sql = "select distinct variant_loc from s_calls where assigned = 0;"
        self.dc.execute(sql)
        B1 = [itm[0] for itm in self.dc.fetchall()]
        B0 = list(set(A)-set(B1))
        print("B0 - variants that don't have negs")
        print(B0) #B0 - variants that don't have negs
        print("---")
        print("B1 - variants that have negs")
        print(B1) #B1 - variants that have negs 
        #sys.exit()

        print("===")
        print("FILTER - step: C")
        print("===")
        sql = "select variant_loc,count(*) as cnt from s_calls where assigned = 1 group by variant_loc;"
        self.dc.execute(sql)
        F = self.dc.fetchall()
        Fa = list(filter(lambda x: x[1]==(len(A)-1), F))
        Fb = [(a) for a,b in Fa] #strip out 2nd element, the count
        C1 = list(set(B1) & set(Fb)) #intersection
        C0 = list(set(B1)-set(C1))
        print("list of *all* one person +ve's")
        print(Fa)
        print("---")
        print("not singletons")
        print(C0) #C0 - not singletons
        print("---")
        print("singletons") #C1 - singletons
        print(C1)
        #sys.exit()

        print("===")
        print("FILTER - step: D")
        print("===")
        sql = "select distinct variant_loc from s_calls where assigned is null group by variant_loc;"
        self.dc.execute(sql)
        F = self.dc.fetchall()
        D0 = list(set(C1)-set(F))
        D1 = list(set(F)-set(D0))
        print("list of variants that are sometimes not called")
        print(F)
        print("---")
        print("imperfect variants")
        print(D0) #D0 - imperfect variants
        print("---")
        print("calls of perfect share variants - these go through the next PROCESS, SORT")
        print(D1) #D1 - perfect share variants
        #sys.exit()

        #NOTE: a type study {{{
        #------------------------------

        #I'm thinking byK is what we're looking to work with

        #byV = [
        #    { v1:
        #        ({pos:[k4,k3]},{cnt:2},{neg:[]},{unk:[]})
        #        },
        #    { v2:
        #        ({pos:[k1,k2,k3]},{cnt:3},{neg:[]},{unk:[]})
        #        },
        #    { v3:
        #        ({pos:[k4,k3]},{cnt:2},{neg:[]},{unk:[]})
        #        },
        #    { v4:
        #        ({pos:[k2,k1]},(cnt:2},{neg:[]},{unk:[]})
        #        }
        #    ]

        #byK = [
        #    { k1:
        #        ({pos:[v2,v4]},{cnt:2},{neg:[]},{unk:[]})
        #        },
        #    { k2:
        #        ({pos,[v2,v4]},{cnt:2},{neg:[]},{unk:[]})
        #        },
        #    { k3:
        #        ({pos,[v1,v2,v3]},{cnt:2},{neg:[]},{unk:[]})
        #        },
        #    { k4:
        #        ({pos,[v1,v3]},{cnt:2},{neg:[]},{unk:[]})
        #        }
        #    ]
            
        #------------------------------ }}}

        print("===")
        print("SORT")
        print("===")
        sql = "select kit_id,variant_loc from s_calls order by kit_id, variant_loc;"
        self.dc.execute(sql)
        F = self.dc.fetchall()
        print("all kits + variants")
        print("+ grabbing kits by variant: F")
        #[(1, 3019783), (1, 6920349), (1, 7378685), (1, 8928037), ... ]
        Fa = list(filter(lambda x: x[1]=="A", F))
        Fb = [(a) for a,b in Fa] #strip out 2nd element, the count
        #[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        print(Fb)
        print("---")
        
        print("===")
        print("SETS")
        print("===")

        #all kits, variant, assignment mixes 
        sql = "select kit_id,variant_loc,assigned from s_calls order by kit_id, variant_loc,assigned;"
        self.dc.execute(sql)
        F = self.dc.fetchall()

        #all unique kits
        print("all unique kits")
        KITS = sorted(list(set([itm[0] for itm in F])))
        print(KITS)

        #all unique variants
        VARIANTS = sorted(list(set([itm[1] for itm in F])))
        VARIANTSp = ['+'+itm[1] for itm in F]
        VARIANTSn = ['-'+itm[1] for itm in F]
        VARIANTSa = sorted(list(set(VARIANTSp+VARIANTSn)))
        print("all unique variants")
        print(VARIANTS)
        print("all unique variants - pos+neg")
        print(VARIANTSa)
        #sys.exit()
        
        #start tree
        VA = {}
        top = Node("top")
        for v in VARIANTS:
            VA[v] = {}
            VA[v]['pos'] = Node('+'+v, parent=top)
            VA[v]['neg'] = Node('-'+v, parent=VA[v]['pos'])
        for pre, fill, node in RenderTree(top):
            print("%s%s" % (pre, node.name))

        ##loop kits looking for relations btw variants

        #kits with positive assignments
        Fp = sorted(list(set([i[0] for i in list(filter(lambda x: x[2]==1, F))])))
        print("kits with positive assignments")
        print(Fp)

        #kits with negative assignments
        Fn = sorted(list(set([i[0] for i in list(filter(lambda x: x[2]==0, F))])))
        print("kits with negative assignments")
        print(Fn)
        #sys.exit()

        print("---")
        KA={}

        #per all the kits with positive variants (build new dict)
        for k in Fp:
            Kp = sorted(list(set(['+'+i[1] for i in list(filter(lambda x: x[0]==k and x[2]==1, F))])))
            #['A+', 'D+', 'F+', 'H+', 'M+']
            print(Kp)
            #sys.exit()
            KA[k] = {'len':len(Kp),'plen':len(Kp),'sort':0,'variants':Kp}

        #per all the kits with negative variants (build new dict)
        for k in Fn:
            Kn = sorted(list(set(['-'+i[1] for i in list(filter(lambda x: x[0]==k and x[2]==0, F))])))
            if k in KA.keys():
                KA[k]['len'] = len(KA[k]['variants'])+len(Kn)
                KA[k]['variants'] = sorted(KA[k]['variants']+Kn)
            else:
                KA[k] = {'len':len(Kn),'plen':0,'sort':0,'variants':Kn}

        #loop dict's variant sets with unique variant types to find relations
        #(step1) A+ {B:(B+,B-)} :: [(A+,B+),(A+,B-)]
        #if x[0] == A+ && x[1] == B+ exists and if x[0] == A+ && x[1] == B- exists:
        #    A+>B-

        #--- so this:
        #    A+>A-
        #    B+>B-
        #--- becomes:
        #    R1+>A-
        #    R1+>B-
        #    unk-relations:
        #    R1:{A1+,B+}

        #(step2) A+ {C:(C+,C-)} :: [(A+,C+),(A+,C-)]
        #if x[0] == A+ && x[1] == C+ exists and if x[0] == A+ && x[1] == C- exists:
        #    A+>C-

        #--- so this:
        #    R1+>A-
        #    R1+>B-
        #    unk-relations:
        #    R1:{A1+,B+}
        #--- becomes:
        #    R1+>A-
        #    R1+>B-
        #    R1+>C-
        #    unk-relations:
        #    R1:{A+,B+,C+}

        #(step3a) B- {C:(C+,C-)} :: [(B-,C+),(B+,C-)]
        #if x[0] == B- && x[1] == C+ exists and if x[0] == B- && x[1] == C- exists:
        #    B+>C-

        #--- so this:
        #    R1+>A-
        #    R1+>B-
        #    R1+>C-
        #    unk-relations:
        #    R1:{A+,B+,C+}
        #--- becomes:
        #    R1+>A-
        #    R1+>B-
        #    B->C-
        #    unk-relations:
        #    R1:{A+,B+,C+}

        #for Vx in VARIANTSa:
        #    for key, value in KA.items():
        #        if key=='variants':
        #            for Vy in value:
                       


        #loop dict to create list
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

        print("---")
        newV3 = {}
        for d in newV2:
            newV3[d['kit']] = d['variants']
            STR = d['kit']+':'+str(d['variants'])
            print(STR.replace("'",""))

        #blocks
        print("---")
        blocks = {}
        #print("blocks")
        for VX in VARIANTS:
            blocks[VX] = {'mix':[],'pos':[],'neg':[]}
            VXP = '+'+VX
            #print(VXP)
            for VY in VARIANTS:
                VYP = '+'+VY
                if VXP == VYP:
                    foo=1
                    #print("-VXP:"+VXP)
                    #print("-VYP:"+VYP)
                    #print("here")
                    #sys.exit()
                else:
                    VYN = '-'+VY
                    #print(VYP)
                    #print(VYN)
                    chk1 = False
                    chk2 = False
                    for d in newV2:
                        if VXP in d['variants']:
                            #print("VXP:"+VXP)
                            #print("VYP:"+VYP)
                            #print("VYN:"+VYN)
                            #print("variants:"+str(d['variants']))
                            if chk1 is False and VYP in d['variants']:
                                #print("here1")
                                chk1 = True
                            if chk2 is False and VYN in d['variants']:
                                #print("here2")
                                chk2 = True
                            #print("finish checks")
                        if chk1 is True and chk2 is True:
                            blocks[VX]['mix'].append(VY)
                            break
                            #return
                    if chk1 is True and chk2 is False:
                        blocks[VX]['pos'].append(VY)
                    if chk2 is True and chk1 is False:
                        blocks[VX]['neg'].append(VY)
                    
        for key, value in blocks.items():
            print(key+'|'+str(value))
            
        #---

        # sample data

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

        # processed raw results:

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

        # rules:

        # mix: (1|2) means 1 is above 2
        # pos: (1|2) means 1 is a direct ancestor or direct descendant or dupe, not a "cousin", "uncle", or 
        #      "sibling" - to 2. (differ to dupe in cases where there's no other clues)
        # neg: (1|2) means 1 is a "cousin" or "uncle" or "sibling" or direct ancestor of 2

        # results when applying these rules manually:

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

        # disputes: 
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

        #---

        sys.exit()

        #and print a json version nof it (debugging)
        print(json.dumps(newV2, indent=4, sort_keys=True))
        sys.exit()

        #{{{ 

        #sql_2b = "select variant_loc,count(*) as pos_v_cnt from s_calls where assigned = 0 group by variant_loc order by count(*) desc;"
        #self.dc.execute(sql_2b)
        #varAn = self.dc.fetchall()
        #print("---")
        #print("variant negative check")
        #print(varAn)

        sql_2b = "select variant_loc,count(*) as pos_v_cnt from s_calls where assigned = 0 group by variant_loc order by count(*) desc;"
        sql_2c = "select variant_loc,count(*) as pos_v_cnt from s_calls where assigned is not null group by variant_loc order by count(*) desc;"
        self.dc.execute(sql_2c)
        varAa = self.dc.fetchall()
        print("---")
        print("variant all check")
        print(varAa)
        #(3) 9 perfectly called variants - execute sort on these
        #(4) 6 imperfectly called variants - do Step A

        sql_3 = "select * from s_calls order by kit_id,assigned;"
        self.dc.execute(sql_3)
        callsA = self.dc.fetchall()
        print("---")
        #[(1, 12060401, None), (1, 6920349, 0), (1, 7378685, 0), (1, 13668461, 0), (1, 19538924, 0), ... ]
        print (callsA);

        #Note: build the default structure with the kits ordered like kitA and the variants ordered like varA
        #[{"k1":[{"v12060401",1)},{"v6920349",1), ... ]}
        #[{"k2":[{"v12060401",1)},{"v6920349",None), ... ]}
        #and display it
        #for call in calls:
        #   for K in kits:
        #    ...
        #   sort_positive_variants(kit_id)

        #}}}

