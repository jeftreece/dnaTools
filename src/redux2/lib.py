# license {{{

# Purpose: Y-DNA NGS analytics
# Git repo: https://github.com/jazdrv/dnaTools
# For free distribution under the terms of the GNU General Public License,
# version 3 (29 June 2007) https://www.gnu.org/licenses/gpl.html

# }}}
# libs {{{

import os,yaml,shutil,glob,re,csv,zipfile,subprocess
from db import *
from sort import *
from collections import defaultdict

# }}}

#debugging {{{

def trace (level, msg):
    print(msg)
    #if level <= config['verbosity']:
    #    print(msg)
    #TODO: below line in clades.py
    #sys.stderr(flush)
    
def debug_chk(TYPE,msg):
    if config[TYPE]:
        print(msg)

#}}}
# conf {{{

try:
    config = yaml.load(open(os.environ['REDUX_CONF']))
except:
    trace(0,"Missing environment variable REDUX_CONF. Aborting.")
    sys.exit()
sys.path.append(config['REDUX_PATH'])

#}}}

# user defined file mapping {{{

names = """
FTDNA345238Newell.zip, 345238, Newell
155941_BigY_RawData_20140911-1.zip, 155941, Unknown
Lee 237414 BigY Raw Data.zip, 237414, Lee
U106_515653_Hogenmiller_BigY_RawData_2016_11_20.zip, 515653, Hogenmiller
bigy-Bettinger57020.zip, 57020, Bettinger
"""

#}}}
# rename dict{{{

rename_dict = {}
for row in csv.reader(names.splitlines()):
    if row and row[0]:
        rename_dict[row[0].strip()] = (row[1].strip(), row[2].strip())

#}}}

# redux2 {{{

# file/dir

def refresh_dir(DIR,cleanFlag=False):
    DIR = config['REDUX_ENV']+'/'+DIR
    #print DIR
    if (os.path.isdir(DIR)):
        files = glob.glob(DIR+'/*')
        if cleanFlag:
            for f in files:
                os.remove(f)
    else:
        os.makedirs(DIR)
    
def delete_file(FILE):
    FILE = config['REDUX_ENV']+'/'+FILE
    if os.path.exists(FILE):
        os.remove(FILE)
    
def touch_file(FILE):
    FILE = config['REDUX_ENV']+'/'+FILE
    if not os.path.exists('merge-ignore.txt'):
        open('merge-ignore.txt','w').close()
    
def cmd_exists(CMD):
    return any(os.access(os.path.join(path, CMD), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

# file/dir - J/H

def setup_dirs():
    shutil.rmtree(config['unzip_dir'],ignore_errors=True)
    os.makedirs(config['unzip_dir'])
    
def extract_zips():

    if not os.path.isdir(config['REDUX_ENV']+'/'+config['zip_dir']):
        trace (0, '   Warn: no directory with zip files: %s' % config['zip_dir'])
        return []

    FILES = os.listdir(config['REDUX_ENV']+'/'+config['zip_dir'])

    # try to parse out at least the kit number by trying a series of regular expressions
    # adding regular expressions at the end of this list is safer than at the beginning
    # order is important - rules at top are matched first

    # constants used in filename regular expressions
    # groupings (?:xxx) are ignored

    ws = r'^[_]?'
    nam1 = r"[a-z]{0,20}|O\&#39;[a-z]{3,20}|O['][a-z]{3,20}"
    cname = r'([\w]{1,20})' #matches unicode chars; also matches digits though
    pnam = r'\('+nam1+r'\)'
    nam2 = r'(?:' +nam1 +'|' +pnam +r')'
    ndate = r'(?:(201[1-8][\d]{4}|201[1-8]-\d\d-\d\d|\d{4}201[1-8]))'
    sep = r'[\-\s\._]'
    seps = r'[\-\s\._]?'
    # sepp = r'[\-\s\._]+'
    sepp = r'_' # only use underscore as field separator
    sept = r'[\-\s\._]{3}'
    bigy = r'(?:big' +seps+ r'y(?:data)?|ydna)'
    rslt = r'(?:results|data|rawdata|vcfdata|raw data|csvexport|raw_data|raw|bigyrawdata)'
    name = r'((?:'+nam2+seps+'){1,3})'
    kit = r'(?:(?:kit|ftdna)?[ #]?)?([enhb1-9][0-9]{3,6})'
    rzip = r'zip(?:.zip)?'
    plac = r'([A-Z]{2})'

    #0 e.g. bigy-Treece-N4826.zip
    #1 e.g. N4826_Treece_US_BigY_RawData_2018-01-03.zip
    #2 e.g. 548872_Lindstrom_Germany_BigY_RawData_2018-01-01.zip

    name_re = [
        (re.compile(ws+sep.join([bigy,name,kit,rzip]), re.I), 'name', 'kit'),
        (re.compile(ws +sepp.join([kit,name,plac,bigy,rslt,ndate])+'.zip', re.I), 'kit', 'name'),
        (re.compile(ws +sepp.join([kit,cname,plac,bigy,rslt,ndate])+'.zip', re.I), 'kit', 'name')
        ]


    trace (25, '   File names mapped, according to which regular expression:')
    # track counts - only for diagnostics
    cnt = defaultdict(int)
    # list of non-matching files
    nomatch=[]
    # all of the file names we could parse
    fname_dict = {}
    
    for line in FILES:
        fname = line.strip()
        if fname in rename_dict:
            # hand-edited filename mappings
            kkit, nname = rename_dict[fname]
            fname_dict[fname] = kkit, nname
            trace(25, '     {3:>2} {0:<50s}{1:<15s}{2:<10s}'.format(fname, nname, kkit, 'd'))
            cnt['d'] += 1
        else:
            if fname[-4:] not in ('.gitignore'):
                if fname[-4:] not in ('.zip'):
                    trace (15, '   Found foreigner hanging out in zip directory: {0}'.format(fname))
                continue
            d = {}
            for ii, (r,k1,k2) in enumerate(name_re):
                s = r.search(line)
                if s:
                    d[k1] = s.groups()[0]
                    if k2:
                        d[k2] = s.groups()[1]
                    else:
                        d['name'] = 'Unknown'
                    try:
                        trace (25, '     {3:>2} {0:<50s}{1:<15s}{2:<10s}'.format(fname,
                                                   d['name'], d['kit'], ii))
                        cnt[ii] += 1
                        fname_dict[fname] = d['kit'], d['name']
                    except:
                        trace (1, '   FAILURE on filename:', fname)
                    break
            else:
                if line not in ('.gitignore'):
                    nomatch.append(line)

    trace (20, '   Number of filenames not matched: {0}'.format(len(nomatch)))
    trace (22, '   Which expressions were matched:')
    for nn,cc in cnt.items():
        trace (22, '     {0:>2}: {1:>4}'.format(nn,cc))

    if len(nomatch) > 0:
        trace (10, '   Files that did not match:')
        for ll in nomatch:
            if ll.strip() not in ('.gitignore'):
                trace (10, '    %s' % ll.strip())
            else:
                nomatch = nomatch - 1

    # keep track of what needs to be cleaned up
    emptydirs = []

    for fname in fname_dict:
        kitnumber, kitname = fname_dict[fname]
        try:
            zf = zipfile.ZipFile(config['REDUX_ENV']+'/'+config['zip_dir']+'/'+fname)
            #zf = zipfile.ZipFile(os.path.join(config['zip_dir'],fname))

        except:
            trace (1, '   ERROR: file %s is not a zip' % fname)
            sys.exit()
        listfiles = zf.namelist()
        bedfile = vcffile = None
        for ff in listfiles:
            dirname, basename = os.path.split(ff)
            if basename == 'regions.bed':
                bedfile = ff
            elif basename == 'variants.vcf':
                vcffile = ff
            if dirname and (dirname not in emptydirs):
                emptydirs.append(dirname)
        if (not bedfile) or (not vcffile):
            trace(1, '   Warn: missing data in '+fname)
            continue
        if (bedfile == None) ^ (vcffile == None):
            trace(1, '   Warn: BED or VCF file is missing for %s' % fname)
            trace(1, '   This is an unexpected error. %s not processed.' % fname)
            continue
        zf.extractall(config['unzip_dir'], [bedfile, vcffile])
        base = '%s-%s' % (kitname, kitnumber)
        try:
            fpath = os.path.join(config['unzip_dir'], '%s')
            trace (40, "      "+fpath % base)
            os.rename(fpath % bedfile, (fpath % base)+'.bed')
            os.rename(fpath % vcffile, (fpath % base)+'.vcf')
        except:
            trace(1, '   Warn: could not identify VCF and/or BED file for '+base)

    # clean up any empty dirs unzip created

    if emptydirs:
        trace (30, '   Trying to remove droppings:')
        for dir in emptydirs:
            try:
                dp = os.path.join(config['unzip_dir'], dir)
                os.removedirs(dp)
                trace (30, '     {0}'.format(dp))
            except FileNotFoundError:
                pass
            except:
                trace (30, '     W! could not remove {0}'.format(dp))
                pass

    # list of file names we unzipped

    files = os.listdir(config['REDUX_ENV']+'/'+config['unzip_dir'])
    return files
    
def unpack():

    # messy problem - messy solution - kit names not consistent - Harald/Jef
    # collect run time statistics

    trace(10,'   Running the unpack-zip script...')
    #setup_dirs()
    refresh_dir(config['unzip_dir'],cleanFlag=False)
    fnames = extract_zips()
    trace (10, '   Number of files: {0}'.format(len(fnames)))
    trace (40, '   Files unpacked:')
    for ff in fnames:
        trace (40, ff)
def skip_to_Hg19(dbo):

    # skip to <= 1 - unpack zips

    if (config['skip_to'] <= 1):
        trace (2, "Unpacking ZIP files...")
        #unpack(REDUX_ENV+'/'+config['zip_dir'],REDUX_ENV+'/'+config['unzip_dir'],config['verbosity'])
        unpack()
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
        trace (5, "Associating unzipped files with kits...")
        
    # skip to <= 10 - associate kits with people

    if (config['skip_to'] <= 10):
        trace (2, "Associating kits with people...")
        
    # skip to <= 11 - generate dictionary of variant positions 

    #   The VCF files are parsed twice. The first pass identifies the list
    #   of variants to be queried. The second pass reads the calls for those
    #   variants. This means we only need to treat positions with a variant,
    #   rather than every position in a chromosome.
    #   Using a dictionary here means only one copy of each variant is
    #   created.

    if (config['skip_to'] <= 11):

        #vcffiles

        trace (2, "Generating database of all variants...")
        vcffiles = [f for f in os.listdir(config['REDUX_ENV']+'/'+config['unzip_dir']) if f.endswith('.vcf')]
        trace (10, "   %i files detected" % len(vcffiles))
        
        #variant_dict

        #print(REDUX_ENV)
        variant_dict = {}
        for file in vcffiles:
            vcf_calls = readHg19Vcf(config['REDUX_ENV']+'/'+config['unzip_dir']+'/'+ file)
            variant_dict.update(vcf_calls)

        trace (10, "   %i variants found" % len(variant_dict))
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)

        # dump variant dict into sorted array

        trace (20, "   Dumping variants into array...")
        variant_array = np.array(list(variant_dict.values()))
        print(variant_array)
        sys.exit()

        # variant_array = np.array([],dtype={'names': ('start', 'anc', 'der'),'formats': ('i4', 'S20', 'S20')})

        trace (30, "      Check variant [0] is %s" % variant_array[0])
        trace (30, "      Check variant [0] position is %s" % variant_array[0][1])
        trace (30, "      Check variant [%s] is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1]))
        trace (30, "      Check variant [%s] position is %s" % (len(variant_dict)-1, variant_array[len(variant_dict)-1][1]))

        #db calls

        trace (20, "   Inserting data into variant array database...")
        dbo.insert_v1_variants(variant_array)
        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
        
    # skip to <= 12 - reading calls for variants
    
    #db calls
    
    if (config['skip_to'] <= 12):
        trace (2, "Generating database of calls...")
        vcffiles = [f for f in os.listdir(config['REDUX_ENV']+'/'+config['unzip_dir']) if f.endswith('.vcf')]
        trace (10, "   %i files detected" % len(vcffiles))
        dbo.insert_v1_calls()

    # skip to <= 13 - name variants and derive ancestral values

    # Some variants are positive in the reference sequence, so we need to
    # look up their ancestral values. We'll get the SNP names while we're
    # at it.

    #db calls

    if (config['skip_to'] <= 13):
        # Read in SNPs from reference lists
        trace (2, "Getting names of variants...")
        trace (10, "   Importing SNP reference lists...")
            
        snp_reference = csv.reader(open(config['REDUX_ENV']+'/'+config['b37_snp_file']))
        #for rec in snp_reference:
        #   print "INSERT INTO v1_hg19(grch37,grch37end,name,anc,der) VALUES (?,?,?,?,?)", (rec[3], rec[4], rec[8], rec[10], rec[11])
        dbo.insert_v1_hg19(snp_reference)
            
        snp_reference = csv.reader(open(config['REDUX_ENV']+'/'+config['b38_snp_file']))
        dbo.insert_v1_hg38(snp_reference)

        # db work - how we doing? {{{

        # Read in SNPs from reference lists
        # Probably doesn't need to be done at this point
        # trace (10, "   Joining reference lists to variant database...")

        # self.dc.execute('''SELECT hg38.grch38, hg38.name
        # FROM hg38
        # INNER JOIN hg19 on hg19.name = hg38.name''')

        # self.dc.execute('''SELECT variants.id, hg38.name
        # FROM variants
        # LEFT OUTER JOIN hg38 on hg38.grch38 = variants.id''')

        # }}}

        t = float((time.clock() - start_time))
        trace(10, '   ...complete after %.3f seconds' % t)
            
    #commit

    dbo.commit()

    # Print final message and exit {{{{

    t = float((time.clock() - start_time))
    trace (1, "Execution finished in: %.3f seconds" % t)

    # }}}

# import vcf hg38 - it looks like clades.py is doing something like this for hg19

def getH38references():
    foo = 1

# }}}
# sample code {{{

def getVCFvariants(FILE):
    cmd = config['REDUX_ENV']+"/getVCFvariants.sh"
    p = subprocess.Popen(cmd, FILE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()

# }}}
# clades {{{

# SNP extraction routines based on original - Harald 
# extracts the SNP calls from the VCF files and
# determines the coverage of SNPs in the BED files of BigY tests.

def analyzeVcf(file):

    #Returns a dict of position -> mutation mappings

    with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
        trace (30, "   Extracting VCF: %s" % vcffile)
        result = {}
        for line in vcffile:
            fields = line.split()
            if (fields[0] == 'chrY' and fields[6] == 'PASS' and fields[3] != '.' and fields[4] != '.'):
                # fix by Jef Treece for fields containing commas:
                result[int(fields[1])] = fields[1] + '.' + fields[3].replace(',', ';') + '.' + fields[4].replace(',', ';')
                # result[int(fields[1])] = fields[1] + '.' + fields[3] + '.' + fields[4]
        return result
    
def analyzeBed(file):

    #Returns an array of path segments.

    with open(os.path.splitext(file)[0] + '.bed') as bedfile:
        trace (30, "   Extracting BED: %s" % bedfile)
        result = []
        for line in bedfile:
            fields = line.split()
            if (fields[0] == 'chrY'):
                result.append((int(fields[1]), int(fields[2])))
        return result
    
def makeCall(pos, index_container, bed_calls):

    #Figure out whether this position is on a segment boundary.
    #Between segments = 'nc'; top of segment = 'cbu'; bottom of segment = 'cbl'.
    #Only call in a single-position segment = 'cblu'.
    #index_container contains first segment to be looked at.
    #This function must only be called for increasing values of pos, and with
    #sorted bed_calls.

    call = ';nc'
    for bed_index in range(index_container[0], len(bed_calls)):
        pos_pair = bed_calls[bed_index]
        index_container[0] = bed_index
        if pos_pair[1] >= pos:
            # Position is before or within this segment.
            if pos_pair[0] <= pos:
                # Position is within this segment.
                if pos_pair[0] == pos_pair[1] and pos_pair[0] == pos:
                    call = ';cblu'
                elif pos_pair[0] == pos:
                    call = ';cbl'
            elif pos_pair[1] == pos:
                call = ';cbu'
            else:
                call = ''
        else:
            # Position is before this segment.
            call = ';nc'
        return call
        # If position is after segment, continue.
    return ';nc' # After end of last segment.
    
def extract(unzip_dir,files,variants):

    d = []
    s = []

    curpath = os.path.abspath(os.curdir)
    with open(os.path.join(curpath, 'variant-list.txt')) as line_headings:
        for line in line_headings:
            d.append(line.rstrip())
            x = line.split(',')
            s.append(int(x[0]))  # s holds the genome position for each line

    for file in files:
        vcf_calls = analyzeVcf(config['unzip_dir'] + file)
        bed_calls = analyzeBed(config['unzip_dir'] + file)
        bed_index = [0]
        for lineno in range(len(d)):
            d[lineno] += ','
            if s[lineno] in vcf_calls:
                d[lineno] += vcf_calls[s[lineno]]
            d[lineno] += makeCall(s[lineno], bed_index, bed_calls)

        for line in d:
            print (line)
    
def file_len(fname):

    #File length, thanks to StackOverflow
    #https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python

    i=-1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# routines - Iain 

def readHg19Vcf(file):

    #Returns a dict of position -> mutation mappings
    #Modified from Harald's analyzeVCF, this version returns every mutation with
    #its derived value, regardless of whether it was ancestral or not

    with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
        trace (30, "   Extracting VCF: %s" % vcffile)
        result = {}
        for line in vcffile:
            fields = line.split()
            if (fields[0] == 'chrY' and int(fields[1]) > 0 and fields[3] != '.' and fields[4] != '.'):
                result[fields[1]] = [int(fields[1]), str(fields[3]), str(fields[4])]
        return result

#}}}

