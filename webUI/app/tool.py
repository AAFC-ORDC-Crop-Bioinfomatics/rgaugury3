from flask import render_template, url_for
import re
import os
from os import makedirs
from time import strftime,localtime
from psutil import Process,wait_procs,pid_exists,cpu_count,cpu_percent
from os.path import isfile, isdir
from app import app, db, models
import zipfile
from config import MOTIF_PATH,BASE_PATH, PRJ_IMG_PATH, PRJ_HOME, GFF_FILE
from shutil import copy

def render(temp, **context):
    temp_name = temp.split('.')[0]
    js = url_for('static', filename='js/'+temp_name+'.js') 
    #if isfile(app.config['PROGRAM_HOME']+js):
    context['js_file'] = js
    css = url_for('static', filename='style/'+temp_name+'.css')
    #if isfile(app.config['PROGRAM_HOME']+css):
    context['css_file'] = css
    return render_template(temp, **context)

# count genen types in a project
def countType(fasta_path, gene_amount):
    types={'nbs':0, 'cnl':0, 'tnl':0, 'cn':0, 'tn':0,\
          'nl':0, 'tx':0, 'other':0,'rlp':0,'rlk':0,'tm-cc':0}
    total = 0
    with open(fasta_path,'r') as f:
        for line in f:
            if re.match(".*\|.*",line):
                genre = line.split('|')[1]
                genre = genre.rstrip(os.linesep).lower()
                types.setdefault(genre,0)
                types[genre] += 1
                total += 1
    
    gene_amount.nbs = types['nbs']
    gene_amount.cnl = types['cnl']
    gene_amount.tnl = types['tnl']
    gene_amount.cn  = types['cn']
    gene_amount.tn = types['tn']
    gene_amount.nl = types['nl']
    gene_amount.tx = types['tx']
    gene_amount.other = types['other']
    gene_amount.rlp = types['rlp']
    gene_amount.rlk = types['rlk']
    gene_amount.tmcc = types['tm-cc']
    gene_amount.total = total

def kill_proc_tree(pid, including_parent=True):
    parent = Process(pid)
    children = parent.children(recursive=True)
    for child in children:
        child.terminate()
    wait_procs(children)
    if including_parent:
        parent.terminate()

def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            subfile = os.path.join(root, file)
            ziph.write(subfile)

def initGeneImage(path):
    query = db.session.query(models.GeneImage)
    with open(path, 'r') as f:
        for line in f:
            list = line.split('\t')
            name = list[2]
            name = 'Image'+name.split('_')[1].zfill(3)+'.png'
            geneImage = query.filter_by(name = name).first()
            if geneImage is None:
                category = list[1]
                geneImage = models.GeneImage(name=name,category=category)
                db.session.add(geneImage)
                db.session.flush()
    db.session.commit()

def getGeneId(prj_id,geneName):
    query = db.session.query(models.Gene)
    genes = query.filter(models.Gene.prj_id == prj_id, models.Gene.name == geneName)
    return genes.first().id

def initGeneSection(path,prj_id):
    with open(path, 'r') as f:
        for line in f:
            list = line.split('\t')
            geneName = list[0]
            geneId = getGeneId(prj_id, geneName)
            category = list[1]
            span = list[3:]
            for s in span:
                if not s.isspace():
                    db.session.add(models.GeneSection(span =s, gene_id = geneId, category=category))
    db.session.commit()

def coutInputAmount(path):
    amount = 0
    with open(path, 'r') as f:
        for line in f:
            if re.match('^>.*',line):
                amount += 1
    return amount

def getValidCPU():
    l = cpu_percent(interval=2, percpu=True)
    count = 0

    for cpu in l:
        if cpu < 20.0:
            count +=1
    return count

def getCPU():
    count = getValidCPU()
    if count < 5:
        return 1
    else:
        return count/2

def getSequence(prj_id, gene):
    path = PRJ_HOME + '/' + prj_id +'/' +prj_id +'.RGA.candidates.fasta'
    with open(path, 'r') as f:
        found = False
        for line in f:
            if found:
                return line
            if re.search(gene, line):
                found = True

#    found = 0
#    result='Can not find the protein sequence'
#    with open(path, 'r') as f:
#        for line in f:
#            if re.search('>', line):
#                if line.strip() == '>'+gene:
#                    result = ''
#                    found = 1
#                elif found == 1:
#                    return result
#            elif found==1:
#                result += line.strip().replace(' ','').replace('\t','')
#    return result

def getGff(prj_id, gene):
    path = PRJ_HOME + '/' + prj_id + '/' + prj_id + GFF_FILE
    gffs = [] 
    with open(path) as f:
        for line in f:
            l = line.split('\t')
            if re.match(r'.*'+gene+'.*', l[8].strip()):
                gffs.append(models.Gff(seqname = l[0], source = l[1], feature =l[2],\
                       start = l[3], end = l[4], score = l[5], \
                       strand  = l[6], frame  = l[7], attribute  =gene))
    return gffs

def get_byte_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size