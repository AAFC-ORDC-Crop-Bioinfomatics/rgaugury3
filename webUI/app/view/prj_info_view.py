from app import app,db,models
from app.tool import render, zipdir
from config import FASTA_EXTENSION, PRJ_HOME,PROJECTS, GENE_FILE,FINAL_STATUS,SUMMARY_FILE,META_FILE
from flask import send_from_directory, jsonify, redirect, url_for
from os.path import isfile, isdir
from zipfile import ZipFile, ZIP_DEFLATED
import os

class td:
  class_name=''
  value=0
  def __init__(self, cn, value):
    self.class_name = cn
    self.value = value
  
@app.route('/prj_info/<prj_id>')
def prj_info(prj_id):
    gff3 = 'No'
    createZip(prj_id,prj_id+FASTA_EXTENSION)
    green = 'light_green_bg'
    pink = 'light_pink_bg'
    project_info = db.session.query(models.Project).get(prj_id)
    if not project_info:
       return redirect('/')
    if project_info.status != FINAL_STATUS:
      return redirect('/status')
    gff3 = project_info.gff3
    gene_amount = db.session.query(models.GeneAmount).get(project_info.gene_amount)
      
    tbs =[td(pink,gene_amount.nbs),td(pink,gene_amount.cnl),td(pink,gene_amount.tnl),\
          td(pink,gene_amount.cn),td(pink,gene_amount.tn),td(pink,gene_amount.nl),\
          td(pink,gene_amount.tx),td(pink,gene_amount.other),td(green,gene_amount.rlp),\
          td(pink,gene_amount.rlk), td(green,gene_amount.tmcc)]
    
    annotation_url =  url_for('index')[:-1]+PROJECTS+'/'+prj_id +'/'+prj_id+FASTA_EXTENSION+'.zip'
    return render('prj_info.html', prj_name=project_info.name, values = tbs, title='Results and Summaries',\
                  annotation_url=annotation_url, prj_id=prj_id, gff3 =gff3, root=url_for('index'))
  
@app.route(PROJECTS+'/<id>/<file>.zip')
def download_prj (id,file):
  directory = PRJ_HOME + '/' + id
  return send_from_directory(directory= directory, filename=file+'.zip')

@app.route('/ds_header/<prj_id>')
def ds_header(prj_id):
  path = PRJ_HOME + '/' + prj_id + '/' + prj_id + GENE_FILE
  with open(path,'r') as f:
    row = f.readline().split('\t')
  return jsonify(columns=row)  


@app.route('/gene_info/<prj_id>/<gene_type>')
def gene_info (prj_id, gene_type):
  gene_type = gene_type.strip()
  prj_id = prj_id.strip()
  genes = db.session.query(models.Gene).filter_by(prj_id = prj_id)
  genes = genes.filter(models.Gene.type.like(gene_type+'%'))

  return jsonify(data=[gene.to_list() for gene in genes])    
  
@app.route('/img/<prj_id>/<path:filename>')
def get_img(prj_id,filename):
  directory = PRJ_HOME + '/' + prj_id + '/'
  path = filename.rsplit('/',1)
  if len(path) == 2:
    directory += path[0]
    filename = path[1]
  return send_from_directory(directory,filename)

def addZip(file, zipf):
    if isfile(file):
        zipf.write(file)

def createZip(id, file):
  directory = PRJ_HOME + '/' + id
  zip_file = file +'.zip'
  if not isfile(zip_file):
    if isdir(directory):
      os.chdir(directory)
      with ZipFile(zip_file,'w', ZIP_DEFLATED) as zipf:
        if isdir('img/'):
          zipdir('img/',zipf)
        addZip(file,zipf)
        addZip(id+META_FILE,zipf)
        addZip(id+GENE_FILE,zipf)
        addZip(id+SUMMARY_FILE,zipf)
        addZip(id+'.CViT.all.txt',zipf)

@app.route('/distribution/<prj_id>')
def distribution(prj_id):
  return render('distribution.html', prj_id=prj_id)

@app.route('/dist_gallery/img/<prj_id>/img/<img>')
def dis_gallery(prj_id,img):
  img_path = '/img/'+prj_id+'/img/'+img
  return render('dist_gallery.html', prj_id=prj_id, img_path=img_path)


