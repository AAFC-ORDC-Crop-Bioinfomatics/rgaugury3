#!/usr/bin/python
import sys
from subprocess import Popen, PIPE
import argparse
from app import db, models
import app
from psutil import Process
from config import PRJ_HOME, META_FILE,FASTA_EXTENSION,GENE_FILE, LOG_FILE
from config import TOTAL_STEPS,DATE_FORMAT, FINAL_STATUS,PENDING
from app.tool import countType, initGeneImage, initGeneSection,coutInputAmount, getCPU
from datetime import datetime
from os.path import exists
import re

if app.app.config['WEB_UI_LOG']:
    from weblog import logging

command = sys.argv[1:14]

if command[8]=='':
	del command[7]
	del command[7]
	gff3 = 'No'
else:
	gff3 = 'Yes'


parser = argparse.ArgumentParser()
parser.add_argument('-pfx', dest="proj_id")
parser.add_argument('-pn', dest="proj_name")
parser.add_argument('-gff', dest="gff3_file") 
parser.add_argument('-e', dest="evalue")
parser.add_argument('-st', dest="start_time")
parser.add_argument('-fp', dest="fingerprint")
parser.add_argument('-c', dest="cpu_num")
parser.add_argument('-p')
args, unknown = parser.parse_known_args()

# add project info into database
prj_id = args.proj_id
project = db.session.query(models.Project).get(prj_id)
project.name=args.proj_name
project.gff3=gff3
project.e_value=args.evalue
project.start_time = args.start_time
gene_amount = models.GeneAmount(total=0)
db.session.add(gene_amount)
db.session.flush()

project.gene_amount=gene_amount.id
project.fingerprint=args.fingerprint
project.status = PENDING
path = PRJ_HOME + '/' + prj_id + '/' + prj_id + '.fasta'
project.input_amount = coutInputAmount(path)
db.session.commit()

if args.cpu_num == '-1':
    command[6] = str(getCPU())
    
print command
if app.app.config['WEB_UI_LOG']:
    logging.info(command)

o,e=Popen(command,stdout=PIPE,stderr=PIPE).communicate()

if app.app.config['WEB_UI_LOG']:
    logging.info(o)
    logging.error(e)

'''
Task has been finished
'''
'''
initialize database
'''
def initDB():
  ### .RGA.info.txt
  path = PRJ_HOME + '/' + prj_id + '/' + prj_id + GENE_FILE
  count = 0
  if exists(path):
    with open(path,'r') as f:
      for line in f:
        # skip head line
        if count != 0:
          if re.match('\s+',line):
            continue
          row = line.split('\t')
          db.session.add(models.Gene(prj_id = prj_id, name=row[0], length=row[1], type=row[2], image=row[3]))
        count = 1
    db.session.commit()
  
  ### .RGA.gene.meta.txt
  if gff3 == 'Yes':
    path = '/'.join([PRJ_HOME,prj_id,prj_id+META_FILE])
    if exists(path):
      initGeneImage(path)
      initGeneSection(path,prj_id)
  
  annotation_url = '/' +prj_id +'/'+prj_id+FASTA_EXTENSION
  path = PRJ_HOME + annotation_url
  if exists(path):
    countType(path, gene_amount)
    db.session.commit()
  
  ### project final status
  prj = db.session.query(models.Project).get(prj_id)
  proj_folder=PRJ_HOME + '/' + prj.id
  log_file = proj_folder + '/' + prj.id + LOG_FILE
  if exists(log_file):
    step = 0
    with open(log_file) as f:
        for line in f:
     		   step += 1
    if step == TOTAL_STEPS:
      start_time = datetime.strptime(prj.start_time, DATE_FORMAT)
      prj.status = FINAL_STATUS
      end_time = datetime.strptime(line.split(' : ')[0], DATE_FORMAT)
      
      delta = end_time - start_time
      list_time = str(delta).split(':')
      list_time[2] = str(int(round(float(list_time[2]))))
      prj.elapsed_time = ':'.join(list_time)
      
      db.session.commit()

  ### gff3 file
  #path = PRJ_HOME + '/' + prj_id + '/' + prj_id + '.gff'
  #if exists(path):
  #  query = db.session.query(models.Gff)
  #  with open(path) as f:
  #      for line in f:
  #          l = line.split('\t')
  #          db.session.merge(models.Gff(seqname = l[0], source = l[1], feature =l[2],\
  #                     start = l[3], end = l[4], score = l[5], \
  #                     strand  = l[6], frame  = l[7], attribute  =l[8].strip()))
  #  db.session.commit()

initDB()