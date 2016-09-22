from app import app, db, models
from flask import request
from time import time, strptime
from config import PRJ_HOME, PERL, RGAUGURY_PL, CPU_TOGGLE, SAMPLE_GFF, WEB_UI_LOG
from config import DATE_FORMAT, SAMPLE_FASTA,START_PIPELINE,ENVIR, BASE_PATH
from os import path, makedirs, chdir
from subprocess import Popen, PIPE
from datetime import datetime
from app.tool import render, get_byte_size
import requests
import re
from psutil import Process
import os 
from random import randint
from shutil import copyfile
if WEB_UI_LOG == True:
    from weblog import logging

InterProScan_PAHT = ''
PANTHER_MIN_SIZE = 10000000900 # byte

## set up environment variables
for key, value in ENVIR.iteritems():
    if key =='PATH':
        InterProScan_PAHT = re.search(':([^:]*?interproscan[^:]*?):',value).group(1)
        os.environ[key] = os.environ[key] +':'+value
    else:
        os.environ[key] = value
    

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        if request.form:
 	        return processForm()
    else:
        cpu_toggle =''
        if CPU_TOGGLE == 0:
            cpu_toggle = 'hidden'
        property_file  = InterProScan_PAHT +'/interproscan.properties'
        panther = 'disabled'
        if path.exists(property_file):
            with open(property_file, 'r') as f:
                for line in f:
                     m = re.search('panther.models.dir\s*=\s*(.+)\s*', line)
                     if m:
                        panther_dir = m.group(1)
                        print panther_dir
            if panther_dir[0] != '/':
                panther_dir = InterProScan_PAHT + '/' +panther_dir
            if path.exists(panther_dir) and get_byte_size(panther_dir) > PANTHER_MIN_SIZE:
                panther = ''
                print "found panther"

        else:
            print "There is no "+property_file
            if WEB_UI_LOG == True:
                logging.error("There is no "+property_file)

        return render('index.html', cpu=cpu_toggle, panther=panther)

@app.route('/latestVersion')
def get_data():
    try:
        content = requests.get('https://www.ebi.ac.uk/interpro/interproscan.html').content
        m = re.search('>InterProScan (.*)<',content)
        latest = None
        now = None
        if m:
            latest = m.group(1)
        else:
            return '-1'

        proc = Popen('interproscan.sh',stdout=PIPE)

        out = proc.stdout.readline()
        
        Process(proc.pid).terminate()
        m = re.search('InterProScan-(.*)',out)
        if m:
            now = m.group(1)
            if now == str(latest):
                return '0'
            else:
                return latest
        else:
            return '-1'
    except:
        connection_error = 'connection error. There might be no Internet connection.'
        print connection_error
        return '-1'
    

def processForm():
    # gent project user input name
    proj_name=request.form.get('proj_name')
    fingerprint = request.form.get('fingerprint')

    if not proj_name:
    	proj_name = ""
    	  
    # generate a project
    # add finger print to avoid confliction between two projects submitted at the same time
    proj_id =fingerprint[-2:] + str(int(time()))

    # generate a project folder
    proj_folder=PRJ_HOME + '/' + proj_id
    if not path.exists(proj_folder):
    	makedirs(proj_folder)
    
    # check if there are sequences in the text area or a loaded file
    # get protein seq file name if loaded
    protein_seq_file = path.join(proj_folder, proj_id + ".fasta")
    
    protein_seq_str = request.form.get('protein_seq')

    if protein_seq_str:
        with open(protein_seq_file, "w") as fh:
            fh.write(protein_seq_str)
            fh.close()
    else:
        seq_f = request.files['seq_file']
        if seq_f:
            seq_f.save(protein_seq_file)

    # get gff3 file name if loaded
    gff3_file = path.join(proj_folder, proj_id + ".gff")
    gff3_f = request.files['gff3_file']
    if gff3_f:
    	gff3_f.save(gff3_file)
    elif request.form.get('gff3') == 'default': 
        copyfile(path.join(BASE_PATH, SAMPLE_GFF), gff3_file)
    else:
        gff3_file=''    

    # get E value	  
    evalue=request.form.get('ev')
    if not evalue:
    	evalue = "1e-5"

    # get CPU number used
    if CPU_TOGGLE == 1:	  
        cpu_num=request.form.get('cpu')

        if not cpu_num:
    	   cpu_num = '2'
    else:
        cpu_num = '-1'

    # get selected databases
    dbs = request.form.getlist('db_hidden')
    dbs_str = ",".join(dbs);
    
    # chnge the current folder to the newly generated project folder to run the pipeline
    chdir(proj_folder)
   
    
    start_time = datetime.fromtimestamp(time()).strftime(DATE_FORMAT)
    command = [START_PIPELINE, RGAUGURY_PL,
               '-p',protein_seq_file,
               '-d',dbs_str,
               '-c', cpu_num,
               '-gff', gff3_file,
               '-pfx',proj_id,
               '-e', evalue,
               '-pn',proj_name,
               '-st',start_time,
               '-fp',fingerprint]

    # Execute the pipeline program
    proc = Popen(command)
    psProcess = Process(proc.pid)

    project = models.Project(id=proj_id, pid = proc.pid, pid_ctime= psProcess.create_time())
    project.database = dbs_str
    
    db.session.add(project)
    db.session.commit()

    status=[proj_id, proj_name, start_time]
    lenlist=range(len(status))
    
    with app.app_context():
      pass
    return render('info.html', status=status, lenlist=lenlist, title='Info')