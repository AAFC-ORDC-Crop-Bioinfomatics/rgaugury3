from app import app, db, models
from app.tool import render, kill_proc_tree
from config import PRJ_HOME, TOTAL_STEPS, FINAL_STATUS, DATE_FORMAT,CANCELED,PENDING
from datetime import datetime
from flask import request,url_for
from psutil import Process, pid_exists
from shutil import rmtree
from os.path import isdir, isfile
from bs4 import BeautifulSoup
from os import listdir

@app.route('/status')
def status():
  projects = db.session.query(models.Project).all()
  syc(projects)

  for prj in projects:  
    proj_folder=PRJ_HOME + '/' + prj.id
    if not isdir(proj_folder):
        cancel(prj.id)
        delete(prj.id)
        projects.remove(prj)

    elif prj.status != FINAL_STATUS and prj.status!=CANCELED:
        proj_folder=PRJ_HOME + '/' + prj.id
        log_file = proj_folder + '/' + prj.id + '.status.log'
        
        if isfile(log_file):
            step = 0
            step_content=''
            with open(log_file) as f:
                for line in f:
                    step_content += (line.strip() + '<br>')
                    step += 1

            #prj.status = 'Step {}'.format(step)
            prj.status = str(int(step)*100/TOTAL_STEPS)+'%'
            prj.step = step_content
            end_time = datetime.now()
            start_time = datetime.strptime(prj.start_time, DATE_FORMAT)
            delta = end_time - start_time
            list_time = str(delta).split(':')
            list_time[2] = str(int(round(float(list_time[2]))))
            prj.elapsed_time = ':'.join(list_time)
        
    db.session.commit()

  return render('status.html', projects=projects, title='Job List Status', root= url_for('index'))

    
@app.route('/cancel',methods=['POST'])
def cancel_prj():
    prj_id = request.form['prj_id'].strip()
    cancel(prj_id)
    return 'success'

## this function only deletes projects and records in the database
## cancel operation must be called first from user
@app.route('/delete',methods=['POST'])
def delete_prj():
    prj_id = request.form['prj_id'].strip()
    delete(prj_id)
    return 'success'

@app.route('/fingerprint/<prj_id>/<fingerprint>')
def checkFingerPrint(prj_id, fingerprint):
  prj_id = prj_id.strip()
  fingerprint = fingerprint.strip()
  prj = db.session.query(models.Project).get(prj_id)
  if fingerprint == str(prj.fingerprint):
    return '0'
  return '1'

@app.route('/all-prj')
def getPrjInfo():
    return getTable(status())

def delete(prj_id):
    # delete records in database
    prj = db.session.query(models.Project).get(prj_id)
    if prj:
        if prj.gene_amount:
            gene_amount = db.session.query(models.GeneAmount).get(prj.gene_amount)
            db.session.delete(gene_amount)
        db.session.delete(prj)
        db.session.commit()
    # delete project folder
    proj_folder=PRJ_HOME + '/' + prj.id
    if isdir(proj_folder):
        rmtree (proj_folder)

def cancel(prj_id):
    prj = db.session.query(models.Project).get(prj_id)

    if prj.status != CANCELED:
        pid = prj.pid
        if pid_exists(pid):
            process = Process(pid)
            # check if the current process is the one spawn by python by
            # comparing the create time.
            if process.create_time() == prj.pid_ctime:
                kill_proc_tree(pid)
        prj.status = 'canceled'
        db.session.commit()

def getTable(status_page):
    soup = BeautifulSoup(status_page, 'html.parser')
    return str(soup.table)

def syc(projects):
    db = []
    for prj in projects:
        db.append(prj.id)
    if not isdir(PRJ_HOME):
        return
    dirs = listdir(PRJ_HOME)
    for folder in dirs:
        if folder not in db:
            rmtree (PRJ_HOME + '/' + folder)
