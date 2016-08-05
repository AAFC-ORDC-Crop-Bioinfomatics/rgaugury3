from app import app
from flask import render_template,request, url_for,redirect
from config import BASE_PATH, APP_HOME, TEMPLATE
from os import path
from config import ADMIN_PASSWORD

pages = ['about', 'help']

@app.route('/admin/<page>')
def admin(page):
    if page in pages:
        return render_template('admin.html', title='Admin')
    else:
        return redirect('/')    

@app.route('/save/<page>',methods=['POST'])
def save(page):
    file = path.join(BASE_PATH,APP_HOME,TEMPLATE,page+'.html')
    with open(file, 'w') as f:
        f.write(request.form['data'].encode('utf8'))
    return 'success'
    

@app.route('/check_pwd/<pwd>')
def check_pwd(pwd):
    if pwd.strip() == ADMIN_PASSWORD:
        return 'success'
    else:
        return 'fail'

