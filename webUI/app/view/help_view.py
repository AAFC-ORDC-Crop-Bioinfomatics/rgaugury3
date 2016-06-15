from app import app
from config import WIKI,USER_NAME, PASSWORD
from app.tool import render
import requests
from bs4 import BeautifulSoup


WIKI_CONTENT = '''
<h3>The help information is synced with bitbucket wiki page.</h3>
<h3>If you are a website administrator, you can add any bitbucket username and password in the config.py file.</h3>
<h3>Or</h3> 
<h3>You can click <a href='https://bitbucket.org/yaanlpc/rgaugury/wiki/Web%20UI%20Help'>Here</a> to view the help page. (You still need a bitbucket account)</h3>

'''
@app.route('/help')
def help():
    content = WIKI_CONTENT
    try:
        wiki_help = requests.get(WIKI, auth=(USER_NAME, PASSWORD)).content
        content = getSoup(wiki_help)
    except:
        pass
    return render('help.html', title='Help', content=content)

def getSoup(content):
    soup = BeautifulSoup(content, 'html.parser')
    return soup('section')[1]