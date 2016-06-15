from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
from os.path import abspath,dirname
from config import ENVIR
import os 

app = Flask(__name__)
app.config['PROGRAM_HOME'] = dirname(abspath(__file__))
app.config.from_object('config')
db = SQLAlchemy(app)

## set up environment variables
for key, value in ENVIR.iteritems():
    if os.environ.get(key, None):
        os.environ[key] = os.environ[key] +':'+value
    else:
        os.environ[key] = value

from app import views, models

