from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
from os.path import abspath,dirname

app = Flask(__name__)
app.config['PROGRAM_HOME'] = dirname(abspath(__file__))
app.config.from_object('config')

db = SQLAlchemy(app)

from app import views, models

