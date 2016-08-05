import sys
from os import path

sys.path.insert(0, path.abspath(path.dirname(__file__)))

from app import app as application
