#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from config import BASE_PATH, APP_HOME, PORT , DEBUG

sys.path.append(BASE_PATH+'/'+APP_HOME)

from app import app

if __name__ == "__main__":
    app.run(debug=DEBUG, host= '0.0.0.0', port=PORT)


