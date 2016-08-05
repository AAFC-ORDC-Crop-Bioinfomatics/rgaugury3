#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from config import BASE_PATH, APP_HOME, PORT

sys.path.append(BASE_PATH+'/'+APP_HOME)

import app

app.app.config['WEB_UI_LOG'] = app.app.debug 

if __name__ == "__main__":
    app.app.run(host= '0.0.0.0', port=PORT)


