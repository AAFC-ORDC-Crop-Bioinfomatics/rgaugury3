#!/usr/bin/python
# -*- coding: utf-8 -*-

from os import path
from sqlite3 import connect, Error
from sys import exit

CTINE = 'CREATE TABLE IF NOT EXISTS '
createTbUsr = CTINE + 'tb_user (id INT, name TEXT)'
createTbPrj = CTINE + 'tb_project ' \
			  '(id INT, user_name TEXT, description TEXT, status TEXT, start_time INT, end_time INT, ' \
			  'GFF3 INT, E_value TEXT)'

def init_db():
    #print ("initializing database...")
    con = None

    try:
      con = connect(path.dirname(__file__)+'/rgaugury.db')
      cur = con.cursor()    
      cur.execute(createTbUsr)       
      cur.execute(createTbPrj)
      
    except Error, e:
      print ('Error:{}'.format(e.args[0]))
      exit(1)

    finally:
      if con:
        con.close()

    #print ("Finished initializing database.")