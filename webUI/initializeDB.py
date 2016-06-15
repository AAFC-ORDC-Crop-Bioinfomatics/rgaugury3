#!/usr/bin/python
# -*- coding: utf-8 -*-

from migrate.versioning import api
from config import SQLALCHEMY_DATABASE_URI, DB_PATH
from config import SQLALCHEMY_MIGRATE_REPO
from app import db
from os import path, remove
from shutil import rmtree

if path.exists(DB_PATH):
    remove(DB_PATH)
if path.exists(SQLALCHEMY_MIGRATE_REPO):
    rmtree(SQLALCHEMY_MIGRATE_REPO)

db.create_all()

if not path.exists(SQLALCHEMY_MIGRATE_REPO):
    api.create(SQLALCHEMY_MIGRATE_REPO, 'database repository')
    api.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
else:
    api.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, api.version(SQLALCHEMY_MIGRATE_REPO))
