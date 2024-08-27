# makedb.py
# aim: allow users to upload their own data to the database

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import re
import os
import sys
import pymongo

def createdb(db_name):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    db = client[db_name]
    return db

def createcollection(collection_name):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    collection = client[collection_name]
    return collection

def insertdoc(collection, doc):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    collection.insert_one(doc)

def insertdata(db, data):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    db.insert_one(data)

def get_db(*args):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    db = client[args]
    return db

def get_collection(*args):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    collection = client[args]
    return collection

def get_data(collection):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    data = client[collection].find()
    return data

def get_doc(collection, doc):
    from pymongo import MongoClient
    client = MongoClient("mongodb://localhost:27017/")
    doc = client[collection].find_one(doc)
    return doc

## bulk insert problem
## check database & what are in the database