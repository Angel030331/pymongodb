import argparse
import sys
import subprocess
from utils.makedb import createdb, createcollection, insertdata, insertdoc, get_db, get_collection, get_data, get_doc
from utils.meta import Genbank_record, Gb_parser
from utils import meta

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import re
import json
import pickle

def main():
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("--db_name", help="The name of the database")
    parser.add_argument("--collection_name", help="The name of the collection")
    parser.add_argument("--input_file", help="The name of the input file")
    parser.add_argument("--list_db", help="List the databases")
    parser.add_argument("--list_collections", help="List the collections")
    parser.add_argument("--search", help="Search the database")
    parser.add_argument("--genbank", help="The name of the genbank file")
    args = parser.parse_args()
    
    db_name = args.db_name
    collection_name = args.collection_name
    input_file = args.input_file
    list_db = args.list_db
    list_collections = args.list_collections
    search = args.search
    genbank = args.genbank

    
    # create the database
    db = createdb(db_name)
    
    # create the collection
    collection = createcollection(collection_name)
    
    # insert the data
    with open(input_file, 'r') as f:
        if input_file:
            if input_file.endswith('.fasta'):
                data = SeqIO.parse(f, 'fasta')
                for record in data:
                    insertdata(db, collection, record)
                print("data inserted to the database")
            elif input_file.endswith('json'):
                data = json.load(f)
                insertdata(db, collection, data)
                print("data inserted to the database")
            elif input_file.endswith('pkl'):
                data = pickle.load(f)
                insertdata(db, collection, data)
                print("data inserted to the database")
            elif input_file.endswith('gb'):
                data = SeqIO.read(f, 'genbank')
                insertdoc(db, collection, data)
                print("data inserted to the collection")
            elif input_file.endswith('csv') or input_file.endswith('tsv') or input_file.endswith('txt'):
                data = f.read()
                insertdoc(db, collection, data)
                print("data inserted to the database")
            else:
                print("The input file is not in the correct format, please use either fasta, json, pkl, genbank, csv, tsv, or txt")
                sys.exit(1)
                # potential add-ons: allow conversion to the required format
        else:
            print("The input file is empty")
            sys.exit(1)

    
    # check the database
    print("The database contains the following collections:")
    print(db.list_collection_names())
    
    # check the collection
    print("The collection contains the following documents:")
    for doc in collection.find():
        print(doc)

    # list the databases
    if list_db:
        print("The databases are:")
        print(db.list_database_names())

    # list the collections
    if list_collections:
        print("The collections are:")
        print(db.list_collection_names())

    # search the database
    if search:
        if search.endswith('.fasta'):
            data = SeqIO.parse(search, 'fasta')
            for record in data:
                print(get_doc(collection, record))
        elif search.endswith('json'):
            data = json.load(search)
            print(get_doc(collection, data))
        elif search.endswith('pkl'):
            data = pickle.load(search)
            print(get_doc(collection, data))
        elif search.endswith('gb'):
            data = SeqIO.read(search, 'genbank')
            print(get_doc(collection, data))
        elif search.endswith('csv') or search.endswith('tsv') or search.endswith('txt'):
            data = search.read()
            print(get_doc(collection, data))
        else:
            print("File not found. Please enter valid file name.")
            sys.exit(1)

    # get metadata from the genbank file
    metadata = Gb_parser(genbank)
    for record in metadata.records:
        insertdoc(db, collection, record)

if __name__ == "__main__":
    main()