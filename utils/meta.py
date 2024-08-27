import re

INPUT_DIR = str(input("Enter the path to the directory containing the genbank files: "))
OUTPUT_DIR = str(input("Enter the path to the directory where the output files will be saved: "))
Metadata_file = str(input("Enter the name of the metadata file output: "))
Fullname_file = str(input("Enter the name of the full name record output: "))

class Genbank_record:
    def __init__(self, block):
        self.accession = self.extract_accession(block) 
        self.locus = self.extract_locus(block)
        self.source = self.extract_source(block)
        self.organism = self.extract_organism(block)
        self.sample_host = self.extract_sample_host(block) 
        self.country = self.extract_country(block) 
        self.collection_data = self.extract_collection_data(block) 
        self.fullname = self.defination(block) 

    def extract_accession(self, block):
        accession = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'ACCESSION\s+(\w+)', line)
            if search:
                accession = search.group(1).replace('"', '')
        return accession
    def extract_locus(self, block):
        locus = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'LOCUS\s+(\w+)', line)
            if search:
                locus = search.group(1).replace('"', '')
        return locus
    
    def extract_source(self, block):
        source = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'SOURCE\s+(.*)', line)
            if search:
                source = search.group(1).replace('"', '')
        return source
    
    def extract_organism(self, block):
        organism = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'ORGANISM\s+(.*)', line)
            if search:
                organism = search.group(1).replace('"', '')
        return organism

    def extract_sample_host(self, block):
        sample_host = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'/host=(.*)', line)
            if search:
                sample_host = search.group(1).replace('"', '')
        return sample_host
    
    def extract_country(self, block):
        country = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'/geo_loc_name=(.*)', line)
            if search:
                country = search.group(1).replace('"', '')
        return country
    
    def extract_collection_data(self, block):
        collection_data = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'/collection_date=(.*)', line)
            if search:
                collection_data = search.group(1).replace('"', '')
        return collection_data
    
    def defination(self, block):
        fullname = None
        for line in block.split('\n'):
            line = line.strip()
            search = re.search(r'DEFINITION\s+(.*)', line)
            if search:
                fullname = search.group(1).replace('"', '')
        return fullname

    def format(self):
        accession = self.accession
        sample_host = self.sample_host
        country = self.country
        collection_data = self.collection_data
        return f"{accession}|{sample_host}|{country}|{collection_data}"
    
    def seq_record(self):
        accession = self.accession
        fullname = self.fullname
        return f"{accession}|{fullname}"
    
    def __str__(self):
        return f"{self.accession}|{self.fullname}|{self.sample_host}|{self.country}|{self.collection_data}"
        


class Gb_parser:
    def __init__(self, genbank):
        self.genbank = genbank
        self.block = self.genbank2block()
        self.records = self.block2genbank_record()

    def genbank2block(self):
        block = []
        in_block = False
        current_block = ''
        with open(self.genbank, 'r') as f:
            lines = f.readlines()

        # each block start from the detection of 'LOCUS' & end by 'ORIGIN'
        for line in lines:
            if line.strip().startswith('LOCUS'):
                in_block = True
            elif line.strip().startswith('ORIGIN'):
                in_block = False
                block.append(current_block)
                current_block = ''
            if in_block:
                current_block += line
        return block

    def block2genbank_record(self):
        gb_ls = []
        for i in self.block:
            gb = Genbank_record(i)
            gb_ls.append(gb)
        print(f'There are {len(gb_ls)} genbank records in {self.genbank}')
        return gb_ls



## How to use example
# rota_B_gb = Gb_parser(os.path.join(INPUT_DIR, Metadata_file))
# # write the metadata & full names to rota_B_metadata.txt & rota_B_fullname.txt
# with open(os.path.join(OUTPUT_DIR, Metadata_file), 'a') as f:
#     for record in rota_B_gb.records:
#         f.write(record.format() + '\n')
# with open(os.path.join(OUTPUT_DIR, Fullname_file), 'a') as f:
#     for record in rota_B_gb.records:
#         f.write(record.seq_record() + '\n')
# num_written = len(rota_B_gb.records)

# print('Total number of records written to {}}: {}'.format(Metadata_file, num_written))
# print('Total number of records written to {}}: {}'.format(Fullname_file, num_written))
