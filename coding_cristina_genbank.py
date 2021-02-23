

from datetime import datetime  # To record time of BLAST query submission
import os  # For management of input and output files
import time  # To enable time delay between query status checks
import xml.etree.ElementTree as ET  # To enable .xml manipulation

# import randfacts  # To alleviate boredom during BLAST searches
import requests  # To query web services
import matplotlib.pyplot as plt # To generate figure of phylogenetic tree

from Bio import Entrez, SeqIO, Phylo  # To get NCBI records and manipulate sequences
from Bio.Blast import NCBIWWW  # To query NCBI BLAST

import os.path
from Bio.Seq import Seq


 # User defines a name for the query (NM_ number, this depends on Carissa's code)
def get_query_name(query_name):
    # User defines a name for the query
    # query_name = input('Please enter a name for this query: ')
    # while not query_name:
        # query_name = input('Please enter a name for this query: ')
    
    
    
    current_dir = os.getcwd() #Directory 
    path_name = current_dir + '/' + query_name

    # Create a new folder using the query name (must be unique)
    sentry = False
    while sentry == False:
        if not os.path.exists(path_name): # os.path.exists is used to check whether the specified path exists or not.
            os.mkdir(path_name) #os.mkdir is used to create a path with the specific numeric mode
            sentry = True
        else:
            print('Folder already exists, please choose another name.')
            query_name = input('Please enter a name for this query: ')
            path_name = current_dir + '/' + query_name
        
    return query_name, path_name

 #Function to extract CDS feature in sequence record
def get_cds_feature_with_qualifier_value(seq_record, name, value):
    """. You can use this for finding features by locus tag, gene ID, or protein ID.
    """
    # Loop over the features
    for feature in genome_record.features:
        if feature.type == "CDS" and value in feature.qualifiers.get(name, []):
            return feature
            
    # Could not find it
    return None

#Function to extract exon and write file
def extract_cds():
    #Extract cds sequence
    gene_sequence = cds_feature.extract(genome_record.seq)
    cds = str(gene_sequence)
    print(cds)

    #Write file.txt, ideally query_name is the name entered by the user
     # Here I use again the function because otherwise it doesn't recognize query_name. how can we avoid that???
    query_name = input('Please enter a name for this query: ')
    while not query_name:
        query_name = input('Please enter a name for this query: ')

    save_path = os.getcwd() + '/' + query_name
    output_file_name = query_name
    completeName = os.path.join(save_path, output_file_name+"_cds.txt")  #Join one or two more path components

    with open(completeName, "a") as file:
        for nucleotide in cds:
             file.write(nucleotide)
           
#Function to translate cdna sequence to protein sequence (this could be maybe combined with the previous function)
def translate_cds_seq():
    #Extract cds sequence (here it is redundant, maybe combine is a good idea...)
    gene_sequence = cds_feature.extract(genome_record.seq)
    cds = str(gene_sequence)
    print(cds)

    cds_object = Seq(cds)
    protein = cds_object.translate(table=1, cds=True, to_stop=True)

    #Write file.txt, ideally query_name is the name entered by the user
     # Here I use again the function because otherwise it doesn't recognize query_name. how can we avoid that???
    query_name = input('Please enter a name for this query: ')
    while not query_name:
        query_name = input('Please enter a name for this query: ')

    save_path = os.getcwd() + '/' + query_name
    output_file_name = query_name
    completeName = os.path.join(save_path, output_file_name+"_protein.txt")  #Join one or two more path components

    with open(completeName, "a") as file:
        for sequence in protein:
             file.write(sequence)


'''get_query_name('test1')'''





# genome_record = SeqIO.read("/Users/carlossilvestreroig/Desktop/Python_API/adjuntos/NM_005816.gb", "genbank") #this should be query_name
# cds_feature = get_cds_feature_with_qualifier_value(genome_record, "protein_id", "NP_005807.1")

#I create a dictionary:
#Dictionary_sequences = {'cdna' = cdna, 
                       # 'protein'}

# extract_cds()
# translate_cds_seq()
