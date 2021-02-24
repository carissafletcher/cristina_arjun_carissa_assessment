

import xml.etree.ElementTree as ET  # To enable .xml manipulation
import requests  # To query web services
import matplotlib.pyplot as plt # To generate figure of phylogenetic tree
import os  # For management of input and output files
import json, re
from Bio import Entrez, SeqIO, Phylo  # To get NCBI records and manipulate sequences
from Bio.Blast import NCBIWWW  # To query NCBI BLAST

Entrez.email = "arjunryatt@hotmail.co.uk"

print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>< Hello! ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

def get_transcript_info(transcript_id, path_name):

    if transcript_id.startswith("NM"):

        handle = Entrez.efetch(db="nucleotide", id=transcript_id, rettype="gb", retmode="text")
        # Parse the handle into a SeqRecord object called record
        record = SeqIO.read(handle, "gb")
        handle.close()
    
        d1 = {}
    
        id = record.id
        cdna = record.seq
        rna = cdna.transcribe()
        d1['id'] = id
        d1['cdna'] = str(cdna)
        d1['rna'] = str(rna).lower()

        element = 0
        for feature in record.features:
            if feature.type == 'CDS':  # Note, a non-coding transcript would have the type ncRNA
                report = 'Printing record.features[%s]' % str(element)
            #print(record.features[element])
                break
            elif feature.type == 'misc_RNA' :  # Note, a non-coding transcript would have the type ncRNA
                report = 'Printing record.features[%s]' % str(element)
                break
            else:
                element = element + 1
    
    
        d1['translation'] = {}
        d1['translation']['id'] = record.features[element].qualifiers['protein_id'][0]
        d1['translation']['prot'] = record.features[element].qualifiers['translation'][0]
    
        return d1
    
    else:
        return 'Protein-coding transcripts only! Please enter a RefSeq accession number that starts with "NM_".'


def blast_search(query_sequence):
    
    # Submit query sequence to NCBI BLASTN and read results
    
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
    blast_results = result_handle.read()
    result_handle.close()
    
    # return blast_results
    
    query_name = 'msa_project'
    path_name = os.getcwd()
    output_file = query_name + '_blast_output.xml'
    output_path = os.path.join(path_name, output_file)
    
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)
    
    return output_path

print(">>>>>>>>>>>>>>>>>>>>>>>>>>< Compiling Transcript Dictionary ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

transcript_id = 'NM_007298.3'
transcript_dictionary = get_transcript_info(transcript_id)
print(transcript_dictionary)


print(">>>>>>>>>>>>>>>>>>< Running Multiple Sequence Alignments through BLASTN ><<<<<<<<<<<<<<<<<<<<<<")

# query_sequence_test = transcript_dictionary['cdna']

query_sequence_test = transcript_dictionary['translation']['prot']

test = blast_search(query_sequence_test)

print(test)
