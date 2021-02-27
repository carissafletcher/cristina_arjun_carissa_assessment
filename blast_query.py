

import xml.etree.ElementTree as ET  # To enable .xml manipulation
import requests  # To query web services
import matplotlib.pyplot as plt # To generate figure of phylogenetic tree
import os, shutil  # For management of input and output files
import json, re
import time  # To enable time delay between query status checks
from Bio import Entrez, SeqIO, Phylo  # To get NCBI records and manipulate sequences
from Bio.Blast import NCBIWWW  # To query NCBI BLAST

print('http://127.0.0.1:5000/')


print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>< Hello! ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

def get_transcript_info(transcript_id):

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

        gene_symbol = record.features[element].qualifiers['gene'][0]
        d1['gene_symbol'] = gene_symbol
    
    
        d1['translation'] = {}
        d1['translation']['id'] = record.features[element].qualifiers['protein_id'][0]
        d1['translation']['prot'] = record.features[element].qualifiers['translation'][0]
    
        return d1
    
    else:
        return 'Protein-coding transcripts only! Please enter a RefSeq accession number that starts with "NM_".'


def blastn_search(query_sequence, transcript_id):
    
    # Submit query sequence to NCBI BLASTN and read results. BLAST returns results in XML format
    
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
    blast_results = result_handle.read()
    result_handle.close()
    
    # return blast_results
    
    query_name = transcript_id
    current_dir = os.getcwd()
    path_name = current_dir + '/' + query_name
    
    if os.path.exists(path_name):
        pass
    else:
        os.mkdir(path_name)
        
    output_file = query_name + '_blast_output.xml'
    output_path = os.path.join(path_name, output_file)
    
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)
        
    
    # Parse the .xml output from BLASTN into a tree
    blast_hits = {}
    tree = ET.parse(output_path)
    root = tree.getroot()

    # Use the tree to get the gene description, transcript accession and
    # sequence identity for each hit
    for hit in root.iter('Hit'):
        blast_hits[hit.find('Hit_id').text] = {}
        gene_description = hit.find('Hit_def').text
        transcript = hit.find('Hit_accession').text

        length = int(hit.find('Hit_len').text)
        for hsp in hit.iter('Hsp'):
            identity = int(hsp.find('Hsp_identity').text)
        percent_identity = 100 * (identity / length)

        blast_hits[hit.find('Hit_id').text]['gene_description'] = gene_description
        
        blast_hits[hit.find('Hit_id').text]['transcript'] = transcript
        
        blast_hits[hit.find('Hit_id').text]['percent_identity'] = percent_identity
    
    # Arrange the hits in order of percent_identity and create a dictionary out of the top 10.
    blast_closest_hits = sorted(blast_hits.values(), key=lambda value: value['percent_identity'], reverse=True)
    
    blast_top10 = {}
    for key, val in blast_hits.items():
        RefSeq_ID = key.split('|')[3]
        if val in blast_closest_hits[:10]:
            blast_top10[RefSeq_ID] = val

    return blast_top10
    
def blastp_search(query_sequence, transcript_id, translation_id):
    
    # Submit query sequence to NCBI BLASTP and read results
    
    result_handle = NCBIWWW.qblast("blastp", "nr", query_sequence)
    blast_results = result_handle.read()
    result_handle.close()
    
    query_name = transcript_id
    current_dir = os.getcwd()
    path_name = current_dir + '/' + query_name
    
    if os.path.exists(path_name):
        if os.path.exists(path_name + '/' + translation_id):
            pass
        else:
            os.mkdir(path_name + '/' + translation_id)
    else:
        os.mkdir(path_name)
        os.mkdir(path_name + '/' + translation_id)
    
    output_file = translation_id + '_blastp_output.xml'
    output_path = os.path.join(path_name, translation_id, output_file)
    
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)

    # Parse the .xml output from BLASTN into a tree
    blast_hits = {}
    tree = ET.parse(output_path)
    root = tree.getroot()

    # Use the tree to get the gene description, transcript accession and
    # sequence identity for each hit
    for hit in root.iter('Hit'):
        blast_hits[hit.find('Hit_id').text] = {}
        gene_description = hit.find('Hit_def').text
        transcript = hit.find('Hit_accession').text

        length = int(hit.find('Hit_len').text)
        for hsp in hit.iter('Hsp'):
            identity = int(hsp.find('Hsp_identity').text)
        percent_identity = 100 * (identity / length)

        blast_hits[hit.find('Hit_id').text]['gene_description'] = gene_description
        
        blast_hits[hit.find('Hit_id').text]['transcript'] = transcript
        
        blast_hits[hit.find('Hit_id').text]['percent_identity'] = percent_identity
        
        
    # Arrange the hits in order of percent_identity and create a dictionary out of the top 10.
    blast_closest_hits = sorted(blast_hits.values(), key=lambda value: value['percent_identity'], reverse=True)
    
    blast_top10 = {}
    for key, val in blast_hits.items():
        RefSeq_ID = key.split('|')[1]
        if val in blast_closest_hits[:10]:
            blast_top10[RefSeq_ID] = val

    return blast_top10
    
def top10_hits_fasta(top_hit_dictionary, transcript_id, sequence_type, translation_id):

    accession_numbers = []
    for key, val in top_hit_dictionary.items():
        accession_numbers.append(key)

    fasta_handle = Entrez.efetch(db = sequence_type, id = accession_numbers, rettype = 'fasta', retmode='text')
    fasta_record = fasta_handle.read()
    fasta_handle.close()

    query_name = transcript_id + '_HOMOLOGUES'
    current_dir = os.getcwd()
    path_name = current_dir + '/' + transcript_id
    
    if sequence_type == 'nucleotide':
        fasta_file = query_name + '.fasta'
        fasta_path = os.path.join(path_name, fasta_file)
    elif sequence_type == 'Protein':
        fasta_file = translation_id + '_HOMOLOGUES.fasta'
        fasta_path = os.path.join(path_name, translation_id, fasta_file)
        
    with open(fasta_path, 'w') as save_file:
        save_file.write(fasta_record)

    return fasta_path
    


def clustal_omega_MSA(email_address, fasta_path, transcript_id, sequence_type, translation_id):
    """Use the formatted FASTA file from Homologene as input to Clustal Omega
    to create a multiple sequence alignment. Output the MSA as a text file.

    Args:
        homologues_path ([file path]): Path to .fasta file containing
        homologous protein sequences.
        gene_symbol ([string]): Relevant human gene symbol for homologous
        proteins.
        path_name ([file path]): Path to directory where output .txt file will
        be stored.
        email_address ([string]): Required to access Clustal Omega API.

    Returns:
        msa [string]: Multiple sequence alignment of homologous protein
        sequences.
    """
    # Open the .fasta file of homologous protein sequences
    with open(fasta_path,'r') as homologues_object:
        fasta_transcripts = homologues_object.read()

    # Submit request to Clustal Omega
    print("\nSending request to Clustal Omega.")
    run = ''
    while not str(run) == '<Response [200]>':
        run = requests.post(
            'https://www.ebi.ac.uk/Tools/services/rest/clustalo/run',
            data = {'email':email_address, 'sequence':fasta_transcripts})
        # Exception Handling to catch if email is invalid
        try:
            run.raise_for_status()
        except requests.exceptions.HTTPError:
            print('Error whilst submitting query.',
            'Please ensure email address is valid.')

    # Decode job ID
    job_id_bytes = run.content
    job_id = job_id_bytes.decode('utf-8')

    # Get status of submitted Clustal Omega job
    status_request = requests.get(
        'https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/'
        + str(job_id))
    status = (status_request.content).decode('utf-8')

    # Check status every 10 seconds
    while status != "FINISHED":
        status_request = requests.get(
            'https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/'
            + str(job_id))
        status = (status_request.content).decode('utf-8')
        print(status)
        if status == "RUNNING":
            time.sleep(10)

    # Get msa results from Clustal Omega once completed
    msa_request = requests.get(
        'https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/'
        + str(job_id) +'/aln-clustal_num')
    msa = (msa_request.content).decode('utf-8')
    
    # Define path to output text file
    job_name = transcript_id
    current_dir = os.getcwd()
    path_name = current_dir + '/' + job_name
    
    if sequence_type == 'Transcript':
        msa_file = job_name + '_msa.aln'
        msa_path = os.path.join(path_name, msa_file)
    elif sequence_type == 'Protein':
        msa_file = translation_id + '_msa.aln'
        msa_path = os.path.join(path_name, translation_id, msa_file)

    # Write output text file
    with open(msa_path, 'w') as msa_file_object:
        msa_file_object.write(msa)
    
    return msa

'''print(">>>>>>>>>>>>>>>>>>>>>>>>>>< Compiling Transcript Dictionary ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

transcript_id = 'NM_007298.3'
transcript_dictionary = get_transcript_info(transcript_id)
print(transcript_dictionary)


print(">>>>>>>>>>>>>>>>>>< Running Multiple Sequence Alignments through BLASTN ><<<<<<<<<<<<<<<<<<<<<<")

query_sequence_test = transcript_dictionary['cdna']
blastn = blastn_search(query_sequence_test, transcript_id)
print(blastn)

print(">>>>>>>>>>>>>>>>>>>>>>>>>>< Getting Protein sequences from BLASTP ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

query_sequence_test = transcript_dictionary['translation']['prot']
blastp = blastp_search(query_sequence_test, transcript_id, transcript_dictionary['translation']['id'])
print(blastp)

print(">>>>>>>>>>>>>>>>>>>>>>>< Preparing FASTA of homologous sequences ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

fasta_path = top10_hits_fasta(blastn, transcript_id, 'nucleotide', '')
print(top_10_transcripts)

print(">>>>>>>>>>>>>>>>>>>>>>>>>>>< Preparing MSA with ClustalOmega ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

email_address = 'arjunryatt@hotmail.co.uk'
print(clustal_omega_MSA(email_address, fasta_path, transcript_id, 'Transcript', ''))

print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>< End of script. Goodbye! ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")'''
