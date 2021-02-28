#!/usr/bin/env python

"""
Introduction to Programming - SPRINT3 Coding Project
Arjun Ryatt :  10770217
Carissa Fletcher-Regez :
Cristina Garcia Moratalla : 10777695

Given a RefSeq NM_ accession number:

1. Generate a dictionary with transcript and protein sequences and IDs. Create directory for outputs.
2. Submit either the nucleotide sequence or amino acid sequence to the BLAST API using Biopython NCBIWWW. Save XML file to respective directory.
3. Extract relevant information from BLAST XML output. Compile dictionary of the top 10 most homologous sequences.
4. Generate .fasta file consisting of homologous sequences in FASTA format by querying RefSeq using the Entrez module from the Bio package. Save output .fasta file to the same directory that contains the BLAST put XML file.
5. Construct a multiple-sequence alignment using Clustal Omega. Save the multiple-sequence alignment in the same directory as the BLAST XML file and .fasta file.
"""


import xml.etree.ElementTree as ET  # To enable .xml manipulation
import requests  # To query web services
import os  # For management of input and output files
import time  # To enable time delay between query status checks
from Bio import Entrez # To get NCBI records and manipulate sequences
from Bio import SeqIO # To get NCBI records and manipulate sequences
from Bio.Blast import NCBIWWW  # To query NCBI BLAST


print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>< Hello! ><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

def get_transcript_info(transcript_id):
    """User enters an NM_ accession number into the Ref_Seq_ID field in the Swagger UI. This should only identify protein-coding transcripts in Genbank. If a corresponding record is found, a dictionary is created which contains the RefSeq transcript ID, cDNA sequence, the RefSeq ID of what is translated by this transcript and the sequence of the corresponding protein.

    Returns:
        d1 [dictionary]: A dictionary containing the relevant information to be queried using BLASTn or BLASTp.
    """
    # The function will only proceed if the user's input begins with 'NM_'.
    if transcript_id.startswith("NM"):

        #Using Entrez.efetch, the nucleotide database in Genbank is queried using the NM_ number.
        handle = Entrez.efetch(db="nucleotide", id=transcript_id, rettype="gb", retmode="text")
        # Parse the handle into a SeqRecord object called 'record'.
        record = SeqIO.read(handle, "gb")
        handle.close()
        
        # Creates a dictionary called 'd1'.
        d1 = {}
        
        # Assigns the the ID found in the 'id' attribute of the record to 'id'.
        id = record.id
        # Assigns the cDNA sequence from the 'seq' attribute of the record to 'cdna'.
        cdna = record.seq
        
        # Creates keys in the dictionary called 'id' and 'cdna' and compiles the corresponding data from the record, to each key.
        d1['id'] = id
        d1['cdna'] = str(cdna)

        # Extracts the contents of the 'features' attribute from the record, including the RefSeq ID of what is translated by this transcript and the sequence of the corresponding protein.
        element = 0
        for feature in record.features:
            if feature.type == 'CDS':  # Note, a non-coding transcript would have the type ncRNA
                report = 'Printing record.features[%s]' % str(element)
                break
            else:
                element = element + 1

        # Creates the key 'translation' in d1 and assigns it an empty dictionary. A dictionary called 'translation' inside of the dictionary called 'd1'
        d1['translation'] = {}
        
        # Creates 'id' key and 'prot' key, in the 'translation' dictionary of d1. Assigns the RefSeq ID of what is translated by this transcript and the sequence of the corresponding protein to the respective key.
        d1['translation']['id'] = record.features[element].qualifiers['protein_id'][0]
        d1['translation']['prot'] = record.features[element].qualifiers['translation'][0]
        # Returns the nucleotide and amino acid sequences, as well as their respective RefSeq accession numbers.
        return d1
    
    # If what is entered by the user does not begin with 'NM_', the script stops functioning.
    else:
        return exit()


def blastn_search(query_sequence, transcript_id):
    """If the user specified 'Transcript' in the Swagger UI, the sequence of the transcript is processed by BLASTN. This returns a list of the top 50 most closely related DNA sequences to the query sequence stored in d1['cdna']. Some of the sequences that are returned are 'PREDICTED' sequences and some are orthologues found in other species.

    The function stores the output XML file from BLAST into a directory named after the RefSeq accession number of the sequence (transcript_id) used to query BLAST. The output is named after the RefSeq accession number followed by '_blast_output.xml'. If the directory already exists in the current directory, the most recent output OVERWRITES the files that share the same name as the output. If the directory does not exist, then the directory is created in the current directory.

    A separate dictionary is returned that is compiled of the top 10 sequences from the BLAST output, that share the highest percentage identity with the query sequence.

    Args:
        query_sequence ([string]): Nucleotide sequence to submit to BLASTN.
        transcript_id ([string]): RefSeq accession number of transcript entered by user.

    Returns:
        blast_top10 ([dictionary]): A dictionary consisting of the top 10 sequences from the BLAST output, that share the highest percentage of identity with the query sequence.
    """


    # Submit query sequence to NCBI BLASTN and read results. BLAST returns results in XML format. Results are assigned to 'blast_results'.
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
    blast_results = result_handle.read()
    result_handle.close()
    
 
    # A path to the directory where the XML file will be stored, is created.
    query_name = transcript_id
    current_dir = os.getcwd()
    path_name = current_dir + '/' + query_name
    
    # if/else statement to check if directory already exists. If it does not exist, the directory is created.
    if os.path.exists(path_name):
        pass
    else:
        os.mkdir(path_name)
    
    # The contents of the 'blast_results' are written into a file named after the transcript RefSeq accession ID followed by '_blast_output.xml'. That file is saved into the directory created previously.
    output_file = query_name + '_blast_output.xml'
    output_path = os.path.join(path_name, output_file)
    
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)
        
    # create an empty dictionary called 'blast_hits'
    blast_hits = {}
    # Parse the .xml output from BLASTN into 'tree'
    tree = ET.parse(output_path)
    root = tree.getroot()

    # For each sequence or 'hit' returned by BLAST...
    for hit in root.iter('Hit'):
        # ...create an empty dictionary within the 'blast_hits' dictionary, named after RefSeq accession of the sequence...
        blast_hits[hit.find('Hit_id').text] = {}
        # ...extract the sequence description...
        gene_description = hit.find('Hit_def').text

        #...extract the length of the sequence...
        length = int(hit.find('Hit_len').text)
        # ...and the NUMBER of times that length shares the same base, at the same position as the query sequence...
        for hsp in hit.iter('Hsp'):
            identity = int(hsp.find('Hsp_identity').text)
        #...and times the ratio of the number of identitical bases:the full length of the sequence, by 100, to get the percentage of the sequence that resembles the query sequence.
        percent_identity = 100 * (identity / length)
        
        # Into the dictionary created at the beginning of this for loop, create the keys 'gene_description' and 'percent_identity'. Assign the related information from each sequence to their respective keys.
        blast_hits[hit.find('Hit_id').text]['gene_description'] = gene_description
        blast_hits[hit.find('Hit_id').text]['percent_identity'] = percent_identity
    
    # Arrange the sequences returned by BLAST in order of percent_identity.
    blast_closest_hits = sorted(blast_hits.values(), key=lambda value: value['percent_identity'], reverse=True)
    
    # Creat ANOTHER empty dictionary.
    blast_top10 = {}
    
    # Stores the RefSeq accession number and accompanying gene_description and percent_identity in the blast_top10 dictionary if the gene_description AND percent_identity are the same as the top 10 most identitcal sequences to the query sequence.
    for key, val in blast_hits.items():
        RefSeq_ID = key.split('|')[3]
        if val in blast_closest_hits[:10]:
            blast_top10[RefSeq_ID] = val

    # returns the blast_top10 dictionary
    return blast_top10
    
def blastp_search(query_sequence, transcript_id, translation_id):
    """If the user specified 'Protein' in the Swagger UI, the sequence of the protein encoded by transcript is processed by BLASTP. This returns a list of the top 50 most closely related amino acid sequences to the query sequence stored in d1['translation']['prot']. Some of the sequences that are returned are 'PREDICTED' sequences and some are orthologues found in other species.

    The function stores the output XML file from BLAST into a directory named after the RefSeq accession number of the sequence (translation_id) used to query BLAST. The output is named after the RefSeq accession number followed by '_blast_output.xml'. The output file will be saved into a directory named after the RefSeq accession number of the query sequence. This directory will be stored in a directory named after the 'NM_' accession number (transcript_id) entered by the user.*
    
    * path = path to current directory / NM_ directory / directory named after protein RefSeq accession number / <protein RefSeq accession number>_blast_output.xml
    
    If the 'NM_' directory or the query sequence directory already exist in their respective locations, the most recent output OVERWRITES the files that share the same name as the output. If the directories do not exist, then they are created in the current directory.

    A separate dictionary is returned that is compiled of the top 10 sequences from the BLAST output, that share the highest percentage of identity with the query sequence.

    Args:
        query_sequence ([string]): Amino acid sequence to submit to BLASTP.
        transcript_id ([string]): RefSeq accession number of transcript entered by user.
        translation_id ([string]): RefSeq accession number of protein encoded by transcript entered by user
    Returns:
        blast_top10 [dictionary]: A dictionary consisting of the top 10 sequences from the BLAST output, that share the highest percentage identity with the query sequence.
    """

    
    # Submit query sequence to NCBI BLASTN and read results. BLAST returns results in XML format. Results are assigned to 'blast_results'.
    result_handle = NCBIWWW.qblast("blastp", "nr", query_sequence)
    blast_results = result_handle.read()
    result_handle.close()

    # A path to the NM_ directory is created.
    query_name = transcript_id
    current_dir = os.getcwd()
    path_name = current_dir + '/' + query_name
    
    # if/else statement to check if NM_ directory and directory named after protein RefSeq accession number already exist. If they do not exist, the directories are created.
    if os.path.exists(path_name):
        if os.path.exists(path_name + '/' + translation_id):
            pass
        else:
            os.mkdir(path_name + '/' + translation_id)
    else:
        os.mkdir(path_name)
        os.mkdir(path_name + '/' + translation_id)

    # The contents of the 'blast_results' are written into a file named after the protein RefSeq accession ID followed by '_blast_output.xml'. That file is saved into the directories created previously.
    output_file = translation_id + '_blastp_output.xml'
    output_path = os.path.join(path_name, translation_id, output_file)
    
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)

    # create an empty dictionary called 'blast_hits'
    blast_hits = {}
    # Parse the .xml output from BLASTN into a tree
    tree = ET.parse(output_path)
    root = tree.getroot()

    # For each sequence or 'hit' returned by BLAST...
    for hit in root.iter('Hit'):
        # ...create an empty dictionary within the 'blast_hits' dictionary, named after RefSeq accession of the sequence...
        blast_hits[hit.find('Hit_id').text] = {}
        # ...extract the sequence description...
        protein_description = hit.find('Hit_def').text
        
        #...extract the length of the sequence...
        length = int(hit.find('Hit_len').text)
        # ...and the NUMBER of times that length shares the same amino acid, at the same position as the query sequence...
        for hsp in hit.iter('Hsp'):
            identity = int(hsp.find('Hsp_identity').text)
        #...and times the ratio of the number of identitical bases:the full length of the sequence, by 100, to get the percentage of the sequence that resembles the query sequence.
        percent_identity = 100 * (identity / length)

        # Into the dictionary created at the beginning of this for loop, create the keys 'protein_description' and 'percent_identity'. Assign the related information from each sequence to their respective keys.
        blast_hits[hit.find('Hit_id').text]['protein_description'] = protein_description
        blast_hits[hit.find('Hit_id').text]['percent_identity'] = percent_identity
        
        
    # Arrange the hits in order of percent_identity and create a dictionary out of the top 10.
    blast_closest_hits = sorted(blast_hits.values(), key=lambda value: value['percent_identity'], reverse=True)

    # Creat ANOTHER empty dictionary.
    blast_top10 = {}
    
    # Stores the RefSeq accession number and accompanying protein_description and percent_identity in the blast_top10 dictionary if the protein_description AND percent_identity are the same as the top 10 most identitcal sequences to the query sequence.
    for key, val in blast_hits.items():
        RefSeq_ID = key.split('|')[1]
        if val in blast_closest_hits[:10]:
            blast_top10[RefSeq_ID] = val

    # returns the blast_top10 dictionary
    return blast_top10
    
def top10_hits_fasta(top_hit_dictionary, transcript_id, sequence_type, translation_id):
    """Identify RefSeq accession number of homologues and associated protein sequences by parsing the keys from the dictionary generated by either the blastn_search or blastp_search functions. Query RefSeq to return the corresponding sequences in FASTA format. Output one  sequence from each species in FASTA format. Save all the sequences returned into .fasta file. Depending on if the user chooses 'Transcript or 'Protein' in the API, the .fasta file is saved into the respective directory.

    Args:
        top_hit_dictionary ([dictionary]): Dictionary containing relevant RefSeq accession IDs.
        transcript_id ([string]): RefSeq accession number of transcript entered by user.
        sequence_type ([string]): 'Transcript' or 'Protein', depending on user preference.
        translation_id ([string]): RefSeq accession number of protein encoded by transcript entered by user
    Returns:
        fasta_path ([fasta_path]): Path to .fasta file consisting of the top 10 sequences from the BLAST output, that share the highest percentage identity with the query sequence.
    """
    # Extract RefSeq accession number from blast_top10 into 'accession_numbers' python list.
    accession_numbers = []
    for key, val in top_hit_dictionary.items():
        accession_numbers.append(key)

    # Query RefSeq to return sequences of each RefSeq accession number in 'accession_numbers' list.
    fasta_handle = Entrez.efetch(db = sequence_type, id = accession_numbers, rettype = 'fasta', retmode='text')
    fasta_record = fasta_handle.read()
    fasta_handle.close()

    # A path to the NM_ directory is created.
    query_name = transcript_id + '_HOMOLOGUES'
    current_dir = os.getcwd()
    path_name = current_dir + '/' + transcript_id

    # if/else statement to check if the user is requesting .fasta file consisting of amino acid or nucleotide sequence homologues. Saves contents from 'fasta_record' into the respective directory.
    if sequence_type == 'nucleotide':
        fasta_file = query_name + '.fasta'
        fasta_path = os.path.join(path_name, fasta_file)
    elif sequence_type == 'Protein':
        fasta_file = translation_id + '_HOMOLOGUES.fasta'
        fasta_path = os.path.join(path_name, translation_id, fasta_file)
        
    with open(fasta_path, 'w') as save_file:
        save_file.write(fasta_record)

    # returns the path to the .fasta file
    return fasta_path
    

def clustal_omega_MSA(email_address, fasta_path, transcript_id, sequence_type, translation_id):
    """Use the formatted FASTA file from the top10_hits_fasta function as input to Clustal Omega
    to create a multiple sequence alignment. Output the MSA as a text file. Ouput is saved to the same location as the input .fasta file.

    Args:
        email_address ([string]): Required to access Clustal Omega API.
        fasta_path ([file path]): Path to .fasta file containing homologous sequences.
        transcript_id ([string]): RefSeq accession number of transcript entered by user.
        sequence_type ([string]): 'Transcript' or 'Protein', depending on user preference.
        translation_id ([string]): RefSeq accession number of protein encoded by transcript entered by user
    Returns:
        msa [string]: Multiple sequence alignment of homologous protein
        sequences.
    """
    # Open the .fasta file of homologous sequences
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
    
    # A path to the NM_ directory is created.
    job_name = transcript_id
    current_dir = os.getcwd()
    path_name = current_dir + '/' + job_name
    
    # if/else statement to check if the user is requesting multi-sequence alignment of amino acid or nucleotide sequence homologues. Saves output MSA text file into the respective directory.
    if sequence_type == 'Transcript':
        msa_file = job_name + '_msa.aln'
        msa_path = os.path.join(path_name, msa_file)
    elif sequence_type == 'Protein':
        msa_file = translation_id + '_msa.aln'
        msa_path = os.path.join(path_name, translation_id, msa_file)

    # Write output text file
    with open(msa_path, 'w') as msa_file_object:
        msa_file_object.write(msa)
    
    # return multi-sequence alignment
    return msa
