#!/usr/bin/env python

"""
The Homologenator.

Introduction to Programming - SPRINT3 Paired Coding Project
Thomas Scott-Adams :  9627185
Jay Miles          : 10806682

Given a query nucleotide sequence in a fasta file:

1. Get user-defined query name, create new directory for output.
2. Get user-defined nucleotide sequence as string or fasta file.
3. Confirm that any .fasta file is of the correct format.
4. Submit sequence to BLASTN using Biopython NCBIWWW.
5. Extract relevant information from BLASTN XML output.
6. Identify closest sequence match and corresponding gene symbol.
7. Identify homologous protein sequences using Homologene API.
8. Format homologous protein sequences as a .fasta file.
9. Construct a multiple sequence alignment using Clustal Omega.
"""

from datetime import datetime  # To record time of BLAST query submission
import os  # For management of input and output files
import time  # To enable time delay between query status checks
import xml.etree.ElementTree as ET  # To enable .xml manipulation

import randfacts  # To alleviate boredom during BLAST searches
import requests  # To query web services
import matplotlib.pyplot as plt # To generate figure of phylogenetic tree

from Bio import Entrez, SeqIO, Phylo  # To get NCBI records and manipulate sequences
from Bio.Blast import NCBIWWW  # To query NCBI BLAST


def get_query_name():
    """User inputs a name for the query. A new directory is created to
    store file outputs of the program. The directory is named based on
    the user input.

    Returns:
        query_name [string]: User-defined name associated with the query.
        path_name [file path]: Path to new directory.
    """
    # User defines a name for the query
    query_name = input('Please enter a name for this query: ')
    while not query_name:
        query_name = input('Please enter a name for this query: ')
    
    current_dir = os.getcwd()
    path_name = current_dir + '/' + query_name

    # Create a new folder using the query name (must be unique)
    sentry = False
    while sentry == False:
        if not os.path.exists(path_name):
            os.mkdir(path_name)
            sentry = True
        else:
            print('Folder already exists, please choose another name.')
            query_name = input('Please enter a name for this query: ')
            path_name = current_dir + '/' + query_name
        
    return query_name, path_name


def fasta_check(filename):
    """Checks whether the specified file is in FASTA format. If the file is
    in FASTA format, SeqIO will be able to parse it as a FASTA file. If not,
    the 'fasta' variable will be empty.

    Args:
        filename ([string]): Name of file to evaluate.

    Returns:
        [bool]: False if the specified file is not in FASTA format.
    """
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            return True
        else:
            return False


def acquire_input():
    """User inputs a query sequence, either by direct entry of a string or by
    specifying an existing FASTA file. If a file is specified, the contents
    are read and stored for subsequent use.

    Returns:
        query_sequence [string]: Sequence to be queried using BLASTN.
    """
    # User specifies format of query sequence
    print('\nIMPORTANT!\nIf your transcript sequence is longer than 1 line,',
        'please ensure you submit it as a .fasta file.\n')
    input_type = input('Transcript format is (1) string (2) fasta file? ')
    while not ((input_type == '1') or (input_type == '2')):
        print('Please enter either 1 or 2.')
        input_type = input(
            'Transcript format is (1) string (2) fasta file? ')

    # Format a directly-input sequence
    if input_type == '1':
        string_input = input('Paste transcript sequence here: ')
        upper_case_input = string_input.upper()
        query_sequence = (upper_case_input.replace(' ', '')).strip()
        print('Sequence to query:\n', query_sequence)

        ok_input = [
            'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R',
            'S', 'T', 'U', 'V', 'W', 'Y', '-'
            ]  # Characters comprising the degenerate nucleotide code
        for character in query_sequence:
            if character not in ok_input:
                print('Warning, unsupported character: ' + character)  
                # Warning only, stills runs query

    # Check correct file extension/format for a sequence in a .fasta file
    elif input_type == '2':
        filename = input('Enter filename or full path to file: ')
        while not os.path.isfile(filename):
            print('File does not exist, please ensure path is correct.')
            filename = input('Enter filename or full path to file: ')
        while not (filename.endswith('.fasta')):
            print('File must be .fasta type.')
            filename = input('Enter filename or full path to file: ')
        while fasta_check(filename) == False:
            print('Data must be in fasta format.')
            filename = input('Enter filename or full path to file: ')
        with open(filename,'r') as file_object:
            query_sequence = file_object.read()

    return query_sequence


def blast_search(query_name, path_name, query_sequence):
    """Submit a query nucleotide sequence to BLASTN, and output the raw data
    as an .xml file.

    Args:
        query_sequence ([string]): Nucleotide sequence to submit to BLASTN.
        query_name ([string]): User-defined name for query.
        path_name ([file path]): Directory to store output file in.

    Returns:
        output_path [file path]: Path to output .xml file.
    """
    # Submit query sequence to NCBI BLASTN and read results
    print('\nQuerying BLAST - This may take a long time during peak times.')
    for i in range(1,4):
        print('Querying BLASTN: Attempt ' + str(i))
        print((datetime.now()).strftime("%Y-%m-%d %H:%M:%S"))
        print('\n...whilst you\'re waiting, here\'s a randomly-generated',
            'interesting fact:')
        print(randfacts.getFact() + "\n")
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence)
            break
        except Exception as ex:
            print('Error occurred: ', ex.__class__.__name__)
            if i == 3:
                print('Unable to access BLASTN. Exiting script.')
                exit()
    print('BLAST query complete.')
    blast_results = result_handle.read()
    result_handle.close()

    # Define file path for output file
    output_file = query_name + '_blast_output.xml'
    output_path = os.path.join(path_name, output_file)

    # Write output file
    with open(output_path, 'w') as save_file:
        save_file.write(blast_results)
    print("\nBLAST output xml file created.")

    return output_path


def get_blast_hits(filepath):
    """Process raw BLASTN output from an .xml file to identify the gene 
    description, transcript accession, and sequence identity associated with 
    each potential sequence match.

    Args:
        filepath ([file path]): .xml file containing raw BLASTN output

    Returns:
        blast_closest_hits [list]: Each list element contains the gene 
        description, transcript accession, and sequence identity associated 
        with a potential sequence match.
    """
    # Parse the .xml output from BLASTN into a tree
    blast_closest_hits = []
    tree = ET.parse(filepath)
    root = tree.getroot()

    # Use the tree to get the gene description, transcript accession and 
    # sequence identity for each hit
    for hit in root.iter('Hit'):
        gene_description = hit.find('Hit_def').text
        transcript = hit.find('Hit_accession').text

        length = int(hit.find('Hit_len').text)
        for hsp in hit.iter('Hsp'):
            identity = int(hsp.find('Hsp_identity').text)
        percent_identity = 100 * (identity / length)

        blast_closest_hits.append([gene_description, transcript, 
            percent_identity])

    print('There are ' + str(len(blast_closest_hits)) + ' sequence matches.')

    return blast_closest_hits


def find_best_match(blast_closest_hits):
    """Identify which list item has the highest sequence identity with the 
    query sequence, and return its associated gene symbol. If multiple 
    sequences share the highest identity, create a list of gene symbols and 
    determine whether they are duplicates.

    Args:
        blast_closest_hits ([list]): Each list element contains the gene 
        description, transcript accession, and sequence identity associated 
        with a potential sequence match.

    Returns:
        gene_symbol [string]: Gene symbol associated with the transcript(s)
        whose sequences most closely match that of the query sequence.
    """
    # First, find the hit with the highest sequence identity.
    # hit[2] contains sequence identity.
    current_best = 0
    best_matches = []
    for hit in blast_closest_hits:
        if hit[2] > current_best:
            current_best = hit[2]
    for hit in blast_closest_hits:
        if hit[2] == current_best:
            best_matches.append(hit)
    
    # Format the gene descriptions of these hits to get their gene symbol
    blast_genes = []                             
    for hit in best_matches:
        gene_desc = hit[0]
        if '(' in gene_desc:
            split1 = gene_desc.split('(')
            split2 = split1[1:]  # take everything after brackets open
            split3 = split2[0]  # in case of multiple pairs of brackets
            split4 = split3.split(')')
            split5 = split4[:1]  # keep only contents before brackets close
            blast_genes.append(split5[0])

    # Remove duplicate gene symbols from the list
    unique_genes = []
    for gene in blast_genes:
        if gene not in unique_genes:
            unique_genes.append(gene)

    # If there's only one unique gene symbol, that's the best match
    if len(unique_genes) == 1:
        gene_symbol = unique_genes[0]
        print('Closest gene match: ' + gene_symbol)

    # If there are no gene symbols identified, 
    elif len(unique_genes) == 0:
        print('No gene symbol was identified from this sequence.')
        exit()

    # Otherwise, the user has to select a gene symbol
    else:
        print('\nIMPORTANT:\nMultiple possible gene matches: ')
        print([gene for gene in unique_genes])
        gene_symbol = input('\nPlease select one gene symbol (case-sensitive): ')
        while gene_symbol not in unique_genes:
            print('Selection must be in list above.')
            gene_symbol = input(
                'Please select one gene symbol (case-sensitive): ')

    return gene_symbol


def find_homologues(email_address, gene_symbol):
    """Identify gene homologues and associated protein sequences in multiple 
    species using Homologene. Output one protein sequence from each species in
    FASTA format.

    Args:
        gene_symbol ([string]): Gene symbol associated with the transcript(s)
        whose sequences most closely match that of the query sequence.

    Returns:
        fasta_record [homologene record]: Sequences of protein homologues in 
        FASTA format.
    """
    # Query Homologene with the gene symbol for the best match
    Entrez.email = email_address
    for i in range(1,4):
        try:
            print('\nGetting Homologene ID: Attempt ' + str(i))
            search_handle = Entrez.esearch(
                db='homologene', 
                term=(gene_symbol + '[Gene Name] AND Homo sapiens[Organism]'))
            break
        except Exception as ex:
            print('Error occurred: ', ex.__class__.__name__)
            if i == 3:
                print('Unable to access Homologene. Exiting script.')
                exit()
    search_record = Entrez.read(search_handle)
    search_handle.close()

    # Extract the Homologene ID
    hg_id = search_record['IdList'][0]
    print('Homologene ID: ' + hg_id)

    # Use the ID to request the Homologene record
    for i in range(1,4):
        try:
            print('Getting Homologene record: Attempt ' + str(i))
            hg_handle = Entrez.efetch(
                db='homologene', id = hg_id, rettype = 'homologene', 
                retmode='text')
            break
        except Exception as ex:
            print('Error occurred: ', ex.__class__.__name__)
            if i == 3:
                print('Unable to access Homologene. Exiting script.')
                exit()
    hg_record = hg_handle.readlines()
    hg_handle.close()

    # Format the record to extract list of protein homologue accession numbers
    accession_list = []
    for line in hg_record:
        line_strip = line.strip()
        if line_strip[-3:] == ' aa':
            accession = (line_strip.split(' '))[-2]
            accession_list.append(accession)
    print('Transcript accession numbers: ', [acs for acs in accession_list])

    # Use the ID again to request the homologous protein sequences as .fasta
    for i in range(1,4):
        try:
            print('Getting sequences as fasta: Attempt ' + str(i))
            fasta_handle = Entrez.efetch(
                db='homologene', id = hg_id, rettype = 'fasta', retmode='text')
            break
        except Exception as ex:
            print('Error occurred: ', ex.__class__.__name__)
            if i == 3:
                print('Unable to access Homologene. Exiting script.')
                exit()
    fasta_record = fasta_handle.read()
    fasta_handle.close()

    return fasta_record


def format_fasta(path_name, gene_symbol, fasta_record):
    """Format the FASTA file from Homologene so that each header contains
    the species name, transcript accession, and protein description.

    Args:
        gene_symbol ([string]): Gene symbol associated with the transcript(s)
        whose sequences most closely match that of the query sequence.
        fasta_record ([Homologene record]): FASTA record of homologous protein
        sequences from Homologene.
        path_name ([file path]): Path to directory where output .fasta file 
        will be stored.

    Returns:
        homologues_path [file path]: Path to the output .fasta file.
    """
    # Identify the start of the first protein sequence header and remove 
    # everything preceding it
    first_seq = int(fasta_record.find('>'))
    clipped_fasta_record = (fasta_record[first_seq:]).strip()
    fasta_file = gene_symbol + '_raw_homologues.fasta'
    
    # Get the fasta file from a single string into a list of lines. This is 
    # a REALLY clunky way of doing this, but it works whilst we think of a 
    # better option.
    og_fasta_path = os.path.join(path_name, fasta_file)
    with open(og_fasta_path, 'w') as output_object:
        output_object.write(clipped_fasta_record)
    with open(og_fasta_path, 'r') as input_object:
        fasta_homologues = input_object.readlines()
    print("Raw homologues fasta file created.")

    # For each sequence in the .fasta file, extract the accession number. Use
    # this to get the relevant organism with efetch. 
    i = 0
    for line in fasta_homologues:
        if line[0] == '>':
            split_line = line.split('|')
            handle = Entrez.efetch(
                db="protein", id=split_line[3], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "gb")
            handle.close()
            species = (record.annotations['organism']).replace(' ', '_')

            # Create a new header for each sequence with species and accession
            line_end = line.split('ref|')
            if '[' in line_end[-1]:
                remove_dup_species = line_end[-1].split('[')
                new_line = '>' + gene_symbol + '|' +  str(species) + '|' + remove_dup_species[0] + '\n'
            else:
                new_line = '>' + gene_symbol + '|' +  str(species) + '|' + line_end[-1]
            fasta_homologues[i] = new_line
        i += 1

    # Define the path to the output .fasta file
    homologues_file = gene_symbol + '_formatted_homologues.fasta'
    homologues_path = os.path.join(path_name, homologues_file)

    # Write the output .fasta file
    with open(homologues_path, 'w') as output_object:
        [output_object.write(line) for line in fasta_homologues]
    print("Formatted homologues fasta file created.")

    return homologues_path


def construct_MSA(email_address, path_name, gene_symbol, homologues_path):
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
    with open(homologues_path,'r') as homologues_object:
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
            email_address = input("Please enter a valid email address: ")

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
    msa_file = gene_symbol + '_msa.aln'
    msa_path = os.path.join(path_name, msa_file)

    # Write output text file
    with open(msa_path, 'w') as msa_file_object:
        msa_file_object.write(msa)
    print("\nMultiple Sequence Alignment text file created.")

    # Get newick results from Clustal Omega once completed
    newick_request = requests.get(
        'https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/result/'
        + str(job_id) +'/tree')
    newick = (newick_request.content).decode('utf-8')
    
    # Define path to output text file
    newick_file = gene_symbol + '_newick.txt'
    newick_path = os.path.join(path_name, newick_file)

    # Write output text file
    with open(newick_path, 'w') as newick_file_object:
        newick_file_object.write(newick)
    print("Newick file created.")

    # Construct phylogenetic tree from alignment and newick values
    tree = Phylo.read(newick_path, "newick")
    tree.ladderize()
    
    # Create figure to display phylogenetic tree in
    figure = plt.figure(figsize=(20,10))
    axes = figure.add_subplot(111)
    Phylo.draw(tree, do_show=False, axes=axes)

    # Define path to phylogenetic tree image
    tree_file = gene_symbol + '_tree.png'
    tree_path = os.path.join(path_name, tree_file)
    
    # Save phylogenetic tree image file
    plt.savefig(tree_path)
    print("Phylogenetic tree image created.")

    return msa 


def main():
    print("""

    You are using...

            -. .-.   .-. .-.   .-. .-.   .  
            ||\|||\ /|||\|||\ /|||\|||\ /|
                *~THE HOMOLOGENATOR~*
            |/ \|||\|||/ \|||\|||/ \|||\||
            ~   `-~ `-`   `-~ `-`   `-~ `-

    """)
    # ASCII art from https://www.asciiart.eu/miscellaneous/dna

    # Email address is required to access web API services
    email_address = input("Please enter a valid email address: ")
    while "@" not in email_address:
        print("Web API services require a valid email address.")
        email_address = input("Please enter a valid email address: ")
    
    query_name, path_name = get_query_name()
    query_sequence = acquire_input()
    blast_output = blast_search(query_name, path_name, query_sequence)
    blast_closest_hits = get_blast_hits(blast_output)
    gene_symbol = find_best_match(blast_closest_hits)
    fasta_record = find_homologues(email_address, gene_symbol)
    homologues_path = format_fasta(path_name, gene_symbol, fasta_record)
    construct_MSA(email_address, path_name, gene_symbol, homologues_path)

# Call to main function
if __name__ == '__main__':
    main()
