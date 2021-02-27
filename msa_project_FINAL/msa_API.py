from Bio import Entrez  # To get NCBI records and manipulate sequences
from Bio import SeqIO  # To get NCBI records and manipulate sequences
from Bio import Phylo  # To get NCBI records and manipulate sequences
from Bio.Blast import NCBIWWW  # To query NCBI BLAST
import xml.etree.ElementTree as ET  # To enable .xml manipulation
from flask import Flask
from flask_restplus import Api
from flask_restplus import Resource
import os # For management of input and output files
import sys
import requests  # To query web services
import json
import re
import socket
from socket import gaierror
import urllib
from urllib.error import URLError
from urllib.error import HTTPError
from urllib.request import urlopen
from requests.exceptions import HTTPError
from requests.exceptions import ConnectionError
# from requests.exceptions import HTTPError
import logging
import logging.handlers as handlers
import time
from msa_package import get_transcript_info
from msa_package import blastn_search
from msa_package import blastp_search
from msa_package import top10_hits_fasta
from msa_package import clustal_omega_MSA

"""
This API queries Genbank using NM_ numbers from RefSeq. Using NCBI's BLAST tool, the API returns a list of the top 10 most closely related orthologues, as well as their corresponding sequences formated into a .FASTA file. The Clustal Omega API is then used to perform a Multi-Sequence Alignment of the sequences in the .FASTA file.

    Args:
        - Email_address ([string]): Required to access Genbank, BLAST and Clustal Omega API.
        - RefSeq Accession number starting with 'NM_'([string]). Required to query Genbank and BLAST.
        - Sequence_type ([string]): The user can interrogate either the base sequence in the NM_ transcript record or the amino acid sequence in the corresponding protein record.

    Returns:
        - A list of the top BLAST hits, including species and percentage sequence identity.
        - A .fasta file containing the sequences (and species of origin) of the top 10 homologous sequences from the BLAST search, used to construct the multiple sequence alignment.
        - A .aln file containing the the multiple sequence alignment.
"""

# The link is provided to a user that is running the API from the command-line.
print('Please visit this site to use the API: http://127.0.0.1:5000/')

# Create an error log for the API
logging.basicConfig(filename='msa_api_error.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s', level=logging.INFO)

# Assigns the log file to 'logger'
logger = logging.getLogger('API')

# Log with a rotating file-handler. This sets the maximum size of the log to 0.5Mb and allows two additional logs
# The logs are then deleted and replaced in rotation
logHandler = handlers.RotatingFileHandler('msa_api_error.log', maxBytes=500000, backupCount=2)

# We want to minimise the amount of information we log to capturing bugs
logHandler.setLevel(logging.ERROR)
logger.addHandler(logHandler)


# Creates the Exception class called 'RefSeqError'- A generic exception used to capture any exception.
class RefSeqError(Exception):
    pass

# Creates the Exception class called 'GatewayTimeout' for when the API cannot connect to the server.
class GatewayTimeout(Exception):
    code = 504


# Define the application as a Flask app with the name defined by __name__ (i.e. the name of the current module)
# Most tutorials define application as "app", but I have had issues with this when it comes to deployment,
# so application is recommended
application = Flask(__name__)
# Define the API as api
api = Api(app=application)
# define the namespace of the API endpoint
MSA_api = api.namespace('Multiple Sequence Alignment API', description='Search Genbank for record, BLAST sequence, Align top 10 homologous sequences using Clustal Omega')
@MSA_api.route("/<string:Email_address>/<string:Sequence_type>/<string:Ref_Seq_ID>")

#define the parameters of the query_strings
@MSA_api.param("Ref_Seq_ID", "Please enter an 'NM_' RefSeq Accession number that exists in Genbank\n"
                             "> e.g. NM_007298.3")
@MSA_api.param("Sequence_type", "Please enter either 'Protein' or 'Transcript'")
@MSA_api.param("Email_address", "Please enter a valid e-mail address")

# Searches Genbank using the RefSeq NM_ accession number
class MsaClass(Resource):
    def get(self, Ref_Seq_ID, Sequence_type, Email_address):
        # uses Entrez to query the various RefSeq APIs. E-mail is required for each query.
        Entrez.email = Email_address
        
        try:
            # creates a dictionary of the query's NM_ number, corresponding NP_ number and their respective sequences.
            record = get_transcript_info(Ref_Seq_ID)
            # If the user chooses 'Transcript' sequence type, the following if statement is executed.
            if Sequence_type == 'Transcript':
            
                # Pulls the nucleotide sequence from the 'record' dictionary.
                query_sequence = record['cdna']
                # Generates a list of other transcripts in RefSeq with similar identity to the query sequence.
                blastn = blastn_search(query_sequence, Ref_Seq_ID)
                # Creates a .fasta file of the top 10 most closely related transcripts
                fasta_path = top10_hits_fasta(blastn, Ref_Seq_ID, 'nucleotide', '')
                # Performs a multiple sequence alignment of the sequences in the .fasta file.
                clustal_omega_MSA(Email_address, fasta_path, Ref_Seq_ID, Sequence_type, '')
            
                # Returns list of top 10 most homolgous sequences to query sequence.
                return blastn


            # If the user chooses 'Transcript' sequence type, the following if statement is executed.
            elif Sequence_type == 'Protein':

                # Pulls the amino acid sequence from the 'record' dictionary.
                query_sequence = record['translation']['prot']
                # Generates a list of other proteins in RefSeq with similar identity to the query sequence.
                blastp = blastp_search(query_sequence, Ref_Seq_ID, record['translation']['id'])
                # Creates a .fasta file of the top 10 most closely related amino acid sequences.
                fasta_path = top10_hits_fasta(blastp, Ref_Seq_ID, Sequence_type, record['translation']['id'])
                # Performs a multiple sequence alignment of the sequences in the .fasta file.
                clustal_omega_MSA(Email_address, fasta_path, Ref_Seq_ID, Sequence_type, record['translation']['id'])
                
                # Returns list of top 10 most homolgous sequences to query sequence.
                return blastp

        # Exceptions seem to be raised by urllib.error. The baseexception is URLerror which sometimes 'inherits' its error message from socket.getaddrinfo. This prevents a better exception from being made for Gateway Errors.
        except Exception as e:
            # The error message contains the exception and the date and time that it occured.
            message = '%s occurred on %s.' % (e, time.ctime())
            # The exception is stored in the msa_api_error.log according to its identity, date and time, and response.
            logger.exception(message, exc_info=True)
            
            # Raises the exception
            raise RefSeqError(e)

            
            
'''Error handlers'''
# A list of error handlers that response to error messages that may arise if the API is not functioning properly. This could be due to poor development code, internal server errors, errors with the resources or with the way that they are  being querried.

# A generic error handler that produces any exception's error message into the response body.
@application.errorhandler(RefSeqError)
def error_response_handler(e):
    return {'Error message': str(e)}
    
# For when the page or response cannot be found.
@application.errorhandler(404)
def not_found_error_handler(e):
    return {'Error message': 'Error code 404: Requested Endpoint could not be found (maybe due to incorrect parameter?).'}

# For when there is an internal server error.
@application.errorhandler(500)
def default_error_handler(e):
    return {'Error message': 'unhandled error: please e-mail thisisafakeemailaddress@i3hs.co.uk'}

# For when the the Gateway times out to a lack of response from the server.
@application.errorhandler(504)
def gate_way_time_out_error_handler(e):
    return {'Error message': 'Error code 504: Gateway Timout. No response from server (maybe due to no internet connection).'}

# Allows app to be run in debug mode
if __name__ == '__main__':
    application.debug = True # Enable debugging mode
    application.run(host="127.0.0.1", port=5000) # Specify a host and port fot the app

