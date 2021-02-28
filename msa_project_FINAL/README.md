# MSA Project
## Sprint_3 Team Coding Project

### User story:
As a bioinformatics trainee, I want to create an API that performs an MSA with as little user input as possible. Clinically valid conservation analyses should be able to use the output files from this API." 
                
### **Introduction**
MSA is a tool that creates Multiple Sequence Alignments (MSA) and analyses protein coding sequences of DNA and their respective amino acid sequences, across different species.

To use this tool the input must be a Refseq accession number that corresponds with a record in Genbank. Accession numbers must begin with the prefix NM_. The NM_ prefix corresponds to protein-coding transcript records in RefSeq.

Users are also required to enter a valid e-mail address and the type of sequences that they wish to align.

Interaction with the MSA tool is through the Swagger UI (user interface). Additional information on the software that the user must have in order use this tool can be found in the requirements.txt file.

### Additional software requirements include:
1. Python 3
2. Biopython v1.78
3. Werkzeug v0.16.1
4. requests v2.25.1

### The tool will perform the following funtions:
1. Generate a dictionary with transcript and protein sequences and IDs. Create directory for outputs.
2. Submit either the nucleotide sequence or amino acid sequence to the BLAST API using Biopython NCBIWWW. Save the .xml file to it's respective directory.
3. Extract relevant information from BLAST XML output. Compile a dictionary of the top 10 most homologous sequences.
4. Generate a .fasta file consisting of homologous sequences in FASTA format by querying RefSeq using the Entrez module from the Biopython package. Save output .fasta file to the same directory that contains the BLAST output .xml file.
5. Construct a multiple-sequence alignment using Clustal Omega. Save the multiple-sequence alignment in the same directory as the BLAST .xml file and .fasta file.

### The ouputs from this tool are:
1. A raw .xml file with MSAs from BLAST
1. A list of the top BLAST hits, including species and percentage sequence identity.
2. A .fasta file containing the sequences (and species of origin) used to construct the multiple sequence alignment.
3. A .aln file containing the multiple sequence alignment.

The output files are saved to the directory in which the MSA is initiated.

### **Limitations of the API**
-BLAST search can take a considerable amount of time depending on server load.\
-Homologous sequences returned by BLAST may include predicted transripts and protein sequences, prefixed with XM_ and XP_.\
-The MSA output by Clustal Omega is not in .fasta format and so may not be appropriate input for downstream tools.\
-This API is not ready for release. Further testing of it's functionality is still required.\ 

*This app was made by Arjun Ryatt, Cristina Garcia Moratalla and Carissa Fletcher-Regez*
