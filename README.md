# cris_assessment
We have to create a tool which creates Multiple sequence alignments.
The input must be a transcript reference sequence starting with NM_ 
Use BioPython to pull back the genbank record for the sequence. Extract the sequence and the Translation ID (NP_) and translation sequence. These MUST come from the genbank record.
Create a blast search for both the transcript sequence and the protein sequence to identify similar sequences in multiple organisms.
Perform a Clustal alignment for both the transcript sequences and the protein sequences.
