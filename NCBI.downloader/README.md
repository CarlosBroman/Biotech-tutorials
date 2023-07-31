# Downloading NCBI sequences in Batch using Biopython

This code can be used to download multiples Fasta sequences from the NCBI using Python.

## Usage

Please remember that there a certain restrictions on downloading data from the NCBI.

In this code we:
 - Provide a txt file with species names, one per line.
 - Create a folder to store our sequences.
 - Download the sequences for a specific gene for all our species.
 - Merge all the files in a single fasta file containing a single sequence per species.

This code can be adapted to retrieve sequences for any species list and any gene.


## Output

The output will be a single fasta file containing a single sequence for each specie.
