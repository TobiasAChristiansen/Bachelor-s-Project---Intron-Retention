#! /usr/bin/python3

from Bio import Entrez, SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
import sys
import os
import time


Entrez.email = "Write_e-mail_here"


#First, we get the gene names and entry id (of the first occurance of each gene):
#We start by opening the file of nucleotide accession numbers for the refseq sequences
infile = open("./LRG_RefSeqGene.txt", 'r')

#Then we create data structures to hold the information.
genenamesfound = set()        #To only take one unique genbank file for each gene (faster lookup than lists, since sets are hashed)
genenames = list()            #To keep track of the order in which they are stored
ids = list()                  #Here we keep the nucleotide database accession numbers

#Now we iterate through the input file
for line in infile:

    #We don't include the first line of the file, since it just explains the layout. The code is hard-coded for this layout
    if not line.startswith("#"):

        #We're only interested in column 3 and 4 (index 2 and 3)
        usefuldata = line.strip().split()[2:4]

        #If if's a genename we haven't seen before, we extract the data
        if usefuldata[0] not in genenamesfound:
            genenamesfound.add(usefuldata[0])
            genenames.append(usefuldata[0])
            ids.append(usefuldata[1])

#Closing the input file after use
infile.close()





#Looping through every genename / accession number we extracted to download the files from NCBI. Print statements are used to track progress since it takes around 5 hours
for i in range(len(genenames)):
    print("Processing", genenames[i])

    #Getting the enty from the internet in the correct format
    handle = Entrez.efetch(db='nuccore', retmode='text', rettype='gb', id=ids[i])
    print("Handle found")

    #There should only be a single record, but just in case, we iterate over it 
    for record in SeqIO.parse(handle, "genbank"):

        #By using SeqIO.write, we can write the genbank entry into a file. This file is put in 00_raw_data with the genename as filename and file extension .bg
        SeqIO.write(record,"./00_raw_data/" + genenames[i] + ".gb", "genbank")
        print("Record written")
        print()
    
    #Sleep is added for ensuring we are not throttled because of too frequent requests.
    time.sleep(1)