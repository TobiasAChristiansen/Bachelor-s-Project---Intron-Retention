#! /usr/bin/python3






############################################      SETUP      ############################################
#importing packages
import sys

#Fetching the filename recieved as an argument
raw_filename = sys.argv[1]
introns_sorted_filename = sys.argv[2]
genename = introns_sorted_filename.split("/")[-1]




############################################      READING      ############################################
#Opening the file and creating variables to keep track of the data and the filereading
infile = open(raw_filename, 'r')
AAseq = ""
nseq = ""
cds_string = ""
cds_flag = False
AAseq_flag = False
origin_flag = False
desc_flag = False
desc = ""
geneconfirmed = False
pname = ""

#Iterating through the file
for line in infile:

    #Finding markers
    if line.strip().startswith("CDS") and cds_string == "":
        cds_flag = True
    elif line.strip().startswith("/translation=") and cds_string != "" and AAseq == "":
        AAseq_flag = True
    elif line.strip().startswith("ORIGIN"):
        origin_flag = True        
    elif line.strip().startswith("DEFINITION"):
        desc_flag = True
    elif cds_string != "" and not geneconfirmed and not cds_flag:
        if line.strip().startswith("/gene="):
            if line.split("\"")[1] in desc:
                geneconfirmed = True
                pname = line.split("\"")[1]
            else:
                cds_string = ""

    #If the first CDS in the file was found, we collect the intervals:
    if cds_flag:
        cds_string += line.strip().split()[-1]

        #If there are more than one interval
        if cds_string.endswith(")"):
            cds_string = cds_string.split("(")[-1].split(")")[0]
            cds_flag = False
        
        #If there's only 1 interval
        elif "(" not in cds_string:
            cds_string = cds_string.split()[-1]
            cds_flag = False
   
    #If the corresponding AA sequence is found to the first cds sequence, save that
    elif AAseq_flag:
        AAseq += line.strip()
        if AAseq.endswith("\""):
            AAseq = AAseq.split("\"")[1]
            AAseq_flag = False

    #Lastly, we save the nucleotide sequence
    elif origin_flag:
        if line.strip() != "//" and line.strip() != "ORIGIN":
            nseq += "".join(line.strip().split()[1:])
    
    #The description is also saved
    elif desc_flag:
        if line.strip().startswith("ACCESSION"):
            desc_flag = False
        else:
            desc += line.strip()

infile.close()


#Proofreading the nseq
for base in nseq:
    if base.upper() not in {"A", "T", "C", "G"}:
        print("INVALID BASE:", base.upper(), pname)
        sys.exit()


#if no cds was found:
if cds_string == "" or geneconfirmed == False:
    print("UNABLE TO CONFIRM GENE:", pname)
    sys.exit()




############################################      DATA FORMATTING      ############################################
#Extracting the coordinates from the cds_string
#cds = [[int(coord) - 1 for coord in coordset.split("..")] for coordset in cds_string.split(",")]
cds = [[coord for coord in coordset.split("..")] for coordset in cds_string.split(",")]
#Iterating through the cds intervals to format them
for i in range(len(cds)):
    for j in range(2):

        #If the interval has a ">", it means up until but not including. Therefore we subtract two (1 extra since we convert to 0-index)
        if ("<" in cds[i][j] and cds[i][j][:-1].isnumeric()) or (">" in cds[i][j] and cds[i][j][1:].isnumeric()):
            cds[i][j] = int(cds[i][j][1:]) - 2
        
        #if it doesn't have > or < which are the only allowed non-numeric characters
        elif cds[i][j].isnumeric() == False:
            print("INVALID CDS:", cds[i][j], pname, cds_string)
            sys.exit()

        #Else, we just convert to 0 index
        else:
            cds[i][j] = int(cds[i][j]) - 1


#Cutting the sequence at the ends (retains introns between cds regions). We retain the stop codon
nseq = nseq[cds[0][0]:cds[-1][1]].upper()

#Matching the cds coordinates to the new sequence
cds = [[coord - cds[0][0] for coord in coords] for coords in cds]

#Preparing a list for the amino acid sequences. We will start by adding the one for the cds without introns
sequences = []

#Fixing file extension. Since we now have fasta files, we use .fa
introns_sorted_filename += ".fa"





############################################      COMBINATIONS, TRANSLATIONS AND WRITING      ############################################

#Here is a codon translation table (stop codon = -):
codontable = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TGC': 'C', 'TGT': 'C', 'TGG': 'W',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R',
    'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A',
    'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TAA': '-', 'TAG': '-', 'TGA': '-' 
}



#First, we iterate over all the introns (number of cds intervals - 1) unless there is 0
if len(cds) - 1 != 0:
    for intron_index in range(len(cds) - 1):

        #Merging the marked cds regions before and after the intron
        all_prev_exons = "".join([nseq[cds[i][0]:cds[i][1] + 1] for i in range(intron_index + 1)])
        all_next_exons = "".join([nseq[cds[i][0]:cds[i][1] + 1] for i in range(len(cds) - 1, intron_index, -1)])


        #Calculating the codon number (whole codon number) of the exons before and after
        AA_len_prev_exons = int(len(all_prev_exons) / 3) #Whole codons before the intron


        #Creating the sequence with intron retention
        current_nseq = all_prev_exons + nseq[cds[intron_index][1] + 1:cds[intron_index + 1][0]] + all_next_exons
        

        #Now that we have the nucleotide sequence for this particular intron retention, we begin translating it until we reach a stop codon
        #We read in codons
        current_AAseq = ""
        for i in range(0, len(current_nseq) - 2, 3):
            codon = current_nseq[i:i+3]
            AA = codontable[codon]

            #Checking if the current codon is a stop codon or not:
            if AA != "-":
                current_AAseq += AA
            else:
                break
        
        #We remove all but 8 whole codons before the intron retention, since it's not relevant information
        #We only do this if there is more than 8 codons before the intron:
        #if AA_len_prev_exons > 8:
        #    current_AAseq = current_AAseq[AA_len_prev_exons - 9:]   #-9 because we have "length", which is not 0-index


        #Checking if above 8 amino acids were found. If not, then we do not include it:
        if len(current_AAseq) >= 9:
            sequences.append([">{}".format(genename + "_i_" + str(intron_index + 1)), current_AAseq])


#If there's no introns to retain        
else:
    sys.exit()




#Terminating the program if no suitable sequences was found so as to not get empty files
if len(sequences) == 0:
    sys.exit()




############################################      WRITING      ############################################
#Writing the fasta file:
outfile = open(introns_sorted_filename, 'w')

#First we iterate through all the sequences:
for sequence in sequences:

    #We print the header and the sequence in fasta format
    print(sequence[0], file = outfile)

    #Printing the sequence data
    for i in range(0, len(sequence[1]), 60):
        print(sequence[1][i:i+60], file = outfile)

#Closing the file after use
outfile.close()

