#! /usr/bin/python3




#Importing packages
import sys

#Taking arguments
inputfile = sys.argv[1]
outputfile = sys.argv[2]


# Reading the results:
alreadyseen = set()
normalcds_binders = set()
binders = list()


#First time reading the file. Here, we extract the normal protein 9mers.
infile = open(inputfile, 'r')
for line in infile:
    #Removing excess whitespaces. Might not be needed
    data = line.split()
    if len(data) >= 10:
        for i in range(len(data) - 1, -1, -1):
            if len(data[i]) == 0:
                del data[i]

        #If it's the normal cds without introns, we save the 9-mer binder
        if "_i_" not in data[10]:
            normalcds_binders.add(data[1] + "_" + data[3])

infile.close()






#Here, we get all the information on intron retained sequences
infile = open(inputfile, 'r')
for line in infile:
    #Removing excess whitespaces. Might not be needed
    data = line.split()
    if len(data) >= 10:
        for i in range(len(data) - 1, -1, -1):
            if len(data[i]) == 0:
                del data[i]
    
        #If we haven't seen the gene name before, we save the intron variant and a set containing the peptide. This set will be added upon in later iterations.
        if data[10] not in alreadyseen:
            if data[1] + "_" + data[3] not in normalcds_binders:
                alreadyseen.add(data[10])
                binders.append([data[10], {data[1] + "_" + data[3]}])
    
        #If we've seen this intron retention variant before, we just add to the set:
        elif data[10] in alreadyseen:
            if data[1] + "_" + data[3] not in normalcds_binders:
                binders[[x[0] for x in binders].index(data[10])][1].add(data[1] + "_" + data[3])         #We place  the sequence into the set that is inside the list
infile.close()


#Now the binders list should contain lists of [intron retention name, set of binders not in normal cds]
#Therefore, we can just print it into an output file
outfile = open(outputfile, 'w')
for data in binders:
    print(data[0], "\t".join(list(data[1])), file = outfile, sep = "\t")
outfile.close()