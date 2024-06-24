#! /usr/bin/python3



#Loading humanatlas data:
humanatlas_gn = set()
infile = open("humanatlas_blood_brain.txt", "r")
for line in infile:
    humanatlas_gn.add(line.strip())
infile.close()


#This code generates a file containing only the fasta entries that are in specified tissues.
tissues = ["blood", "brain"]


#Loading the information on protein in the specified tissue
info_dict = dict()
infile = open("prot_info.txt", "r")
for line in infile:
    GN, AC, ID, OX, OS, RC, PE, NM, SV = line.strip().split("\t")
    for tissue in tissues:
        if tissue.lower() in RC.split(",") or GN in humanatlas_gn:
            info_dict[GN] = [AC, ID, OX, OS, RC, PE, NM, SV]
            break #if one is there, then we don't need to look for the other
infile.close()


#Printing stats
print(len(info_dict), tissues, "proteins were found")


#Reading all the .fa information:
infile = open("collected_decoded.fa", "r")
alldata = "".join(infile.readlines())
alldata = [[line for line in entry.split("\n")] for entry in alldata.split(">")]   #Splitting the information up so every index in the outer list is a list of the lines for a specific entry

#Finding the header coordinates that should be included:
accepted_indexes = list()
for i in range(len(alldata)):
    currentgenename = alldata[i][0].split("_")[0]
    if currentgenename in info_dict:
        accepted_indexes.append(i)



#Writing a new file containing only the entries that fit the tissue criteria
outfile = open("humanatlas_Uniprot_collected_{0}.fa".format("_".join(tissues)), "w")
unique_counter = set()
for index in accepted_indexes:
    #Fetching the genename and the intron retention variant
    current_intron_variant = alldata[index][0]
    currentgenename = current_intron_variant.split("_")[0]
    intronnumber = current_intron_variant.split("_")[-1]
    unique_counter.add(currentgenename)

    #Printing the header in format >XX|YYYYYY|ZZZZ_HUMAN FREE-TEXT OS=organismname OX=taxid GN=genename
    current_info = info_dict[currentgenename]
    print(">ir|{0}|{1} {2} OS={3} OX={4} GN={5} PE={6}, SV={7}".format(current_info[0] + "_i_" + intronnumber, current_info[1], current_info[6], current_info[3], current_info[2], current_intron_variant, current_info[5], current_info[7]), file = outfile)

    #Printing sequence data:
    for seq in alldata[index][1:-1]:
        print(seq, file=outfile)
outfile.close()

print(len(unique_counter), "uniqe genes were included")
print(len(accepted_indexes), "intron retention variants in total")




#Numbers for blood and brain
# Human atlas + PDB definitions:    3106 proteins found, 2909 proteins included, and 37333 intron retentions included
# PDB definitions:                  2652 proteins found, 2493 proteins included, and 33456 intron variants in total
# Human atlas definitions:          1052 proteins found, 964  proteins included, and 10916 intron variants in total
#
#
#
#
#
#
#