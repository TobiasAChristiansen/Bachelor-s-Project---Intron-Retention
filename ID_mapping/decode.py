#! /usr/bin/python3

#Loading decoding table
decoder = dict()
infile = open("decoding_table.txt", "r")
for line in infile:
    encoded, decoded = line.strip().split()
    decoder[encoded] = decoded
infile.close()


#Decoding the file:
infile = open("collected.fa", "r")
outfile = open("collected_decoded.fa", "w")
for line in infile:
    if line.startswith(">"):
        print(">" + decoder[line[1:5]] + line[5:-1], file = outfile)
    else:
        print(line.strip(), file = outfile)
infile.close()
outfile.close()


