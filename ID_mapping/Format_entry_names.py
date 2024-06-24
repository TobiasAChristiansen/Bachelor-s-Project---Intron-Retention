#! /usr/bin/python3


#Loading all the data in a single string (only 180MB, so it's managable):
infile = open("all_info.txt", "r")
entries = "".join(infile.readlines())

#Data formatting
entries = entries.split("//\n")[:-1]    #Separating the different entries
entries = [[line for line in entry.split("\n")] for entry in entries]   #Separating the different lines into list indexes



# RC ...TISSUE=TISSUE_NAME;           note: if "{" not in line or "}" not in line. There can be more than one in the first entry, but they will appear separately in later lines.
# OS ...Organism Name (alias)
# OX ...NCBI_TaxID=TAXID;
# GN ...Name=GENENAME;
# PE ...PE_NUMBER: useless description


#Creating a list for all the info
entryinfo_list = list()

#Data extraction
for entry in entries:

    #Creating the data types
    RC = set() #Since there can be multiple tissues, and we only want unique tissues
    OS = ""
    OX = ""
    GN = ""
    ID = ""
    PE = ""
    AC = ""
    SV = ""

    for line in entry:

        #Tissue information
        if line.startswith("RC"):
            if ("{" not in line or "}" not in line) and "TISSUE=" in line:

                #If there is more than one:
                if "," in line:

                    #If the last one contains "and"
                    if ", and " in line:
                        [RC.add(a) for a in line.split("=")[1].split(",")[:-1]]

                        #If the line ends with a semicolon
                        if ";" in line:
                            RC.add(line.split(", and ")[1][:-1])
                        
                        else:
                            RC.add(line.split(", and ")[1])
                    else:
                        if ";" in line:
                            [RC.add(a) for a in line.split("=")[1].split(";")[0].split(", ")]
                        else:
                            [RC.add(a) for a in line.split("=")[1].split(", ")]
                
                #If there is only one tissue
                else:
                    #If it ends in a semicolon
                    if ";" in line.split("=")[1]:
                        RC.add(line.split("=")[1].split(";")[0])    #Also removing the semicolon at the end
                    
                    else:
                        RC.add(line.split("=")[1])
        




        #Organism information. We're only interested in the latin name
        elif line.startswith("OS") and OS == "":     #Only want one per entry
            if " (" in line:
                OS = " ".join(line.split(" (")[0].split()[1:])
            else:
                OS = " ".join(line.split()[1:])





        #Tax-id
        elif line.startswith("OX") and OX == "":     #Only want one per entry
            OX = line.split("=")[1][:-1]
            if " " in OX:
                OX = OX.split()[0]





        #Gene name
        elif line.startswith("GN") and GN == "":     #Only want one per entry
            GN = line.split("=")[1]

            #Accounting for different formats that can appear
            if ";" in GN:
                GN = GN.split(";")[0]
                if " " in GN:
                    GN = GN.split()[0]
            else:
                GN = GN.split()[0]
        




        #Evidence score
        elif line.startswith("PE") and PE == "":     #Only want one per entry
            PE = line.split(":")[0].split()[-1]
        

        #Full name of entry;
        elif line.startswith("DE   RecName: Full="):
            if " {" in line:
                NM = line.split("=")[1].split(" {")[0]
            else:
                NM = line.split("=")[1].split(";")[0]

        elif line.startswith("ID   ") and ID == "":
            ID = line[5:].split()[0]

        elif line.startswith("AC   ") and AC == "":
            AC = line[5:].split()[0][:-1]
        
        elif line.startswith("DT   ") and SV == "":
            if "sequence version" in line:
                SV = line.strip().split()[-1][:-1] #We remove the "." after the sequence version number

    
    #Adding the collected information to the entryinfo list:
    entryinfo_list.append([GN, AC, ID, OX, OS, RC, PE, NM, SV])


#Converting RC to a comma seperated string:
entryinfo_list = [[data[0], data[1], data[2], data[3], data[4], ",".join(list(data[5])).lower(), data[6], data[7], data[8]] for data in entryinfo_list]


outfile = open("prot_info.txt", "w")
outfile.write("\n".join(["\t".join(data) for data in entryinfo_list]))
outfile.close()

