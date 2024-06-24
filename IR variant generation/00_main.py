#! /usr/bin/python3

import subprocess
import sys



############################################      SETUP      ############################################

#Defining the paths we will use
target_dir = "./"   #Directory containing the scripts
raw_location = target_dir + "00_raw_data/"
introns_sorted_location = target_dir + "01_formatted_data/"
NETMHC_res_location = target_dir + "02_NETMHC_res/"
SB_WB_MHC_location = target_dir + "03_SB_WB_NETMHC/"
new_binders_location = target_dir + "04_new_binders/"
collected_data_location = target_dir + "05_collected_data/"
final_results_location = target_dir + "06_final_results/"



#Creating temporary folders
for directory_path in [raw_location, introns_sorted_location, NETMHC_res_location, SB_WB_MHC_location, new_binders_location, collected_data_location]:
    subprocess.run(["mkdir", directory_path[:-1]])










#Finding all the filenames of the raw data files
filenames = subprocess.run(["ls", raw_location], capture_output = True, text = True)
filenames = filenames.stdout.split()

#Dictionary with shortened names:
filenames_dict = dict()
for i in range(len(filenames)):
    filenames_dict[filenames[i]] = "{0:04}".format(i)



outfile = open("./decoding_table.txt", "w")
for key in filenames_dict:
    outfile.write(filenames_dict[key] + "\t" + key[:-3] + "\n")
outfile.close()








####### Script 01
#Executing a worker job for every file to create fastas containing possible intron retentions for all files
for filename in filenames:
    subprocess.run(["./01_formatting.py", raw_location + filename, introns_sorted_location + filenames_dict[filename]])

#Changing the filenames we're working with to the files in the 01_introns_sorted directory
filenames = subprocess.run(["ls", introns_sorted_location], capture_output = True, text = True)
filenames = filenames.stdout.split()










sys.exit()



####### Script 02
#Running NETMHCpan and NETMHCIIpan
#subprocess.run(["./02_NETMHC.sh", NETMHCpan_res_location])
#
#











####### Script 03
#Fetching strong and weak binders
subprocess.run([target_dir + "03_fetch_data.sh", NETMHC_res_location[:-1], SB_WB_MHC_location[:-1]])

#Finding new filenames by using the bash function ls:
filenames = subprocess.run(["ls", SB_WB_MHC_location], capture_output = True, text = True)
filenames = filenames.stdout.split()

















####### Script 04
#Finding MHC binders created by intron retention
for filename in filenames:
    subprocess.run([target_dir + "04_new_binders.py", SB_WB_MHC_location + filename, new_binders_location + "new_" + filename])

#Getting data for file collection
filenames = subprocess.run(["ls", new_binders_location[:-1]], capture_output = True, text = True)
filenames = filenames.stdout.split()









####### Script 05
#Collecting the data
subprocess.run(["./05_collect.sh", new_binders_location[:-1], collected_data_location[:-1]])










############################################     DECODING     ############################################
#Getting data for file collection
filenames = subprocess.run(["ls", collected_data_location[:-1]], capture_output = True, text = True)
filenames = filenames.stdout.split()

#Flipping the filename encoder dictionary to function as a decoder
filenames_dict = {value:key[:-3] for key, value in filenames_dict.items()}        #the [:-3] is to cut off the file extension


#Looping through the filenames
for file in filenames:

    #Opening the input and output files
    infile = open(collected_data_location + file, 'r')
    outfile = open(final_results_location + file.split("_")[0] + "_" + file.split("_")[2], 'w')

    #Iterating through the lines in the input file
    for line in infile:
        id, sequence = line.strip().split()


        #If the current line has data from intron retention (it should in all cases, but this is to make sure)
        if "_" in id:
            id = id.split("_")   #Splitting it so we have a list containing [id, intron number retained]
            id[0] = filenames_dict[id[0]]    #decrypting the id
            id = "_".join(id)   #Mering the id string again so we have a string containing the correct genename and the intron number it has retained.
        
        #If the current line is by some miracle a complete CDS without introns
        else:
            id = filenames_dict[id] #Just decrypting the ID since we dont have intron retention information in it
            
        print(id, sequence, file = outfile)
    
    #Closing the files after use
    infile.close()
    outfile.close()











############################################     Cleanup     ############################################
#The program that is run below will delete all the temporary folders (every folder except for 06_final_results).
#The 00_raw_data is also removed since the actual starting point is the .gb file in the current directory (./)
#This line can be commented out to keep all the results from the other codes as well
#subprocess.run([target_dir + "07_cleanup.sh", target_dir])







