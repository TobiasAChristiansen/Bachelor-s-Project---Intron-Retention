#!/usr/bin/env python3

import sys
import re
import time

#Getting the read file name:

if not len(sys.argv) >= 2: #You can pass 'data_short-NETMHCPan.txt' as the second argument. 
	userInput = input("Write the name of the file containing the NetMHCPan data - e.g. 'data_res_collected.txt': ")
	fileName = userInput
else:
	fileName = sys.argv[1]

#Opening file:

try:
	openedFile = open(fileName,"r")
except IOError as errText:
	print(errText)
	sys.exit()


#Opening the file where we write the results.

resultWriteFileName = "result_MHCbinderResults.txt"

try:
	writeFile = open(resultWriteFileName,"w")
except IOError as errText:
	print("Could not open/ create a file where the results should be store ('result_MHCbinderResults.txt'). ",errText)
	sys.exit()


#Opening the file with the decoding table for the encrypted names:

tableFileName = "data_decoding_table.txt"

decodingEncryptionDict = dict()

try:
	openedTranslationTableFile = open(tableFileName,"r")
except IOError as errText:
	print(errText)
	sys.exit()

for line in openedTranslationTableFile:
	encryptedName, geneName = line.strip().split(sep="\t")

	decodingEncryptionDict[encryptedName] = geneName #Saving the gene name under the encrypted name, in a dict. 

openedTranslationTableFile.close()


def deencryptGeneName(name): #This function will be for translating encrypted gene names, to their actual gene name, using the dict we just established above. 

	#We don't know if the input has '0001_i_12' or not - so the code should be able to handle both.
	intronInformation = ""

	intronName_regex = re.search(r"(\d{4})(_i_\d+)",name)

	if intronName_regex:
		name = intronName_regex.group(1) #The input will only be the gene name number, and not the e.g. '_i_18' part. This way we can search in the decode table file.  
		intronInformation = intronName_regex.group(2)

	if name in decodingEncryptionDict.keys():
		
		if not intronInformation: #If the name is of this format: 0012
			return [decodingEncryptionDict[name]] #Return ['AATK']

		else: #If the name is of this format: 0012_i_18
			return [decodingEncryptionDict[name], intronInformation] #Return ['AATK','_i_18']

	else:
		raise ValueError(f"{name} not in decoding file")
		sys.exit()


#Reading file: 

dataDict = dict() #We are interested in storing: Peptid sequence (9mer), protein name (remember that it is encrypted), HLA allele and strong- / weak binder. 

print("\nThe script will start reading the defined NetMHCPan result file - takes approximately 20 seconds")

for line in openedFile:
	
	if line.strip().endswith("<= WB") or line.strip().endswith("<= SB"): #We are only interested in lines with strong- and weak binders

		#Important:
		#The file is not tab-seperated, so we need to use regex to extract everything...  

		#All protein letters are: ACDEFGHIKLMNPQRSTVWYX.
		peptide_regex = re.search(r"\s([ACDEFGHIKLMNPQRSTVWYX]{9})\s",line.strip()) #We need to e.g. find string: ' AMARKLLPG '.  
		#Note, since it is a 9mer, we want to search for a length of 9. 
		#Note, we don't want it to match '-', since NetMHCPan might sometimes fill that into a sequence.

		if peptide_regex:
			if not "-" in peptide_regex.group(1) and len(peptide_regex.group(1)) == 9:
				peptideSequence = peptide_regex.group(1) #Note, this should only fetch the content from the 4'th column. Since the data isn't tab-seperated, we can't use split(), to get to it. 
			else:
				print("\tBad peptide input on this line! - Skipping (err type 1)...", line,"\n")
				peptide_regex = None
				time.sleep(5)
				continue #Something is faulty with this 9mer sequence. 
		else:
				print("\tBad peptide input on this line! - Skipping (err type 2)...", line,"\n")
				peptide_regex = None
				time.sleep(5)
				continue #Something is faulty with this 9mer sequence. 


		geneName_regex = re.search(r"(\d{4}_i_\d+)",line.strip()) #We need to extract e.g. "0010_i_13"

		if geneName_regex:
			encryptedGeneName = geneName_regex.group(1)

		HLA_regex = re.search(r"(HLA-\w+\*\w+\:\w+)",line.strip()) #Extracting the HLA name, e.g. of the format "HLA-A*02:01"

		if HLA_regex:
			HLA_name = HLA_regex.group(1)


		WB_regex = re.search(r"<= WB",line.strip())
		SB_regex = re.search(r"<= SB",line.strip())


		if peptideSequence and encryptedGeneName and HLA_name and len(peptideSequence) == 9:
			#Check if we extracted all the data. 

			if WB_regex:
				bindingAffinity = "WB"
			elif SB_regex:
				bindingAffinity = "SB"


			#Storing the data in the dict. 

			if not peptideSequence in dataDict.keys(): #If result is not already stored in the dict. 
				dataDict[peptideSequence] = [encryptedGeneName,bindingAffinity,HLA_name] #Store peptideSequence (key) with list (value). 

			else: #If we already have encountered this sequence before.
				#In that case, we need to compare and see if the gene name, binding affinity and HLA name are still the same. 


				if not [encryptedGeneName,bindingAffinity,HLA_name] in dataDict[peptideSequence]: #'dataDict[peptideSequence]' will give us the list of all the currently stored gene names (encrypted), binding affinities and HLA types. 
					
					#If not either the gene name, binding affinity or HLA type has been previously observed for this specific type of sequence, we want to note it down.
					#Since this could e.g. mean that a certain amino acid 9mer sequence binds to different types of HLA, or that the sequence occurs in multiple genes - and this might be useful info. 

					#print(dataDict[peptideSequence])

					dataDict[peptideSequence].extend([encryptedGeneName,bindingAffinity,HLA_name]) #The extend function allows us to append the 3 values the pre-existing list, but by appending them as individual strings and not a list of strings. 

					#print(dataDict[peptideSequence])

					#time.sleep(0.5)



		else:
			raise ValueError("Line is missing either: core sequence (or doesn't have a length of 9), (encrypted) gene name or HLA name:\n",line.strip())
			print(peptideSequence)
			sys.exit()


		#Resetting variables for next iteration.

		geneName_regex = None
		peptide_regex = None
		HLA_regex = None
		WB_regex = None
		SB_regex = None
		encryptedGeneName = ""
		peptideSequence = ""
		HLA_name = ""
		bindingAffinity = ""



openedFile.close()


print(f"Amount of predicted MHC binders: {len(dataDict.keys())}")

seeMoreOfDictPrompt = input("\n\nDo you want to see entire dict + write the results in a file? (type 'enter'): ")


for sequence, resultList in dataDict.items():

	#print(resultList)

	#time.sleep(0.5)

	#print(f"9mer sequence: {sequence}")

	for i in range(0,len(resultList),3): #Format of result list might be: [gene name 1, WB, HLA type 1, gene name 2, SB, HLA type 2, gene name 3, SB, HLA type 3]. So we need to iterate over every 3'rd element. 

		encrypGeneName = resultList[i]
		bindaff = resultList[i+1]
		HLA = resultList[i+2]

		#Translating encrypted gene name. 

		geneNameList = deencryptGeneName(encrypGeneName) #The translator function will either return ['AATK'] or ['AATK','_i_18']
		#Although it should only return of the format: ['AATK','_i_18']. It is still good to take some safety measures. 

		if len(geneNameList) == 2: #Intron-retained variant
			geneName, intronNameInfo = geneNameList[0], geneNameList[1]
			outputGeneName = geneName+intronNameInfo

		elif len(geneNameList) == 1: #Normal variant
			geneName = geneNameList[0]
			outputGeneName = geneName


		#print(f"\n\tGene name with intron retained variation: {outputGeneName}\n\tBinding affinity: {bindaff}\n\tHLA allele: {HLA}")



	#Writing 'important' results in a file, so it can later be compared to e.g. 9mers in the human proteome. 


	#V3: 'extraResultString' was a stupid concept which I removed. 

	print(outputGeneName,sequence,HLA,bindaff,sep="\t",file=writeFile) #Will be stored in a file named 'result_MHCbinderResults.txt'. 

	#print(f"{outputGeneName}\t{sequence}\t{extraResultString}")


writeFile.close()


print(f"Done writing results in '{resultWriteFileName}'...")

