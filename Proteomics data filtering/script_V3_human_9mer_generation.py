#!/usr/bin/env python3


#In this script, we are going to loop through the entire human proteome. Then we are going to generate 9mers from every protein. 
#Then these 9mers will be inserted in a 'set', where only unique objects can exist. 
#If the script has previously been used (if you are reading this, then that is the case), 
#then the proteome 9mer data will already have been generated, and this will save some processing time. 

#After establishing the set of unique 9mer, the set will be compared to the predicted 9mer MHC binders.
#These predicted 9mer MHC binders are generated/ predicted based on translated protein from intron retained regions. 


import sys
import re
import time
import math

print("Started the script. Please wait a bit...")

peptideFragmentSet = set() #The set where all the unique human 9mer peptides will be store. It is important to remember this name. 

def peptideFragmentGenerator(AAseq, peptideLen):

	#AA is abreviation for 'amino acid'. The peptideLen argument decides the length of the fragments we generate. 

	if len(AAseq) < peptideLen: #If the amino acid sequence is shorter than desired peptide length (i.e. 9, since we want 9mers), we can't do anything with it. 
		pass #Not anything useful can be done. 

	elif len(AAseq) == peptideLen: #Length of 9, means it is already a peptide fragment.
		peptideFragmentSet.add(AAseq)

	else: #If longer than 9 residues - we can generate multiple 9mers. 
		for i in range(0,len(AAseq),1):
			
			if len(AAseq[i:i+peptideLen]) == peptideLen: #If it is 9 residues long (and not shorter). 

				peptideFragmentSet.add(AAseq[i:i+peptideLen]) #Save the 9mer to the set. 


###

#Opening the human proteome data file:

proteomeFile = "data_big_proteome"

try:
	openedFile = open(proteomeFile,"r")
except:
	raise IOError(f"File does not exist:{proteomeFile}")
	sys.exit()


#Opening the file with the predicted MHC 9mers and gene names + intron retention number. 

treatedMHCdataFileName = "result_MHCbinderResults.txt"

try:
	openedMHCresultFile = open(treatedMHCdataFileName,"r")
except:
	raise IOError(f"File does not exist: {treatedMHCdataFileName}")
	sys.exit()


###

#Check if we have previously generated saved all the unique human 9mer peptide fragments - else, we need to generate it from scratch: 

human9merFileName = "result_uniqueHuman9mers.txt"

try:
	openedHuman9merFile = open(human9merFileName,"r") #If we have already done this before, this file should be present. 

	if openedHuman9merFile:

		print("Since you already have 'result_uniqueHuman9mers.txt' in your directory, this is not the first time you run this script, and the peptide fragment (9mer) set will be made from this data... - approx 19 sec")

		for line in openedHuman9merFile:
			uniquePeptideFragment = line.strip()

			peptideFragmentSet.add(uniquePeptideFragment)

	openedHuman9merFile.close()

except:
	print(f"Bad file name: {human9merFileName} \nThe file containing unique human proteome 9mers has not been generated before - it will be generated now...\n")


	continueInput = input("\t\tYou are about to generate a lot of unique human proteome 9mers. Press 'enter' to continue: ") #This doesn't actually do anything. It is only to have some checkpoints in the code. 

	print("Processing will take approximately 40 seconds.")


	peptideLength = 9 #We want 9mers

	totalEntries = 146566
	#Bonus: There are approximately 146566 entries.


	aminoAcidSequence = ""
	count = 0
	printFlag = False

	for line in openedFile:

		headerDetected = line.strip().startswith(">")


		#Option 1: Encounter accession number.
		#Option 2: Else we encounter amino acid sequence.

		if headerDetected: #Option 1 - Accession number.

			count += 1

			#This is to keep track of the processing progress - how many percent have been covered: 

			if math.floor(count/totalEntries*100) in [10,20,30,40,50,60,70,80,90] and printFlag == False: #Print the percent progress. 
				print(f"{math.floor(count/totalEntries*100)}% - Count: {count} out of {totalEntries}") #146,566 total entries
				printFlag = True
			elif not math.floor(count/totalEntries*100) in [10,20,30,40,50,60,70,80,90]:
				printFlag = False


			if aminoAcidSequence: #The very first time we encounter an accession number, we don't have an amino acid sequence from previous iteration to save from. Therefore, we need to make this check.

				#Save and reset results from previous iteration: 

				peptideFragmentGenerator(aminoAcidSequence,peptideLength) #Generate the 9mer peptide fragments from the amino acid sequence. And adds it to our set 'peptideFragmentSet'.

				aminoAcidSequence = ""


				#Getting the gene name:
				#Important, "[gene=CENPS-CORT]" will give problems, if we use that search string! 
				#We also get trouble at TRBD1 or TRBJ1-1. 

				regularExpression_geneName = re.search(r"\[gene=(\w+.*?)\]",line.strip()) #Regular expression to extract the gene name. 


				if regularExpression_geneName: #If the count equals the count limit, we have not gotten the amino acid. 
					geneName = regularExpression_geneName.group(1)

				else:
					raise ValueError(f"!!!Regular expression could not find gene name on this line:\n{line}")
					sys.exit()
					


		elif not headerDetected: #Option 2 - When encounter amino acid sequence after header. 
			aminoAcidSequence += line.strip() #Simply append the amino acid sequence on this line to the existing amino acid string. 


	#As soon as we reached end of file-reading loop: 
	peptideFragmentGenerator(aminoAcidSequence,peptideLength) #Generate 9mer for the very last amnio acid sequence in the file (since it won't have any subsequent headers). 
	aminoAcidSequence = ""


	### Saving the results in a file - for next time you run the file (to save time). 

	print("Now saving these results in a file named 'result_uniqueHuman9mers.txt', to save time, next time running this code...")

	human9merResultFile = open("result_uniqueHuman9mers.txt","w")

	for peptideFragment in peptideFragmentSet: #Iterating all the 9mers, and writing them into a file.
	
		print(peptideFragment,file=human9merResultFile) #This way, each line in the file should contain a peptide fragment. 

	proteomeFile.close()
	human9merResultFile.close()

###

#Printing results: 

print(f"\nLenght of the set containing all unique 9mers in the human proteome: {len(peptideFragmentSet)}")
print(f"Theoretical amount of 9mer combinations (20^9): {20**9}\n\t\tThat means a percent coverage of: {round(len(peptideFragmentSet)/(20**9)*100,4)}%")
print(f"Size of the list: {sys.getsizeof(peptideFragmentSet)} bytes = {sys.getsizeof(peptideFragmentSet)/10**6} MB")




###########

#Finally! Comparing the 9mers in human proteome with predicted MHC binding 9mers, created by intron retention. 


continueInput = input("\n\t\tDo you want to compare these human proteome sequences with predicted MHC binders caused by intron retention? (press 'enter'): ")
#This doesn't actually do anything. It is only to have some checkpoints in the code.  


print("This processing might take some time...\n")


#Saving the predicted 9mer MHC binders in a dict.   

mhcbinderDict = dict()

for line in openedMHCresultFile:
	lineList = line.strip().split("\t") #This file was made in another script called 'collect_NetMHC.py'. It contains gene name and 9mers that are predicted to bind to MHC. 
	
	geneName, sequence, HLA, bindingAffinity = lineList
	
	HLA = HLA.replace("HLA-","") #Remove the irrelevant part. 

	if not sequence in mhcbinderDict.keys():
		
		mhcbinderDict[sequence] = [geneName,HLA,bindingAffinity]
	else:
		raise ValueError("Sequence was encountered twice in the file!", sequence) 
		sys.exit()

		######################Important notice!!!!!! For some reason, all of the sequences are unique!!! No repeats are present! (Maybe NETMHCPan does some filtering??)

openedMHCresultFile.close()


print("Finished gathering all the predicted MHC binders - The script will now compare the predicted MHC binders with the 9mers found in the human proteome")



noneSelfSimilar_MHCBindingers_Outfile = open("result_None-self_MHC-binders.txt","w")



uniqueMHCbinders_set = set()

count_uniqueIntronRetainedSequences = 0 #This is the important number. If a predicted MHC binder is not already found in the human proteome, then there is a higher chance that it could cause autoimmunity. 

for mhcBindingSequence in mhcbinderDict.keys():

	if len(mhcBindingSequence) == 9:

		if not mhcBindingSequence in peptideFragmentSet:

			noneSelfSimilar_MHCBindingSequence = mhcBindingSequence

			#Write the sequence, gene name, HLA name, and binding affinity:

			geneName,HLA,bindingAffinity = mhcbinderDict[mhcBindingSequence]

			print(noneSelfSimilar_MHCBindingSequence,geneName,HLA,bindingAffinity,sep="\t",file = noneSelfSimilar_MHCBindingers_Outfile) #Saving these none-self-similar, MHC-binding sequences, intron retain variants. They can potentially lead to autoimmunity. 

			count_uniqueIntronRetainedSequences += 1

			uniqueMHCbinders_set.add(mhcBindingSequence) #New - getting a set of all unique MHC binding 9mers, for appendix. 

	else:
		raise ValueError("MHC binder did not have a lenght of 9!", mhcBindingSequence, len(mhcBindingSequence))
		sys.exit()




print(f"\nAmount of unique sequences created by intron retention: {count_uniqueIntronRetainedSequences} out of the {len(mhcbinderDict)} total predicted MHC binders.")
print(f"\tThat equals {round(count_uniqueIntronRetainedSequences/len(mhcbinderDict)*100,4)}%\n") #Printing the percentage. 
print(f"Had the intron retention caused creating of 9mers at complete random, we would expect 100-{round(len(peptideFragmentSet)/(20**9)*100,4)} = {100-round(len(peptideFragmentSet)/(20**9)*100,4)}%\n")

print(f"Amount of unique 9mers in human proteome data: {len(peptideFragmentSet)}")


#print(f"\nUnique 9mers:\n{uniqueMHCbinders_set}\n")
#print(len(uniqueMHCbinders_set))

noneSelfSimilar_MHCBindingers_Outfile.close()


print("\nSee results in 'result_None-self_MHC-binders'...")

