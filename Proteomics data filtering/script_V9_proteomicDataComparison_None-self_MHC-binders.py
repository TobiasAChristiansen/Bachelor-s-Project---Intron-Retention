#!/usr/bin/env python3

#Date: 20/04-2024

import sys
import re
import time


print("Reading results from previous analysis - Getting none self-similar MHC binding 9mer sequences...")


noneSelfMHCbinderDict = dict() #If any of the MHC binding sequences (which have undergone a self-similarity check, to ensure they are not found in the human proteome) are found in the proteomics data, it will be HUGE, since that means that intron-retained variant will be present in actual proteomic data! 

noneSelfMHCbinderFile = open("result_None-self_MHC-binders.txt", "r")

for line in noneSelfMHCbinderFile:

	lineList = line.strip().split(sep="\t")

	noneSelfMHCbinder_9mer, geneName, HLA, bindingAffinity = lineList

	try:
		noneSelfMHCbinderDict[noneSelfMHCbinder_9mer] #Look-up.

		#We expect the look-up to fail. If it doesn't we want to raise an error: 
		print("Error! Did not expect to find the same MHC 9mer sequence in the file!")
		sys.exit()

	except: #If look-up fails.
		
		noneSelfMHCbinderDict[noneSelfMHCbinder_9mer] = (geneName, HLA, bindingAffinity)


	#noneSelfMHCbinderSet.add(noneSelfMHCbinder_9mer)

noneSelfMHCbinderFile.close()



#Correction - will later be converted to sets: 
proteomic_Intrapolated_MHCbinderMatchList = list() #Will be used for not repeating the same results in the output, since some fragments will be repeat a bunch of times in the same gene. 

proteomic_Extrapolated_MHCbinderMatchList = list() #When we have fragments shorter than 9mers, we need to extrapulate sequence data based on the expected sequence. 




print("Lenght of the set containing the none self-similar MHC binding 9mer sequences:", len(noneSelfMHCbinderDict))

print("\nPreparing to read proteomic data...\n")


#proteomicFileName = "data_new_only_ir_2018.tsv"
#proteomicFileName = "data_new_only_ir_2022.tsv"

#V9:
proteomicFileName = "data_new_only_ir_2022.tsv" #"BREG22_data_full.tsv"


proteomicDataFile = open(proteomicFileName,"r")



#noneSelfMHCbinder_proteomicFormat_resultOutfile = open("result_FilteredNoneSelfMHCbinder_presentInProteomicData_proteomicFormat_2018.tsv","w") #Saving only the lines where we get a hit! 

noneSelfMHCbinder_proteomicFormat_resultOutfile = open("result_FilteredNoneSelfMHCbinder_presentInProteomicData_proteomicFormat_2022.tsv","w") #Saving only the lines where we get a hit! 


#Outdated - I no longer create 2 files to analyse essentially the same thing: 
#tooShort_proteomicFormat_resultOutfile = open("result_tooShort_proteomicFormat_2022.tsv","w") #Saving the lines where the sequences are shorter than 9 residues. 





###################################################################### MAKING A GENENAME TO SEQUENCE DICT: ##################################################

#This will be usefull for when we encounter sequences below a length of 9 - in that case we need to compare the fragment to peptide sequence it is expected to have. 
#In other words, we make a look-up/ extrapolate to get the missing pieces of peptide. 

geneToSequenceDict = dict()

intronFile = open("data_humanatlas_Uniprot_collected_blood_brain.fa","r") #This is the file Tobias made with all the intron variants - brain only. 

header = None
sequence = ""

headerCount = 0

for line in intronFile:


	if line.strip().startswith(">"): 

		headerCount += 1

		#Expect header of the format:
		#>ir|O00429_i_20|DNM1L_HUMAN Dynamin-1-like protein OS=Homo sapiens OX=9606 GN=DNM1L_i_20 PE=1, SV=2
		
		if header == None:
			header = line.strip()

		elif header: #This is not the first time encountering a header. We need to save the sequence from the previous header. 

			headerList = header.split(sep="|")

			geneName = headerList[1]

			if len(geneName) <= len("Q13685"): #Just a random gene name I picked, missing the intron part. All of the names should have the same length. Should have been: Q13685_i_2.
				raise ValueError("Too short gene name!", geneName)
				sys.exit()


			try:
				geneToSequenceDict[geneName] #There should NOT already be an entry in the dict for this gene name! It could still be present in other genes tho.  

				raise ValueError("The same gene name should not be 2 times in the file!", geneName, geneToSequenceDict[geneName]) #If the previous don't trigger an exception, the we want to raise an error. 
				sys.exit()

			except:

				#print(geneName, sequence)
				#time.sleep(1)


				if sequence:
					geneToSequenceDict[geneName] = sequence #Save the sequence in the dict. 

					sequence = "" #remember to clear it, after we saved the sequence string. 

					header = line.strip() #New header is the one found on the current line. 

				else:
					raise ValueError("No sequence!", line)
					sys.exit()

	else:
		if line == "":
			print("Empty line...")
		else:
			sequence += line.strip() #Otherwise append to the 'sequence' string


#Repeating the code above, for the last line. 

####
if header:
	headerList = header.split(sep="|")
	geneName = headerList[1]
	if len(geneName) <= len("Q13685"): 
		raise ValueError("Too short gene name!", geneName)
		sys.exit()
	try:
		geneToSequenceDict[geneName]
		raise ValueError("The same gene name should not be 2 times in the file!", geneName, geneToSequenceDict[geneName]) 
		sys.exit()
	except:
		if sequence:
			geneToSequenceDict[geneName] = sequence #Save the very last sequence and header to the dict. 
		else:
			raise ValueError("No sequence!", line)
			sys.exit()
####

intronFile.close()




#Printing results

print("Amount of gene to sequence dict entries (amount of headers):", len(geneToSequenceDict))
print("Should be equals:", headerCount)

if len(geneToSequenceDict) != headerCount:
	raise ValueError("The check did not match - insufficient dict was established!")
	sys.exit()


###################################################################### Getting patient info - survial time and age: #################################

patientInfoFile = open("List_Brains_24.04.24.csv", "r")

#Patient ID to diagnosis:
patientID_age_dict = dict()

ID_to_diagnosisDict = dict()

id_to_patientSex_dict = dict() #V8

for line in patientInfoFile:

	lineList = line.strip().split(sep=";")
	#This should give:
	#No.	K-number	Diagnosis	Age	Sex	Disease duration

	if lineList[0] == "No.": #First line can be skipped
		continue

	patientID = lineList[1].lstrip('kK') #Patient id (k-number) without the 'k' or 'K'.

	diagnosis = lineList[2]


	age = lineList[3]

	sex = lineList[4]

	diaseaseDuration = lineList[5]

	
	try:
		ID_to_diagnosisDict[patientID] #look-up, which should fail. 

		print("Error! Same ID is present more than 1 time!")
		sys.exit()

	except: #When look-up fails. 
		ID_to_diagnosisDict[patientID] = diagnosis

	#V7: 
	if not patientID in patientID_age_dict:
		patientID_age_dict[patientID] = age
	else:
		raise ValueError("Dict error - patient already present", patientID)
		sys.exit()

	if not patientID in id_to_patientSex_dict:
		id_to_patientSex_dict[patientID] = sex
	else:
		raise ValueError("Dict error - patient already present", patientID)
		sys.exit()



###################################################################### Getting statistic HLA frequency data: #################################

HLAfreqFile = open("HLA_data/result_HLAfrequencies_1field.txt","r") #Note: My file is in a directory. 

HLAfreqStatistic_Dict = dict()
#HLA -> Diagnosis -> (freq, std)

for line in HLAfreqFile:

	if line.strip().startswith("#"): #Skip first line. 
		continue


	lineList = line.strip().split(sep="\t")

	HLA, diagnosis, frequency, standardDeviation, EUfreq, EUstd = lineList

	diagnosis = diagnosis.split(sep=";")
	frequency = [float(x) for x in frequency.split(sep=";")]
	standardDeviation = [float(x) for x in standardDeviation.split(sep=";")]


	#A line could e.g. look like this:

	#A*11	MSA;PSP;PD;CTRL	0;0.09091;0.05;0.0625	0.0;0.08428;0.06092;0.06052

	for i, diagnosis_i in enumerate(diagnosis):

		try:
			HLAfreqStatistic_Dict[HLA][diagnosis_i] #Look-up. We expect it to fail.

			print("Error! This entry does already exist in the dict!", HLA, diagnose_i, HLAfreqStatistic_Dict[HLA][diagnosis_i])

		except:

			try: #If an entry already exists this HLA with other diagnoses
				HLAfreqStatistic_Dict[HLA][diagnosis_i] = (frequency[i],standardDeviation[i])

			except: #Otherwise, create a new entry - a dict in a dict. 
				HLAfreqStatistic_Dict[HLA] = {diagnosis_i: (frequency[i],standardDeviation[i])}



#print(HLAfreqStatistic_Dict)


input("Are you ready to read proteomic data? - [Press 'enter']: ")

####################################################################################################################################

#Patient specific HLA subtypes:

#This will be important, when we later examine how many of the patients actually have the HLA subtype, 
#Which bind the matched MHC binding 9mer.


patientHLAoverviewFile = open("HLA_data/result_HLA_patient_overview.txt","r") #Note: In another folder. As part of the HLA freq. study. 

patientSpecificHLAdict = dict()

for line in patientHLAoverviewFile:

	lineList = line.strip().split(sep="\t")

	print(lineList)

	brainSampleID = lineList[0]

	diagnosis = lineList[1]

	HLAsubtypeList = lineList[2].split(sep=";") #All the subtypes are semicolon seperated. 


	if brainSampleID in patientSpecificHLAdict:
		raise ValueError("Sample ID already in dict!")
		sys.exit()

	else:
		patientSpecificHLAdict[brainSampleID] = (diagnosis, HLAsubtypeList)




####################################################################################################################################
#
####################################################################################################################################








print("\nNow reading proteomic data...\n")

###################################################################### READING THE PROTEOMIC DATA: ##################################################

outfile_significantResults = open("result_HLAsignificantMatches.txt","w")

print("#MHC binding non-self 9mer sequence","Gene name", "HLA","Diagnosis","Brain sample ID","Weak- or strong binder","HLA freq. (diagnosis group)","HLA std (diagnosis group)","Result quality","Age","Brain region","Sex",sep="\t",file = outfile_significantResults)



significantMatchCount = 0
totalMatchCount = 0
lineCount = 0

peptide_minLength = None
peptide_maxLength = None

#V6:
modifiedPeptideDict = dict()
modifiedPeptide_flag = False

#V7:
modifiedSequenceFile = open("result_modifiedIntronSequence.txt","w")
print("#Gene name","Patient ID","Diagnosis","Sequence (only modifications)",sep="\t",file = modifiedSequenceFile)

intronSequenceLength_diagnosisDict = dict() #Total amount of sequence length pr. patient group (not divided)
diagnosis_size_dict = dict() #E.g. size of CTRL group is 8 patients. 

patientID_set = set()

#V8:
diagnosis_intronFragmentCount_dict = dict()
geneName_MHCnonSelfMatchCount_dict = dict() #Collecting all of the gene names' associated with the 'relevant MHC non-self match'. 

allUniqueIntronGenes_set = set()
nonselfMHCmatch_intronGene_set = set()



for line in proteomicDataFile:

	lineCount += 1

	#The format of a line will be:
	#Not Defined	20230103_EXPL2_Evo0_EaS_OOE_DIA42m_BREG18-RR073-CBX_022	1	O00429_i_20	Dynamin-1-like protein	DNM1L_HUMAN	79.6%	False	0	5730140.5	_SSVLESLVGR_.2	3.6386417358168286E-07	2244988

	#It is tab-seperated (\t). 

	lineList = line.strip().split(sep="\t")

	patientID = lineList[1].split(sep="-")[1]
	geneIntronInfo = lineList[3] #Aka. gene name. UniProt format. 
	sequence = lineList[10].lstrip("_")

	brainRegion = lineList[1].split(sep="-")[-1].split(sep="_")[0] #V7




	patientID_set.add(patientID)

	#An example of how the sequence can look like: '_YYSGLIYTYSGLFC[Carbamidomethyl (C)]VVINPYK_.3'

	#V6 - Detect modifictaion in sequence: 
	if "[" in sequence:
		modifiedPeptide_flag = True

	#V5 - Remove irrelevant stuff on right side: 
	regex_sequence = re.sub(r"\_\.\d+","",sequence)


	#V7 - Diagnosis look-up from patient ID: 
	if patientID == "BLK": #Blank patient entry.
		diagnosis = "N/A"
	else:
		try:
			diagnosis = ID_to_diagnosisDict[patientID]
		except:
			diagnosis = "unclear diagnosis" #Some of the patients such as '1400', are difficult diagnose, and were therefore not included in our provided patientlists. 

	if patientID in patientID_age_dict:
		age = patientID_age_dict[patientID]
	else:
		age = "N/A"

	#V8:
	if patientID in id_to_patientSex_dict:
		patient_sex = id_to_patientSex_dict[patientID]
	else:
		patient_sex = "N/A"

	#V6 - Saving the modified residues:  
	if modifiedPeptide_flag == True:
		modifiedSequence = regex_sequence

		#Quick statistics - Save the modified intron residue, under the given diagnosis:
		if diagnosis in modifiedPeptideDict:
			modifiedPeptideDict[diagnosis].append(modifiedSequence)
		else:
			modifiedPeptideDict[diagnosis] = [modifiedSequence]
		
		#Saving the modified peptides in a seperate file:

		print(geneName,patientID,diagnosis,modifiedSequence,sep="\t",file = modifiedSequenceFile)

		modifiedPeptide_flag = False


	regex_sequence = re.sub(r"\[[^\]]*\]","",regex_sequence)



	#Converting sequences to non-modified versions: 
	if regex_sequence:
		sequence = regex_sequence
	else:
		print("Bad RE match:", sequence)
		sys.exit()


	if not diagnosis in intronSequenceLength_diagnosisDict:
		intronSequenceLength_diagnosisDict[diagnosis] = len(sequence)
	else:
		intronSequenceLength_diagnosisDict[diagnosis] += len(sequence) #Calculating the amount of intron sequence, for each group (not divided by group size). 


	if not diagnosis in diagnosis_size_dict:
		diagnosis_size_dict[diagnosis] = {patientID}
	else:
		diagnosis_size_dict[diagnosis].add(patientID)



	if peptide_maxLength == None:
		peptide_maxLength = len(sequence)
	if peptide_minLength == None:
		peptide_minLength = len(sequence)

	if len(sequence) > peptide_maxLength:
		peptide_maxLength = len(sequence)
	if len(sequence) < peptide_minLength:
		peptide_minLength = len(sequence)



	#V8:

	if not diagnosis in diagnosis_intronFragmentCount_dict:
		diagnosis_intronFragmentCount_dict[diagnosis] = 1
	else:
		diagnosis_intronFragmentCount_dict[diagnosis] += 1 #Add +1 to the fragment count. 


	#V8 - get a set of intron variants:
	allUniqueIntronGenes_set.add(geneIntronInfo)




	if len(sequence) >= 9:

		for i in range(0,len(sequence)):

			movingWindow = sequence[i:i+9] #Since we are looking for 9mers! 


			if len(movingWindow) != 9: #The moving window should not have a length below 9. 
				break

			else:


				if movingWindow in noneSelfMHCbinderDict: #The dict contains the 9mer sequences as keys. 

					proteomic_Intrapolated_MHCbinderMatchList.append(movingWindow) #Save to this set, so we don't get repeated results. 
					
					#V3 - include HLA type and binding affinity:
					geneName, HLA, bindingAffinity = noneSelfMHCbinderDict[movingWindow] #Get the important information from the found MHC-binding non-self 9mer. 

					#Important: 'geneName' is of the GenBank accession number format (!)
					#And 'geneIntronInfo' is UniProt accession number! 


					HLA_1field = re.sub(r":\d{2}$","",HLA) #Remove the last 2 digits and :. A*01:02 -> A*01

					#HLAfreqStatistic_Dict:
					#HLA -> Diagnosis -> (freq, std)
					if diagnosis == "N/A" or diagnosis == "unclear diagnosis": #If the diagnosis is unclear, we can't assign it a statistical value
						HLA_frequency, HLA_std = ("N/A", "N/A")

					else: #Else, we can make a look-up to get the statistical values. 
						
						try:
							HLA_frequency, HLA_std = HLAfreqStatistic_Dict[HLA_1field][diagnosis] #V9 - The look-up can now sometimes fuck up. E.g. in the case of 'B*58'. 
						except:
							HLA_frequency, HLA_std = ("N/A", "N/A")


					print(f"{line.rstrip()}\t{movingWindow}\t{geneName}\t{HLA_1field}\t{bindingAffinity}\t{diagnosis}\t{HLA_frequency}\t{HLA_std}\tIntrapolated", file = noneSelfMHCbinder_proteomicFormat_resultOutfile)
					

					#V8:
					if not geneIntronInfo in geneName_MHCnonSelfMatchCount_dict:
						geneName_MHCnonSelfMatchCount_dict[geneIntronInfo] = 1
					else:
						geneName_MHCnonSelfMatchCount_dict[geneIntronInfo] += 1

					nonselfMHCmatch_intronGene_set.add(geneIntronInfo)

					###
					
					totalMatchCount += 1

					#Important new update - See if any of the patient have the HLA type matching that which can bind the matched 9mer.

					if patientID in patientSpecificHLAdict: #Blanks and mixed-diagnosis patients should not be present. 

						diagnosis_temp, HLA_specificToThisPatient_List = patientSpecificHLAdict[patientID]

						if diagnosis_temp != diagnosis: 
							raise ValueError("Diagnosis expected to be matching, did not match!", diagnosis_temp, diagnosis)
							sys.exit()

						if HLA_1field in HLA_specificToThisPatient_List: 	#If the HLA subtype, which binds to this 9mer, 
																			#is actually present in this patient,
																			#Then this is a really significant match!  

							#print("Significant match!", movingWindow, geneName, HLA_1field, diagnosis,patientID, f"| HLA statistics: {HLA_frequency} +- {HLA_std}")

							significantMatchCount += 1

							#Saving the significant result in another file:
							print(f"{movingWindow}\t{geneName}\t{HLA_1field}\t{diagnosis}\t{patientID}\t{bindingAffinity}\t{HLA_frequency}\t{HLA_std}\tIntrapolated\t{age}\t{brainRegion}\t{patient_sex}", file = outfile_significantResults)



	###################### Fragment is shorter than 9mer:  

	else: 
		#If the fragment is too short to compare to MHC 9mer, we want to find the expected 9mer version of the fragment, and based on this data, compare it all of the none-self MHC binding 9mers.
		
		entireGeneSequence = geneToSequenceDict[geneIntronInfo] #Look up the entire gene sequence, so the rest can be extrapolated. 

		seqCount = entireGeneSequence.count(sequence) #Amount of times this sequence (7- or 8mer) is present in its intron-retained gene variant. 

		if seqCount >= 2: #More than 1 repeat of the sequence is unlikely... but not impossible - deal with that part, if we actually encounter problems with it. 
			raise ValueError("Did not expect more than 1 match for a 7-8mer sequence in a single gene sequence", sequence)
			sys.exit()

		elif seqCount == 1:

			lenDif = 9 - len(sequence) #e.g. 9 - 7, would give a difference of 2. This means that we should take 2 extra nucleotides from the big gene sequence, which we can use for moving window. 

			indexValue = entireGeneSequence.index(sequence)

			fragmentSequence = sequence
			
			sequence = entireGeneSequence[indexValue-lenDif:indexValue+len(sequence)+lenDif] #If we have a fragment of lenght 7, we want to include 2 nucleotides up- and downstream, so we can create all potential moving windows. 
			#So if we e.g. get a match at index 18, and the fragment is 7 long, we need to take [16:27], to get the 2 flanking nucleotides as well. 



			#### The following part might appear like redundant/ repeated code - but the results produced by this, are more error-prone since the data is extrapolated! 

			for i in range(0,len(sequence)):

				movingWindow = sequence[i:i+9] #Since we are looking for 9mers! 


				if len(movingWindow) != 9: #The moving window should not have a length below 9. 
					break

				if not fragmentSequence in movingWindow:
					continue
					#Note: This is kinda lazy coding on my part. 
					#But if the sequence of 7- or 8-mer fragment is not 
					#present in the moving window, 
					#we need to skip to the next i-iteration.  
				else:
					if movingWindow in noneSelfMHCbinderDict:
						proteomic_Extrapolated_MHCbinderMatchList.append(movingWindow) #Save to this set, so we don't get repeated results. 

						#New, V3 - include HLA type and binding affinity:
						geneName, HLA, bindingAffinity = noneSelfMHCbinderDict[movingWindow] #Get the important information from the found MHC-binding non-self 9mer. 

						HLA_1field = re.sub(r":\d{2}$","",HLA) #Remove the last 2 digits and :. A*01:02 -> A*01

						#HLAfreqStatistic_Dict:
						#HLA -> Diagnosis -> (freq, std)

						if diagnosis == "N/A" or diagnosis == "unclear diagnosis": #If the diagnosis is unclear, we can't assign it a statistical value
							HLA_frequency, HLA_std = ("N/A", "N/A")

						else: #Else, we can make a look-up to get the statistical values. 
												
							try:
								HLA_frequency, HLA_std = HLAfreqStatistic_Dict[HLA_1field][diagnosis] #V9 - The look-up can now sometimes fuck up. E.g. in the case of 'B*58'. 
							except:
								HLA_frequency, HLA_std = ("N/A", "N/A")

						print(f"{line.rstrip()}\t{movingWindow}\t{geneName}\t{HLA_1field}\t{bindingAffinity}\t{diagnosis}\t{HLA_frequency}\t{HLA_std}\tExtrapolated", file = noneSelfMHCbinder_proteomicFormat_resultOutfile)
						

						#V8:
						if not geneIntronInfo in geneName_MHCnonSelfMatchCount_dict:
							geneName_MHCnonSelfMatchCount_dict[geneIntronInfo] = 1
						else:
							geneName_MHCnonSelfMatchCount_dict[geneIntronInfo] += 1

						nonselfMHCmatch_intronGene_set.add(geneIntronInfo)

						#####	

						totalMatchCount += 1

						#Important new update - See if any of the patient have the HLA type matching that which can bind the matched 9mer.

						if patientID in patientSpecificHLAdict: #Blanks and mixed-diagnosis patients should not be present. 

							diagnosis_temp, HLA_specificToThisPatient_List = patientSpecificHLAdict[patientID]

							if diagnosis_temp != diagnosis: 
								raise ValueError("Diagnosis expected to be matching, did not match!", diagnosis_temp, diagnosis)
								sys.exit()

							if HLA_1field in HLA_specificToThisPatient_List: 	#If the HLA subtype, which binds to this 9mer, 
																				#is actually present in this patient,
																				#Then this is a really significant match!  

								#print("(*)Significant match!", movingWindow, geneName, HLA_1field, diagnosis,patientID, f"| HLA statistics: {HLA_frequency} +- {HLA_std}")

								significantMatchCount += 1

								#Saving the significant result in another file:
								print(f"{movingWindow}\t{geneName}\t{HLA_1field}\t{diagnosis}\t{patientID}\t{bindingAffinity}\t{HLA_frequency}\t{HLA_std}\tExtrapolated\t{age}\t{brainRegion}\t{patient_sex}", file = outfile_significantResults)


						





#Total amount of results:

print("\n\tTotal amount of results (with repeats):", len(proteomic_Intrapolated_MHCbinderMatchList)+len(proteomic_Extrapolated_MHCbinderMatchList))
print("\tTotal amount of results (unique):", len(set(proteomic_Intrapolated_MHCbinderMatchList))+len(set(proteomic_Extrapolated_MHCbinderMatchList))) #Elements in a set can only be unique. 

print("\n\tIntrapolated results (with repeats):", len(proteomic_Intrapolated_MHCbinderMatchList))
print("\tIntrapolated results (unique):", len(set(proteomic_Intrapolated_MHCbinderMatchList)))

print("\n\tExtrapolated results (with repeats):", len(proteomic_Extrapolated_MHCbinderMatchList))
print("\tExtrapolated results (unique):", len(set(proteomic_Extrapolated_MHCbinderMatchList)))

###


total_MHCbinders_set = set(proteomic_Intrapolated_MHCbinderMatchList).union(set(proteomic_Extrapolated_MHCbinderMatchList))
print("\nBonus check:")
#print(len(set(set(proteomic_Intrapolated_MHCbinderMatchList) + set(proteomic_Extrapolated_MHCbinderMatchList))))
print(len(total_MHCbinders_set))



print("\n\nTotal matches:", totalMatchCount)
print("Significant matches:", significantMatchCount)

print(f"Significant matches, ratio (%): {round(significantMatchCount/totalMatchCount*100,2)}%")

###

print(f"\nMax length of peptide fragment: {peptide_maxLength}\nMin length of peptide fragment: {peptide_minLength}")

###

modificationCounts = [string.count("[") for string in modifiedPeptideDict[diagnosis] for diagnosis in modifiedPeptideDict.keys()]

uniqueSeqModificationCounts = [string.count("[") for string in set(modifiedPeptideDict[diagnosis]) for diagnosis in modifiedPeptideDict.keys()]

print(f"\nInfo - Modified residues in intron regions (might be faulty)\nDiagnosis: {modifiedPeptideDict.keys()}\nTotal amount of modifications: {sum(modificationCounts)}\nUnique results, total: {sum(uniqueSeqModificationCounts)}\n")

for diagnosis in modifiedPeptideDict:
	print(diagnosis,sum([seq.count("[") for seq in modifiedPeptideDict[diagnosis]]),sep="\t")


###

print("\nPatient group sizes:")
for diagnosis in diagnosis_size_dict:
	print(f"{diagnosis}:\t{len(diagnosis_size_dict[diagnosis])}")

###

print(f"\nTotal intron sequence length for each diagnosis group (pr. patient):")

for diagnosis in intronSequenceLength_diagnosisDict:
	print(f"\n--{diagnosis}--\nTotal:\t{intronSequenceLength_diagnosisDict[diagnosis]}")
	print(f"Pr. patient:\t{intronSequenceLength_diagnosisDict[diagnosis]/len(diagnosis_size_dict[diagnosis])}") #Getting the number pr. patient. 


print("\nBonus:", diagnosis_size_dict["PD"])

###

print(f"\nAmount of patients in proteomic data: {len(patientID_set)}\n")


###

#V8:

print(f"\nIntron fragment count for different diagnoses: \n{diagnosis_intronFragmentCount_dict}\n")


sortedGeneNames = sorted(geneName_MHCnonSelfMatchCount_dict.keys(), key = geneName_MHCnonSelfMatchCount_dict.get)

print([str(sortedGene)+":"+str(geneName_MHCnonSelfMatchCount_dict[sortedGene]) for sortedGene in sortedGeneNames])

print(", ".join([str(sortedGene) for sortedGene in sortedGeneNames]), end = "\n\n")

###

print("\nUnique intron variants observed in all of the patients:\n", ", ".join(allUniqueIntronGenes_set), sep="")
print("\nAmount of unique genes in the proteomic data: \n",len(allUniqueIntronGenes_set), sep="")

print("\nUnique intron variants observed in nonself MHC binding matches:\n", ", ".join(nonselfMHCmatch_intronGene_set), sep="")
print("\nAmount of unique genes associated with nonself MHC binding matches: \n",len(nonselfMHCmatch_intronGene_set), sep="")





#print(f"Another dict:\n{[f"{gene} : {geneName_MHCnonSelfMatchCount_dict[gene]}" for gene in sorted(geneName_MHCnonSelfMatchCount_dict.keys(), key = geneName_MHCnonSelfMatchCount_dict.get)]}")
print(f"\nAnother dict, length:\n{len(geneName_MHCnonSelfMatchCount_dict)}")

#print("\nSaved results - all of the important lines are retained, but irrelevant ones are removed - results in 'result_FilteredNoneSelfMHCbinder_presentInProteomicData_proteomicFormat_2018.tsv'")

print("\nSaved results - all of the important lines are retained, but irrelevant ones are removed - results in 'result_FilteredNoneSelfMHCbinder_presentInProteomicData_proteomicFormat_2022.tsv'")

print("Significant matches (where HLA of patient and HLA corresponding to predicted 9mer MHC-binder, are of the same type) is saved in 'result_HLAsignificantMatches.txt'")


#time.sleep(2)
#print("Unique results observed in proteomic data:")
#print(total_MHCbinders_set)



