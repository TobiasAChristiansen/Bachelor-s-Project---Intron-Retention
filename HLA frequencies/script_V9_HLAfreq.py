#!/usr/bin/env python3

import sys
import time
import math
import re


with open("data_header_info.txt","r",encoding = "ISO-8859-1") as infile:
	for line in infile:
		headerInfoList = line.strip().split(sep="\t")



################################### Shortening the HLA names

print("Before - Relevant HLA info:",headerInfoList) #[2:-14] is the actual important interval! 

for i, header in enumerate(headerInfoList):
	headerInfoList[i] = header[:-2] #This way we change 'HLA-DRB1_1' to 'HLA-DRB1'. This will be relevant when comparing with the reference frequencies, from another dataset. 

	if headerInfoList[i].startswith("HLA-"):
		headerInfoList[i] = headerInfoList[i][4:] #Removing the 'HLA-' part, from the name. 

print("\nAfter - Relevant HLA info:",headerInfoList)


#####################################


#field2dataFile = "DKMSLSL_X_DKBBHNF_RES_3090_39_20240308_customized_2_field.txt"

#V9:
superTypeHLAPatient_dataFileName = "HLA_A_and_B_alleles_translated.txt" #This files contains the HLA types converted to their supertypes, so we can compare them with the more of the predicted types of MHC binders. 
#Tobias made the file based on supertype HLA info. 



try:
	#openedDataFile = open(field2dataFile,"r")
	openedDataFile = open(superTypeHLAPatient_dataFileName,"r") #V9
except IOError as errText:
	print(errText)
	sys.exit()

#There should be 39 brain samples.
#The data is from patients with MSA, PSP and PD - this info can be viewed in the 'Barcoding order_DNA HLA-Genotyping_with diagnose' file. 

barcodeFile = "Barcoding order_DNA HLA-Genotyping_with diagnose.txt"

try: 
	openedBarcodeFile = open(barcodeFile,"r")
except IOError as errText:
	print(errText)
	sys.exit()

N = 0 #Amount of patients. 
barcodeToDiagnoseDict = dict() #Dict will be used to translate data barcode to patient diagnosis - very important feature. 
diagnoseList = [] #This list will contain all the diagnoses

diagnose_patientSizeDict = dict() #This dict will store the amount of patients with diagnosis X. 
#Information: Diagnosis > Count


barcodeToSampleID_dict = dict() #V7


for line in openedBarcodeFile:
	lineList = line.strip().split(sep="\t")

	if lineList[2].startswith("IML"): #We don't want to e.g. include the very first line - so we need a control check of the line's content. 
		
		N += 1 #Should be 39, at the end, since there are 39 patients in all. 
		barcode = lineList[2]
		diagnose = lineList[1]

		sampleID = lineList[0] #V7

		if not diagnose in diagnoseList: #If it is not already stored in the list, store it!
			diagnoseList.append(diagnose)
			diagnose_patientSizeDict[diagnose] = 1

		else: #Else, add to the couter associated with the diagnosis type. 
			diagnose_patientSizeDict[diagnose] = diagnose_patientSizeDict[diagnose] + 1

		if not barcode in barcodeToDiagnoseDict.keys():
			barcodeToDiagnoseDict[barcode] = diagnose #Add diagnosis to look-up dict. 

		else:
			raise ValueError("Did not expect to encounter a barcode which has already been assigned a diagnose!", line)



		#V7:

		if barcode in barcodeToSampleID_dict:
			raise ValueError("Barcode already in the dict!")
			sys.exit()

		else:
			barcodeToSampleID_dict[barcode] = sampleID #This way we can convert barcode to brain sample ID. 





#Print results: 
print("\nTotal amount of patients:", N)
print("All diagnoses:", diagnoseList)
print("Distrubtion of diagnoses:", diagnose_patientSizeDict,"\n")



############################################### Important change - double the size of N. Since we have both HLA-A_1 and HLA-A_2 (the HLA-A inherited from both parents)

N = N*2 #If e.g. a 100% of patients have an allele, there would be 39*2 = 78 of this allele (a HLA from each parent, 2 HLAs). Therefore, we would also need to divide by 78, to get a frequency of 1.0. 

###############################################

#Now that we have established the dict, we can easily make a simple function to translate barcode to diagnose. 
def barcodeToDiagnose(barcode):
	
	if barcode in barcodeToDiagnoseDict.keys():
		return barcodeToDiagnoseDict[barcode] #Returns the disease. 

	else:
		raise ValueError("This barcode does not have an assigned disease!", barcode) 

### Reference frequences: 

#Before we get started, we also need to be able to compare our frequencies with some European average value. The data is from 'AlleleFrequencies.net'.

averageFrequencyFile = "result_averageHLAfreq_EU_field1.txt" #V8 - new data file. 

#These results are produced by the Python script named: 'script_HLA_webdataFreq'.

try:
	openedAverageFrequencyFile = open(averageFrequencyFile,"r")
except IOError as errText:
	print(errText)
	sys.exit()





HLAaverageFrequencyDict = dict() 
		#This dict will be used to store the average frequency and its respective standard deviation, 
		#for the reference EU data. 
	
	#Information: A*01 > Reference frequency and standard deviation

	#Note: Now the data is field 1! 

for line in openedAverageFrequencyFile:

	lineList = line.strip().split(sep="\t")

	if len(lineList) == 3: #Column 1 contains HLA subtype names (2 fields) and column 2 contains the respective frequency average. Column 3 contains the amount of studies which included data of the frequency, for this HLA.
		HLAname = lineList[0]
		frequencyAverage = float(lineList[1])
		standardDeviation = float(lineList[2])

		if not HLAname in HLAaverageFrequencyDict.keys():
			HLAaverageFrequencyDict[HLAname] = (frequencyAverage, standardDeviation)

		else:
			raise ValueError(f"The same HLA subtype name was found twice in the {averageFrequencyFile} file! {HLAname}, {HLAaverageFrequencyDict[HLAname]}")

########################################################################## Question: if an HLA allele is not present in a given dataset, can we then assume a frequency of 0? This would significantly decrease our average frequencies.


#V7 - Saving information about HLA subtypes specific to a patient in a new file

outfile_HLAoverview = open("result_HLA_patient_overview.txt","w")



##########################################################################



###
#A line in the file could e.g. look like this:
#IML007,ID22218783,02:01/02:01L/02:01Q/02:09,NNNN,NNNN,,,,,,,,,,,3.50.0,23.08.23,DNA,DKMS-LSL
#NNNN means that no HLA exist for this HLA type. This is since some of the HLAs are mutual exclusive. 
#',,,,,,,,,,,3.50.0,23.08.23,DNA,DKMS-LSL' refers to 'HLA-E_1	HLA-E_2	ABO	RHD	CCR5_1	CCR5_2	KIR	MICA	MICB	CMV	HLA_Version	lab_date	material	lab_name'. 
#This last part will be ignored, since it doesn't appear to have a relevance to the study at hand. 

######################## Important part: ########################






#HLAcountDict = dict() #Information is stored in the dict like this... 
#Information: HLA name (A*01:02) > Diagnosis (MSA) > Count

HLA_1fieldDict = dict() 
#Information: HLA 1-field name (A*01) > Diagnosis (MSA) > Count

N_sizeDict = dict() #This dict will be used to tell how much HLA data is stored for e.g. HLA-A_1 and HLA-A_2. 

#Information: HLA gene (B) > Amount of HLA data for the 2x HLA-B columns. 


#V5 change - New dict. For calculating frequency differently: 

HLA_1fieldDict_singleCount = dict() 
#Information: HLA 1-field name (A*01) > Diagnosis (MSA) > Count (only +1 or +2 for each patient), 


	#and id for the last list, where it was observed. 
	
		#This dict will only count 1-2 times for each patient! 
		#This has been added as of version 5, 
		#to ensure that the final frequency is calculated properly. 




lineCount = 0



for line in openedDataFile:

	print(line)
	time.sleep(0.1)

	lineCount += 1

	lineList = line.strip().split(sep=",")

	if len(lineList) == 38: #Note, this is the length that the format should have. 

		if lineList[0].startswith("IML"): #If the line starts with a barcode
			barcode = lineList[0]
			diagnose = barcodeToDiagnose(barcode)

			brainSampleID = barcodeToSampleID_dict[barcode] #V7


		if lineList[1].startswith("ID"): #If the subsequenct column contains an ID. 
			dataID = lineList[1]
		
						#By looking in the Excel, you get a more clear picture. 
						#[barcode,id,HLA-A_1,HLA-A_2,HLA-B_1,HLA-B_2,HLA-C_1,HLA-C_2,HLA-DRB1_1,HLA-DRB1_2,DRB3_1,DRB3_2,DRB4_1,DRB4_2,DRB5_1,DRB5_2,DQA1_1, etc.]
						#A quick look at the data and format would help understand what I do here. 
		

		#A new line means a new patient, with different data:

		relevantHLAdataList = lineList[2:-14] 				#Contains all the HLA specific to the HLA A, B, etc.
		relevantHeaderInfoList = headerInfoList[2:-14] 		#A, A, B, B, etc.





		patientSpecificHLASet = set() #V7 - Set containing all unique HLA subtypes for this patient. 


		for i in range(0,len(relevantHLAdataList)): 	#Iterate over all subtypes i.e. HLA-A_1, HLA-A_2, etc. 

			HLAgeneName = relevantHeaderInfoList[i] 	#e.g., A, A, B, B, etc.

			subTypeList_i = relevantHLAdataList[i].split(sep="/") 	

						#A patient might have, e.g. 'HLA-A_2' with 40 different subtypes. 
						#These subtypes are contained in this subtype_list.  

			
			#dataID = id(subTypeList_i) #Using dataID gave me an error, since it only generated 5 unique IDs.... :// 


			#Adding to N_sizeDict:

			if not HLAgeneName in N_sizeDict.keys():
				N_sizeDict[HLAgeneName] = len(subTypeList_i)
			
			else:
				N_sizeDict[HLAgeneName] = N_sizeDict[HLAgeneName] + len(subTypeList_i)



			for HLAsubtype in subTypeList_i:

				if not HLAsubtype: #Ignore empty elements
					continue
				

				############################################################################ Important: If we uncomment this, our summed frequencies won't give 1. Hence the check we do when printing the most frequent 1-field HLA types! 
				#if HLAsubtype == "NNNN":
				#	continue #We don't want the 'NNNN' as part of our final results. 


				#fullName_HLAsubtype = HLAgeneName + '*' + HLAsubtype #e.g.: A + * + 01:02

				HLAsubtype_number = HLAsubtype.split(sep=':')[0] #V9

				print("Check:",HLAsubtype_number)
				if not HLAsubtype_number[0].isalpha():
					field1Only_HLAsubtype = HLAgeneName + '*' + HLAsubtype_number #e.g. A + * + 01
				else:
					field1Only_HLAsubtype = HLAsubtype_number #The HLA is already of the correct format. 



				#Note: As of V8, this should be the same. 


				##### V7 - Make a set containing all the unique HLA subtypes (field 1), specific to this patient: 

				patientSpecificHLASet.add(field1Only_HLAsubtype)




				#######
				#Same thing, but now with only 1 field!  

				if not field1Only_HLAsubtype in HLA_1fieldDict.keys(): #If this is the first time encountering this HLA subtype - create diagnose and count dict. 
					HLA_1fieldDict[field1Only_HLAsubtype] = {diagnose : 1} #Create dict entry, if it is not already there. 


				elif field1Only_HLAsubtype in HLA_1fieldDict.keys(): #If this is not the first time encountering this HLA subtype:

					try: 
						if HLA_1fieldDict[field1Only_HLAsubtype][diagnose]: #Check if this diagnose has already been associated with this HLA subtype...

							HLA_1fieldDict[field1Only_HLAsubtype][diagnose] = HLA_1fieldDict[field1Only_HLAsubtype][diagnose] + 1 #... if so, update the count.


					except: #If the HLA subtype hasn't appeared with this diagnose before - make a new count entry in the dict, for this diagnose.  

						HLA_1fieldDict[field1Only_HLAsubtype][diagnose] = 1

						#Note, we don't make a new dict entry using '{diagnose : 1}'. 



				######

				#V5	-	New population number (N) for new frequency calculations:


				if not field1Only_HLAsubtype in HLA_1fieldDict_singleCount.keys(): #If not e.g. 'A*01' is found.
					
					#Information: HLA 1-field name (A*01) > Diagnosis (MSA) > Count (only +1 or +2 for each patient), and id for the last list, where it was observed. 


					HLA_1fieldDict_singleCount[field1Only_HLAsubtype] = {diagnose : (1, lineCount, i)} #Tuple with count and ID specific to the current data list. 


				elif field1Only_HLAsubtype in HLA_1fieldDict_singleCount.keys():

					try: 
						if HLA_1fieldDict_singleCount[field1Only_HLAsubtype][diagnose]: #Check if this diagnosis has already been associated with this HLA subtype...

							#If this HLA type and diagnosis has been encountered before
							
							countValue, prev_lineCount, prev_index = HLA_1fieldDict_singleCount[field1Only_HLAsubtype][diagnose] 

										#Extract count value and previous ID (so we can compare it with the current one 
										#- if they are the same, we are looking at the same data list. 
										#We don't want this, since we only want to count the HLA occurence for each patient 
										#in each diagnosis group (although there can max be 2 occurences pr. patient))


							#If the data is not from the same list, it will comply with these criteria: 
							if prev_lineCount == lineCount and prev_index != i: #This is a check to ensure that the same HLA-A*01 is not counted multiple times in the same patient. 
								HLA_1fieldDict_singleCount[field1Only_HLAsubtype][diagnose] = (countValue+1, lineCount, i) #Should max have a count value around 20, for each diagnosis (patient group size multiplied by 2). 
							elif prev_lineCount != lineCount and prev_index == i:
								HLA_1fieldDict_singleCount[field1Only_HLAsubtype][diagnose] = (countValue+1, lineCount, i) #Should max have a count value around 20, for each diagnosis (patient group size multiplied by 2). 
							elif prev_lineCount != lineCount and prev_index != i:
								HLA_1fieldDict_singleCount[field1Only_HLAsubtype][diagnose] = (countValue+1, lineCount, i) #Should max have a count value around 20, for each diagnosis (patient group size multiplied by 2). 



					except: #If the HLA subtype hasn't appeared with this diagnose before - make a new count entry in the dict, for this diagnose.  

						HLA_1fieldDict_singleCount[field1Only_HLAsubtype][diagnose] = (1, lineCount, i) #Tuple with count and index (last index where match was encountered).
						#Making it for this diagnosis, for the first time. 



		#V7 - Unique HLA subtypes for each patient:
		
		if patientSpecificHLASet: #V9
		#If it is not empty. 
		#Idk, it is a lazy fix, since the last two lines would otherwise contain 1495 with empty HLA list. 

			sorted_patientSpecificHLA_List = sorted(list(patientSpecificHLASet))
			
			print(brainSampleID,diagnose,";".join(sorted_patientSpecificHLA_List),sep="\t",file = outfile_HLAoverview)



	else:
		raise ValueError("Wrong type of file/ line format - should have 38 commas:", line)
		sys.exit()



#Printing results: 
#print("\nAmount of different subtypes of HLA in data:",len(HLAcountDict)) #2440

print("\nAmount of different field 1 HLA types in data:",len(HLA_1fieldDict)) #153



continueInput = input("\n\tDo you want to see frequencies of HLA data from all patients? (Press 'enter'): ")



###


#Add an entry to the dict, which contains the total amount of counts/ HLA subtypes, i.e. the sum from all diagnoses:


###
#Now for only 1-field info:


total_1field_HLAcountDict = dict()
#Information: HLA-A*01 > Count

for HLA, diagnoseCountDict in HLA_1fieldDict.items():

	totalCount = sum(diagnoseCountDict.values())


	#regex_HLAtype = re.search(r"(\w+)\*",HLA) #Getting e.g. 'DPB1' from 'DPB1*02'
	#if regex_HLAtype:
		#HLAgene = regex_HLAtype.group(1) #This will e.g. be 'DPB1'. 
		#N_HLA = N_sizeDict[HLAgene]


	if not HLA in total_1field_HLAcountDict.keys():
		total_1field_HLAcountDict[HLA] = totalCount #Storing the total count in a new dict, for only 1-field HLA, such as amounts of HLA-A*01 data. 
	
	else:
		raise ValueError("2 of identical HLA in the dict:", HLA, total_1field_HLAcountDict[HLA])
		sys.exit()



#Making a sorted list for the results. 

#sortedHLAcountList = sorted(totalHLAcountDict,key=totalHLAcountDict.get) #Not relevant as of V8. 

###
#1 field sorting a list:

sorted_1field_HLAcountList = sorted(total_1field_HLAcountDict,key=total_1field_HLAcountDict.get)




#Printing sorted data, for the total frequencies:


#############

#Now doing the same thing! But only for 1 FIELD of the HLA:


#Printing sorted data, for the total frequencies:

HLA_1field_specificFrequencyDict = dict() # A dict which will be used to store HLA information in this order: DPB1 > DPB1*02 > Frequency. 



for HLA in sorted_1field_HLAcountList: #Looping through the subtype HLA names that have been sorted based on their count value. 


	regex_HLAtype = re.search(r"(\w+)\*", HLA) #Getting e.g. 'DPB1' from 'DPB1*02'

	if regex_HLAtype:

		HLAgene = regex_HLAtype.group(1) #This will e.g. be 'DPB1'. 

		N_HLA = N_sizeDict[HLAgene] #Would e.g. get the total count of 'DPB1'. 


	else:
		raise ValueError("RE did not succed on HLA name: ", HLA)
		sys.exit()



	frequency = float(total_1field_HLAcountDict[HLA]/N_HLA) #The probability of getting having an allele, aka. the frequency. 

	if frequency > 1:
		print("Too high frequency:", frequency)
		raise ValueError(f"Frequency should not be higher than 1! Count of {HLA} is {total_1field_HLAcountDict[HLA]} out of {HLAgene} in total {N_HLA}")
		sys.exit()

	#! Note the different use of 'N_HLA' instead of 'N' (39*2)

	standardDeviation = math.sqrt(N_HLA*frequency*(1-frequency)) #Binomial std formula. 

	if HLA in HLAaverageFrequencyDict.keys(): #We only print HLAs with a correspondning reference. 
		
		reference_frequencyAverage = HLAaverageFrequencyDict[HLA][0]
		reference_standardDeviation = HLAaverageFrequencyDict[HLA][1]

		print(f"\nHLA subtype: {HLA}\nOccurrence: {total_1field_HLAcountDict[HLA]} +/- {round(standardDeviation,3)} \nFrequency: ({total_1field_HLAcountDict[HLA]}/{N_HLA}) = {round(frequency,3)} +/- {round(standardDeviation/N_HLA,3)}")
		print(f"(!)\tEU average HLA frequency: {round(reference_frequencyAverage,4)} +/- {round(reference_standardDeviation,4)}\n")

	else:
		
		
		#try:
		#	standardDeviation = math.sqrt(N_HLA*frequency*(1-frequency)) #Binomial std formula. 
		#except:
		#	raise ValueError("Could not calculate the standard deviation:", (N_HLA*frequency*(1-frequency)), frequency, N_HLA)

		print(f"\nHLA subtype: {HLA}\nOccurrence: {total_1field_HLAcountDict[HLA]} +/- {round(standardDeviation,3)} \nFrequency: ({total_1field_HLAcountDict[HLA]}/{N_HLA}) = {round(frequency,3)} +/- {round(standardDeviation/N_HLA,3)}")
		print(f"\tEU average HLA frequency: N/A\n")





	#Now making a dict which will later be used to only get the X amount of most frequent HLA subtype of each HLA type. 

	#Dict should work like this: DPB1 > DPB1*02 > Frequency. 

	regex_HLAtype = re.search(r"(\w+)\*",HLA) #Getting e.g. 'DPB1' from 'DPB1*02'


	if regex_HLAtype:

		if HLAgene not in HLA_1field_specificFrequencyDict.keys(): #If e.g. 'DPB1' not in dict containing 'DPB1*01', 'DPB1*02', etc.

			HLA_1field_specificFrequencyDict[HLAgene] = {HLA : frequency}

		else:
			HLA_1field_specificFrequencyDict[HLAgene][HLA] = frequency




continueInput = input("\tDo you want to continue to the most common field-1 HLA gene types? [Press 'enter']: ")




print("\n#################################### Most frequent HLA subtypes ####################################\n")

#Getting the most frequent HLA subtypes: 

numberOfMaxToDisplay = 10 #We want to display the 10 highest frequencies for e.g. each of the HLA type 'DPB1' allele subtypes. 


#Not relevant as of V8: 
"""

# For 2-fields: 
for HLAtype, HLAfrequencyDict in HLAspecificFrequencyDict.items(): #For e.g. 'DPB1' and dict containing all 'DPB1 subtypes' and their associated 'frequency'.
	#print(f"\nMost frequent HLA of class: {HLAtype}") #E.g. 'DPB1'

	sortedFrequencyList = sorted(HLAfrequencyDict,key=HLAfrequencyDict.get, reverse = True)

	frequencyList = [HLAspecificFrequencyDict[HLAtype][subtype] for subtype in sortedFrequencyList[0:numberOfMaxToDisplay]] #Making a list with the 10 first elements in the list.

	#print(f"\tTop {numberOfMaxToDisplay} most frequent subtypes...\nSubclass: \n\t{[subtype for subtype in sortedFrequencyList[0:numberOfMaxToDisplay]]}\nFrequency: \n\t{[round(frequency,4) for frequency in frequencyList]}\n")

"""


###

# For 1-field HLA data:

for HLAtype, HLA_1field_frequencyDict in HLA_1field_specificFrequencyDict.items(): #For e.g. 'DPB1' and dict containing all 'DPB1 subtypes' (DPB1*01, DPB1*02, etc.) and their associated 'frequencies'.
	print(f"\nMost frequent HLA of class: {HLAtype}") #E.g. 'DPB1'

	sorted_1field_FrequencyList = sorted(HLA_1field_frequencyDict,key=HLA_1field_frequencyDict.get, reverse = True)


	#Making a list with the 10 first elements in the list/ aka. the top 10 highest frequencies.
	frequency_1field_List = [HLA_1field_specificFrequencyDict[HLAtype][subtype] for subtype in sorted_1field_FrequencyList[0:numberOfMaxToDisplay]] 

	print(f"\tTop {numberOfMaxToDisplay} most frequent subtypes...\nSubclass: \n\t{[subtype for subtype in sorted_1field_FrequencyList[0:numberOfMaxToDisplay]]}\nFrequency: \n\t{[round(frequency,5) for frequency in frequency_1field_List]}")



#### Bonus check - The sum of the frequencies, should give 1. 

	print(f"The sum of all the {HLAtype}-frequencies should be 1: {sum(HLA_1field_specificFrequencyDict[HLAtype].values())}\n")





continueInput = input("\tDo you also wanna see the occurrence specific to the different diagnoses? (press 'enter'): ")


#Binomial standard deviation: 

#Calculating the average for each HLA out of total amount of patients (N):
#This is a binomial distribution - since a HLA allele is either present or not. 
#The average aka. frequency aka. probability (p) is calculated using:
#p = v/N
#Where N is amount of patients (39 patients). v is the probability of success in N amount of trials (amount of people with a given allele). 
#p could e.g. be 17 out of 39 patients, resulting in a probability of = 0.44
#The standard deviation in binomial distribution is calculated this way:
#std = sqrt(N*p(1-p))
#Following previous example: std = sqrt(39*0.44*(1-0.44)) = 3.10


### 2-field information: 


#Not relevant as of V8: 
"""

#outfile = open("result_HLAfrequencies.txt","w")


for HLA in sortedHLAcountList: #Does not guarantee that the sorting from last time is accurate. But it is a place to start. 

	temp_outputList = [[],[],[],[]] #For each new HLA subtype, clear this list.


	for diagnose in diagnoseList: 

		if diagnose in HLAcountDict[HLA].keys(): #Check whether information (count, aka. occurrence) regarding diagnose X is has been observed in any of the HLA subtypes

			HLAoccurrence_thisDiagnoseOnly = HLAcountDict[HLA][diagnose] #Aka. 'v'
			
			N_diagnosePatientSize = diagnose_patientSizeDict[diagnose]*2 #It is wrong to use N. We need to calculate using a N-value specific to the diagnose population size. 
			
			#Note, again we multiply N with 2, since there are 2 alleles - one from both parents.

			frequency = round(HLAoccurrence_thisDiagnoseOnly/N_diagnosePatientSize,5) #Aka. 'p'

			if frequency > 1:
				raise ValueError("A frequency can not be higher than 1!", HLAoccurrence_thisDiagnoseOnly, N_diagnosePatientSize, frequency)
				sys.exit()


			try:
				standardDeviation = round(math.sqrt(N_diagnosePatientSize*frequency*(1-frequency)),3) #The standard deviation for HLAoccurrence_thisDiagnoseOnly
			except:
				raise ValueError("Could not take the sqrt. Relevant numbers:", N_diagnosePatientSize, frequency)

			temp_outputList[0].append(diagnose)
			temp_outputList[1].append(HLAoccurrence_thisDiagnoseOnly)
			temp_outputList[2].append(frequency)
			temp_outputList[3].append(standardDeviation)

	#After we have iterated over all the diagnoses for this HLA subtype, we just want to print all the info for this HLA subtype:

	temp_diagnoseList = temp_outputList[0]
	temp_occurrenceList = temp_outputList[1]
	temp_frequencyList = temp_outputList[2]
	temp_standardDeviationList = temp_outputList[3] #This standard deviation is for the occurrences

	frequencyStandardDeviationList = [std/N_diagnosePatientSize for std in temp_standardDeviationList]


	#print(f"HLA: {HLA}\nDiagnosis:\t{temp_diagnoseList}\nOccurrences: \t{temp_occurrenceList} +/- {temp_standardDeviationList}\nFrequencies: \t{temp_frequencyList} +/- {frequencyStandardDeviationList}")

	#print(f"HLA: {HLA}\nDiagnosis:\t{temp_diagnoseList}\nOccurrences: \t{temp_occurrenceList} +/- {temp_standardDeviationList}\nFrequencies: \t{temp_frequencyList} +/- {frequencyStandardDeviationList}", file=outfile)

	if HLA in HLAaverageFrequencyDict.keys(): #If reference data actually exists for this subtype, we are interested in printing it. 
		reference_frequencyAverage = HLAaverageFrequencyDict[HLA][0]
		reference_standardDeviation = HLAaverageFrequencyDict[HLA][1]
		
		#print(f"(!)\tReference EU frequency: {round(reference_frequencyAverage,5)} +/- {round(reference_standardDeviation,5)}\n")
		#print(f"(!)\tReference EU frequency: {round(reference_frequencyAverage,5)} +/- {round(reference_standardDeviation,5)}\n", file=outfile)

	else: 
		#print("\tReference EU frequency: N/A\n", file=outfile)
		pass

#print("Done writing the results into the file name, 'result_HLAfrequencies.txt'")

#outfile.close()

"""


##########################

#We need to get all count of e.g. gene HLA-A counts specific to a certain diagnosis group.

#Remember: HLA_1fieldDict 
#Information: HLA 1-field name (A*01) > Diagnosis (MSA) > Count

HLAgeneDiagnosisSpecificCountDict = dict()
#Information: HLA gene (A) > Diagnosis (MSA) > count


for HLA_field1, diagnosis_countDict in HLA_1fieldDict.items():

	regex_HLAtype = re.search(r"(\w+)\*",HLA_field1) #Getting e.g. 'DPB1' from 'DPB1*02'

	if regex_HLAtype:

		HLAgene = regex_HLAtype.group(1) #This will e.g. be 'DPB1'. 

	for diagnosis, count in diagnosis_countDict.items(): 

		try:
			if HLAgeneDiagnosisSpecificCountDict[HLAgene][diagnosis]:
				HLAgeneDiagnosisSpecificCountDict[HLAgene][diagnosis] = HLAgeneDiagnosisSpecificCountDict[HLAgene][diagnosis] + count #If this diagnosis already has an entry for this HLA gene (A) and diagnosis (MSA), then add to the count. 
		
		except:
			
			try:
				if HLAgeneDiagnosisSpecificCountDict[HLAgene]: #If entries have previously been made for this HLA gene, but not for this specific diagnosis. 
					
					HLAgeneDiagnosisSpecificCountDict[HLAgene][diagnosis] = count

			except:
				HLAgeneDiagnosisSpecificCountDict[HLAgene] = {diagnosis : count} #If no entries have been made for this HLA gene previously, make one for this specific diagnosis, now. 




#print(HLAgeneDiagnosisSpecificCountDict["A"])
#print("Check:",sum(HLAgeneDiagnosisSpecificCountDict["A"].values()),"=", N_sizeDict["A"]) #Should give the same number. 


##########################



### 1-field information, diagnosis specific: 

outfile_1field = open("result_HLAfrequencies_1field.txt","w")

print('#HLA subtype','Diagnosis','Frequency','Standard Deviation',"EU Freq", "EU STD",sep="\t",file = outfile_1field)


HLA_1field_statisticsDiagnosisDict = dict()
#Information: HLA (A*01) > Diagnosis (MSA) > (Frequency, Standard deviation)

#print(HLA_1fieldDict_singleCount["DRB3*NNNN"])


#for HLA in sorted_1field_HLAcountList: #Does not guarantee that the sorting from last time is accurate. But it is a place to start. 


for HLA in HLA_1fieldDict_singleCount:
#Information: HLA 1-field name (A*01) > Diagnosis (MSA) > Count (only +1 or +2 for each patient), and id for the last list, where it was observed. 


	temp_outputList = [[],[],[],[]] #For each new HLA subtype, clear this list.

	for diagnose in diagnoseList: 

			if diagnose in HLA_1fieldDict_singleCount[HLA].keys(): #Note that the diagnosis has an entry. 


				regex_HLAtype = re.search(r"(\w+)\*",HLA) #Getting e.g. 'DPB1' from 'DPB1*02'

				if regex_HLAtype:

					HLAgene = regex_HLAtype.group(1) #This will e.g. be 'DPB1'. 



					##### V5 - New HLAoccurrence_thisDiagnoseOnly and N_HLA value!! 

					#Remember: HLA_1fieldDict_singleCount
					#Information: HLA 1-field name (A*01) > Diagnosis (MSA) > Count (only +1 or +2 for each patient), and iteration number, where last iteration was found. 

					HLAoccurrence_thisDiagnoseOnly = HLA_1fieldDict_singleCount[HLA][diagnose][0]


					N_HLA = diagnose_patientSizeDict[diagnose]*2 #N*2, since there can be a HLA from both parents. 
					

					#Remember: diagnose_patientSizeDict
					#Information: Diagnosis > Count


				else:
					raise ValueError("RE did not succed on HLA name: ", HLA)
					sys.exit()

				

				frequency = round(float(HLAoccurrence_thisDiagnoseOnly/N_HLA),5) #The probability of getting having an allele, aka. the frequency. 


				if frequency > 1:
					print("Too high frequency:", frequency)
					raise ValueError(f"Frequency should not be higher than 1! Count of {HLA} is {total_1field_HLAcountDict[HLA]} out of {HLAgene} in total {N_HLA}")
					sys.exit()


				try:
					standardDeviation = round(math.sqrt(N_HLA*frequency*(1-frequency)),5) #The standard deviation for HLAoccurrence_thisDiagnoseOnly
				except:
					raise ValueError("Could not take the sqrt. Relevant numbers:", N_HLA, frequency)

				temp_outputList[0].append(diagnose)
				temp_outputList[1].append(HLAoccurrence_thisDiagnoseOnly)
				temp_outputList[2].append(frequency)
				temp_outputList[3].append(standardDeviation)


			else:
				temp_outputList[0].append(diagnose)
				temp_outputList[1].append(0) #New as of V6: Just fill in blank output.
				temp_outputList[2].append(0)
				temp_outputList[3].append(0)



	#After we have iterated over all the diagnoses for this HLA subtype, we just want to print all the info for this HLA subtype:

	temp_diagnoseList = temp_outputList[0]
	temp_occurrenceList = temp_outputList[1]
	temp_frequencyList = temp_outputList[2]
	temp_standardDeviationList = temp_outputList[3] #This standard deviation is for the occurrences

	frequencyStandardDeviationList = [round(std/N_HLA,5) for std in temp_standardDeviationList]


	#V8 - I added back the part where you print EU reference data:

	if HLA in HLAaverageFrequencyDict:
		EU_freq, EU_std = HLAaverageFrequencyDict[HLA]

	else:
		EU_freq, EU_std = ("N/A", "N/A")


	print(f"HLA: {HLA}\nDiagnosis:\t{temp_diagnoseList}\nOccurrences: \t{temp_occurrenceList} +/- {temp_standardDeviationList}\nFrequencies: \t{temp_frequencyList} +/- {frequencyStandardDeviationList}\nEU freq: \t{EU_freq}\nEU std: \t{EU_std}")



	temp_frequencyList = [str(num) for num in temp_frequencyList] #Convert to string, so we can use join method on the numbers. 

	frequencyStandardDeviationList = [str(num) for num in frequencyStandardDeviationList]

	print(f"{HLA}\t{';'.join(temp_diagnoseList)}\t{';'.join(temp_frequencyList)}\t{';'.join(frequencyStandardDeviationList)}\t{EU_freq}\t{EU_std}", file = outfile_1field)




	###################################### 


	#Write the results into the dict: 
	if not HLA in HLA_1field_statisticsDiagnosisDict.keys():

		for i in range(len(temp_diagnoseList)): #Iterate every diagnosis for this HLA, and write its info into the dict. 

			try:
				HLA_1field_statisticsDiagnosisDict[HLA][temp_diagnoseList[i]] = (temp_frequencyList[i],frequencyStandardDeviationList[i]) #Saving the frequency and standard deviation under: A*01 > MSA > (frequency, standard deviation)
			except:
				HLA_1field_statisticsDiagnosisDict[HLA] = {temp_diagnoseList[i] : (temp_frequencyList[i],frequencyStandardDeviationList[i])} #The first time, we need to make an entry for the HLA. 

	else:
		raise ValueError("Did not expect to encounter the same HLA 2 times!", HLA, HLA_1field_statisticsDiagnosisDict[HLA])



#Sort the data into most frequent e.g. HLA type A, for each diagnosis. 

outfile_1field.close()

print("\nDone writing the results into the file name, 'result_HLAfrequencies_1field.txt'\n")





################ Make a script which picks the highest scoring HLA types for each diagnosis and each HLA gene type (A, B, etc.)

topHLAgeneDiagnosisDict = dict() 
#Information: HLA gene ('A') > diagnosis (MSA) > (Top scoring HLA 1-field subtype (A*01), frequency, standard deviation)


for HLA, infoDict in HLA_1field_statisticsDiagnosisDict.items():

	#Reminder: HLA (A*01) > Diagnosis (MSA) > (Frequency, Standard deviation)


	#Extracting the HLA gene name ('A'): 
	regex_HLAtype = re.search(r"(\w+)\*",HLA) #Getting e.g. 'DPB1' from 'DPB1*02'

	if regex_HLAtype:
		HLAgene = regex_HLAtype.group(1) #This will e.g. be 'DPB1'. 
		

	for diagnosis, frequency_standardDeviationTuple in infoDict.items():

		frequency, standardDeviation = frequency_standardDeviationTuple[0], frequency_standardDeviationTuple[1]

		try:
			EU_freq, EU_std = HLAaverageFrequencyDict[HLA] #V8
		except:
			EU_freq, EU_std = ("N/A","N/A") #V8


		try: 
			currentHighestFrequency = topHLAgeneDiagnosisDict[HLAgene][diagnosis][1] #Calling this dict entry - will raise exception, if not there. 
			
			if currentHighestFrequency < frequency: #If the frequency of the current HLA subtype (e.g. A*02 instead of A*01, which we are comparing to), we need to replace the highest scoring entry. 

				#topHLAgeneDiagnosisDict[HLAgene][diagnosis] = (HLA, frequency, standardDeviation)
				topHLAgeneDiagnosisDict[HLAgene][diagnosis] = (HLA, frequency, standardDeviation, EU_freq, EU_std) #V8

		except: #If this is the first time we try looking up e.g.: 'A' > 'MSA'. Then establish a dict entry (for now, this will be the highest scoring entry). 

			try:
				if len(topHLAgeneDiagnosisDict[HLAgene]) >= 1: #If we already have an entry for another diagnosis, just insert an entry specific to this current diagnose.
					#topHLAgeneDiagnosisDict[HLAgene][diagnosis] = (HLA, frequency, standardDeviation)

					topHLAgeneDiagnosisDict[HLAgene][diagnosis] = (HLA, frequency, standardDeviation, EU_freq, EU_std) #V8.

			except: #Else, create an entry - this is the first time we encounter this diagnosis (MSA) for this type of HLA gene (A). 
				topHLAgeneDiagnosisDict[HLAgene] = {diagnosis : (HLA, frequency, standardDeviation, EU_freq, EU_std)} #V8.



outfile_plotData = open("result_HLAFrequencyPlotData.txt","w")




print("Getting the most frequent HLA gene type for each diagnosis:")

for j in range(0,len(relevantHeaderInfoList),2):

	HLAgene = relevantHeaderInfoList[j]

	mostFrequentHLAgene_thisDiagnoseOnly_dict = topHLAgeneDiagnosisDict[relevantHeaderInfoList[j]]

	#Information: HLA gene ('A') > diagnosis (MSA) > (Top scoring HLA 1-field subtype (A*01), frequency, standard deviation, EU freq, EU std)



	print(f"\nTop score - {HLAgene} - {mostFrequentHLAgene_thisDiagnoseOnly_dict}")


	for diagnosis_selected, statisticInfoTuple in mostFrequentHLAgene_thisDiagnoseOnly_dict.items():

		HLA = statisticInfoTuple[0] #Getting the most frequent HLA subtype (A*01) for each diagnose (MSA) of this HLA gene type (A, B, C, etc.)

		EU_freq, EU_std = statisticInfoTuple[-2],statisticInfoTuple[-1]

		print(f"Associated values for other diagnoses:")
		print(f"\t*{diagnosis_selected}\t{HLA}")

		print(f"\t{diagnosis_selected}\t{HLA}",file=outfile_plotData)

		for diagnosis in diagnoseList:
			
			try:
				frequency = HLA_1field_statisticsDiagnosisDict[HLA][diagnosis][0]
			except:
				frequency = 0

			try:
				standardDeviation = HLA_1field_statisticsDiagnosisDict[HLA][diagnosis][1]
			except:
				standardDeviation = 0

			print(HLA,diagnosis,frequency,standardDeviation,EU_freq, EU_std,sep="\t")

			print(HLA,diagnosis,frequency,standardDeviation,sep="\t",file=outfile_plotData)
			#print(HLA,diagnosis,frequency,standardDeviation,EU_freq, EU_std,sep="\t",file=outfile_plotData)
		

		print(HLA,"EU ref.",EU_freq,EU_std,sep="\t",file=outfile_plotData) #We want the EU ref. data to be its own group. 

print("\n\nDone saving the plot data in 'result_HLAFrequencyPlotData.txt'...")


outfile_plotData.close()




