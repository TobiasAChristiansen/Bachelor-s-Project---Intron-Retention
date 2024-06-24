#!/usr/bin/env python3

import sys
import math
import time

try:
	infile = open("data_HLAfreq_EU_field1.txt","r") #New file, as of V4. 
except IOError as errText:
	print(errText)
	sys.exit()


HLAdict = dict()

for line in infile:
	lineList = line.strip().split(sep="\t")

	if lineList[0].replace(",","").isdigit(): #Remove the comma. Since some of the data looks like this '1,605'. 

		lineNum = lineList[0].replace(",","")
		HLAname = lineList[1]
		populationInfo = lineList[2]
		irrelevant_alleleFrequencyInWholeStudy = lineList[3]
		alleleFrequency = lineList[4]

		if alleleFrequency:
			alleleFrequency = float(alleleFrequency)
		else:
			alleleFrequency = 0 #see entry '3,320' for 'England Lancaster' - it is missing its frequency, but it has sample size! In such case, we assume that the frequency is 0. 

		sampleSize = float(lineList[5].replace(",","")) #This will be important for getting weighted average. 
		location = lineList[6]

		if HLAname:
			
			if not HLAname in HLAdict.keys():
				HLAcount = 1 #How many times a study has been performed searching for this HLA subtype. If there e.g. are 30 studies for HLA*XX:YY, we will need '30' to calculate the average frequency, by dividing that number with the summed frequency. 
				
				N_studies = 1 #Amount of studies where this HLA is present. 

				sum_freq2 = alleleFrequency**2 #Aka. (freq)^2

				sampleSum = sampleSize
				sizeFreqProductSum = sampleSize*alleleFrequency #Will be needed for weighted average

				HLAdict[HLAname] = [alleleFrequency, HLAcount, N_studies,sum_freq2, sampleSum, sizeFreqProductSum]

			else:
				HLAcount = HLAdict[HLAname][1] + 1
				summedAlleleFrequency = HLAdict[HLAname][0] + alleleFrequency
				N_studies = HLAdict[HLAname][2] + 1
				sum_freq2 = HLAdict[HLAname][3] + alleleFrequency**2
				sampleSum = HLAdict[HLAname][4] + sampleSize
				sizeFreqProductSum = HLAdict[HLAname][5] + sampleSize*alleleFrequency

				HLAdict[HLAname] = [summedAlleleFrequency, HLAcount, N_studies, sum_freq2, sampleSum, sizeFreqProductSum]

infile.close()

#Calculating the average HLA frequency and standard deviation for all the HLA alleles. 

for HLAname, valueList in HLAdict.items():

	summedAlleleFrequency, HLAcount, N_studies, sum_freq2, sampleSum, sizeFreqProductSum = valueList[0], valueList[1], valueList[2], valueList[3], valueList[4], valueList[5]

	try:
		#Normal average:
		#averageFrequency = float(summedAlleleFrequency)/float(HLAcount)

		#Weighted average:
		averageFrequency = float(sizeFreqProductSum)/float(sampleSum)

	except:
		raise ValueError(f"Improper values for calculating average frequency: {sizeFreqProductSum}/{sampleSum}")

	if float(averageFrequency) < 0: #If the average is negative - something went wrong!
		raise ValueError("Error! Negative average:", HLAname, averageFrequency)


	HLAdict[HLAname].append(averageFrequency) #Add the average to the list in the dict which stores all the HLA subtypes. 


#############################################################################################################

# Repeat code - we need to calculate the standard deviation:

infile = open("data_HLAfreq_EU_field1.txt","r") #V4. 

for line in infile:
	lineList = line.strip().split(sep="\t")

	if lineList[0].replace(",","").isdigit(): #Remove the comma. Since some of the data looks like this '1,605'. 

		lineNum = lineList[0].replace(",","")
		HLAname = lineList[1]
		populationInfo = lineList[2]
		irrelevant_alleleFrequencyInWholeStudy = lineList[3]
		alleleFrequency = lineList[4]

		if alleleFrequency:
			alleleFrequency = float(alleleFrequency)
		else:
			alleleFrequency = 0 #see entry '3,320' for 'England Lancaster' - it is missing its frequency, but it has sample size! In such case, we assume that the frequency is 0. 

		sampleSize = float(lineList[5].replace(",","")) #This will be important for getting weighted average. 
		location = lineList[6]

		if HLAname:
			
			if HLAname in HLAdict.keys():
				
				averageFrequency = HLAdict[HLAname][6]

			else:
				raise ValueError("Second time reading the file, should not lead to discovery of new HLA in the file!", HLAname)
			

			if len(HLAdict[HLAname]) == 8:

				sumValue = (alleleFrequency-averageFrequency)**2 #The sum value that will be used for calculating the standard deviation.
				

				HLAdict[HLAname][-1] += sumValue

				#if HLAname == "A*01:01":
					#print("(!)Sum value:", sumValue)
					#print((alleleFrequency-averageFrequency))
					#print("freq:", alleleFrequency)
					#print("Avrg:", averageFrequency)
					#print("Sum, so far:", HLAdict[HLAname][-1])

			else:
				sumValue = (alleleFrequency-averageFrequency)**2
				
				HLAdict[HLAname].append(sumValue)


				#if HLAname == "A*01:01":
					#print("(!)First - Sum value:", sumValue)
					#print("freq:", alleleFrequency)
					#print("Avrg:", averageFrequency)
					#print("Sum, so far:", HLAdict[HLAname][-1])


infile.close()



#Calculating the standard deviation: 

averageHLAdict = dict()


for HLAname, valueList in HLAdict.items():

	summedAlleleFrequency, HLAcount, N_studies, sum_freq2, sampleSum, sizeFreqProductSum, averageFrequency, sumValue = valueList[0], valueList[1], valueList[2], valueList[3], valueList[4], valueList[5], valueList[6], valueList[7]

	if HLAname == "DPB1*04":
		print("Sum value:",sumValue)
		print("N_studies:",N_studies)
		time.sleep(2)
		

	if (N_studies-1) <= 0: #If it is 0 or negative, something is wrong. 
		standardDeviation = 0
	else:
		standardDeviation = math.sqrt((sumValue)/(N_studies-1))

		if standardDeviation < 0: #Standard deviation should not be negative!
			raise ValueError("Negative standard deviation", HLAname, standardDeviation)
			sys.exit()



	if not HLAname in averageHLAdict.keys():
		averageHLAdict[HLAname] = [averageFrequency,standardDeviation]


	else:
		raise ValueError("Already have this HLA type in the dict - can't calculate the average for 2 seperate entries. Try to merge them. Something went wrong", HLAname, averageHLAdict[HLAname])
		sys.exit()




#Finally, save the average frequencies in a file:

outfile = open("result_averageHLAfreq_EU_field1.txt","w")

for HLAname, averageAlleleFrequency_std_List in averageHLAdict.items():
	averageAlleleFrequency, standardDeviation = averageAlleleFrequency_std_List[0], averageAlleleFrequency_std_List[1]

	#print(f"{HLAname}\t{averageAlleleFrequency}\t{standardDeviation}")
	print(f"{HLAname}\t{averageAlleleFrequency}\t{standardDeviation}",file=outfile)


outfile.close()

print("Done - Wrote the results in a file named: 'result_averageHLAfreq_EU_field1.txt'")






