
# coding: utf-8

# #Part 1.1
# 
# read1 = R1_test_file.fastq.gz (sequence)
# index1 = R2_test_file.fastq.gz (barcode)
# 
# index2 = R3_test_file.fastq.gz (barcode)
# read2 = R4_test_file.fastq.gz (sequence)

# In[53]:


#Part 1.2 Generate per base call distribution of quality scores for each read

import numpy as np
import gzip 

read1 = "/projects/bgmp/ahowells/multiplexing/R1_test_file.fastq.gz"

#create a function to convert phred scores in files
def convert_phred(letter):
    """Converts a single character into a phred score"""
    phred_scores = ord(letter) - 33 
    return phred_scores

#initialize an array 
read1_array = np.zeros((101),dtype=float)
read1_array[0]
    

with gzip.open(read1,"rt") as fh:
    LN = 0
    #loop through every line and select
    for line in fh: 
        line = line.strip('\n')
            #print(line)
            # this loop evaluates every 4th line
        if LN % 4 == 3: 
            #print(line)
            index = 0  
            #for every letter in the 4th line of the FASTQ file
            for letter in line:
                # convert every letter to a numeric qual. score by calling the convert_phred function created
                qscore = (convert_phred(letter)) 
                #print(qscore)
                # for every index in the array function, add ongoing sum while looping
                read1_array[index] += qscore 
                index +=1
        LN +=1
            

mean = np.zeros(((101)), dtype=float)
mean = read1_array/(LN/4)
#print(mean)

print("# Base Pair\tMean Quality Score")
for j in range(101):
    print (j,mean[j], sep ='\t')


# In[55]:



plt.plot(read1_array)
plt.ylabel('Mean Quality Scores')
plt.xlabel('Base pairs')
plt.title ('Mean Quality Scores per Base Pairs')
plt.savefig('Read1_qualityscores')

