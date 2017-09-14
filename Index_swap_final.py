#!/usr/bin/python


import argparse


parser = argparse.ArgumentParser(description="index swapping assessment")
parser.add_argument("-q", "--min_qscore", help='min q score allowed',required=True,type= int)
parser.add_argument("-f1", "--file1", help='sequence file 1',required=True,type= str)
parser.add_argument("-f2", "--file2", help='index file 2',required=True,type= str)
parser.add_argument("-i", "--index_file", help='file with possible indexes',required=True,type= str)
parser.add_argument("-o", "--output_directory", help='output directory for unknown dictionary',required=False,type= str)
parser.add_argument("-o2", "--output_directory2", help='output directory for known dictionary',required=True,type= str)

args = parser.parse_args()

file1=args.file1
file2=args.file2

q=args.min_qscore
i=args.index_file
out=args.output_directory 
out2=args.output_directory2

indexs=[]
index_comb=[]
index_dict={}
unk_index_dict={}
N_index_dict={}
print('beginning')
with open(file1,'r') as fh, open (file2, 'r') as fh2, open(i, "r") as fh3:
    def rc (my_sequence): #define a function to generate the reverse compliment of the R3 reads
    	
    	my_dictionary = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N'}
    	
    	return "".join([my_dictionary[base] for base in reversed(my_sequence)])
    
    for a in fh3: #generate a dictionary with the known index pairs that we expect to have 
        
        a=a.strip('\n').split('\t')
      
        indexs.append(a[4])
 
        
    for b in indexs[1:len(indexs)]:
        
        index_pair=(b+'_'+b)
        index_dict[index_pair]=0

    i=0
    for z,x in zip(fh,fh2):

        i+=1
        x=x.strip('\n')
        z=z.strip('\n') 
        if i % 5000000 == 0: # print if the line is divisible by 5,000,000 as a sanity check
            print('working on' + str(i), flush=True)
                         #INDEX READ 2
        if i % 4==2:
            save_me=z #grab index sequence from file 1
            save_me2=x  #grab index sequence from file 2
        if i % 4 ==0:
            
            
            sum_phred=0
            sum_phred2=0
#       
             
            for c in z:
                    
                phred=ord(c)-33 #generate the phred score for each character
                if phred >= q: #check if phred score is above the quality cut off
              
                    sum_phred+=1 #add one to the counter if it is
                else:
                    sum_phred=0 #set the counter to 0 if it isn't
                
       
            
            for ch in x:
                
                phred2=ord(ch)-33 #generate the phred score for each character
                if phred2 >= q: #check if phred score is above the quality cut off
                   
                    sum_phred2+=1 #add one to the counter if it is
                else:
                    
                    sum_phred2=0 #set the counter to 0 if it isn't

            if sum_phred2 == 8 and sum_phred == 8: #if the two counters from above are equal to 8, which means all 
            #8 nucleotides passed, then continue
                
                thing=(save_me+'_'+rc(save_me2)) #save the index sequence from R2 and reverse compliment R3
                
                
                if thing in index_dict: #if the index pair is in the known index dictionary, add one
                    
                    index_dict[thing]+=1
           
                    
                else: #if it isn't, check if its already in the unknown index dictionary and if it isn't add it or
                #increment the counter if it is
                	if thing in unk_index_dict: 
                		unk_index_dict[thing]+=1
                	else:
                		unk_index_dict[thing]=1
                	
                	
 
print('writing unknown pairs to file') #write the unknown pairs to an output
for k,v in unk_index_dict.items():
    o=open(out, 'a+')
    o.write(str(k)+'\t'+str(v)+'\n')
#print(unk_index_dict)
print('writing known pairs to file') #write the known pairs to an output
for k,v in index_dict.items():
    o=open(out2, 'a+')
    o.write(str(k)+'\t'+str(v)+'\n')
              
        