import Bio
from Bio import SeqIO
import pandas as pd 
import numpy as np 
import re
import sys
import os
import os.path
import datetime

"""
function that identifies the position of stretches of N in a sequence

find the start and end position of a stretch of N for a sequence
parameters are;
motif (number of consecutive N from the user input)
sequence (sequence to look into)
"""
def find_stretch_N(motif, sequence):
    df_pos_N=pd.DataFrame(columns=['start', 'end'])
    #find the occurence of the motif in the sequence, 
    #once a motif is found, the algorythm starts at the end of the motif
    str_sequence=str(sequence.seq)
    end_last_stretch=-1
    for match in re.finditer(motif, str_sequence):
        start_pos = match.start() 
        
        #start at the end of last stretch
        if (start_pos>end_last_stretch):
            end_pos = match.end() 
            #more than X consecutive N
            while(end_pos<len(sequence) and sequence[end_pos]=="N" ):
                end_pos+=1
            end_last_stretch=end_pos
            #append the dataframe of N stretches
            df_pos_N=df_pos_N.append({'start':start_pos, 'end':end_pos},ignore_index=True)

    return df_pos_N


        
""" 
function that return the list of flagged SNP whithin a list of SNP position or whithin the whole sequence

find the SNP around the stretches of N that are in the interval
parameters are;
nuc_interval (number of snp around stretch)
df_pos_N (data frame of start-end of stretch of N for a sequence)
list_SNP (list of position of SNP from a VCF in bed format)
sequence (sequence to look into)
"""
def find_mutation_not_SNP(nuc_interval, df_pos_N, list_SNP, reference, sequence):
        
    #list of flagged position to return
    df_flag=pd.DataFrame(columns=['sequence id', 'pos', 'dist from N', 'protected', 'nucleotide flag'])
    #creates a df for sequence and their respective number of flagged sites
    df_nb_flag=pd.DataFrame(columns=['sequence id', 'number of flagged sites'])   
    nb_flag=0
    
    #sequence to look at in list format
    print(sequence.id)
    list_seq=list(sequence.seq)
    
    #list of position to eliminate duplicates
    pos_all=list()
    
    #loop through the rows of the stretch of N
    for index, row in df_pos_N.iterrows():
               
        start_N=row['start']
        end_N=row['end']
        i=1
        
        #get the Y number nucleotides before and after
        while (i<=nuc_interval):
            
            pos_before=start_N - i
            pos_after=end_N + i
            
            #if the pos_before is before the beginning of the sequence
            if(pos_before>=0 and (pos_before not in pos_all)):
                #nuc at the position
                nuc_before=list_seq[pos_before]
            
                #compare to the reference sequence to see if this position is mutated
                #if the pos_end is after the end of the sequence or pos_before is before the start
                #add the pos in the list of flag SNP if they aren't a SNP position
                if(nuc_before!=reference[pos_before] and nuc_before!="N"):
                    if(pos_before in list_SNP):
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_before+1), 'dist from N':i, 'protected':0, 'nucleotide flag':nuc_before},ignore_index=True)
                        pos_all.append(pos_before)
                        nb_flag+=1
                        
            #if the pos_end is after the end of the sequence
            if(pos_after<len(sequence)  and (pos_after not in pos_all)):
                
                #nuc at the position
                nuc_after=list_seq[pos_after]
                
                #if the pos_end is after the end of the sequence or pos_before is before the start
                #add the pos in the list of flag SNP if they aren't a SNP position
                if(nuc_after!=reference[pos_after] and nuc_after!="N"): 
                    if(pos_after in list_SNP):
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_after+1), 'dist from N':i, 'protected':0, 'nucleotide flag':nuc_after},ignore_index=True)
                        pos_all.append(pos_after)
                        nb_flag+=1
            
            #iterate
            i+=1
        
        #return a list of the sorted position of the SNP flagged and remove duplicates
        sorted_flag=df_flag.sort_values(by=['pos'])
        #rows_before = sorted_flag.count
        #sorted_flag.drop_duplicates(subset =['sequence id', 'pos'], keep = 'first', inplace = True) 
        #rows_after = sorted_flag.count
        #deduce the number of duplicates removed for the number of flag
        #nb_flag=nb_flag -(rows_before - rows_after)
        
    df_nb_flag=df_nb_flag.append({'sequence id':sequence.id, 'number of flagged sites':nb_flag}, ignore_index=True)
    
    return sorted_flag, df_nb_flag



"""
function return the list of flagged SNP that aren't in the list of protected positions
it also returns the list of flagged SNP that are in the list of protected positions for statistics

parameters are;
nuc_interval (number of snp around stretch)
df_pos_N (data frame of start-end of stretch of N for a sequence)
list_protected_SNP (list of position of SNP from a VCF in bed format that are correct/protected)
sequence (sequence to look into)

"""
def find_mutation_not_protected_SNP(nuc_interval, df_pos_N, df_protected_SNP, reference, sequence):
    
    list_protected_SNP=df_protected_SNP['START'].tolist()
    
    #list of flagged position to return
    df_flag=pd.DataFrame(columns=['sequence id', 'pos', 'dist from N', 'protected', 'nucleotide flag'])
    #list of protected SNP that are close to the stretch of N
    protected_flag=[]
    print(sequence.id)
    #creates a df for sequence and their respective number of flagged sites
    df_nb_flag=pd.DataFrame(columns=['sequence id', 'number of flagged sites'])
    #count the number of flagged position for this sequence
    nb_flag=0
    
     #list of position to eliminate duplicates
    pos_all=list()
    
    #sequence to look at in list format
    list_seq=list(sequence.seq)
    
    #loop through the rows of the stretch of N
    for index, row in df_pos_N.iterrows():
        
        start_N=row['start']
        end_N=row['end']
        i=1
        
        #get the Y number nucleotides before and after
        while (i<=nuc_interval):
            
            pos_before=start_N - i
            pos_after=end_N + i
            
            #if the pos_before is before the beginning of the sequence
            if(pos_before>=0 and (pos_before not in pos_all)):
                #nuc at the position
                nuc_before=list_seq[pos_before]
                pos_all.append(pos_before)
                #compare to the reference sequence to see if this position is mutated
                #if the pos_end is after the end of the sequence or pos_before is before the start
                #add the pos in the list of flag SNP if they aren't a SNP position
                if(nuc_before!=reference[pos_before] and nuc_before!="N"):
                    if(pos_before in list_protected_SNP):
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_before+1), 'dist from N':i, 'protected':1, 'nucleotide flag':nuc_before},ignore_index=True)
                        protected_flag.append(pos_before+1)
                    else:
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_before+1), 'dist from N':i, 'protected':0, 'nucleotide flag':nuc_before},ignore_index=True)
                        nb_flag+=1
            
            
            #if the pos_end is after the end of the sequence
            if(pos_after<len(sequence)  and (pos_after not in pos_all)):
                
                #nuc at the position
                nuc_after=list_seq[pos_after]
                pos_all.append(pos_after)

                #if the pos_end is after the end of the sequence or pos_before is before the start
                #add the pos in the list of flag SNP if they aren't a SNP position
                if(nuc_after!=reference[pos_after] and nuc_after!="N"): 
                    if(pos_after in list_protected_SNP):
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_after+1), 'dist from N':i, 'protected':1, 'nucleotide flag':nuc_after},ignore_index=True)
                        protected_flag.append(pos_after+1)
                    else:
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_after+1), 'dist from N':i, 'protected':0, 'nucleotide flag':nuc_after},ignore_index=True)
                        nb_flag+=1

            
            
            #iterate
            i+=1
        
        #return a list of the sorted position of the SNP flagged and protected, also remove duplicates
        sorted_flag=df_flag.sort_values(by=['pos'])
        #sorted_flag.drop_duplicates(subset =['sequence id', 'pos'], keep = 'first', inplace = True) 
        protected_flag=list(set(protected_flag))
        protected_flag=sorted(protected_flag)
            
    df_nb_flag=df_nb_flag.append({'sequence id':sequence.id, 'number of flagged sites':nb_flag},ignore_index=True)
            
    return sorted_flag, protected_flag, df_nb_flag
    
    
    
"""
function return the list of flagged SNP that aren't in the list of protected positions and are different than the specified allele
it also returns the list of flagged SNP that are in the list of protected positions for statistics

parameters are;
nuc_interval (number of snp around stretch)
df_pos_N (data frame of start-end of stretch of N for a sequence)
list_protected_SNP (list of position of SNP from a VCF in bed format that are correct/protected)
sequence (sequence to look into)
"""   
def find_mutation_not_protected_allele_SNP(nuc_interval, df_pos_N, df_protected_SNP, reference, sequence):
        
    list_protected_SNP=df_protected_SNP['START'].tolist()
    #list of flagged position to return
    df_flag=pd.DataFrame(columns=['sequence id', 'pos', 'dist from N', 'protected', 'nucleotide flag'])
    #list of protected SNP that are close to the stretch of N
    protected_flag=[]
    
    #creates a df for sequence and their respective number of flagged sites
    df_nb_flag=pd.DataFrame(columns=['sequence id', 'number of flagged sites'])
    #count the number of flagged position for this sequence
    nb_flag=0
    
    #sequence to look at in list format
    list_seq=list(sequence.seq)
    
     #list of position to eliminate duplicates
    pos_all=list()
    
    #loop through the rows of the stretch of N
    for index, row in df_pos_N.iterrows():
        
        start_N=row['start']
        end_N=row['end']
        i=1
        
        #get the Y number nucleotides before and after
        while (i<=nuc_interval):
            
            pos_before=start_N - i
            pos_after=end_N + i
            
            #if the pos_before is before the beginning of the sequence and isn't a duplicate
            if(pos_before>=0 and (pos_before not in pos_all)):
                #nuc at the position
                nuc_before=list_seq[pos_before]
            
                #compare to the reference sequence to see if this position is mutated
                #if the pos_end is after the end of the sequence or pos_before is before the start
                #add the pos in the list of flag SNP if they aren't a SNP position
                if(nuc_before!=reference[pos_before] and nuc_before!="N"):
                    
                    pos_all.append(pos_before)
                    
                    #allele specific protection                          
                    if(pos_before in list_protected_SNP):
                        row_pos = df_protected_SNP.loc[df_protected_SNP['START'] == pos_before]
                        row_pos_allele = row_pos['MAJ_ALLELE']
                        tmp_row_pos_allele=list(row_pos_allele)
                        str_row_pos_allele=str(tmp_row_pos_allele[0])
                        list_row_pos_allele=str_row_pos_allele.split(',')
                        #if it isn't corresponding to the alternate allele 
                        if(nuc_before not in list_row_pos_allele):
                            df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_before+1), 'dist from N':i, 'protected':1, 'nucleotide flag':nuc_before},ignore_index=True)
                            protected_flag.append(pos_before+1)
                            nb_flag+=1
                        #if the nuc is one of the alternative alleles
                        if(nuc_before in list_row_pos_allele):
                            df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_before+1), 'dist from N':i, 'protected':2, 'nucleotide flag':nuc_before},ignore_index=True)
                            protected_flag.append(pos_before+1)
                            #nb_flag+=1
                    else:
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_before+1), 'dist from N':i, 'protected':0, 'nucleotide flag':nuc_before},ignore_index=True)
                        nb_flag+=1
            
            
            #if the pos_end is after the end of the sequence
            if(pos_after<len(sequence)  and (pos_after not in pos_all)):
                
                #nuc at the position
                nuc_after=list_seq[pos_after]
        
                #if the pos_end is after the end of the sequence or pos_before is before the start
                #add the pos in the list of flag SNP if they aren't a SNP position
                if(nuc_after!=reference[pos_after] and nuc_after!="N"): 
                    
                    pos_all.append(pos_after)
                    
                    #if(pos_after not in list_protected_SNP):
                    if(pos_after in list_protected_SNP):
                        row_pos = df_protected_SNP.loc[df_protected_SNP['START'] == pos_after]
                        row_pos_allele = row_pos['MAJ_ALLELE']
                        tmp_row_pos_allele=list(row_pos_allele)
                        str_row_pos_allele=str(tmp_row_pos_allele[0])
                        list_row_pos_allele=str_row_pos_allele.split(',')
                        
                        #if it isn't corresponding to the alternate allele 
                        if(nuc_after not in list_row_pos_allele):
                            df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_after+1), 'dist from N':i, 'protected':1, 'nucleotide flag':nuc_after},ignore_index=True)
                            protected_flag.append(pos_after+1)
                            nb_flag+=1
                        #if the nuc is one of the alternative alleles
                        if(nuc_after in list_row_pos_allele):
                            df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_after+1), 'dist from N':i, 'protected':2, 'nucleotide flag':nuc_after},ignore_index=True)
                            protected_flag.append(pos_after+1)
                            #nb_flag+=1
                    else:
                        df_flag=df_flag.append({'sequence id':sequence.id, 'pos':(pos_after+1), 'dist from N':i, 'protected':0, 'nucleotide flag':nuc_after},ignore_index=True)
                        nb_flag+=1
            #iterate
            i+=1
        
        #return a list of the sorted position of the SNP flagged and protected, also remove duplicates
        sorted_flag=df_flag.sort_values(by=['pos'])
        #sorted_flag.drop_duplicates(subset =['sequence id', 'pos'], keep = 'first', inplace = True) 
        protected_flag=list(set(protected_flag))
        protected_flag=sorted(protected_flag)
            
    df_nb_flag=df_nb_flag.append({'sequence id':sequence.id, 'number of flagged sites':nb_flag},ignore_index=True)
            
    return sorted_flag, protected_flag, df_nb_flag
    
    
    
"""
function returns the current sequence with the flagged SNP masked (replace by N) in the sequence

parameters are;
flag_SNP list of SNP that are flagged
sequence current fasta sequence in analysis
"""
def out_fasta(flag_SNP, sequence):
    
    str_sequence=str(sequence.seq)
    #return to a 0 based list
    flag_SNP.pos=(flag_SNP.pos-1)
    #replace the flaged SNP by N
    for index, row in flag_SNP.iterrows():
        pos=row['pos']
        #do not replace the protected SNP 
        if(row['protected']==0):
            str_sequence=str_sequence[:pos]+"N"+str_sequence[pos+1:]
            
    return str_sequence
    
    
"""
function that writes the dataframe of SNP that were flag and if they are protected or not

parameters are;
flag_SNP dataframe and output file name
"""
def write_flag(output_file, flagged_SNP):
    
    #return to a 1 based sequence
    flagged_SNP.pos=(flagged_SNP.pos+1)
    form='%s', '%d', '%d', '%d', '%s'
    f=open(output_file+'.txt','ab')
    np.savetxt(f, flagged_SNP, fmt=form, delimiter="\t")
    f.close()
    
    
"""
function write in a file the current sequence with the flagged SNP masked (replace by N) in the sequence

parameters are;
flag_SNP list of SNP that are flagged
Sequence that is currently beeing analyse
sequence with replacement from function out_fasta
"""
def write_fasta(output_file, sequence, masked_sequence):
    
    #write the fasta file with the masked position
    with open(output_file+".fasta", "a") as out_fasta_file:
        out_fasta_file.write(">"+sequence.id + "\n" + masked_sequence + "\n")
        
"""
function that writes the stats for the position that were protected in an output file

parameters are;
flagged_protected_SNP list of SNP that are flagged and protected
sequence current fasta sequence in analysis
name of the output file
"""
def write_protected_stats(output_file, sequence, flagged_protected_SNP):
    
    with open(output_file+"protected.stats", "a") as stats_file:
        stats_file.write(sequence.id+ " as "+ str(len(flagged_protected_SNP))+ "protected SNP within the interval to the stretches of N.\n")
        for pos in flagged_protected_SNP:
            stats_file.write(sequence.id + "\t")
            stats_file.write("%s\n" % pos)

"""
function that writes the number of SNP flag in each sequence in a file

parameters are;
df_nb_flag number of SNP flag per sequence (sequence id, number of snp flag)
name of the output file
"""            
def write_nb_flag_per_seq(output_file, df_nb_flag):

    form='%s', '%d'
    f=open(output_file+'_nb_flag_per_seq.txt','ab')
    np.savetxt(f, df_nb_flag, fmt=form, delimiter="\t")
    f.close()
    
    
#MAIN
def main():
    
    number_N=int(input("Minimum number of consecutive N: "))
    nuc_interval=int(input("Interval of nucleotide to look at around the stretch of N: "))

    motif="N" * number_N
    print("motif to find "+motif)

    #fasta of the alignement
    fasta_file=input("path and filename for the alignement file in this format: USER/EXEMPLE/file.fasta ->")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    #Reference fasta file
    reference_file=input("path and filename for the reference sequence file in this format: USER/EXEMPLE/reference.fasta ->")
    reference_fasta = list(SeqIO.parse(reference_file, "fasta"))
    reference_seq=list(reference_fasta[0].seq)
    
    #bed file of the positions of the SNPs
    bed_file=input("path and filename for the bed file of all SNP in this format: USER/EXEMPLE/file.bed \n leave empty if you want the program to look at every positions ->")
    if(bed_file!=""):
        df_SNP=pd.read_table(bed_file, sep='\t', header=None)
        df_SNP.columns = ['CHROM', 'START', 'END', 'MAJ_ALLELE']
        list_SNP=df_SNP['START'].tolist()
        fct="SNP"
    else:
        fct=""
        
    protected_file=input("path and filename for the bed of protected SNP file in this format: USER/EXEMPLE/protected.bed \n leave empty if you don't have a list of SNP that you want to protect during the analyses ->")
    if(protected_file!=""):
        df_protected_SNP=pd.read_table(protected_file, sep='\t', header=None)        
        if(len(df_protected_SNP.columns)==3):
            df_protected_SNP.columns = ['CHROM', 'START', 'END']
            fct="protected"
        if(len(df_protected_SNP.columns)==4):
            df_protected_SNP.columns = ['CHROM', 'START', 'END', 'MAJ_ALLELE']
            fct="protected_allele"
        if(len(df_protected_SNP.columns)>4):
            print ("Invalid number of columns in protected bed file. Expected 3 or 4 column ('CHR', 'START', 'END', *'Major allele'*)")
            sys.exit(0)
        
        #check to see if any columns are duplicated in the bed file
        duplicates_row=df_protected_SNP.duplicated(subset=None, keep=False)
        list_duplicates_row=list(duplicates_row)
        #if any line is duplicated (whole line)
        if(True in list_duplicates_row):
            print("*************Your bed file contains duplicated lines. multi-allelic positions should be listed in one line as Allele1,Allele2*************")

    if(bed_file=="" and protected_file==""):
        fct=""
    if(bed_file!="" and protected_file!=""):
        print("*************You can not give a list of SNP and a list of SNP to protect*************")
        sys.exit(0)
    
    
    #name the output file
    output_file=input("path and filename for the output prefix  in this format: USER/EXEMPLE/output ->")
    
    
    #time 
    begin_time = datetime.datetime.now()
    
    
    #check if a file doesn't already exists and if it does, it deletes it
    if(os.path.isfile(output_file + ".txt")):
        print("\n")
        print("*************File "+ output_file + ".txt" " already exist. File was deleted*************")
        print("\n")
        #sys.exit(0)
        os.remove(output_file+".txt")
    if(os.path.isfile(output_file + ".fasta")):
        print("\n")
        print("*************File "+ output_file + ".fasta" + " already exist. File was deleted*************")
        print("\n")
        #sys.exit(0)
        os.remove(output_file+".fasta")
    if(os.path.isfile(output_file + ".stats")):
        print("\n")
        print("*************File "+ output_file + ".stats"+ " already exist. File was deleted*************")
        print("\n")
        #sys.exit(0)
        os.remove(output_file + ".stats")
    if(os.path.isfile(output_file + "_nb_flag_per_seq.txt")):
        print("\n")
        print("*************File "+ output_file + "_nb_flag_per_seq.txt"+ " already exist. File was deleted*************")
        print("\n")
        #sys.exit(0)
        os.remove(output_file + "_nb_flag_per_seq.txt")
    
    #loop through all the sequences in the fasta file 
    for sequence in records:
        
        #check if the sequence has the right length or stop
        
        if(len(sequence.seq)!= 29903):
            print("*************               ERROR : sequence ", sequence.id, " has the wrong lenght ( reference = 29903).          ******************")
            sys.exit(0)
        
        
        #find the stretch of at least X (parameter) N.
        pos_stretch=find_stretch_N(motif, sequence)
        #make sure the dataframe isn't empty
        if(pos_stretch.empty):
            flagged_SNP=pd.DataFrame(columns=['sequence id', 'pos', 'dist from N', 'protected', 'nucleotide flag'])
            df_nb_flag=pd.DataFrame(columns=['sequence id', 'number of flagged sites'])
            df_nb_flag=df_nb_flag.append({'sequence id':sequence.id , 'number of flagged sites': 0},ignore_index=True)
            masked_sequence=out_fasta(flagged_SNP, sequence)
            write_fasta(output_file, sequence, masked_sequence)
            write_nb_flag_per_seq(output_file, df_nb_flag)
            
        else:
            #find the SNP that are flagged for this sequence
            if(fct=="protected"):
                flagged_SNP, flagged_protected_SNP, df_nb_flag=find_mutation_not_protected_SNP(nuc_interval, pos_stretch, df_protected_SNP, reference_seq, sequence)
                #write_protected_stats(output_file, sequence, flagged_protected_SNP)
            if(fct=="protected_allele"):
                flagged_SNP, flagged_protected_SNP, df_nb_flag=find_mutation_not_protected_allele_SNP(nuc_interval, pos_stretch, df_protected_SNP, reference_seq, sequence)
                #write_protected_stats(output_file, sequence, flagged_protected_SNP)
            if(fct=="SNP"):
                flagged_SNP, df_nb_flag=find_mutation_not_SNP(nuc_interval, pos_stretch, list_SNP, reference_seq, sequence)
            if(fct==""):
                list_SNP_all=range(len(reference_seq))
                flagged_SNP, df_nb_flag=find_mutation_not_SNP(nuc_interval, pos_stretch, list_SNP_all, reference_seq, sequence)
            

            masked_sequence=out_fasta(flagged_SNP, sequence)
            #write the fasta sequences with masked SNP
            write_fasta(output_file, sequence , masked_sequence)
            #write the results into a file
            write_flag(output_file, flagged_SNP) 
            write_nb_flag_per_seq(output_file, df_nb_flag)
            
        
    #time it has taken
    total_time=(datetime.datetime.now() - begin_time)
    print("execution time: " + str(total_time))
        

if __name__ == "__main__":
    main()
    
    
