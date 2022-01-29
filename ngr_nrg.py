from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import regex
import sys
import pandas as pd
import os

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])
def base_change(seq,variation1,position1,position2,strand):
    seq_list=list(seq)
    if strand == "+":
     seq_list[position1]=variation1   
     if seq_list[position2]  == 'T':
       seq_list[position2] ='A'
       sgrna="".join(seq_list)
       return sgrna
     elif seq_list[position2] == 'C':
       seq_list[position2] ='G'
       sgrna="".join(seq_list)
       return sgrna
     elif seq_list[position2] == 'A':
       seq_list[position2] ='T'
       sgrna="".join(seq_list)
       return sgrna
     elif seq_list[position2] == 'G':
       seq_list[position2] ='C'
       sgrna="".join(seq_list)
       return sgrna        
    elif strand == "-":
        seq_list[position1]=reverse_complement(variation1)
        if seq_list[position2]  == 'T':
         seq_list[position2] ='A'
         sgrna="".join(seq_list)
         return sgrna
        elif seq_list[position2] == 'C':
         seq_list[position2] ='G'
         sgrna="".join(seq_list)
         return sgrna
        elif seq_list[position2] == 'A':
         seq_list[position2] ='T'
         sgrna="".join(seq_list)
         return sgrna
        elif seq_list[position2] == 'G':
         seq_list[position2] ='C'
         sgrna="".join(seq_list)
         return sgrna        

    


# In[ ]:


def ontarget(genome,vcf_ls):

   lst=[]
   vcf_sgrna_ls=[]
   records = SeqIO.to_dict(SeqIO.parse(open(genome), 'fasta'))
   for line1 in vcf_ls:
    #print (line1[0])
    long_seq_record = records[line1[0]]
    long_seq = long_seq_record.seq
    #alphabet = long_seq.alphabet
    short_seq1=str(long_seq)[int(line1[1])-1]
    #short_seq1_tmp=str(long_seq)[int(line1[1])]
    if short_seq1 == line1[3]:
     #print line1[3], "match", short_seq1
     short_seq2_6 = str(long_seq)[int(line1[1])-1:int(line1[1])+4] ####mutation at 2nd position and random at 6
     short_seq6_2 = str(long_seq)[int(line1[1])-1:int(line1[1])+8] ####mutation at 6th position and random at 2  
     short_seq16_19 = str(long_seq)[int(line1[1])-1:int(line1[1])+18] #####mutation at 16th position and random at 19
     short_seq19_16 = str(long_seq)[int(line1[1])-1:int(line1[1])+21] #####mutation at 19th position and random at 16    
     short_seq18 = str(long_seq)[int(line1[1])-1:int(line1[1])+20] #####mutation at 18th position
     short_seq17 = str(long_seq)[int(line1[1])-1:int(line1[1])+19] #####mutation at 17th position
     short_seq3_4 = str(long_seq)[int(line1[1])-3:int(line1[1])] ####mutation at 3rd position and random at 4/5 (LwCas13a)
     short_seq11_12 = str(long_seq)[int(line1[1])-11:int(line1[1])] ####mutation at 11th position and random at 12 Cas14a/ssdna  
     short_seq12_11 = str(long_seq)[int(line1[1])-12:int(line1[1])] ####mutation at 11th position and random at 12  Cas14a/ssdna  
     short_seq11_12_pam = str(long_seq)[int(line1[1])-15:int(line1[1])] ####mutation at 11th position and random at 12 Cas14a
     short_seq12_11_pam = str(long_seq)[int(line1[1])-16:int(line1[1])] ####mutation at 11th position and random at 12 Cas14a   
     short_seq1_5 = str(long_seq)[int(line1[1])-4:int(line1[1])] ####mutation at 1nd position and random at 5 (1 and 4)
     short_seq5_1 = str(long_seq)[int(line1[1])-8:int(line1[1])] ####mutation at 5th position and random at 1  (5 and 11)
     short_seq4_1 = str(long_seq)[int(line1[1])-7:int(line1[1])] ####mutation at 4th position and random at 1 (4 and 11)
     short_seq16_19_cas12b = str(long_seq)[int(line1[1])-19:int(line1[1])] ####mutation at 16th position and random at 19
     short_seq19_16_cas12b = str(long_seq)[int(line1[1])-22:int(line1[1])] ####mutation at 19th position and random at 16   
     short_seq11_5 = str(long_seq)[int(line1[1])-14:int(line1[1])] ####mutation at 11th position and random at 5 
     short_seq1 = str(long_seq)[int(line1[1])-5:int(line1[1])]  
     short_seq2 = str(long_seq)[int(line1[1])-6:int(line1[1])]
     short_seq3 = str(long_seq)[int(line1[1])-7:int(line1[1])]
     short_seq4 = str(long_seq)[int(line1[1])-8:int(line1[1])]
     short_seq5 = str(long_seq)[int(line1[1])-9:int(line1[1])]
     short_seq6 = str(long_seq)[int(line1[1])-10:int(line1[1])]
     short_seq7 = str(long_seq)[int(line1[1])-11:int(line1[1])]
     for matchr in regex.finditer(r'(TTT[ACGT]{2})',str(short_seq1),regex.BESTMATCH,overlapped=True):
      #print (matchr)
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TTT[ACGT]{3})',str(short_seq2),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(TTT[ACGT]{4})',str(short_seq3),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(TTT[ACGT]{5})',str(short_seq4),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)       
     for matchr in regex.finditer(r'(TTT[ACGT]{6})',str(short_seq5),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)     
     for matchr in regex.finditer(r'(TTT[ACGT]{7})',str(short_seq6),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)     
     for matchr in regex.finditer(r'(TTT[ACGT]{8})',str(short_seq7),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)     
   
   
     for matchr in regex.finditer(r'(TT[ACGT]{2})',str(short_seq1_5),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{6})',str(short_seq5_1),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{5})',str(short_seq4_1),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{17})',str(short_seq16_19_cas12b),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['AaCas12b']
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{20})',str(short_seq19_16_cas12b),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['AaCas12b']
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{12})',str(short_seq11_5),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
   
     for matchr in regex.finditer(r'(TTTA[ACGT]{12})',str(short_seq12_11_pam),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TTTA[ACGT]{11})',str(short_seq11_12_pam),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)   
        
     for matchr in regex.finditer(r'([ACGT]{12})',str(short_seq12_11),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'([ACGT]{11})',str(short_seq11_12),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)   
   
     for matchr in regex.finditer(r'([ACGT]{3})',str(short_seq3_4),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)

     for matchr in regex.finditer(r'(?:([ACGT]{19}[G][AG]))',str(short_seq18),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{18}[G][AG]))',str(short_seq17),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)

     for matchr in regex.finditer(r'(?:([ACGT]{3}[G][AG]))',str(short_seq2_6),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{7}[G][AG]))',str(short_seq6_2),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
   
     for matchr in regex.finditer(r'(?:([ACGT]{17}[G][AG]))',str(short_seq16_19),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      #print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(?:([ACGT]{20}[G][AG]))',str(short_seq19_16),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      #print (final) 
      lst.append(final)   
     
     ####NAG
     for matchr in regex.finditer(r'(?:([ACGT]{19}AG))',str(short_seq18),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{18}AG))',str(short_seq17),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)

     for matchr in regex.finditer(r'(?:([ACGT]{3}AG))',str(short_seq2_6),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{7}AG))',str(short_seq6_2),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
   
     for matchr in regex.finditer(r'(?:([ACGT]{17}AG))',str(short_seq16_19),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      #print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(?:([ACGT]{20}AG))',str(short_seq19_16),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[1]) + matchr.start() 
      endr=int(line1[1])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      #print (final) 
      lst.append(final)   
 
     
                    
    
   for line2 in lst:
       print (line2)
     #  jk1="\t".join(jk)
      # print (jk[0])
       #jk=jk.rstrip()
       #line2=re.split("\t",jk)
       #line2=re.split("\t",jk)
       #print (line2[0])
       long_seq_record1 = records[line2[0]]
       long_seq1 = long_seq_record1.seq
       if len(line2[5]) == 5 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[1])-5:int(line2[1])+17]
        #print (short_seq18)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         #print (mr1)
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-4# + matchr1.start() 
         endr1=int(line2[1])+17 #+ matchr1.end()
        
         seq_list=list(mseqr1)
         seq_list[4]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',4,'None' ]
        
            
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']

         #print (final1)  
         vcf_sgrna_ls.append(final1)
       
       elif len(line2[5]) == 6 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[1])-6:int(line2[1])+16]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-5# + matchr1.start() 
         endr1=int(line2[1])+16 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[5]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',5,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']

        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1) 
       
       elif len(line2[5]) == 7 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[1])-7:int(line2[1])+15]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-6# + matchr1.start() 
         endr1=int(line2[1])+15 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[6]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',6,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']

        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1) 
       elif len(line2[5]) == 8 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[1])-8:int(line2[1])+14]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-7# + matchr1.start() 
         endr1=int(line2[1])+14 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[7]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',7,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']

        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
       elif len(line2[5]) == 9 and 'LbCas12a' in line2 :
            
        short_seq18=str(long_seq1)[int(line2[1])-9:int(line2[1])+13]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-8# + matchr1.start() 
         endr1=int(line2[1])+13 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[8]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',8,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
       elif len(line2[5]) == 10 and 'LbCas12a' in line2 :
            
        short_seq18=str(long_seq1)[int(line2[1])-10:int(line2[1])+12]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-9# + matchr1.start() 
         endr1=int(line2[1])+12 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[9]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',9,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
       elif len(line2[5]) == 11 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[1])-11:int(line2[1])+11]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-10# + matchr1.start() 
         endr1=int(line2[1])+11 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[10]=line2[4]   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1,jatayu1r,'+',10,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
 
       elif len(line2[5]) == 4 : 
        short_seq23=str(long_seq1)[int(line2[1])-4:int(line2[1])+19]
         ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-3# + matchr1.start() 
         endr1=int(line2[1])+19 #+ matchr1.end()
         ##print (line2[2],line2[3],line2[4])
         jatayu1_5 = base_change(mseqr1,line2[4],3,7,'+')
         #print (line2[0],line2[1],line2[2],line2[3],line2[4],jatayu1_5)   
         jatayu1_5r=  str( Seq(jatayu1_5).reverse_complement())
         jatayu1_4 = base_change(mseqr1,line2[4],3,6,'+')
         jatayu1_4r=  str( Seq(jatayu1_4).reverse_complement())
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1_5,jatayu1_5r,'+',3,7 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu1_4,jatayu1_4r,'+',3,6 ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)   
       elif len(line2[5]) == 8 :
            
        short_seq23=str(long_seq1)[int(line2[1])-8:int(line2[1])+15]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-7# + matchr1.start() 
         endr1=int(line2[1])+15 #+ matchr1.end()
        
         jatayu5_1 = base_change(mseqr1,line2[4],7,3,'+')
         jatayu5_1r=  str( Seq(jatayu5_1).reverse_complement())
         jatayu5_8 = base_change(mseqr1,line2[4],7,10,'+')
         jatayu5_8r=  str( Seq(jatayu5_1).reverse_complement())
         jatayu5_11 = base_change(mseqr1,line2[4],7,13,'+')
         jatayu5_11r=  str( Seq(jatayu5_1).reverse_complement())
         
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu5_1,jatayu5_1r,'+',7,3 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu5_8,jatayu5_8r,'+',7,10 ]
         match_info3 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu5_11,jatayu5_11r,'+',7,13 ]

         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
         final3 =line2 + match_info3 +oligo+oligo1+['AaCas12b']   
            
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
         vcf_sgrna_ls.append(final3)
       
       elif len(line2[5]) == 7 :
            
        short_seq23=str(long_seq1)[int(line2[1])-7:int(line2[1])+16]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-6# + matchr1.start() 
         endr1=int(line2[1])+16 #+ matchr1.end()
        
         jatayu4_1 = base_change(mseqr1,line2[4],6,3,'+')
         jatayu4_1r=  str( Seq(jatayu4_1).reverse_complement())
         jatayu4_11 = base_change(mseqr1,line2[4],6,13,'+')
         jatayu4_11r=  str( Seq(jatayu4_11).reverse_complement())
         jatayu4_16 = base_change(mseqr1,line2[4],6,18,'+')
         jatayu4_16r=  str( Seq(jatayu4_16).reverse_complement())
         
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu4_1,jatayu4_1r,'+',6,3 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu4_11,jatayu4_11r,'+',6,13 ]
         match_info3 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu4_16,jatayu4_16r,'+',6,18 ]   
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
         final3 =line2 + match_info3 +oligo+oligo1   +['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
         vcf_sgrna_ls.append(final3)
        
       elif len(line2[5]) == 19 and 'AaCas12b' in line2 :
            
        short_seq23=str(long_seq1)[int(line2[1])-19:int(line2[1])+4]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-18# + matchr1.start() 
         endr1=int(line2[1])+4 #+ matchr1.end()
        
         jatayu16_19 = base_change(mseqr1,line2[4],18,21,'+')
         jatayu16_19r=  str( Seq(jatayu16_19).reverse_complement())
         jatayu16_4 = base_change(mseqr1,line2[4],18,6,'+')
         jatayu16_4r=  str( Seq(jatayu16_4).reverse_complement())
         
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu16_19,jatayu16_19r,'+',18,21 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu16_4,jatayu16_4r,'+',18,6 ]
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2[:-1] + match_info2 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
       
       elif len(line2[5]) == 22 and 'AaCas12b' in line2:
            
        short_seq23=str(long_seq1)[int(line2[1])-22:int(line2[1])+1]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-21# + matchr1.start() 
         endr1=int(line2[1])+1 #+ matchr1.end()
        
         jatayu19_16 = base_change(mseqr1,line2[4],21,18,'+')
         jatayu19_16r=  str( Seq(jatayu19_16).reverse_complement())
         
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu19_16,jatayu19_16r,'+',21,18 ]
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
             
       elif len(line2[5]) == 14 :
            
        short_seq23=str(long_seq1)[int(line2[1])-14:int(line2[1])+9]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-13# + matchr1.start() 
         endr1=int(line2[1])+9 #+ matchr1.end()
        
         jatayu11_5 = base_change(mseqr1,line2[4],13,7,'+')
         jatayu11_5r=  str( Seq(jatayu11_5).reverse_complement())
         jatayu11_4 = base_change(mseqr1,line2[4],13,6,'+')
         jatayu11_4r=  str( Seq(jatayu11_4).reverse_complement())
         
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu11_5,jatayu11_5r,'+',13,7 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu11_4,jatayu11_4r,'+',13,6 ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
 
       elif len(line2[5]) == 16 :
            
        short_seq24=str(long_seq1)[int(line2[1])-16:int(line2[1])+8]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTTA[ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-15# + matchr1.start() 
         endr1=int(line2[1])+8 #+ matchr1.end()
        
         jatayu12_11 = base_change(mseqr1,line2[4],15,14,'+')
         jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[15]=line2[4]   
         jatayu12="".join(seq_list)
           
         jatayu12r=  str( Seq(jatayu12).reverse_complement())
         #print (jatayu12,jatayu12r)   
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu12_11,jatayu12_11r,'+',15,14 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu12,jatayu12r,'+',15,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['Cas14a']
         final2 =line2 + match_info2 +oligo+oligo1+['Cas14a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
        
       elif len(line2[5]) == 15 :
            
        short_seq24=str(long_seq1)[int(line2[1])-15:int(line2[1])+9]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTTA[ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-14# + matchr1.start() 
         endr1=int(line2[1])+9 #+ matchr1.end()
        
         jatayu11_12 = base_change(mseqr1,line2[4],14,15,'+')
         jatayu11_12r=  str( Seq(jatayu11_12).reverse_complement())
        
         
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu11_12,jatayu11_12r,'+',14,15]
         
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['Cas14a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
 
       elif len(line2[5]) == 12 :
            
        short_seq24=str(long_seq1)[int(line2[1])-12:int(line2[1])+8]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'([ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-11# + matchr1.start() 
         endr1=int(line2[1])+8 #+ matchr1.end()
        
         jatayu12_11 = base_change(mseqr1,line2[4],11,10,'+')
         jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[11]=line2[4]   
         jatayu12="".join(seq_list)
           
         jatayu12r=  str( Seq(jatayu12).reverse_complement())
         #print (jatayu12,jatayu12r)   
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu12_11,jatayu12_11r,'+',11,10 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu12,jatayu12r,'+',11,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['Cas14a/ssDNA']
         final2 =line2 + match_info2 +oligo+oligo1+['Cas14a/ssDNA']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
        
       elif len(line2[5]) == 11 :
            
        short_seq24=str(long_seq1)[int(line2[1])-11:int(line2[1])+9]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'([ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-10# + matchr1.start() 
         endr1=int(line2[1])+9 #+ matchr1.end()
        
         jatayu11_12 = base_change(mseqr1,line2[4],10,11,'+')
         jatayu11_12r=  str( Seq(jatayu11_12).reverse_complement())
        
         
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu11_12,jatayu11_12r,'+',10,11]
         
        
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['Cas14a/ssDNA']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
 
       elif len(line2[5]) == 3 :
            
        short_seq28=str(long_seq1)[int(line2[1])-3:int(line2[1])+25]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'([ACGT]{28})',str(short_seq28),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-2# + matchr1.start() 
         endr1=int(line2[1])+25 #+ matchr1.end()
         ##print (line2[2],line2[3],line2[4])
         jatayu3_4 = base_change(mseqr1,line2[4],2,3,'+')
         
         jatayu3_4r=  str( Seq(jatayu3_4).reverse_complement())
         jatayu3_5 = base_change(mseqr1,line2[4],2,4,'+')
         jatayu3_5r=  str( Seq(jatayu3_5).reverse_complement())
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu3_4,jatayu3_4r,'+',2,3 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu3_5,jatayu3_5r,'+',2,4 ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['LwCas13a']
         final2 =line2 + match_info2 +oligo+oligo1+['LwCas13a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)   
 
       #alphabet1 = long_seq1.alphabet
       elif len(line2[5]) == 5 : ###2 and 6
        short_seq23=str(long_seq1)[int(line2[1])-19:int(line2[1])+4]
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-18# + matchr1.start() 
         endr1=int(line2[1])+4 #+ matchr1.end()
        
         jatayu2_6 = base_change(mseqr1,line2[4],18,14,'+')
         jatayu2_6r=  str( Seq(jatayu2_6).reverse_complement())
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         #print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu2_6,jatayu2_6r,'+',18,14 ]
         #match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu2_16,jatayu2_16r,'+',18,4 ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         #print (oligo)
         Cas=[cas_system(mseqr1)]   
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         #final2 =line2 + match_info2 +oligo+oligo1
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-18# + matchr1.start() 
         endr1=int(line2[1])+4 #+ matchr1.end()
        
         jatayu2_6 = base_change(mseqr1,line2[4],18,14,'+')
         jatayu2_6r=  str( Seq(jatayu2_6).reverse_complement())
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'+')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         #print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu2_6,jatayu2_6r,'+',18,14 ]
         #match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu2_16,jatayu2_16r,'+',18,4 ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         #print (oligo)
         Cas=[cas_system(mseqr1)]    
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         #final2 =line2 + match_info2 +oligo+oligo1
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)

       elif len(line2[5]) == 9 : ###6 and 2
        short_seq23=str(long_seq1)[int(line2[1])-15:int(line2[1])+8]
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-14# + matchr1.start() 
         endr1=int(line2[1])+8 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         jatayu6_2 = base_change(mseqr1,line2[4],14,18,'+')
         jatayu6_2r=  str (Seq(jatayu6_2).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu6_2,jatayu6_2r,'+',14,18 ]
   
         #match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1) ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [str(long_seq1) [int(strtr1) -400    :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)] 
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-14# + matchr1.start() 
         endr1=int(line2[1])+8 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         jatayu6_2 = base_change(mseqr1,line2[4],14,18,'+')
         jatayu6_2r=  str (Seq(jatayu6_2).reverse_complement())
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu6_2,jatayu6_2r,'+',14,18 ]
   
         #match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1) ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [str(long_seq1) [int(strtr1) -400    :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)] 
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)

       elif len(line2[5]) == 19:###16/19 and 16
        short_seq23=str(long_seq1)[int(line2[1])-5:int(line2[1])+18]
        #print (short_seq23,line2[0],int(line2[1])-19,int(line2[1])+4)
        #for matchr1 in regex.finditer(r'(?:(G\S{20}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True):
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-4# + matchr1.start() 
         endr1=int(line2[1])+18 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         jatayu16_19 = base_change(mseqr1,line2[4],4,1,'+')
         jatayu16_19r= str( Seq(jatayu16_19).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[4]=line2[4]   
         jatayu16="".join(seq_list)
         jatayu16r=  str (Seq(jatayu16).reverse_complement() )
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu16_19,jatayu16_19r,'+',4,1 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu16,jatayu16r,'+',4,'None' ]
   

   
         #match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1) ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         final2 =line2 + match_info2 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    ) 
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-4# + matchr1.start() 
         endr1=int(line2[1])+18 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         jatayu16_19 = base_change(mseqr1,line2[4],4,1,'+')
         jatayu16_19r= str( Seq(jatayu16_19).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[4]=line2[4]   
         jatayu16="".join(seq_list)
         jatayu16r=  str (Seq(jatayu16).reverse_complement() )
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu16_19,jatayu16_19r,'+',4,1 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu16,jatayu16r,'+',4,'None' ]
   

   
         #match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1) ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]   
   
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         final2 =line2 + match_info2 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    ) 

       elif len(line2[5]) == 20: ###17
        short_seq23=str(long_seq1)[int(line2[1])-4:int(line2[1])+19]
        #print (short_seq23,line2[0],int(line2[1])-19,int(line2[1])+4)

        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-3# + matchr1.start() 
         endr1=int(line2[1])+19 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         
         seq_list=list(mseqr1)
         seq_list[3]=line2[4]   
         jatayu17="".join(seq_list)
         jatayu17r=  str (Seq(jatayu17).reverse_complement() )
         
         
   
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu17,jatayu17r,'+',3,'None' ]
         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]   
         final1 =line2 + match_info1 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
        
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-3# + matchr1.start() 
         endr1=int(line2[1])+19 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         
         seq_list=list(mseqr1)
         seq_list[3]=line2[4]   
         jatayu17="".join(seq_list)
         jatayu17r=  str (Seq(jatayu17).reverse_complement() )
         
         
   
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu17,jatayu17r,'+',3,'None' ]
         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]   

         final1 =line2 + match_info1 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )

       elif len(line2[5]) == 21: ###18
        short_seq23=str(long_seq1)[int(line2[1])-3:int(line2[1])+20]
        #print (short_seq23,line2[0],int(line2[1])-19,int(line2[1])+4)
        #for matchr1 in regex.finditer(r'(?:(G\S{20}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True):
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-2# + matchr1.start() 
         endr1=int(line2[1])+20 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         
         seq_list=list(mseqr1)
         seq_list[2]=line2[4]   
         jatayu18="".join(seq_list)
         jatayu18r=  str (Seq(jatayu18).reverse_complement() )
         
         
   
         
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu18,jatayu18r,'+',2,'None' ]

         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo+oligo1+Cas

         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-2# + matchr1.start() 
         endr1=int(line2[1])+20 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         
         seq_list=list(mseqr1)
         seq_list[2]=line2[4]   
         jatayu18="".join(seq_list)
         jatayu18r=  str (Seq(jatayu18).reverse_complement() )
         
         
   
         
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu18,jatayu18r,'+',2,'None' ]

         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]

         final1 =line2 + match_info1 +oligo+oligo1+Cas

         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )

       elif len(line2[5]) == 22: ###19 and 16 
        short_seq23=str(long_seq1)[int(line2[1])-2:int(line2[1])+21]
        #print (short_seq23,line2[0],int(line2[1])-19,int(line2[1])+4)
        #for matchr1 in regex.finditer(r'(?:(G\S{20}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True):
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-1# + matchr1.start() 
         endr1=int(line2[1])+21 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         jatayu19_16 = base_change(mseqr1,line2[4],1,4,'+')
         jatayu19_16r=  str(Seq(jatayu19_16).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[1]=line2[4]   
         jatayu19="".join(seq_list)
         jatayu19r=  str (Seq(jatayu19).reverse_complement() )
         
         
   
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu19_16,jatayu19_16r,'+',1,4 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu19,jatayu19r,'+',1,'None' ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]

         final1 =line2 + match_info1 +oligo+oligo1+Cas
         final2 =line2 + match_info2 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    )   
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[1])-1# + matchr1.start() 
         endr1=int(line2[1])+21 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         jatayu19_16 = base_change(mseqr1,line2[4],1,4,'+')
         jatayu19_16r=  str(Seq(jatayu19_16).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[1]=line2[4]   
         jatayu19="".join(seq_list)
         jatayu19r=  str (Seq(jatayu19).reverse_complement() )
         
         
   
         match_info1 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu19_16,jatayu19_16r,'+',1,4 ]
         match_info2 =  [str(mseqr1) ,str(strtr1) , str(endr1),jatayu19,jatayu19r,'+',1,'None' ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [str(long_seq1) [int(strtr1) -400  :int(strtr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         Cas=[cas_system(mseqr1)]

         final1 =line2 + match_info1 +oligo+oligo1+Cas
         final2 =line2 + match_info2 +oligo+oligo1+Cas
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    )   

     
   return vcf_sgrna_ls


# In[22]:


def ontarget1(genome,vcf_ls):

   lst=[]
   vcf_sgrna_ls=[]
   records = SeqIO.to_dict(SeqIO.parse(open(genome), 'fasta'))
   for line1 in vcf_ls:
    ##print (line1[5])
    long_seq_record = records[line1[0]]
    num=len(long_seq_record.seq) - (line1[1]-1)
    line1.insert(5, num)
    #print (line1,len(long_seq_record.seq))
    long_seq = long_seq_record.seq
    #alphabet = long_seq.alphabet
    short_seq1=str(long_seq)[int(line1[5])-1]
    #short_seq1_tmp=str(long_seq)[int(line1[5])]
    ##print(line1[3])
    #line1[3]=reverse_complement(line1[3])
    if short_seq1 == reverse_complement(line1[3]):
     #print (reverse_complement(line1[3]), "match", short_seq1)
     short_seq2_6 = str(long_seq)[int(line1[5])-1:int(line1[5])+4] ####mutation at 2nd position and random at 6
     short_seq6_2 = str(long_seq)[int(line1[5])-1:int(line1[5])+8] ####mutation at 6th position and random at 2  
     short_seq16_19 = str(long_seq)[int(line1[5])-1:int(line1[5])+18] #####mutation at 16th position and random at 19
     short_seq19_16 = str(long_seq)[int(line1[5])-1:int(line1[5])+21] #####mutation at 19th position and random at 16 
     short_seq18 = str(long_seq)[int(line1[5])-1:int(line1[5])+20] #####mutation at 18th position
     short_seq17 = str(long_seq)[int(line1[5])-1:int(line1[5])+19] #####mutation at 17th position
     short_seq3_4 = str(long_seq)[int(line1[5])-3:int(line1[5])] ####mutation at 3rd position and random at 5 (3 and 4)
     short_seq11_12 = str(long_seq)[int(line1[5])-11:int(line1[5])] ####mutation at 11th position and random at 12 Cas14a/ssDNA
     short_seq12_11 = str(long_seq)[int(line1[5])-12:int(line1[5])] ####mutation at 12th position and random at 11 Cas14a/ssDNA   
     short_seq11_12_pam = str(long_seq)[int(line1[5])-15:int(line1[5])] ####mutation at 11th position and random at 12 
     short_seq12_11_pam = str(long_seq)[int(line1[5])-16:int(line1[5])] ####mutation at 12th position and random at 11    
     short_seq1_5 = str(long_seq)[int(line1[5])-4:int(line1[5])] ####mutation at 1nd position and random at 5
     short_seq5_1 = str(long_seq)[int(line1[5])-8:int(line1[5])] ####mutation at 5th position and random at 1  
     short_seq4_1 = str(long_seq)[int(line1[5])-7:int(line1[5])] ####mutation at 4th position and random at 1
     short_seq16_19_cas12b = str(long_seq)[int(line1[5])-19:int(line1[5])] ####mutation at 4th position and random at 1
     short_seq19_16_cas12b = str(long_seq)[int(line1[5])-22:int(line1[5])] ####mutation at 4th position and random at 1   
     short_seq11_5 = str(long_seq)[int(line1[5])-14:int(line1[5])] ####mutation at 4th position and random at 1   
     short_seq1 = str(long_seq)[int(line1[5])-5:int(line1[5])]  
     short_seq2 = str(long_seq)[int(line1[5])-6:int(line1[5])]
     short_seq3 = str(long_seq)[int(line1[5])-7:int(line1[5])]
     short_seq4 = str(long_seq)[int(line1[5])-8:int(line1[5])]
     short_seq5 = str(long_seq)[int(line1[5])-9:int(line1[5])]
     short_seq6 = str(long_seq)[int(line1[5])-10:int(line1[5])]
     short_seq7 = str(long_seq)[int(line1[5])-11:int(line1[5])]
     for matchr in regex.finditer(r'(TTT[ACGT]{2})',str(short_seq1),regex.BESTMATCH,overlapped=True):
      #print (matchr)
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TTT[ACGT]{3})',str(short_seq2),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(TTT[ACGT]{4})',str(short_seq3),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(TTT[ACGT]{5})',str(short_seq4),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)       
     for matchr in regex.finditer(r'(TTT[ACGT]{6})',str(short_seq5),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)     
     for matchr in regex.finditer(r'(TTT[ACGT]{7})',str(short_seq6),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)     
     for matchr in regex.finditer(r'(TTT[ACGT]{8})',str(short_seq7),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['LbCas12a']
      ##print (final) 
      lst.append(final)     
   
   
         
     for matchr in regex.finditer(r'(TT[ACGT]{2})',str(short_seq1_5),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{6})',str(short_seq5_1),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{5})',str(short_seq4_1),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{17})',str(short_seq16_19_cas12b),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['AaCas12b']
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{20})',str(short_seq19_16_cas12b),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info +['AaCas12b']
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TT[ACGT]{12})',str(short_seq11_5),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
        
       
     for matchr in regex.finditer(r'(TTTA[ACGT]{12})',str(short_seq12_11_pam),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(TTTA[ACGT]{11})',str(short_seq11_12_pam),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)   
   
     for matchr in regex.finditer(r'([ACGT]{12})',str(short_seq12_11),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'([ACGT]{11})',str(short_seq11_12),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)   
   
     for matchr in regex.finditer(r'([ACGT]{3})',str(short_seq3_4),regex.BESTMATCH,overlapped=True):
      ###print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      ##print (final) 
      lst.append(final)

     for matchr in regex.finditer(r'(?:([ACGT]{19}[G][AG]))',str(short_seq18),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{18}[G][AG]))',str(short_seq17),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(?:([ACGT]{3}[G][AG]))',str(short_seq2_6),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{7}[G][AG]))',str(short_seq6_2),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)
   
     for matchr in regex.finditer(r'(?:([ACGT]{17}[G][AG]))',str(short_seq16_19),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(?:([ACGT]{20}[G][AG]))',str(short_seq19_16),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)   
     ###NAG               
     for matchr in regex.finditer(r'(?:([ACGT]{19}AG))',str(short_seq18),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{18}AG))',str(short_seq17),regex.BESTMATCH,overlapped=True):
      #print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info 
      #print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(?:([ACGT]{3}AG))',str(short_seq2_6),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)
     for matchr in regex.finditer(r'(?:([ACGT]{7}AG))',str(short_seq6_2),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)
   
     for matchr in regex.finditer(r'(?:([ACGT]{17}AG))',str(short_seq16_19),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)   
     for matchr in regex.finditer(r'(?:([ACGT]{20}AG))',str(short_seq19_16),regex.BESTMATCH,overlapped=True):
      ##print matchr
      mr=matchr.groups()
      mseqr=re.sub('[,\')(\s]','',str(mr).rstrip())
      strtr=int(line1[5]) + matchr.start() 
      endr=int(line1[5])-1 + matchr.end()   
      #listToStr = '\t'.join([str(elem) for elem in line1]) 
      match_info =  [str(mseqr) ,str(strtr) , str(endr) ] 
      final = line1 + match_info
      ##print (final) 
      lst.append(final)   
 
   for line2 in lst:
       ##print (line2[6],line2[5])
     #  jk1="\t".join(jk)
      # #print (jk[0])
       #jk=jk.rstrip()
       #line2=re.split("\t",jk)
       #line2=re.split("\t",jk)
       #print (line2)
       #print (line2[6]) 
       long_seq_record1 = records[line2[0]]
       
       long_seq1 = long_seq_record1.seq
       if len(line2[6]) == 5 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[5])-5:int(line2[5])+17]
        #print (short_seq18)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         #print (mr1)
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-4# + matchr1.start() 
         endr1=int(line2[5])+17 #+ matchr1.end()
        
         seq_list=list(mseqr1)
         seq_list[4]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-4)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+17-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',4,'None' ]
        
            
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
         #print (final1)  
         vcf_sgrna_ls.append(final1)
       
       elif len(line2[6]) == 6 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[5])-6:int(line2[5])+16]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-5# + matchr1.start() 
         endr1=int(line2[5])+16 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[5]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-5)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+16-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',5,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1) 
       
       elif len(line2[6]) == 7 and 'LbCas12a' in line2 :
            
        short_seq18=str(long_seq1)[int(line2[5])-7:int(line2[5])+15]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-6# + matchr1.start() 
         endr1=int(line2[5])+15 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[6]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-6)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+15-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',6,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']

        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1) 
       elif len(line2[6]) == 8 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[5])-8:int(line2[5])+14]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-7# + matchr1.start() 
         endr1=int(line2[5])+14 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[7]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-7)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+14-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',7,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
       elif len(line2[6]) == 9 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[5])-9:int(line2[5])+13]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-8# + matchr1.start() 
         endr1=int(line2[5])+13 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[8]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-8)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+13-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',8,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']

        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
       elif len(line2[6]) == 10 and 'LbCas12a' in line2 :
            
        short_seq18=str(long_seq1)[int(line2[5])-10:int(line2[5])+12]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-9# + matchr1.start() 
         endr1=int(line2[5])+12 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[9]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-9)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+12-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',9,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   
       elif len(line2[6]) == 11 and 'LbCas12a' in line2:
            
        short_seq18=str(long_seq1)[int(line2[5])-11:int(line2[5])+11]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTT[ACGT]{19})',str(short_seq18),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-10# + matchr1.start() 
         endr1=int(line2[5])+11 #+ matchr1.end()
        
         #jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         #jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[10]=reverse_complement(line2[4])   
         jatayu1="".join(seq_list)
           
         jatayu1r=  str( Seq(jatayu1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-10)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+11-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu1,jatayu1r,'-',10,'None' ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str(long_seq1) [int(endr1) -400    :int(endr1)+400   ]]
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['LbCas12a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)   

       elif len(line2[6]) == 4 :
            
        short_seq23=str(long_seq1)[int(line2[5])-4:int(line2[5])+19]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-3# + matchr1.start() 
         endr1=int(line2[5])+19 #+ matchr1.end()
         ##print (line2[2],line2[3],line2[4])
         jatayu1_5 = base_change(mseqr1,line2[4],3,7,'-')
         #print (line2[0],line2[5],line2[2],line2[3],line2[4],jatayu1_5)   
         jatayu1_5r=  str( Seq(jatayu1_5).reverse_complement())
         jatayu1_4 = base_change(mseqr1,line2[4],3,6,'-')
         jatayu1_4r=  str( Seq(jatayu1_4).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-3)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+19-1) #+ matchr1.end()
   
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu1_5,jatayu1_5r,'-',3,7 ]
         match_info2 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu1_4,jatayu1_4r,'-',3,6 ]
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = ['NA']#.reverse_complement()
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)   
       elif len(line2[6]) == 8 :
            
        short_seq23=str(long_seq1)[int(line2[5])-8:int(line2[5])+15]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-7# + matchr1.start() 
         endr1=int(line2[5])+15 #+ matchr1.end()
        
         jatayu5_1 = base_change(mseqr1,line2[4],7,3,'-')
         jatayu5_1r=  str( Seq(jatayu5_1).reverse_complement())
         jatayu5_8 = base_change(mseqr1,line2[4],7,10,'-')
         jatayu5_8r=  str( Seq(jatayu5_1).reverse_complement())
         jatayu5_11 = base_change(mseqr1,line2[4],7,13,'-')
         jatayu5_11r=  str( Seq(jatayu5_1).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-7)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+15-1) #+ matchr1.end()
         
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu5_1,jatayu5_1r,'-',7,3 ]
         match_info2 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu5_8,jatayu5_8r,'-',7,10 ]
         match_info3 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu5_11,jatayu5_11r,'-',7,13 ]

         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] 
         oligo1 = ['NA']
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
         final3 =line2 + match_info3 +oligo+oligo1 +['AaCas12b']  
            
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
         vcf_sgrna_ls.append(final3)
       
       elif len(line2[6]) == 7 :
            
        short_seq23=str(long_seq1)[int(line2[5])-7:int(line2[5])+16]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-6# + matchr1.start() 
         endr1=int(line2[5])+16 #+ matchr1.end()
        
         jatayu4_1 = base_change(mseqr1,line2[4],6,3,'-')
         jatayu4_1r=  str( Seq(jatayu4_1).reverse_complement())
         jatayu4_11 = base_change(mseqr1,line2[4],6,13,'-')
         jatayu4_11r=  str( Seq(jatayu4_11).reverse_complement())
         jatayu4_16 = base_change(mseqr1,line2[4],6,18,'-')
         jatayu4_16r=  str( Seq(jatayu4_16).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-6)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+16-1) #+ matchr1.end()
   
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu4_1,jatayu4_1r,'-',6,3 ]
         match_info2 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu4_11,jatayu4_11r,'-',6,13 ]
         match_info3 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu4_16,jatayu4_16r,'-',6,18 ]   
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] 
         oligo1 = ['NA']         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
         final3 =line2 + match_info3 +oligo+oligo1+['AaCas12b']   
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
         vcf_sgrna_ls.append(final3)
        
       elif len(line2[6]) == 19 and 'AaCas12b' in line2:
            
        short_seq23=str(long_seq1)[int(line2[5])-19:int(line2[5])+4]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-18# + matchr1.start() 
         endr1=int(line2[5])+4 #+ matchr1.end()
        
         jatayu16_19 = base_change(mseqr1,line2[4],18,21,'-')
         jatayu16_19r=  str( Seq(jatayu16_19).reverse_complement())
         jatayu16_4 = base_change(mseqr1,line2[4],18,6,'-')
         jatayu16_4r=  str( Seq(jatayu16_4).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-18)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+4-1) #+ matchr1.end()
   
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu16_19,jatayu16_19r,'-',18,21 ]
         match_info2 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu16_4,jatayu16_4r,'-',18,6 ]
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] 
         oligo1 = ['NA']         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2[:-1] + match_info2 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
       
       elif len(line2[6]) == 22 and 'AaCas12b' in line2 :
            
        short_seq23=str(long_seq1)[int(line2[5])-22:int(line2[5])+1]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-21# + matchr1.start() 
         endr1=int(line2[5])+1 #+ matchr1.end()
        
         jatayu19_16 = base_change(mseqr1,line2[4],21,18,'-')
         jatayu19_16r=  str( Seq(jatayu19_16).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-21)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+1-1) #+ matchr1.end()
   
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu19_16,jatayu19_16r,'-',21,18 ]
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] 
         oligo1 = ['NA']         ###print (oligo)
         final1 =line2[:-1] + match_info1 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
             
       elif len(line2[6]) == 14 :
            
        short_seq23=str(long_seq1)[int(line2[5])-14:int(line2[5])+9]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TT[ACGT]{21})',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-13# + matchr1.start() 
         endr1=int(line2[5])+9 #+ matchr1.end()
        
         jatayu11_5 = base_change(mseqr1,line2[4],13,7,'-')
         jatayu11_5r=  str( Seq(jatayu11_5).reverse_complement())
         jatayu11_4 = base_change(mseqr1,line2[4],13,6,'-')
         jatayu11_4r=  str( Seq(jatayu11_4).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-13)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+9-1) #+ matchr1.end()
   
              #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r=  str( Seq(str(jatayu2_16)).reverse_complement())
         ###print (type(jatayu2_16),type(jatayu2_16r))
         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu11_5,jatayu11_5r,'-',13,7 ]
         match_info2 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu11_4,jatayu11_4r,'-',13,6 ]
        
            
         #oligo = [str(long_seq1) [int(line2[5]) -252   :int(line2[5]) +400   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] 
         oligo1 = ['NA']         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['AaCas12b']
         final2 =line2 + match_info2 +oligo+oligo1+['AaCas12b']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)

       #alphabet1 = long_seq1.alphabet
       elif len(line2[6]) == 16 :
            
        short_seq24=str(long_seq1)[int(line2[5])-16:int(line2[5])+8]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTTA[ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-15# + matchr1.start() 
         endr1=int(line2[5])+8 #+ matchr1.end()
         
         jatayu12_11 = base_change(mseqr1,line2[4],15,14,'-')
         jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[15]=reverse_complement(line2[4])   
         jatayu12="".join(seq_list)
           
         jatayu12r=  str( Seq(jatayu12).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-15)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+8-1) #+ matchr1.end()   
         #print (jatayu12,jatayu12r)   
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu12_11,jatayu12_11r,'-',15,14 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu12,jatayu12r,'-',15,'None' ]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] #.reverse_complement()
         oligo1 = ['NA']#.reverse_complement()
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         ###print (oligo)
         final1 =line2 + match_info1 +oligo+oligo1+['Cas14a']
         final2 =line2 + match_info2 +oligo+oligo1+['Cas14a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
        
       elif len(line2[6]) == 15 :
            
        short_seq24=str(long_seq1)[int(line2[5])-15:int(line2[5])+9]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'(TTTA[ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-14# + matchr1.start() 
         endr1=int(line2[5])+9 #+ matchr1.end()
        
         jatayu11_12 = base_change(mseqr1,line2[4],14,15,'-')
         jatayu11_12r=  str( Seq(jatayu11_12).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-14)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+9-1) #+ matchr1.end()   

         
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu11_12,jatayu11_12r,'-',14,15]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())] #.reverse_complement()
         oligo1 = ['NA']#.reverse_complement()
         final1 =line2 + match_info1 +oligo+oligo1+['Cas14a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         
 
       elif len(line2[6]) == 12 :
            
        short_seq24=str(long_seq1)[int(line2[5])-12:int(line2[5])+8]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'([ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-11# + matchr1.start() 
         endr1=int(line2[5])+8 #+ matchr1.end()
         
         jatayu12_11 = base_change(mseqr1,line2[4],11,10,'-')
         jatayu12_11r=  str( Seq(jatayu12_11).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[11]=reverse_complement(line2[4])   
         jatayu12="".join(seq_list)
           
         jatayu12r=  str( Seq(jatayu12).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-11)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+8-1) #+ matchr1.end()   
         #print (jatayu12,jatayu12r)   
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu12_11,jatayu12_11r,'-',11,10 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu12,jatayu12r,'-',11,'None' ]
         oligo = [ 'NA'] #.reverse_complement()
         oligo1 = [str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
            
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         ###print (oligo)
         final1 =line2 + match_info1 +oligo1+oligo+['Cas14a/ssDNA']
         final2 =line2 + match_info2 +oligo1+oligo+['Cas14a/ssDNA']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)
        
       elif len(line2[6]) == 11 :
            
        short_seq24=str(long_seq1)[int(line2[5])-11:int(line2[5])+9]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'([ACGT]{20})',str(short_seq24),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-10# + matchr1.start() 
         endr1=int(line2[5])+9 #+ matchr1.end()
        
         jatayu11_12 = base_change(mseqr1,line2[4],10,11,'-')
         jatayu11_12r=  str( Seq(jatayu11_12).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-10)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+9-1) #+ matchr1.end()   

         
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu11_12,jatayu11_12r,'-',10,11]
         oligo = [ 'NA'] #.reverse_complement()
         oligo1 = [str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         final1 =line2 + match_info1 +oligo1+oligo+['Cas14a/ssDNA']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
 
       elif len(line2[6]) == 3 :
            
        short_seq28=str(long_seq1)[int(line2[5])-3:int(line2[5])+25]
        ##print (short_seq23)
        for matchr1 in regex.finditer(r'([ACGT]{28})',str(short_seq28),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-2# + matchr1.start() 
         endr1=int(line2[5])+25 #+ matchr1.end()
         ##print (line2[2],line2[3],line2[4])
         jatayu3_4 = base_change(mseqr1,line2[4],2,3,'-')
         
         jatayu3_4r=  str( Seq(jatayu3_4).reverse_complement())
         jatayu3_5 = base_change(mseqr1,line2[4],2,4,'-')
         jatayu3_5r=  str( Seq(jatayu3_5).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-2)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+25-1) #+ matchr1.end()   
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu3_4,jatayu3_4r,'-',2,3 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu3_5,jatayu3_5r,'-',2,4 ]
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +400   ]]
         #oligo = [ str(long_seq1) [int(strtr1) -400    :int(strtr1)+1000   ]]
         #oligo1 = [str(long_seq1) [int(strtr1) -1000    :int(strtr1)+400   ]]
         oligo = ['NA'] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -400    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         ###print (oligo)
         final1 =line2 + match_info1 +oligo1+oligo+['LwCas13a']
         final2 =line2 + match_info2 +oligo1+oligo+['LwCas13a']
        #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
         vcf_sgrna_ls.append(final2)   
 
       elif len(line2[6]) == 5 :
        short_seq23=str(long_seq1)[int(line2[5])-19:int(line2[5])+4]
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1= (int(line2[5])-18)# + matchr1.start() 
         endr1= (int(line2[5])+4) #+ matchr1.end()
         jatayu2_6 = base_change(mseqr1,line2[4],18,14,'-')
         jatayu2_6r=  str(Seq(str(jatayu2_6)).reverse_complement())
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r= str(Seq(str(jatayu2_16)).reverse_complement() )
         strtr11=len(long_seq1) - (int(line2[5])-18)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+4-1) #+ matchr1.end()

         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu2_6,jatayu2_6r,'-',18,14 ]
         
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]

         final1 =line2 + match_info1 +oligo1+oligo+Cas
         
         vcf_sgrna_ls.append(final1)
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1= (int(line2[5])-18)# + matchr1.start() 
         endr1= (int(line2[5])+4) #+ matchr1.end()
         jatayu2_6 = base_change(mseqr1,line2[4],18,14,'-')
         jatayu2_6r=  str(Seq(str(jatayu2_6)).reverse_complement())
         #jatayu2_16 = base_change(mseqr1,line2[4],18,4,'-')
         #jatayu2_16r= str(Seq(str(jatayu2_16)).reverse_complement() )
         strtr11=len(long_seq1) - (int(line2[5])-18)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+4-1) #+ matchr1.end()

         match_info1 =  [str(mseqr1) ,str(endr11),str(strtr11),jatayu2_6,jatayu2_6r,'-',18,14 ]
         
         #oligo = [str(long_seq1) [int(line2[1]) -252   :int(line2[1]) +250   ]]
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]

         final1 =line2 + match_info1 +oligo1+oligo+Cas
         
         vcf_sgrna_ls.append(final1)
        
   
       elif len(line2[6]) == 9 :
        short_seq23=str(long_seq1)[int(line2[5])-15:int(line2[5])+8]
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-14# + matchr1.start() 
         endr1=int(line2[5])+8 #+ matchr1.end()
         jatayu6_2 = base_change(mseqr1,line2[4],14,18,'-')
         jatayu6_2r=str(Seq(str(jatayu6_2)).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-14)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+8-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu6_2,jatayu6_2r,'-',14,18 ]

         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]

         final1 =line2 + match_info1 +oligo1+oligo+Cas
                  #oligo1r=oligo1.reverse_complement()
         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1)
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-14# + matchr1.start() 
         endr1=int(line2[5])+8 #+ matchr1.end()
         jatayu6_2 = base_change(mseqr1,line2[4],14,18,'-')
         jatayu6_2r=str(Seq(str(jatayu6_2)).reverse_complement())
         strtr11=len(long_seq1) - (int(line2[5])-14)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+8-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu6_2,jatayu6_2r,'-',14,18 ]

         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas
         vcf_sgrna_ls.append(final1)

       elif len(line2[6]) == 19:
        short_seq23=str(long_seq1)[int(line2[5])-5:int(line2[5])+18]
        ##print (short_seq23,line2[0],int(line2[5])-19,int(line2[5])+4)
        #for matchr1 in regex.finditer(r'(?:(G\S{20}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True):
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-4# + matchr1.start() 
         endr1=int(line2[5])+18 #+ matchr1.end()
         
         jatayu16_19 = base_change(mseqr1,line2[4],4,1,'-')
         jatayu16_19r=  str(Seq(str(jatayu16_19)).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[4]=reverse_complement(line2[4])   
         
         jatayu16="".join(seq_list)
         jatayu16r=  str (Seq(jatayu16).reverse_complement() )

         
         strtr11=len(long_seq1) - (int(line2[5])-4)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+18-1) #+ matchr1.end()
   
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu16_19,jatayu16_19r,'-',4,1 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu16,jatayu16r,'-',4,'None' ]  
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas
         final2 =line2 + match_info2 +oligo1+oligo+Cas
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    )   
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-4# + matchr1.start() 
         endr1=int(line2[5])+18 #+ matchr1.end()
         
         jatayu16_19 = base_change(mseqr1,line2[4],4,1,'-')
         jatayu16_19r=  str(Seq(str(jatayu16_19)).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[4]=reverse_complement(line2[4])   
         
         jatayu16="".join(seq_list)
         jatayu16r=  str (Seq(jatayu16).reverse_complement() )

         
         strtr11=len(long_seq1) - (int(line2[5])-4)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+18-1) #+ matchr1.end()
   
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu16_19,jatayu16_19r,'-',4,1 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu16,jatayu16r,'-',4,'None' ]  
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas
   
         final2 =line2 + match_info2 +oligo1+oligo+Cas
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    )   

       elif len(line2[6]) == 20: ###17
        short_seq23=str(long_seq1)[int(line2[5])-4:int(line2[5])+19]
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-3# + matchr1.start() 
         endr1=int(line2[5])+19 #+ matchr1.end()
        
         seq_list=list(mseqr1)
         seq_list[3]=reverse_complement(line2[4])   
         
         jatayu17="".join(seq_list)
         jatayu17r=  str (Seq(jatayu17).reverse_complement() )
        
         strtr11=len(long_seq1) - (int(line2[5])-3)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+19-1) #+ matchr1.end() 
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu17,jatayu17r,'-',3,'None' ]
         oligo = [str(long_seq1) [int(endr1) -400  :int(endr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(endr1) -1000    :int(endr1)+400   ]]
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas   
         

         vcf_sgrna_ls.append(final1    )
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-3# + matchr1.start() 
         endr1=int(line2[5])+19 #+ matchr1.end()
        
         seq_list=list(mseqr1)
         seq_list[3]=reverse_complement(line2[4])   
         
         jatayu17="".join(seq_list)
         jatayu17r=  str (Seq(jatayu17).reverse_complement() )
        
         strtr11=len(long_seq1) - (int(line2[5])-3)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+19-1) #+ matchr1.end() 
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu17,jatayu17r,'-',3,'None' ]
         oligo = [str(long_seq1) [int(endr1) -400  :int(endr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(endr1) -1000    :int(endr1)+400   ]]
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas   
         

         vcf_sgrna_ls.append(final1    )

       elif len(line2[6]) == 21: ###18
        short_seq23=str(long_seq1)[int(line2[5])-3:int(line2[5])+20]
        #print (short_seq23,line2[0],int(line2[1])-19,int(line2[1])+4)
        #for matchr1 in regex.finditer(r'(?:(G\S{20}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True):
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-2# + matchr1.start() 
         endr1=int(line2[5])+20 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         seq_list=list(mseqr1)

         seq_list[2]=reverse_complement(line2[4])   
         
         jatayu18="".join(seq_list)
         jatayu18r=  str (Seq(jatayu18).reverse_complement() )
         strtr11=len(long_seq1) - (int(line2[5])-2)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+20-1) #+ matchr1.end()
         
   
         
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu18,jatayu18r,'-',2,'None' ]

         oligo = [str(long_seq1) [int(endr1) -400  :int(endr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(endr1) -1000    :int(endr1)+400   ]]
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas
         

         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-2# + matchr1.start() 
         endr1=int(line2[5])+20 #+ matchr1.end()
         #mseqr1=list(mseqr1)
         seq_list=list(mseqr1)

         seq_list[2]=reverse_complement(line2[4])   
         
         jatayu18="".join(seq_list)
         jatayu18r=  str (Seq(jatayu18).reverse_complement() )
         strtr11=len(long_seq1) - (int(line2[5])-2)+1# + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+20-1) #+ matchr1.end()
         
   
         
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu18,jatayu18r,'-',2,'None' ]

         oligo = [str(long_seq1) [int(endr1) -400  :int(endr1)+1000   ]]
         oligo1 = [str(long_seq1) [int(endr1) -1000    :int(endr1)+400   ]]
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas
         

         #line22="\t".join(line2)
         vcf_sgrna_ls.append(final1    )

         
       elif len(line2[6]) == 22:
        short_seq23=str(long_seq1)[int(line2[5])-2:int(line2[5])+21]
        ##print (short_seq23,line2[0],int(line2[5])-19,int(line2[5])+4)
        #for matchr1 in regex.finditer(r'(?:(G\S{20}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True):
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}[G][AG]))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-1# + matchr1.start() 
         endr1=int(line2[5])+21 #+ matchr1.end()
         jatayu19_16 = base_change(mseqr1,line2[4],1,4,'-')
         jatayu19_16r=  str(Seq(str(jatayu19_16)).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[1]=reverse_complement(line2[4])   
         
         jatayu19="".join(seq_list)
         jatayu19r=  str (Seq(jatayu19).reverse_complement() )

         
         
         strtr11=len(long_seq1) - (int(line2[5])-1) +1 # + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+21-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu19_16,jatayu19_16r,'-',1,4 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu19,jatayu19r,'-',1,'None']
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas   
         
         final2 =line2 + match_info2 +oligo1+oligo+Cas
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    )   
        for matchr1 in regex.finditer(r'(?:([ACGT]{21}AG))',str(short_seq23),regex.BESTMATCH,overlapped=True): 
         mr1=matchr1.groups()
         mseqr1=re.sub('[,\')(\s]','',str(mr1).rstrip())
         strtr1=int(line2[5])-1# + matchr1.start() 
         endr1=int(line2[5])+21 #+ matchr1.end()
         jatayu19_16 = base_change(mseqr1,line2[4],1,4,'-')
         jatayu19_16r=  str(Seq(str(jatayu19_16)).reverse_complement())
         seq_list=list(mseqr1)
         seq_list[1]=reverse_complement(line2[4])   
         
         jatayu19="".join(seq_list)
         jatayu19r=  str (Seq(jatayu19).reverse_complement() )

         
         
         strtr11=len(long_seq1) - (int(line2[5])-1) +1 # + matchr1.start() 
         endr11=len(long_seq1) - (int(line2[5])+21-1) #+ matchr1.end()
         match_info1 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu19_16,jatayu19_16r,'-',1,4 ]
         match_info2 =  [str(mseqr1) ,str(endr11) , str(strtr11),jatayu19,jatayu19r,'-',1,'None']
         oligo = [ str((long_seq1) [int(endr1) -400    :int(endr1)+1000   ].reverse_complement())] #.reverse_complement()
         #oligor=oligo.reverse_complement()   
         oligo1 = [str((long_seq1) [int(endr1) -1000    :int(endr1)+400   ].reverse_complement())]#.reverse_complement()
         Cas=[cas_system(mseqr1)]
         final1 =line2 + match_info1 +oligo1+oligo+Cas   
         
         final2 =line2 + match_info2 +oligo1+oligo+Cas   
         
         vcf_sgrna_ls.append(final1    )
         vcf_sgrna_ls.append(final2    )   
     
   return vcf_sgrna_ls
    
    


# In[ ]:


def cas_system (row):
   #print (row['gene_biotype'])
   if row.endswith('GG') == True  :
      #print ('protein_coding1')  
      return 'Fn/enFnCas9'
   else: #row['2bpam'].endswith('AG') ==  :
      return 'enFnCas9'
    
def sg_df(variant_summ_fil):
    chrom_ref=pd.read_csv("chromosome_refseq.txt",sep="\t")
    variant_summ_coord= variant_summ_fil[["Chr", "Position" ,"rsID", "Ref", "Alt"]]
    variant_summ_snp_coord1= variant_summ_coord.copy()
    variant_summ_snp_ls = variant_summ_snp_coord1.values.tolist()
    variant_summ_snp_ls1=variant_summ_snp_ls
    vcf_sgrna_feluda1=ontarget("/home/asgar/asgar_work/asgar_crispr/wdout_aug_all_cells_epigenetics_features/nat_comm_data/hg38/GRCh38.primary_assembly.genome.fa",variant_summ_snp_ls)
    vcf_sgrna_feluda2=ontarget1("/home/asgar/asgar_work/asgar_crispr/wdout_aug_all_cells_epigenetics_features/nat_comm_data/hg38/GRCh38.primary_assembly.genome_reverse.fa",variant_summ_snp_ls1)
    vcf_sgrna_feluda_df1 = pd.DataFrame(vcf_sgrna_feluda1)
    vcf_sgrna_feluda_df2 = pd.DataFrame(vcf_sgrna_feluda2)
    vcf_sgrna_feluda_df2.drop(vcf_sgrna_feluda_df2.columns[[5,7,8]], axis = 1, inplace = True)
    vcf_sgrna_feluda_df1.drop(vcf_sgrna_feluda_df1.columns[[6, 7]], axis = 1, inplace = True)
    vcf_sgrna_feluda_df1.columns=  ["Chr", "Position" ,"rsID", "Ref", "Alt","2bpam","presgrna","prestart","preend","jatayu","jatayur" ,"Strand","mut_loc_cr","ran_loc_cr"    ,"oligo","oligo1","Cas_system"]
    vcf_sgrna_feluda_df2.columns=  ["Chr", "Position" ,"rsID", "Ref", "Alt","2bpam","presgrna","prestart","preend","jatayu","jatayur" ,"Strand","mut_loc_cr","ran_loc_cr"    ,"oligo","oligo1","Cas_system"]
    vcf_sgrna_feluda_df=pd.concat([vcf_sgrna_feluda_df1,vcf_sgrna_feluda_df2])
    vcf_sgrna_feluda_df_full=pd.merge(vcf_sgrna_feluda_df,variant_summ_fil,on=["Chr", "Position" ,"rsID", "Ref", "Alt"])
    vcf_sgrna_feluda_df_full=pd.merge(vcf_sgrna_feluda_df_full,chrom_ref,on=["Chr"])
    vcf_sgrna_feluda_df_full['Strand']=vcf_sgrna_feluda_df_full['Strand'].str.replace(r'+', 'fw')
    vcf_sgrna_feluda_df_full['Strand']=vcf_sgrna_feluda_df_full['Strand'].str.replace(r'-', 'rev')
    return vcf_sgrna_feluda_df_full







    






    







    







    


# In[6]:


#cols_to_use=[0,4,10,5,6]


# In[47]:


#sgene=pd.read_csv("/home/asgar/flask/crispsnps/sgrna_results/cr_ed977456ee9f4b85b281a4eba71e15a5_gene",sep="\t", header=None)


# In[48]:


#sgene=sgene.iloc[:,[0,4,10,5,6]]


# In[49]:


#sgene.columns=["CHROMOSOME", "POS" ,"RS", "REF", "ALT"]


# In[50]:


#sgene


# In[53]:


#time sg_df(sgene)

