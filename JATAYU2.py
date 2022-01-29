from flask import Flask, render_template,request,g,flash,redirect,url_for,jsonify,session
from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField,  PasswordField, TextAreaField, RadioField, SelectField
from wtforms.validators import InputRequired,Length,AnyOf,ValidationError,NoneOf,Regexp
from matplotlib import pyplot as plt              
from subprocess	import	call
import subprocess
import sys
import re
import numpy as np
from io import BytesIO
import base64
import os
import time
import ntpath
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
#from Bio.Alphabet import generic_dna
import pandas as pd
import sqlite3
import string
import random
import primer3
import uuid
from matplotlib.pyplot import figure
import ngr_nrg
#import pyodbc
app = Flask(__name__)
app.config['SECRET_KEY'] = ''

def connect_db():
    sql=sqlite3.connect('/home/asgar/JATAYU2/db/ngr_nrg_snp_wd_freq.db')
    sql.row_factory=sqlite3.Row
    return sql
def get_db():
    if not hasattr (g,'sqlite3'):
        g.sqlite_db=connect_db()
    return g.sqlite_db        

@app.teardown_appcontext
def close_db(error):
    if hasattr (g,'sqlite_db'):
        g.sqlite_db.close()


def Convert(string):
    li = list(string.split("|"))
    return li

"""
def base_check(form, field):
    pattern = r"^rs[\d]+$"
    #print (pattern)
    if re.findall(pattern, field.data):
               raise ValidationError('Not a Correct RSID')
"""




class jatayu(FlaskForm):

    variation = SelectField('Type of Mutation', choices=[('A', 'A'), ('C', 'C'),('G', 'G'), ('T', 'T') ]   )
    organism = SelectField('Organism', choices=[('Homo Sapiens', 'Sars-CoV-2') ]   )
def base_check(form, field):
    
    pattern = "[^ACTG]+"
    if re.findall(pattern, field.data):
               #flash('Provide Valid DNA sequence','error')

               raise ValidationError('Provide Valid DNA sequence')
     
class seqcrispr(jatayu):
    sequence = TextAreaField('Sequence', validators=[InputRequired('Please provide the genomic sequence'),base_check,Length(min=20,max=30,message='Must be between 20 and 30') ] ) 
    position = IntegerField('Position', validators=[InputRequired('Position of the mutation')])


class MyForm(FlaskForm):
    RSID = TextAreaField('RSID', validators=[ InputRequired()  ] ) 

"""
def base_check1(form, field):
    pattern = 'ORF1ab|S|ORF3a|E|M|ORF6|ORF7a|ORF8|N|ORF10' #r"^rs[\d]+$"''
    #print (pattern)
    if not field.data.contains(pattern):
        raise ValidationError('Not a Correct Mutation')
"""



    #if re.findall(pattern, field.data):
    #           raise ValidationError('Not a Correct RSID')
#ORF1ab|S|ORF3a|E|M|ORF6|ORF7a|ORF8|N|ORF10|    
#counter = 1    

@app.route('/')
def index():
    return render_template('front1.html')

#@app.route('/covid', methods=['GET', 'POST'])
#def form():

def cas_offinder(input,output):
    cmdArgs2=['cas-offinder' ,input,'C', output ] 
    print(cmdArgs2)
    call(cmdArgs2)
    #output.flush()
    columns=['crRNA','Chromosome','Position','Off-targets','Strand','No of Mismatches']
    ots=pd.read_csv(output,sep="\t",header=None,names=columns)
    #ots
    ots['Chromosome'] = ots['Chromosome'].str.split().str[0]

    ots=ots.sort_values(by='No of Mismatches')
    return (ots)
    #return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])

@app.route('/offtargets',methods=['GET','POST'])
def offtargets():
        my_var = request.args.get('my_var', None)
        my_var1 = request.args.get('my_var1', None)
        my_var2 = request.args.get('my_var2', None)

        print ("off1-target:",my_var,my_var1)
        print ("off1-target:",my_var)
        path1="offtargets/"+"cr"+ uuid.uuid4().hex  
        f1= path1+ ".txt"
        #file_only1=ntpath.basename(f1)
        file1 = open(f1,"w")
        if my_var2 == 'Human':
            ref='/home/asgar/jatayufeluda/GRCh38.p13.genome.fa'
        else:
            ref='/home/asgar/JATAYU2/sequence.fasta'    
        if my_var1 =="Fn/enFnCas9":
            NGG1='NNNNNNNNNNNNNNNNNNNNNGG'
            #my_varl=list(myvar)

            #my_var[21]='N'
            j1=ref+ "\n"+NGG1 + "\n" + my_var[0:20]+'NGG' + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])
        elif my_var1 =="enFnCas9" and my_var.endswith("GA") :
            NGA='NNNNNNNNNNNNNNNNNNNNNGA'
            j1=ref+ "\n"+NGA + "\n" + my_var[0:20]+'NGA' + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])
        elif my_var1 =="enFnCas9" and my_var.endswith("AG") :
            NAG='NNNNNNNNNNNNNNNNNNNNNAG'
            j1=ref+ "\n"+NAG + "\n" + my_var[0:20]+'NAG' + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])
        elif my_var1 =="AaCas12b"  :
            TTN='TTNNNNNNNNNNNNNNNNNNNNN'
            j1=ref+ "\n"+TTN + "\n" + 'TTN' + my_var[3:] + "\t" +str(4)

            #j1=ref+ "\n"+TTN + "\n" + 'TTN' + my_var[3:] + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])
        elif my_var1 =="LbCas12a"  :
            TTTN='TTTNNNNNNNNNNNNNNNNNNN'
            j1=ref+ "\n"+TTTN + "\n" + 'TTTN'+ my_var[4:] + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])    
        elif my_var1 =="Cas14a"  :
            TTTA='TTTANNNNNNNNNNNNNNNNNNNN'
            j1=ref+ "\n"+TTTA + "\n" + 'TTTA' + my_var[4:] + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])
        elif my_var1 =="Cas14a/ssDNA"  :
            ssdna='NNNNNNNNNNNNNNNNNNNN'
            j1=ref+ "\n"+ssdna + "\n" + my_var + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])    
        elif my_var1 =="LwCas13a"  :
            H='NNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            j1=ref+ "\n"+H + "\n" + my_var + "\t" +str(4)
            print (j1)
            file1.write(j1)
            file1.flush()
            f2=path1+"_ots"+".txt"
            file_only2=ntpath.basename(f2)
            file2= "offtargets/"+file_only2 
            print(file2)
            ots=cas_offinder(f1,file2)
            return render_template('ots.html',tables=[ots.to_html(classes='primer',index=None)],titles = [ 'na' ,'Off-targets'])

def align(wtype,len1,position,variation,ref,organism):
  db=get_db()
  path1="sgrna_results/"+"cr_"+ uuid.uuid4().hex     
  my_seqs = SeqRecord(Seq(wtype), id = "randomsequence")
  filename = path1 + ".fasta"
  SeqIO.write(my_seqs, filename, "fasta")
  filename_sam = path1 + ".sam"
  filename_bed = path1 + ".bed"
  filename_bed_gene = path1 + "_gene"
  #filename_bed_gene1 = path1 + "_gene1"
  #filename_bed_gene1_sgrna = path1 + "_sgrna"
  if organism =='Homo Sapiens':
    cmdArgs=['bwa' ,'mem'  ,'-t' ,'8', '-k', str(len1) , '-T' , '20', '-a' ,'/home/asgar/jatayufeluda/GRCh38.p13.genome.fa',filename, '-o' ,filename_sam] # 'bwa', 'mem' ,'a','-t', '4'  ,'GRCh38.p13.genome.fa' , 'example.fasta',  '-o', 'example.sam'  
    call(cmdArgs)
    df=pd.read_csv(filename_sam,sep='\t',header=None,comment="@")
    df=df.iloc[:,0:12]
    df.columns=['QNAME','Strand','CHROMOSOME','Start','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','NM']
    print(df)
    if df['CHROMOSOME'][0] == '*':
        print ("NOT HUMAN SEQ")
        return "Not a Homo Sapiens Sequence"
    else :    
     df['CIGAR'] = df['CIGAR'].str.replace(r'\D', '').astype(int)    
     df['End']=df['Start'] + df['CIGAR']   
     df=df.loc[df['NM']=='NM:i:0', ['CHROMOSOME','Start', 'End','Strand'] ]
     df.loc[df['Strand'] == 0, 'Strand'] = '+'
     df.loc[df['Strand'] == 256, 'Strand'] = '+'
     df.loc[df['Strand'] == 16, 'Strand'] = '-'
     df.loc[df['Strand'] == 272, 'Strand'] = '-'
     df['POS']=df['Start']-1+position
     df['REF']=ref
     df['ALT']=variation
     genome_pos=str(df['POS'].values[0])
     df.to_csv(filename_bed,sep='\t',index=None,header=None)
     print(df)
     cur=db.execute('select  jatayu,mut_loc_cr,ran_loc_cr, strand,CHROMOSOME,CHROM,POS,RS,REF,ALT,presgrna,prestart,preend,GENEINFO,CLNDN,oligo,oligo1,FREQ,Cas_system from ngr_nrg_snp_wd_freq where POS=?',[genome_pos])
     results=cur.fetchall()
     print(genome_pos,results)
     if results:
            df = pd.DataFrame(results, columns =[  'jatayu',   'mut_loc_cr' ,  'ran_loc_cr' ,'Strand','Chr','chrom','Position','rsID','Ref','Alt','presgrna','prestart','preend','Gene','Disease','oligo','oligo1','FREQ','Cas_system']) 
            df["Variation"] = df["Ref"].astype(str) + ">" +df["Alt"].astype(str)
            df1=df
            print(df1)
            alt_list=df1.Alt.unique()
            if df['FREQ'].values[0]=="Not available":
                plot_url=""
            else:    
             figfile = BytesIO() ####for plotting variation frequency####
             list1=Convert(df1['FREQ'][0]) 
             df_freq=pd.DataFrame(list1)
             df_freq.columns=['freq']
             df_freq['freq'] = df_freq['freq'].str.replace(r':', ',')
             #df_freq = df_freq.replace(dict.fromkeys(['.'],'0')) 
             print("df_Freq:",df_freq)
             df_freq1=df_freq.freq.apply(lambda x: pd.Series(str(x).split(",")))
             if len(df_freq1.columns) ==2:
              df_freq1.columns = ['Dataset',"Variant"]
              df_freq1.set_index('Dataset',inplace=True)
              df_freq1.rename(columns={'Variant':alt_list[0] + " Allele" }  , inplace=True) 
             elif len(df_freq1.columns) ==3:
              df_freq1.columns = ['Dataset',"Variant1","Variant2"]
              df_freq1.set_index('Dataset',inplace=True)
              df_freq1.rename(columns={'Variant1':alt_list[0] + " Allele",'Variant2':alt_list[1] + " Allele" }  , inplace=True) 
             elif len(df_freq1.columns) ==4:
              df_freq1.columns = ['Dataset',"Variant1","Variant2","Variant3"]
              df_freq1.rename(columns={'Variant1':alt_list[0] + " Allele",'Variant2':alt_list[1] + " Allele",'Variant3':alt_list[2] + " Allele" }  , inplace=True)
              df_freq1.set_index('Dataset',inplace=True)
             
             df_freq1 = df_freq1.replace(dict.fromkeys(['.'],'0'))
             print("df_Freq1:",df_freq1)
             df_freq1 = df_freq1.astype(float)
             pal = ["#00BFC4"  ,"#F8766D", "#00BA38", "#619CFF"]
             
            #plt.figure(figsize=(15,20))
             df_freq1.plot(kind="bar",color=pal,width=0.3,linewidth=0.4 ,  stacked=True ,edgecolor="black" ,figsize=(10, 5))
             plt.title("Variant allele frequency in different datasets")
             plt.xlabel("")
             plt.xticks(rotation=30,ha='right')
             plt.ylabel("Allele Frequency")
            #plt.tight_layout()
             ax = plt.subplot(111)
            #ax.xticks( rotation=90, ha='right')
             box = ax.get_position()
             ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
             ax.legend(loc='upper center', bbox_to_anchor=(0.9, 1.175), fancybox=True, shadow=True, ncol=3)
            #plt.legend( bbox_to_anchor=(0.5, 0., 0.5, 0.5),loc='best', borderaxespad=0.)
            
             plt.savefig(figfile, format='png')
             figfile.seek(0)
             plot_url = base64.b64encode(figfile.getvalue()).decode('utf8')

            chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=datafetch(df1,'Human')
            chrom=df1['chrom'][0]
            disease=df1['Disease'][0]#.str.replace("_"," ", regex=True)
            disease=disease.replace('_',' ')
            variation=df1["Ref"][0]+"/"+"/".join(alt_list)
            alleles =df1["Ref"][0]+" > " +"/".join(alt_list)
            comments1=rsID + " " +variation
            return chr_num,chrom,gene,rsID,var_pos,disease,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13
     else:
            cmdArgs3=['bedtools', 'intersect', '-a', str(filename_bed),'-b', '/home/asgar/jatayufeluda/gencode.v32.annotation_genes_filtered.gtf' , '-wb'   ]
            call(cmdArgs3,stdout=open(filename_bed_gene,'w'))
            bedgene=pd.read_csv(filename_bed_gene,header=None,sep="\t")
            bedgene=bedgene.iloc[:,[0,4, 5,6,10]]
            
            

            bedgene.columns=["Chr", "Position" , "Ref", "Alt","Gene"]
            bedgene.insert(2, "rsID", bedgene['Ref'].astype(str)+bedgene['Position'].astype(str)+bedgene['Alt'].astype(str), True)
            df=ngr_nrg.sg_df(bedgene)
            df["Variation"] = df["Ref"].astype(str) + ">" +df["Alt"].astype(str)

            chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=datafetch(df,'Human')
            print (chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13)
            chrom=df['CHROM'][0]
            disease="Not available"
            variation=df["Variation"][0]
            print (variation)
            alleles =df["Variation"][0]
            comments1=rsID + " " +variation
            plot_url=""
            print(var_pos)
            return chr_num,chrom,gene,rsID,var_pos,disease,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13
  elif organism =='Sars-CoV-2':
    cmdArgs=['bwa' ,'mem'  ,'-t' ,'8', '-k', str(len1) , '-T' , '20', '-a' ,'sequence.fasta',filename, '-o' ,filename_sam] # 'bwa', 'mem' ,'a','-t', '4'  ,'GRCh38.p13.genome.fa' , 'example.fasta',  '-o', 'example.sam'  
    call(cmdArgs)
    df=pd.read_csv(filename_sam,sep='\t',header=None,comment="@")
    df=df.iloc[:,0:12]
    df.columns=['QNAME','Strand','CHROMOSOME','Start','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','NM']
    print(df)
    if df['CHROMOSOME'][0] == '*':
        print ("NOT Sars-CoV-2")
        return "Not a Sars-CoV-2 Sequence"
    else :    
     df['CIGAR'] = df['CIGAR'].str.replace(r'\D', '').astype(int)    
     df['End']=df['Start'] + df['CIGAR']   
     df=df.loc[df['NM']=='NM:i:0', ['CHROMOSOME','Start', 'End','Strand'] ]
     df.loc[df['Strand'] == 0, 'Strand'] = '+'
     df.loc[df['Strand'] == 256, 'Strand'] = '+'
     df.loc[df['Strand'] == 16, 'Strand'] = '-'
     df.loc[df['Strand'] == 272, 'Strand'] = '-'
     df['POS']=df['Start']-1+position
     df['REF']=ref
     df['ALT']=variation
     genome_pos=str(df['POS'].values[0])
     df.to_csv(filename_bed,sep='\t',index=None,header=None)
     print(df)
     cur=db.execute('select  jatayu,mut_loc_cr,ran_loc_cr, strand,CHROMOSOME,POS,RS,REF,ALT,presgrna,prestart,preend,Gene,oligo,oligo1,Virus_num,Annotation_Type,Gene_Pos_Codons,Cas_system from ngr_nrg_snp_wd_freq_covid where POS=?',[genome_pos])

     #cur=db.execute('select  jatayu,mut_loc_cr,ran_loc_cr, strand,CHROMOSOME,CHROM,POS,RS,REF,ALT,presgrna,prestart,preend,GENEINFO,CLNDN,oligo,oligo1,FREQ,Cas_system from ngr_nrg_snp_wd_freq_covid where POS=?',[genome_pos])
     results=cur.fetchall()
     print(genome_pos,results)
     if results:
            df = pd.DataFrame(results, columns =[  'jatayu',   'mut_loc_cr' ,  'ran_loc_cr' ,'Strand','Chr','Position','rsID','Ref','Alt','presgrna','prestart','preend','Gene','oligo','oligo1','FREQ','Annotation_Type','Gene_Pos_Codons','Cas_system']) 
            df = df[df['Alt'].str.contains(variation)]
            #df=df['Alt'].str.contains(variation)
            print ("variation:", df['Alt'],variation)
            df["Variation"] = df["Ref"].astype(str) + ">" +df["Alt"].astype(str)
            df1=df.reset_index()
            #print (df1)
            alt_list=df1.Alt.unique()
            #print (alt_list)
            figfile = BytesIO() ####for plotting variation frequency####
            print (df1['FREQ']) 
            #list1=Convert(df1['FREQ'][0]) 
            df_freq=df1[['FREQ','Alt']]

            df_freq.columns=['freq','alt']
            df_freq.sort_values("alt", inplace = True)

            df_freq.drop_duplicates(keep = 'first', inplace = True)
 
            print("df_Freq:",df_freq)

            #df_freq['freq'] = df_freq['freq'].str.replace(r':', ',')
            
            #df_freq = df_freq.replace(dict.fromkeys(['.'],'0')) 
            print(type(df_freq))

            df_freq1=df_freq
             
            print("df_Freq1:",df_freq1)
            pal = ["#00BFC4"  ,"#F8766D", "#00BA38", "#619CFF"]
            df_freq1.plot(kind="bar",color=pal,width=0.2,linewidth=0.4 ,x='alt',y='freq',legend=False, edgecolor="black" ,figsize=(8, 4))
            
            plt.title("No. of Virus genomes with mutation (Total 3567389 Sequences)")
            plt.xlabel("")
            
            plt.xticks(rotation=0,ha='right')
            plt.ylabel("No. of genomes")
            #plt.tight_layout()
            ax = plt.subplot(111)
            #ax.xticks( rotation=90, ha='right')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
            
            plt.savefig(figfile, format='png')
            figfile.seek(0)
            plot_url = base64.b64encode(figfile.getvalue()).decode('utf8')
              
            chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=datafetch(df1,'Sars-CoV-2')

            variation=df1["Ref"][0]+"/"+"/".join(alt_list)
            alleles =df1["Ref"][0]+" > " +"/".join(alt_list)
            var_pos=var_pos+1
            #rsID=
            comments1=rsID + " " +variation
            return chr_num,gene,rsID,var_pos,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13
     else:
            cmdArgs3=['bedtools', 'intersect', '-a', str(filename_bed),'-b', 'sars_cov2_filtered.gtf' , '-wb'   ]
            call(cmdArgs3,stdout=open(filename_bed_gene,'w'))
            bedgene=pd.read_csv(filename_bed_gene,header=None,sep="\t")
            bedgene=bedgene.iloc[:,[0,4, 5,6,10]]
            
            

            bedgene.columns=["Chr", "Position" , "Ref", "Alt","Gene"]
            bedgene.insert(2, "rsID", bedgene['Ref'].astype(str)+bedgene['Position'].astype(str)+bedgene['Alt'].astype(str), True)
            df=ngr_nrg.sg_df(bedgene)
            df["Variation"] = df["Ref"].astype(str) + ">" +df["Alt"].astype(str)

            chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=datafetch(df,'Human')
            print (chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13)
            chrom=df['CHROM'][0]
            disease="Not available"
            variation=df["Variation"][0]
            print (variation)
            alleles =df["Variation"][0]
            comments1=rsID + " " +variation
            plot_url=""
            print(var_pos)
            return chr_num,chrom,gene,rsID,var_pos,disease,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13
          


@app.route('/seq', methods=['GET', 'POST'])
def form2():
    form=seqcrispr()
    if request.method == 'POST' :

     #print(form.sequence.data,form.position.data)
     s =list(form.sequence.data)
     #pattern=['GG','GA','AG']
     #mut_pam = ("".join((s[form.position.data-1:form.position.data+4 ] )) )
     print(form.sequence.data,form.position.data,form.variation.data)
     if (form.variation.data == s[form.position.data-1]) :
         flash('Mutation is same as wild type','error')
         return redirect(url_for('form2'))
     else:
        #mutation_type=  s[form.position.data-1] + ">" + form.variation.data  
        ref=s[form.position.data-1] 
        s[form.position.data-1]=form.variation.data
        mutation_sequence="".join(s)
        sq_len=len(form.sequence.data)
        print(ref,mutation_sequence,form.organism.data)
        organism=form.organism.data
        try:
         if organism == "Sars-CoV-2":
            chr_num,gene,rsID,var_pos,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13 =align(form.sequence.data,sq_len,form.position.data,form.variation.data,ref,organism)
              #print (chr_num,gene,rsID,var_pos,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13)
            if plot_url== '':
                html2='seqcrisp_results_novel_covid.html'
            else:
                html2='seqcrisp_results_covid.html'     
            return render_template(html2,chr=chr_num,gene=gene, comments1=comments1,rsid=rsID,allele=alleles, primer1_fw=primer1_fw_fn ,primer1_rw=primer1_rw_fn,primer2_fw=primer2_fw_fn,primer2_rw=primer2_rw_fn,var_pos=var_pos ,coo=coo,coo1=coo1_fn,coo1_cas12b=coo1_cas12b,  plot=plot_url ,data1=Fn,cas12b=cas12b,cas12a=cas12a,cas14=cas14,cas14ss =cas14ss,cas13=cas13, tables=[primer_fn.to_html(classes='primer',index=None),primer_cas12b.to_html(classes='primer',index=None),primer_cas12a.to_html(classes='primer',index=None),primer_cas14.to_html(classes='primer',index=None),primer_cas14ss.to_html(classes='primer',index=None),primer_cas13.to_html(classes='primer',index=None)  ],
          titles = [ 'na' ,'Primers','Primers',  'Primers'  , 'Primers' ,'Primers','Primers'])
         elif organism == "Homo Sapiens":
            chr_num,chrom,gene,rsID,var_pos,disease,variation,alleles,comments1,plot_url, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=align(form.sequence.data,sq_len,form.position.data,form.variation.data,ref,organism)
            if rsID.startswith('rs') and plot_url :
                html1='seqcrisp_results.html'     
            elif rsID.startswith('rs') and plot_url=='':
                html1='seqcrisp_results_nofreq.html'    

                
            else:
                html1='seqcrisp_results_novel.html'
            return render_template(html1,chr=chr_num,gene=gene,disease= disease,  chrom=chrom,comments1=comments1,rsid=rsID,allele=alleles, primer1_fw=primer1_fw_fn ,primer1_rw=primer1_rw_fn,primer2_fw=primer2_fw_fn,primer2_rw=primer2_rw_fn,var_pos=var_pos ,coo=coo,coo1=coo1_fn,coo1_cas12b=coo1_cas12b,  plot=plot_url ,data1=Fn,cas12b=cas12b,cas12a=cas12a,cas14=cas14,cas14ss =cas14ss,cas13=cas13, tables=[primer_fn.to_html(classes='primer',index=None),primer_cas12b.to_html(classes='primer',index=None),primer_cas12a.to_html(classes='primer',index=None),primer_cas14.to_html(classes='primer',index=None),primer_cas14ss.to_html(classes='primer',index=None),primer_cas13.to_html(classes='primer',index=None)  ],
          titles = [ 'na' ,'Primers','Primers',  'Primers'  , 'Primers' ,'Primers','Primers']) 
        except:
             if align(form.sequence.data,sq_len,form.position.data,form.variation.data,ref,organism) == "Not a Sars-CoV-2 Sequence" :
                  flash('Not A Sars-CoV-2 Sequence, Try A Valid Sequence','error')
                  return redirect(url_for('form2'))
             elif align(form.sequence.data,sq_len,form.position.data,form.variation.data,ref,organism) == "Not a Homo Sapiens Sequence" :
                  flash('Not A Homo Sapiens Sequence, Try A Valid Sequence','error')
                  return redirect(url_for('form2'))
 

        #if form.validate_on_submit():
        # try:
          
          
            #align(form.sequence.data,sq_len,form.position.data,form.variation.data)
            #return render_template(html1,chr=chr_num,gene=gene,comments1=comments1,rsid=rsID,allele=alleles, primer1_fw=primer1_fw_fn ,primer1_rw=primer1_rw_fn,primer2_fw=primer2_fw_fn,primer2_rw=primer2_rw_fn,var_pos=var_pos ,coo=coo,coo1=coo1_fn,coo1_cas12b=coo1_cas12b,  plot=plot_url ,data1=Fn,cas12b=cas12b,cas12a=cas12a,cas14=cas14,cas14ss =cas14ss,cas13=cas13, tables=[primer_fn.to_html(classes='primer',index=None),primer_cas12b.to_html(classes='primer',index=None),primer_cas12a.to_html(classes='primer',index=None),primer_cas14.to_html(classes='primer',index=None),primer_cas14ss.to_html(classes='primer',index=None),primer_cas13.to_html(classes='primer',index=None)  ],
          #titles = [ 'na' ,'Primers','Primers',  'Primers'  , 'Primers' ,'Primers','Primers'])
          
         #except:
          #    if align(form.sequence.data,sq_len,form.position.data,form.variation.data,ref,organism) == "Not a Homo Sapiens Sequence" :
           #       flash('Not A Homo Sapiens Sequence, Try A Valid Sequence')
            #      return redirect(url_for('form2'))
             # elif align(form.sequence.data,sq_len,form.position.data,form.variation.data,ref,organism) == "Not a Sars-CoV-2 Sequence" :
              #    flash('Not A Sars-CoV-2 Sequence, Try A Valid Sequence')
               #   return redirect(url_for('form2'))
        
         #if organism == "Sars-CoV-2":
         # print (organism)   
          #try:
              
          #except:
           
          
 
         

    return render_template('seqcrisp.html',form=form)#,visitors=counter)

def datafetch(df1,organism):
            if organism =='Human':
                list_fnh=[ [200,300],  [300,400],[400,500], [500,600], [600,700],[700,800],[800,900],[900,1000],[1000,1100],[1100,1200]   ]
                cut_pos1=418
                cut_pos2=1018
                list1=[ [100,200],[200,300] , [300,400],[400,500],[500,600],[600,700],[700,800]   ]
            elif organism =='Sars-CoV-2':
                list_fnh=[ [100,200],  [200,300],  [300,400],[400,500], [500,600], [600,700],[700,800],[800,900]  ]
                cut_pos1=318
                cut_pos2=618
                list1=[ [100,200],[200,300] , [300,400],[400,500 ],[500,600 ],[600,700 ]   ]
            df1['crRNA Position']=df1['Chr'] +":" +df1['prestart']+"-"+df1['preend']
            df1.rename(columns={"jatayu": "crRNA Sequence","presgrna":"Wt Sequence",'Cas_system' : "Cas system"   },inplace=True     )
            df1 = df1.astype({"prestart": int}) 
            Fn = df1.loc[df1['Cas system'].str.contains ("Fn/enFnCas9|enFnCas9")]
            cas13 = df1[df1['Cas system'] == "LwCas13a"]
            cas14 = df1[df1['Cas system'] == "Cas14a"]
            cas14ss = df1[df1['Cas system'] == "Cas14a/ssDNA"]
            cas12b = df1[df1['Cas system'] == "AaCas12b"]
            cas12a =df1[df1['Cas system'] == "LbCas12a"]
            if Fn.empty : 
             Fn = pd.DataFrame({'crRNA Sequence':'No crRNA','Wt Sequence' :'-','crRNA Position':"-",  'Strand':"-",'Variation':"-",'Cas system':"-",'mut_loc_cr':[1],'ran_loc_cr':[4] })   
             primer_fn=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
             primer1_fw_fn=primer1_rw_fn=primer2_fw_fn=primer2_rw_fn= coo=coo1_fn=""
            elif Fn['oligo'].values[0] == 'Not available' :
              Fn=Fn.reset_index()  
              Fn['Off-Targets']="Click for Off-targets"
              #'jatayu',   'mut_loc_cr' ,  'ran_loc_cr' ,'Strand','Chr','Position','rsID','Ref','Alt','presgrna','prestart','preend','Gene','oligo','oligo1','FREQ','Annotation_Type','Gene_Pos_Codons','Cas_system'
              Fn=Fn[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation',  'Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]                  
              primer_fn=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
              primer1_fw_fn=primer1_rw_fn=primer2_fw_fn=primer2_rw_fn= coo=coo1_fn=""
            else:   
             Fn=Fn.reset_index()  
             #session['crRNA']=Fn['crRNA Sequence'].to_list()    
             crr1_fn=Fn['prestart'][0]
             crr2_fn=Fn['preend'][0]
             coo1_fn=str(crr1_fn) + ":" +str(crr2_fn)
             #list_fnh=[ [200,300],  [300,400],[400,500], [500,600], [600,700],[700,800],[800,900],[900,1000],[1000,1100],[1100,1200]   ]

             primer1_fn=desg_primer_wd_pam(Fn['oligo'][0],Fn.iloc[0]['prestart'] ,cut_pos1,list_fnh)
             primer2_fn=desg_primer_wd_pam(Fn['oligo1'][0], Fn.iloc[0]['prestart'],cut_pos2,list_fnh)
             #plist= [primer1,primer2]
             primerpre_fn=pd.concat([primer1_fn,primer2_fn])
             #primer=primerpre.drop_duplicates()
             #primerpre['combined_primer']=primerpre['Forward Primer']   + primerpre['Reverse Primer']
             primer_fn = primerpre_fn.drop_duplicates( subset = [  'Forward Primer','Reverse Primer'], keep = 'first').reset_index(drop = True) 
             #primer = primerpre.groupby(['Forward Primer', 'Reverse Primer']).first() 
             primer_fn=primer_fn.sort_values(by='Amplicon (bp)').head(5)
             primer_fn.insert(0, 'Pair No.', range(1, 1 + len(primer_fn)))
             ##print (primer)
             ##print (primer.iloc[0]['primer_locus'])
             primer1_fw_fn = str( primer_fn.iloc[0]['primer_locus']) + ":" + str(  primer_fn.iloc[0]['primer_locus'] + primer_fn.iloc[0]['Length']  -1    )
             primer1_rw_fn = str( primer_fn.iloc[0]['primer_locus1']) + ":" + str(  primer_fn.iloc[0]['primer_locus1'] - primer_fn.iloc[0]['Length'] -1     )
             primer2_fw_fn = str( primer_fn.iloc[1]['primer_locus']) + ":" + str(  primer_fn.iloc[1]['primer_locus'] + primer_fn.iloc[1]['Length']   -1 )
             primer2_rw_fn = str( primer_fn.iloc[1]['primer_locus1']) + ":" + str(  primer_fn.iloc[1]['primer_locus1'] - primer_fn.iloc[1]['Length'] -1     ) 
             coo=str(primer_fn.iloc[0:2]['primer_locus'].min()-50) + ":" +str(primer_fn.iloc[0:2]['primer_locus1'].max()+50)
             ##print (coo,primer.iloc[0:1]['primer_locus'])
             primer_fn=primer_fn.drop(['Start','End','primer_locus','primer_locus1'],axis=1)
             #Fn['crRNA Position']=Fn['Chr'] +":" +Fn['prestart'].astype(str)+"-"+Fn['preend'].astype(str)
             #link_frmt = lambda x: '''<a href="" target="popup" onclick="window.open('offtargets/','crRNA Sequence','width=600,height=400')">{0}</a>'''.format(x)
             #link_frmt = lambda x:'''<a href="{{'url_for('templates/', filename='ots.html')}}">{0}</a>'''.format(x)
             #link_frmt= lambda x:'''<a href="{{ url_for('tempalates', 'offtargets', my_var='TGACTGTCC') }}">{0}</a>'''.format(x)
             
             #<a href="{{ url_for('offtargets',my_var= data.iloc[k,0]) }}" target="_blank"> {{data.iloc[k,5]}} </a>
             Fn['Off-Targets']="Click for Off-targets"
             #link_frmt= lambda x:'''<a href="{{ url_for('offtargets',my_var= {0}) }}" target="_blank"> {0} </a>'''.format(x)
             #Fn['Off-Targets']=Fn['crRNA Sequence']
             print (Fn.columns,Fn)
             Fn=Fn[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]
            if cas12b.empty: 
             cas12b = pd.DataFrame({'crRNA Sequence':'No crRNA','Wt Sequence' :'-','crRNA Position':"-",  'Strand':"-",'Variation':"-",'Cas system':"-",'mut_loc_cr':[1],'ran_loc_cr':[4] })   
             primer_cas12b=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
             primer1_fw_cas12b=primer1_rw_cas12b=primer2_fw_cas12b=primer2_rw_ca12=coo_cas12b=coo1_cas12b=""
            elif cas12b['oligo'].values[0] == 'Not available' :
              cas12b=cas12b.reset_index()  
              cas12b['Off-Targets']="Click for Off-targets"
              cas12b=cas12b[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]                  
              primer_cas12b=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
              primer1_fw_cas12b=primer1_rw_cas12b=primer2_fw_cas12b=primer2_rw_ca12=coo_cas12b=coo1_cas12b="" 
            else:   
             cas12b=cas12b.reset_index()
             ##print (cas12b)
             crr1_cas12b=cas12b['prestart'][0]
             crr2_cas12b=cas12b['preend'][0]
             coo1_cas12b=str(crr1_cas12b) + ":" +str(crr2_cas12b)
             #print ("cas12b:" ,cas12b['oligo1'][0],cas12b['preend'][0])
             #chrom=df['Position'][0]
             ##print (chrom)
             ##print (coo)
             ##print(df1['oligo'][0],df1.iloc[0]['prestart'])
             ##print(df1['oligo1'][0],df1.iloc[0]['prestart'])
             #list1=[ [100,200],[200,300] , [300,400],[400,500],[500,600],[600,700],[700,800]   ]

             primer1_cas12b=desg_primer_cas(cas12b['oligo'][0],cas12b.iloc[0]['prestart'],20,list1)
             #primer2_cas12b=desg_primer_cas(cas12b['oligo1'][0], cas12b.iloc[0]['prestart'],268)
             #plist= [primer1,primer2]
             #primerpre_cas12b=pd.concat([primer1_cas12b,primer2_cas12b])
             primer_cas12b = primer1_cas12b.drop_duplicates( subset = [  'Forward Primer','Reverse Primer'], keep = 'first').reset_index(drop = True) 
             #primer = primerpre.groupby(['Forward Primer', 'Reverse Primer']).first() 
             #print ('primer12: ',primer_cas12b)
             primer_cas12b=primer_cas12b.sort_values(by='Amplicon (bp)').head(5)
             primer_cas12b.insert(0, 'Pair No.', range(1, 1 + len(primer_cas12b)))
             primer1_fw_cas12b = str( primer_cas12b.iloc[0]['primer_locus']) + ":" + str(  primer_cas12b.iloc[0]['primer_locus'] + primer_cas12b.iloc[0]['Length']  -1    )
             primer1_rw_ca12 = str( primer_cas12b.iloc[0]['primer_locus1']) + ":" + str(  primer_cas12b.iloc[0]['primer_locus1'] - primer_cas12b.iloc[0]['Length'] -1     )
             primer2_fw_cas12b = str( primer_cas12b.iloc[1]['primer_locus']) + ":" + str(  primer_cas12b.iloc[1]['primer_locus'] + primer_cas12b.iloc[1]['Length']   -1 )
             primer2_rw_cas12b = str( primer_cas12b.iloc[1]['primer_locus1']) + ":" + str(  primer_cas12b.iloc[1]['primer_locus1'] - primer_cas12b.iloc[1]['Length'] -1     ) 
             coo_cas12b=str(primer_cas12b.iloc[0:2]['primer_locus'].min()-50) + ":" +str(primer_cas12b.iloc[0:2]['primer_locus1'].max()+50)
             primer_cas12b=primer_cas12b.drop(['Start','End','primer_locus','primer_locus1','IVC Product (bp)'],axis=1)
             ##print ('primer12':primer_cas12b)
             cas12b['Off-Targets']="Click for Off-targets"

             cas12b=cas12b[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]
            if cas12a.empty: 
             cas12a = pd.DataFrame({'crRNA Sequence':'No crRNA','Wt Sequence' :'-','crRNA Position':"-",  'Strand':"-",'Variation':"-",'Cas system':"-",'mut_loc_cr':[1],'ran_loc_cr':[4] })   
             primer_cas12a=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
             primer1_fw_cas12a=primer1_rw_cas12a=primer2_fw_cas12a=primer2_rw_ca12=coo_cas12a=coo1_cas12a=""
            elif cas12a['oligo'].values[0] == 'Not available' :
              cas12a=cas12a.reset_index()  
              cas12a['Off-Targets']="Click for Off-targets"
              cas12a=cas12a[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]                  
              primer_cas12a=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
              primer1_fw_cas12a=primer1_rw_cas12a=primer2_fw_cas12a=primer2_rw_ca12=coo_cas12a=coo1_cas12a=""
 
            else:   
             cas12a=cas12a.reset_index()
             ##print (cas12a)
             crr1_cas12a=cas12a['prestart'][0]
             crr2_cas12a=cas12a['preend'][0]
             coo1_cas12a=str(crr1_cas12a) + ":" +str(crr2_cas12a)
             #print ("cas12a:" ,cas12a['oligo1'][0],cas12a['preend'][0])
             #chrom=df['Position'][0]
             ##print (chrom)
             ##print (coo)
             ##print(df1['oligo'][0],df1.iloc[0]['prestart'])
             ##print(df1['oligo1'][0],df1.iloc[0]['prestart'])
             #list1=[ [100,200],[200,300] , [300,400],[400,500],[500,600],[600,700],[700,800]   ]

             primer1_cas12a=desg_primer_cas(cas12a['oligo'][0],cas12a.iloc[0]['prestart'],20,list1)
             #primer2_cas12a=desg_primer_cas(cas12a['oligo1'][0], cas12a.iloc[0]['prestart'],268)
             #plist= [primer1,primer2]
             #primerpre_cas12a=pd.concat([primer1_cas12a,primer2_cas12a])
             primer_cas12a = primer1_cas12a.drop_duplicates( subset = [  'Forward Primer','Reverse Primer'], keep = 'first').reset_index(drop = True) 
             #primer = primerpre.groupby(['Forward Primer', 'Reverse Primer']).first() 
             #print ('primer12: ',primer_cas12a)
             primer_cas12a=primer_cas12a.sort_values(by='Amplicon (bp)').head(5)
             primer_cas12a.insert(0, 'Pair No.', range(1, 1 + len(primer_cas12a)))
             primer1_fw_cas12a = str( primer_cas12a.iloc[0]['primer_locus']) + ":" + str(  primer_cas12a.iloc[0]['primer_locus'] + primer_cas12a.iloc[0]['Length']  -1    )
             primer1_rw_ca12 = str( primer_cas12a.iloc[0]['primer_locus1']) + ":" + str(  primer_cas12a.iloc[0]['primer_locus1'] - primer_cas12a.iloc[0]['Length'] -1     )
             primer2_fw_cas12a = str( primer_cas12a.iloc[1]['primer_locus']) + ":" + str(  primer_cas12a.iloc[1]['primer_locus'] + primer_cas12a.iloc[1]['Length']   -1 )
             primer2_rw_cas12a = str( primer_cas12a.iloc[1]['primer_locus1']) + ":" + str(  primer_cas12a.iloc[1]['primer_locus1'] - primer_cas12a.iloc[1]['Length'] -1     ) 
             coo_cas12a=str(primer_cas12a.iloc[0:2]['primer_locus'].min()-50) + ":" +str(primer_cas12a.iloc[0:2]['primer_locus1'].max()+50)
             primer_cas12a=primer_cas12a.drop(['Start','End','primer_locus','primer_locus1','IVC Product (bp)'],axis=1)
             ##print ('primer12':primer_cas12a)
             cas12a['Off-Targets']="Click for Off-targets"

             cas12a=cas12a[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']] 
            if cas14.empty:
             cas14 = pd.DataFrame({'crRNA Sequence':'No crRNA','Wt Sequence' :'-','crRNA Position':"-",  'Strand':"-",'Variation':"-",'Cas system':"-",'mut_loc_cr':[1],'ran_loc_cr':[4] })   
             primer_cas14=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
             primer1_fw_cas14 = primer1_rw_cas14 = primer2_fw_cas14 = primer2_rw_cas14=coo_cas14=coo1_cas14= ""
            elif cas14['oligo'].values[0] == 'Not available' :
              cas14=cas14.reset_index()  
              cas14['Off-Targets']="Click for Off-targets"
              cas14=cas14[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]                  
              primer_cas14=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
              primer1_fw_cas14 = primer1_rw_cas14 = primer2_fw_cas14 = primer2_rw_cas14=coo_cas14=coo1_cas14= ""
 
            else:
             cas14=cas14.reset_index()
             #print (cas14)
             crr1_cas14=cas14['prestart'][0]
             crr2_cas14=cas14['preend'][0]
             coo1_cas14=str(crr1_cas14) + ":" +str(crr2_cas14)
             #print (cas14['preend'][0])
             #chrom=df['Position'][0]
             ##print (chrom)
             ##print (coo)
             ##print(df1['oligo'][0],df1.iloc[0]['prestart'])
             ##print(df1['oligo1'][0],df1.iloc[0]['prestart'])
             #list1=[ [100,200],[200,300] , [300,400],[400,500],[500,600],[600,700],[700,800]   ]

             primer1_cas14=desg_primer_cas(cas14['oligo'][0],cas14.iloc[0]['prestart'],20,list1)
             #primer2_cas14=desg_primer_cas(cas14['oligo1'][0], cas14.iloc[0]['prestart'],268)
             #plist= [primer1,primer2]
             #primerpre_cas14=pd.concat([primer1_cas14,primer2_cas14])
             primer_cas14 = primer1_cas14.drop_duplicates( subset = [  'Forward Primer','Reverse Primer'], keep = 'first').reset_index(drop = True) 
             #primer = primerpre.groupby(['Forward Primer', 'Reverse Primer']).first() 
             primer_cas14=primer_cas14.sort_values(by='Amplicon (bp)').head(5)
             primer_cas14.insert(0, 'Pair No.', range(1, 1 + len(primer_cas14)))
             primer1_fw_cas14 = str( primer_cas14.iloc[0]['primer_locus']) + ":" + str(  primer_cas14.iloc[0]['primer_locus'] + primer_cas14.iloc[0]['Length']  -1    )
             primer1_rw_cas14 = str( primer_cas14.iloc[0]['primer_locus1']) + ":" + str(  primer_cas14.iloc[0]['primer_locus1'] - primer_cas14.iloc[0]['Length'] -1     )
             primer2_fw_cas14 = str( primer_cas14.iloc[1]['primer_locus']) + ":" + str(  primer_cas14.iloc[1]['primer_locus'] + primer_cas14.iloc[1]['Length']   -1 )
             primer2_rw_cas14 = str( primer_cas14.iloc[1]['primer_locus1']) + ":" + str(  primer_cas14.iloc[1]['primer_locus1'] - primer_cas14.iloc[1]['Length'] -1     ) 
             coo_cas14=str(primer_cas14.iloc[0:2]['primer_locus'].min()-50) + ":" +str(primer_cas14.iloc[0:2]['primer_locus1'].max()+50)
             primer_cas14=primer_cas14.drop(['Start','End','primer_locus','primer_locus1','IVC Product (bp)'],axis=1)
             cas14['Off-Targets']="Click for Off-targets"

             cas14=cas14[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']] 
            if cas14ss.empty:
             cas14ss = pd.DataFrame({'crRNA Sequence':'No crRNA','Wt Sequence' :'-','crRNA Position':"-",  'Strand':"-",'Variation':"-",'Cas system':"-",'mut_loc_cr':[1],'ran_loc_cr':[4] })   
             primer_cas14ss=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
             primer1_fw_cas14ss = primer1_rw_cas14ss = primer2_fw_cas14ss = primer2_rw_cas14ss=coo_cas14ss=coo1_cas14ss= ""
            elif cas14ss['oligo'].values[0] == 'Not available' :
              cas14ss=cas14ss.reset_index()  
              cas14ss['Off-Targets']="Click for Off-targets"
              cas14ss=cas14ss[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]                  
              primer_cas14ss=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
              primer1_fw_cas14ss = primer1_rw_cas14ss = primer2_fw_cas14ss = primer2_rw_cas14ss=coo_cas14ss=coo1_cas14ss= ""
 
            else:
             cas14ss=cas14ss.reset_index()
             #print (cas14ss)
             crr1_cas14ss=cas14ss['prestart'][0]
             crr2_cas14ss=cas14ss['preend'][0]
             coo1_cas14ss=str(crr1_cas14ss) + ":" +str(crr2_cas14ss)
             #print (cas14ss['preend'][0])
             #chrom=df['Position'][0]
             ##print (chrom)
             ##print (coo)
             ##print(df1['oligo'][0],df1.iloc[0]['prestart'])
             ##print(df1['oligo1'][0],df1.iloc[0]['prestart'])
             #list1=[ [100,200],[200,300] , [300,400],[400,500],[500,600],[600,700],[700,800]   ]
             primer1_cas14ss=desg_primer_cas(cas14ss['oligo'][0],cas14ss.iloc[0]['prestart'],20,list1)
             #primer2_cas14ss=desg_primer_cas(cas14ss['oligo1'][0], cas14ss.iloc[0]['prestart'],268)
             #plist= [primer1,primer2]
             #primerpre_cas14ss=pd.concat([primer1_cas14ss,primer2_cas14ss])
             primer_cas14ss = primer1_cas14ss.drop_duplicates( subset = [  'Forward Primer','Reverse Primer'], keep = 'first').reset_index(drop = True) 
             #primer = primerpre.groupby(['Forward Primer', 'Reverse Primer']).first() 
             primer_cas14ss=primer_cas14ss.sort_values(by='Amplicon (bp)').head(5)
             primer_cas14ss.insert(0, 'Pair No.', range(1, 1 + len(primer_cas14ss)))
             primer1_fw_cas14ss = str( primer_cas14ss.iloc[0]['primer_locus']) + ":" + str(  primer_cas14ss.iloc[0]['primer_locus'] + primer_cas14ss.iloc[0]['Length']  -1    )
             primer1_rw_cas14ss = str( primer_cas14ss.iloc[0]['primer_locus1']) + ":" + str(  primer_cas14ss.iloc[0]['primer_locus1'] - primer_cas14ss.iloc[0]['Length'] -1     )
             primer2_fw_cas14ss = str( primer_cas14ss.iloc[1]['primer_locus']) + ":" + str(  primer_cas14ss.iloc[1]['primer_locus'] + primer_cas14ss.iloc[1]['Length']   -1 )
             primer2_rw_cas14ss = str( primer_cas14ss.iloc[1]['primer_locus1']) + ":" + str(  primer_cas14ss.iloc[1]['primer_locus1'] - primer_cas14ss.iloc[1]['Length'] -1     ) 
             coo_cas14ss=str(primer_cas14ss.iloc[0:2]['primer_locus'].min()-50) + ":" +str(primer_cas14ss.iloc[0:2]['primer_locus1'].max()+50)
             primer_cas14ss=primer_cas14ss.drop(['Start','End','primer_locus','primer_locus1','IVC Product (bp)'],axis=1)
             cas14ss['Off-Targets']="Click for Off-targets"

             cas14ss=cas14ss[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]  
            if cas13.empty:
             cas13 = pd.DataFrame({'crRNA Sequence':'No crRNA','Wt Sequence' :'-','crRNA Position':"-",  'Strand':"-",'Variation':"-",'Cas system':"-",'mut_loc_cr':[1],'ran_loc_cr':[4] })   
             primer_cas13=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
             primer1_fw_cas13 = primer1_rw_cas13 = primer2_fw_cas13 = primer2_rw_cas13=coo_cas13=coo1_cas13= ""
            elif cas13['oligo'].values[0] == 'Not available' :
              cas13=cas13.reset_index()  
              cas13['Off-Targets']="Click for Off-targets"
              cas13=cas13[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system','Off-Targets', 'mut_loc_cr','ran_loc_cr']]                  
              primer_cas13=pd.DataFrame({'Pair No.':'1','Forward Primer':['-'], 	'Length':"-", 	'Tm (Cel.)':'-', 'GC%':'-','Self-Com. (TH)':'-','Reverse Primer':'-', 'Length ':'-','Tm (Cel.) ': '-', 'GC% ':'-','Self-Com. (TH)':'-','Amplicon (bp)':'-','IVC Product (bp)':'-' })
              primer1_fw_cas13 = primer1_rw_cas13 = primer2_fw_cas13 = primer2_rw_cas13=coo_cas13=coo1_cas13= ""
 
            else:
             cas13=cas13.reset_index()
             #print (cas13)
             crr1_cas13=cas13['prestart'][0]
             crr2_cas13=cas13['preend'][0]
             coo1_cas13=str(crr1_cas13) + ":" +str(crr2_cas13)
             #print (cas13['preend'][0])
             #chrom=df['Position'][0]
             ##print (chrom)
             ##print (coo)
             ##print(df1['oligo'][0],df1.iloc[0]['prestart'])
             ##print(df1['oligo1'][0],df1.iloc[0]['prestart'])
             #list1=[ [100,200],[200,300] , [300,400],[400,500],[500,600],[600,700],[700,800]   ]

             primer1_cas13=desg_primer_cas(cas13['oligo'][0],cas13.iloc[0]['prestart'],28,list1)
             #primer2_cas13=desg_primer_cas(cas13['oligo1'][0], cas13.iloc[0]['prestart'],268)
             #plist= [primer1,primer2]
             #primerpre_cas13=pd.concat([primer1_cas13,primer2_cas13])
             primer_cas13 = primer1_cas13.drop_duplicates( subset = [  'Forward Primer','Reverse Primer'], keep = 'first').reset_index(drop = True) 
             #primer = primerpre.groupby(['Forward Primer', 'Reverse Primer']).first() 
             primer_cas13=primer_cas13.sort_values(by='Amplicon (bp)').head(5)
             primer_cas13.insert(0, 'Pair No.', range(1, 1 + len(primer_cas13)))
             print (primer_cas13)
             primer1_fw_cas13 = str( primer_cas13.iloc[0]['primer_locus']) + ":" + str(  primer_cas13.iloc[0]['primer_locus'] + primer_cas13.iloc[0]['Length']  -1    )
             primer1_rw_cas13 = str( primer_cas13.iloc[0]['primer_locus1']) + ":" + str(  primer_cas13.iloc[0]['primer_locus1'] - primer_cas13.iloc[0]['Length'] -1     )
             primer2_fw_cas13 = str( primer_cas13.iloc[1]['primer_locus']) + ":" + str(  primer_cas13.iloc[1]['primer_locus'] + primer_cas13.iloc[1]['Length']   -1 )
             primer2_rw_cas13 = str( primer_cas13.iloc[1]['primer_locus1']) + ":" + str(  primer_cas13.iloc[1]['primer_locus1'] - primer_cas13.iloc[1]['Length'] -1     ) 
             coo_cas13=str(primer_cas13.iloc[0:2]['primer_locus'].min()-50) + ":" +str(primer_cas13.iloc[0:2]['primer_locus1'].max()+50)
             primer_cas13=primer_cas13.drop(['Start','End','primer_locus','primer_locus1','IVC Product (bp)'],axis=1)
             cas13['Off-Targets']="Click for Off-targets"
             cas13=cas13[['crRNA Sequence','Wt Sequence','crRNA Position',  'Strand','Variation','Cas system', 'Off-Targets','mut_loc_cr','ran_loc_cr']] 
             #print (cas13)
              
            chr_num=df1['Chr'][0]
            gene=df1['Gene'][0]
            rsID=df1['rsID'][0]
            var_pos=df1['Position'][0]-1
            return chr_num,gene,rsID, var_pos,primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13
            
@app.route('/snvs', methods=['GET', 'POST'])
def form():
    #global counter
    #counter += 1
    db=get_db()
    form = MyForm()
    #form = MyForm()
    if request.method == 'POST' :
        #rsid1=request.form['RSID']
        rsid1=form.RSID.data
        #return render_template()
        #return '<h1>RSID:{}'.format(rsid1)
        #print(rsid1)
        #POS=81737298

        cur=db.execute('select  jatayu,mut_loc_cr,ran_loc_cr, strand,CHROMOSOME,CHROM,POS,RS,REF,ALT,presgrna,prestart,preend,GENEINFO,CLNDN,oligo,oligo1,FREQ,Cas_system from ngr_nrg_snp_wd_freq where RS=?',[rsid1])
        #cur=db.execute('select  * from jatayu_dbsnps_updated_full where RS=?',[rsid1])

        results=cur.fetchall()
        #print (results)
        if not results:
            flash ('Something went wrong!','error')
            flash ('The entered rsID is not valid or cannot be targeted by CRISPRDx','info')

            return redirect(url_for('form'))
        else:
            df = pd.DataFrame(results, columns =[  'jatayu',   'mut_loc_cr' ,  'ran_loc_cr' ,'Strand','Chr','chrom','Position','rsID','Ref','Alt','presgrna','prestart','preend','Gene','Disease','oligo','oligo1','FREQ','Cas_system']) 
            #print (df)
            #df['jatayu'],df['jatayur']=np.where(df['Strand']=='rev'  ,(df['jatayur'],df['jatayu']),(df['jatayu'],df['jatayur']))
            df["Variation"] = df["Ref"].astype(str) + ">" +df["Alt"].astype(str)
            df1=df
            #print (df1)
            alt_list=df1.Alt.unique()
            #print (alt_list)
            figfile = BytesIO() ####for plotting variation frequency####
            list1=Convert(df1['FREQ'][0]) 
            df_freq=pd.DataFrame(list1)
            df_freq.columns=['freq']
            df_freq['freq'] = df_freq['freq'].str.replace(r':', ',')
            
            #df_freq = df_freq.replace(dict.fromkeys(['.'],'0')) 
            print("df_Freq:",df_freq)


            df_freq1=df_freq.freq.apply(lambda x: pd.Series(str(x).split(",")))
            if len(df_freq1.columns) ==2:
             df_freq1.columns = ['Dataset',"Variant"]
             df_freq1.set_index('Dataset',inplace=True)
             df_freq1.rename(columns={'Variant':alt_list[0] + " Allele" }  , inplace=True) 
            elif len(df_freq1.columns) ==3:
             df_freq1.columns = ['Dataset',"Variant1","Variant2"]
             df_freq1.set_index('Dataset',inplace=True)
             df_freq1.rename(columns={'Variant1':alt_list[0] + " Allele",'Variant2':alt_list[1] + " Allele" }  , inplace=True) 
             
            elif len(df_freq1.columns) ==4:
             df_freq1.columns = ['Dataset',"Variant1","Variant2","Variant3"]
             df_freq1.rename(columns={'Variant1':alt_list[0] + " Allele",'Variant2':alt_list[1] + " Allele",'Variant3':alt_list[2] + " Allele" }  , inplace=True)
             df_freq1.set_index('Dataset',inplace=True)
             
            df_freq1 = df_freq1.replace(dict.fromkeys(['.'],'0'))
            print("df_Freq1:",df_freq1)

            df_freq1 = df_freq1.astype(float)
            pal = ["#00BFC4"  ,"#F8766D", "#00BA38", "#619CFF"]
             
            #plt.figure(figsize=(15,20))
            df_freq1.plot(kind="bar",color=pal,width=0.3,linewidth=0.4 ,  stacked=True ,edgecolor="black" ,figsize=(10, 5))
            
            plt.title("Variant allele frequency in different datasets")
            plt.xlabel("")
            
            plt.xticks(rotation=30,ha='right')
            plt.ylabel("Allele Frequency")
            #plt.tight_layout()
            ax = plt.subplot(111)
            #ax.xticks( rotation=90, ha='right')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
            ax.legend(loc='upper center', bbox_to_anchor=(0.9, 1.175), fancybox=True, shadow=True, ncol=3)
            #plt.legend( bbox_to_anchor=(0.5, 0., 0.5, 0.5),loc='best', borderaxespad=0.)
            
            plt.savefig(figfile, format='png')
            figfile.seek(0)
            plot_url = base64.b64encode(figfile.getvalue()).decode('utf8')

            chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=datafetch(df1,'Human')
            chrom=df1['chrom'][0]
            disease=df1['Disease'][0]
            disease=disease.replace('_',' ')
            variation=df1["Ref"][0]+"/"+"/".join(alt_list)
            alleles =df1["Ref"][0]+" > " +"/".join(alt_list)
            comments1=rsID + " " +variation

            
            print(var_pos)
            if form.validate_on_submit():
             return render_template('jatayu2_results.html',chr=chr_num,gene=gene,disease= disease,  chrom=chrom,comments1=comments1,rsid=rsID,allele=alleles, primer1_fw=primer1_fw_fn ,primer1_rw=primer1_rw_fn,primer2_fw=primer2_fw_fn,primer2_rw=primer2_rw_fn,var_pos=var_pos ,coo=coo,coo1=coo1_fn,coo1_cas12b=coo1_cas12b,  plot=plot_url ,data1=Fn,cas12b=cas12b,cas12a=cas12a,cas14=cas14,cas14ss =cas14ss,cas13=cas13, tables=[primer_fn.to_html(classes='primer',index=None),primer_cas12b.to_html(classes='primer',index=None),primer_cas12a.to_html(classes='primer',index=None),primer_cas14.to_html(classes='primer',index=None),primer_cas14ss.to_html(classes='primer',index=None),primer_cas13.to_html(classes='primer',index=None)  ],
          titles = [ 'na' ,'Primers','Primers',  'Primers'  , 'Primers' ,'Primers','Primers'] )

        #return '<h1>chromosome {}.pos {}.rs {}.ref{}.</h1>'.format(results[0]['CHROMOSOME'],results[0]['POS'],results[0]['RS'],results[0]['REF'])

    #return '<h1>RSID:{}'.format(rsid1)
    #return redirect(url_for('form'))
    return render_template('index.html',form=form)#,visitors=counter)



@app.route('/covid_snvs', methods=['GET', 'POST'])
def form1():
    #global counter
    #counter += 1
    db=get_db()
    form = MyForm()
    #form = MyForm()
    if request.method == 'POST' :
        #rsid1=request.form['RSID']
        rsid1=form.RSID.data
        #return render_template()
        #return '<h1>RSID:{}'.format(rsid1)
        #print(rsid1)

        cur=db.execute('select  jatayu,mut_loc_cr,ran_loc_cr, strand,CHROMOSOME,POS,RS,REF,ALT,presgrna,prestart,preend,Gene,oligo,oligo1,Virus_num,Annotation_Type,Gene_Pos_Codons,Cas_system from ngr_nrg_snp_wd_freq_covid where RS=?',[rsid1])
        #cur=db.execute('select  * from jatayu_dbsnps_updated_full where RS=?',[rsid1])

        results=cur.fetchall()
        #print (results)
        if not results:
            flash ('Something went wrong!','error')
            flash ('The entered SARS-CoV-2 mutations is not valid or cannot be targeted by CRISPRDx','info')

            return redirect(url_for('form1'))
        else:
            df = pd.DataFrame(results, columns =[  'jatayu',   'mut_loc_cr' ,  'ran_loc_cr' ,'Strand','Chr','Position','rsID','Ref','Alt','presgrna','prestart','preend','Gene','oligo','oligo1','FREQ','Annotation_Type','Gene_Pos_Codons','Cas_system']) 
            df["Variation"] = df["Ref"].astype(str) + ">" +df["Alt"].astype(str)
            df1=df
            #print (df1)
            alt_list=df1.Alt.unique()
            #print (alt_list)
            figfile = BytesIO() ####for plotting variation frequency####
            print (df1['FREQ']) 
            #list1=Convert(df1['FREQ'][0]) 
            df_freq=df1[['FREQ','Alt']]

            df_freq.columns=['freq','alt']
            df_freq.sort_values("alt", inplace = True)

            df_freq.drop_duplicates(keep = 'first', inplace = True)
 
            print("df_Freq:",df_freq)

            #df_freq['freq'] = df_freq['freq'].str.replace(r':', ',')
            
            #df_freq = df_freq.replace(dict.fromkeys(['.'],'0')) 
            print(type(df_freq))

            df_freq1=df_freq
            #if len(df_freq1.columns) ==2:
            # df_freq1.columns = ['Dataset',"Variant"]
            # df_freq1.set_index('Dataset',inplace=True)
            # df_freq1.rename(columns={'Variant':alt_list[0] + " Allele" }  , inplace=True) 
            """ 
            elif len(df_freq1.columns) ==3:
             df_freq1.columns = ['Dataset',"Variant1","Variant2"]
             df_freq1.set_index('Dataset',inplace=True)
             df_freq1.rename(columns={'Variant1':alt_list[0] + " Allele",'Variant2':alt_list[1] + " Allele" }  , inplace=True) 
             
            elif len(df_freq1.columns) ==4:
             df_freq1.columns = ['Dataset',"Variant1","Variant2","Variant3"]
             df_freq1.rename(columns={'Variant1':alt_list[0] + " Allele",'Variant2':alt_list[1] + " Allele",'Variant3':alt_list[2] + " Allele" }  , inplace=True)
             df_freq1.set_index('Dataset',inplace=True)
            """ 
            print("df_Freq1:",df_freq1)

            #df_freq1 = df_freq1.astype(float)
            pal = ["#00BFC4"  ,"#F8766D", "#00BA38", "#619CFF"]
             
            #plt.figure(figsize=(15,20))
            df_freq1.plot(kind="bar",color=pal,width=0.2,linewidth=0.4 ,x='alt',y='freq',legend=False, edgecolor="black" ,figsize=(8, 4))
            
            plt.title("No. of Virus genomes with mutation (Total 3567389 Sequences)")
            plt.xlabel("")
            
            plt.xticks(rotation=0,ha='right')
            plt.ylabel("No. of genomes")
            #plt.tight_layout()
            ax = plt.subplot(111)
            #ax.xticks( rotation=90, ha='right')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
            #ax.legend(loc='upper center', bbox_to_anchor=(0.9, 1.175), fancybox=True, shadow=True, ncol=3)
            #plt.legend( bbox_to_anchor=(0.5, 0., 0.5, 0.5),loc='best', borderaxespad=0.)
            
            plt.savefig(figfile, format='png')
            figfile.seek(0)
            plot_url = base64.b64encode(figfile.getvalue()).decode('utf8')
            chr_num,gene,rsID,var_pos, primer1_fw_fn , primer1_rw_fn,primer2_fw_fn,primer2_rw_fn ,coo,coo1_fn,coo1_cas12b, Fn,cas12b,cas12a,cas14,cas14ss,cas13, primer_fn,primer_cas12b,primer_cas12a,primer_cas14,primer_cas14ss,primer_cas13=datafetch(df1,'Sars-CoV-2')

            variation=df1["Ref"][0]+"/"+"/".join(alt_list)
            alleles =df1["Ref"][0]+" > " +"/".join(alt_list)
            var_pos=var_pos+1
            comments1=rsID + " " +variation
                       
            if form.validate_on_submit():
             return render_template('jatayu2_results_covid.html',chr=chr_num,gene=gene, comments1=comments1,rsid=rsID,allele=alleles, primer1_fw=primer1_fw_fn ,primer1_rw=primer1_rw_fn,primer2_fw=primer2_fw_fn,primer2_rw=primer2_rw_fn,var_pos=var_pos ,coo=coo,coo1=coo1_fn,coo1_cas12b=coo1_cas12b,  plot=plot_url ,data1=Fn,cas12b=cas12b,cas12a=cas12a,cas14=cas14,cas14ss =cas14ss,cas13=cas13, tables=[primer_fn.to_html(classes='primer',index=None),primer_cas12b.to_html(classes='primer',index=None),primer_cas12a.to_html(classes='primer',index=None),primer_cas14.to_html(classes='primer',index=None),primer_cas14ss.to_html(classes='primer',index=None),primer_cas13.to_html(classes='primer',index=None)  ],
          titles = [ 'na' ,'Primers','Primers',  'Primers'  , 'Primers' ,'Primers','Primers'] )

        #return '<h1>chromosome {}.pos {}.rs {}.ref{}.</h1>'.format(results[0]['CHROMOSOME'],results[0]['POS'],results[0]['RS'],results[0]['REF'])

    #return '<h1>RSID:{}'.format(rsid1)
    #return redirect(url_for('form'))
    return render_template('covid_snv.html',form=form)#,visitors=counter)



"""
@app.route('/viewresules', methods=['GET', 'POST'])
def viewresules():
    db=get_db()
    if request.method == 'POST' :
        rsid1=request.form['RSID']
    cur=db.execute('select  CHROMOSOME,POS,RS,REF from jatayu_dbsnps where RS=?',[rsid1])
    results=cur.fetchall()
    return '<h1>chromosome {}.pos {}.rs {}.ref{}.</h1>'.format(results[0]['CHROMOSOME'],results[0]['POS'],results[0]['RS'],results[0]['REF'])

"""

def desg_primer_cas(oligo,prestart,sgrna_len,list1):
    #def desg_primer_wd_pam(oligo):
    #primer_list=[]
    ls=[]
    #cutpos1 =cutpos -18
    #print (len(oligo))
    #list1=[[300,500],[500,700], [700,900], [900,1100],[1100,1300],[1300,1400]   ]
    #list1=[ [100,200] [200,300],  [300,400]   ]
    for rg in list1:
     if (len(oligo) > rg[1] or len(oligo) > rg[0])  :    
        
      primer =primer3.bindings.designPrimers(
     {
        'SEQUENCE_ID': 'random',
        'SEQUENCE_TEMPLATE': oligo,
        'SEQUENCE_INCLUDED_REGION': [0,len(oligo)]
     },
     {
        'PRIMER_OPT_SIZE': 24,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 75.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_LEFT_NUM_RETURNED':5,
        'PRIMER_RIGHT_NUM_RETURNED':5,
        'PRIMER_PRODUCT_SIZE_RANGE': rg
     })
     
    
     ##print (primer)
     #cutpos=418#418#268  #318 #168
     #
     #ls1=[]
     ##print (len(primer_list))
     #for primer in primer_list:
     ##print (primer)
      for i in range(len(list1)):
        #print (oligo)
        ##print (len(primer))
        ##print (primer)
        ##print (type (primer['PRIMER_LEFT_{0}'.format(i)][0] ))
        if primer['PRIMER_PAIR_NUM_RETURNED'] ==0:
         continue   
         #ls.extend(['NA']*12)
        elif 'PRIMER_LEFT_{0}'.format(i) not in primer:
         continue   
         #ls.extend(['NA']*12)   
        else:
         #print (primer['PRIMER_RIGHT_{0}'.format(i)][0])   
         ivc1 = 400 -  (int (primer['PRIMER_LEFT_{0}'.format(i)][0]  )-1 )
         ivc2 = (int (primer['PRIMER_RIGHT_{0}'.format(i)][0])) - (400+sgrna_len)
         ##print (ivc1,ivc2)
         #if i
         #ivc1=ivc1+1
         #ivc2=ivc2+1   
         #ivcdiff1= (ivc2/ivc1)
         #ivcdiff2= (ivc1/ivc2)
         ##print (ivc)
         if  ivc1 <= 0 or ivc2 <= 0 :
             continue
             #ls.extend(['NA'] * 12)
         
          
         #elif ((ivc2/ivc1 >= 1.7 or ivc1/ivc2 >= 1.7 ) and (ivc1 >=90 and ivc2 >=90) ):
         elif (ivc1 >=10 and ivc2 >=sgrna_len ):    
          #print (ivc1,ivc2)
          data= {
                'Forward Primer': primer['PRIMER_LEFT_{0}_SEQUENCE'.format(i)] ,
                'Start': primer['PRIMER_LEFT_{0}'.format(i)][0] ,
                 'primer_locus' : (prestart - 400  + (primer['PRIMER_LEFT_{0}'.format(i)][0])    ),
                 'Length': primer['PRIMER_LEFT_{0}'.format(i)][1] ,
                  'Tm (Cel.)' : int(primer['PRIMER_LEFT_{0}_TM'.format(i)]),
                  'GC%': int (primer['PRIMER_LEFT_{0}_GC_PERCENT'.format(i)]),
                     'Self-Com. (TH)' : int (primer['PRIMER_LEFT_{0}_SELF_ANY_TH'.format(i)]),
                    'Reverse Primer': primer['PRIMER_RIGHT_{0}_SEQUENCE'.format(i)],
                                    'End': primer['PRIMER_RIGHT_{0}'.format(i)][0] ,
                 'primer_locus1' : (prestart - 400  + (primer['PRIMER_RIGHT_{0}'.format(i)][0]) +2   ),

                 'Length ': primer['PRIMER_RIGHT_{0}'.format(i)][1],
                 'Tm (Cel.) ' : int(primer['PRIMER_RIGHT_{0}_TM'.format(i)]) ,
                 'GC% ': int (primer['PRIMER_RIGHT_{0}_GC_PERCENT'.format(i)]),
                 'Self-Com. (TH) ' : int (primer['PRIMER_RIGHT_{0}_SELF_ANY_TH'.format(i)]),
                     'Amplicon (bp)'         : primer['PRIMER_PAIR_{0}_PRODUCT_SIZE'.format(i)],
            'IVC Product (bp)':    ( (400 -  (int (primer['PRIMER_LEFT_{0}'.format(i)][0]  )-1 )) , ( int (primer['PRIMER_RIGHT_{0}'.format(i)][0]) -400   ) ) 

               }  
          ls.append(data)     
    
    primer_df = pd.DataFrame(ls) 
    ##print (primer_df)
    return primer_df   
def desg_primer_wd_pam(oligo,prestart,cutpos,list1):
    #def desg_primer_wd_pam(oligo):
    #primer_list=[]
    ls=[]
    cutpos1 =cutpos -18
    #list1=[[300,500],[500,700], [700,900], [900,1100],[1100,1300],[1300,1400]   ]
    #print(list1[-1][1],len(oligo))

    for rg in list1:
     if (len(oligo) > rg[1] or len(oligo) > rg[0])  :
        print(list1[-1][1],len(oligo))     
        primer =primer3.bindings.designPrimers(
     {
        'SEQUENCE_ID': 'random',
        'SEQUENCE_TEMPLATE': oligo,
        'SEQUENCE_INCLUDED_REGION': [0,len(oligo)]
     },
     {
        'PRIMER_OPT_SIZE': 24,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 75.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_LEFT_NUM_RETURNED':5,
        'PRIMER_RIGHT_NUM_RETURNED':5,
        'PRIMER_PRODUCT_SIZE_RANGE': rg
     })
     
    
     ##print (primer)
     #cutpos=418#418#268  #318 #168
     #
     #ls1=[]
     ##print (len(primer_list))
     #for primer in primer_list:
     ##print (primer)
        for i in range(len(list1)):
        ##print (oligo)
        ##print (lenprimer))
        ##print (primer)
        ##print (type (primer['PRIMER_LEFT_{0}'.format(i)][0] ))
         if primer['PRIMER_PAIR_NUM_RETURNED'] ==0:
          continue   
         #ls.extend(['NA']*12)
         elif 'PRIMER_LEFT_{0}'.format(i) not in primer:
          continue   
         #ls.extend(['NA']*12)   
         else:
          ivc1 = cutpos -  (int (primer['PRIMER_LEFT_{0}'.format(i)][0]  )-1 )
          ivc2 = int (primer['PRIMER_RIGHT_{0}'.format(i)][0]) -cutpos   
         ##print (ivc1,ivc2)
         #if i
         #ivc1=ivc1+1
         #ivc2=ivc2+1   
         #ivcdiff1= (ivc2/ivc1)
         #ivcdiff2= (ivc1/ivc2)
         ##print (ivc)
          if  ivc1 <= 0 or ivc2 <= 0 :
             continue
             #ls.extend(['NA'] * 12)
         
          
          elif ((ivc2/ivc1 >= 1.7 or ivc1/ivc2 >= 1.7 ) and (ivc1 >=90 and ivc2 >=90) ):
           print(primer['PRIMER_LEFT_{0}_SEQUENCE'.format(i)])   
          #print (ivc1,ivc2)
           data= {
                'Forward Primer': primer['PRIMER_LEFT_{0}_SEQUENCE'.format(i)] ,
                'Start': primer['PRIMER_LEFT_{0}'.format(i)][0] ,
                 'primer_locus' : (prestart - cutpos1  + (primer['PRIMER_LEFT_{0}'.format(i)][0])    ),
                 'Length': primer['PRIMER_LEFT_{0}'.format(i)][1] ,
                  'Tm (Cel.)' : int(primer['PRIMER_LEFT_{0}_TM'.format(i)]),
                  'GC%': int (primer['PRIMER_LEFT_{0}_GC_PERCENT'.format(i)]),
                     'Self-Com. (TH)' : int (primer['PRIMER_LEFT_{0}_SELF_ANY_TH'.format(i)]),
                    'Reverse Primer': primer['PRIMER_RIGHT_{0}_SEQUENCE'.format(i)],
                                    'End': primer['PRIMER_RIGHT_{0}'.format(i)][0] ,
                 'primer_locus1' : (prestart - cutpos1  + (primer['PRIMER_RIGHT_{0}'.format(i)][0]) +2   ),

                 'Length ': primer['PRIMER_RIGHT_{0}'.format(i)][1],
                 'Tm (Cel.) ' : int(primer['PRIMER_RIGHT_{0}_TM'.format(i)]) ,
                 'GC% ': int (primer['PRIMER_RIGHT_{0}_GC_PERCENT'.format(i)]),
                 'Self-Com. (TH) ' : int (primer['PRIMER_RIGHT_{0}_SELF_ANY_TH'.format(i)]),
                     'Amplicon (bp)'         : primer['PRIMER_PAIR_{0}_PRODUCT_SIZE'.format(i)],
            'IVC Product (bp)':    ( (cutpos -  (int (primer['PRIMER_LEFT_{0}'.format(i)][0]  )-1 )) , ( int (primer['PRIMER_RIGHT_{0}'.format(i)][0]) -cutpos   ) ) 

               }  
           ls.append(data)     
    
    primer_df = pd.DataFrame(ls) 
    ##print (primer_df)
    return primer_df     

    


    

if __name__ == '__main__':
    app.run(debug=True)
