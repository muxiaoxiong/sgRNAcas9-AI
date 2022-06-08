#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2022-02-21 11:13:33
# @Author  : Muxiaoxiong 
# @Email   : xiongweinie@foxmail.com

# sgRNAcas9-AI Online Tool
# sgRNACas9-AI Offline Tool Version V2.2.1
# http://123.57.239.141:8080/home

import logging
import argparse
import datetime
import sys
import os
import regex as re
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from tensorflow.keras.models import load_model
import numpy as np
import pandas as pd


__version__ = "2.2.1"

logging.basicConfig(
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
error   = logger.critical
warn    = logger.warning
debug   = logger.debug
info    = logger.info


dd={'A':'T','G':'C','T':'A','C':'G'}
cas_name=['eSpCas9','evoCas9','HypaCas9','Sniper-Cas9','SpCas9','SpCas9-HF1','SpCas9-NG','VRQR','xCas9']

def get_logo():
    return (r'''
 ########################### WELCOME  ##############################
#                                                                   #
#                          sgRNACas9-AI                             #
# ---a program for prediction of sgRNA activity using deep learning #
#                                                                   #
# Homepage: http://123.57.239.141:8080/home                         #
# Huazhong Agricultural University                                  #
# Offline Tool Version V2.2.1                                       #
# LAST REVISED: 2022.06.06                                          #
# Support:                                                          #
#   - SpCas9                                                        #
#   - evoCas9(1:1)                                                  #
#   - HypaCas9                                                      #
#   - Sniper-Cas9                                                   #
#   - eSpCas9                                                       #
#   - SpCas9-HF1                                                    #
#   - SpCas9-NG                                                     #
#   - VRQR                                                          #
#   - xCas9                                                         #
 ###################################################################
''')

def get_header():
    """
    Creates the header string with the header_str
    """
    term_width = 65
    logo = get_logo()
    description_str = logo+ "\n"
    description_str += ('[sgRNAcas9-AI version ' + __version__ + ']').center(term_width)+'\n'
    return description_str

def check_src():
    if not os.path.exists('./bin/crisflash'):
        if os.path.exists('src') and os.path.exists('Makefile'):
            os.system('make')
            return True
        else:
            warn('src folder or Makefile is missing')
            sys.exit(1)
    else:
        return True

def check_model():
    if os.path.exists('model'):
        return True
    else:
        warn('The model file is missing')
        sys.exit(1)

def load_gene(file):
    info('Load %s'%file)
    gene_dic={}
    with open(file) as ff:
        for line in ff:
            line=line.strip()
            if line.startswith('>'):
                name=line.replace('>','')
                gene_dic[name]=[]
            else:
                line=line.upper()
                gene_dic[name].append(line)
    return gene_dic


def Fasta_reverse(sequence):
    #sequence reverse complement
    sequence=sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.upper()
    return sequence[::-1]

def encode(lines):
    data=np.zeros((len(lines),23,4),dtype=int)
    for index,seq in enumerate(lines):
        for i,j in enumerate(seq):
            if j =='A':
                data[index,i,0]=1
            elif j=='T':
                data[index,i,1]=1
            elif j=='C':
                data[index,i,2]=1
            elif j=='G':
                data[index,i,3]=1
    return data

def count_mismatch(mismatch,file):
    """
    count sgRNA N20 mismatch number
    ruturn dict
    """
    if not os.path.exists('Temp.txt'):
        #If there is no temporary file, the reason may be that the input target gene and reference genome sequences are not in standard FASTA format
        error('File error, please check input file')
        exit(1)
    info('count %s off-target'%file)
    seqdic={}
    mismatch=int(mismatch)
    with open('Temp.txt') as ff:
        for line in ff:
            line=line.strip()
            test=line.split()
            seq=test[0][:20]
            num=int(test[-1])
            if seqdic.get(seq):
                seqdic[seq][num]+=1
            else:
                seqdic[seq]=[0 for i in range(mismatch+1)]
    os.remove('Temp.txt')
    return seqdic
############################################################################
############################################################################

def cal_off(genefile,genomefile,pam,mismatch,process):
    """
    run crisflash
    """
    info('cal %s off-target'%genefile)
    if os.path.exists('Temp.txt'):
        os.remove('Temp.txt')
    #update file

    os.system('./bin/crisflash -g %s -s %s -o Temp.txt -m %s -p %s -C -t %s >log.txt'%(genomefile,genefile,mismatch,pam,process))
    os.remove('log.txt')

def run(args):
    """
    seqdic={N20:[0,0,0]}
    genedic={'gene':[v1,v2]}
    """
    seqdic=count_mismatch(args.mismatch,args.input)
    genedic=load_gene(args.input)
    #write title
    out=open(args.output,'a+')
    out.write('sgR_ID,sgR_seq+PAM,sgR_seq,PAM,strand,start,end,GC%,')
    for i in range(int(args.mismatch)+1):
        out.write('%sM,'%(i))
    out.write('Total,')
    if  args.active !='':
        out.write('4Ts_motif,sgR_efficiency\n')
    else:
        out.write('4Ts_motif\n')

    pam=args.pam
    if pam[0]=='N':
        flag1=pam[1:]
        flag2=[dd[i] for i in flag1]
        flag2=''.join(flag2)
        pattern1 = re.compile('.{21}%s'%flag1)
        pattern2 = re.compile('%s.{21}'%flag2)
    else:
        flag1=pam
        flag2=[dd[i] for i in flag1]
        flag2=''.join(flag2)
        pattern1 = re.compile('.{20}%s'%flag1)
        pattern2 = re.compile('%s.{20}'%flag2)

    info('output %s off-target'%args.input)
    for key,value in genedic.items():
        value=''.join(value)
        sgRNA_L=pattern1.finditer(value,overlapped=True)
        sgRNA_F=pattern2.finditer(value,overlapped=True)
        count=0
        for sgRNA in sgRNA_L:
            count+=1
            name=key+'_s_'+str(count)
            start=sgRNA.start()
            end=sgRNA.end()
            seq=sgRNA.group()
            N20=seq[:20]
            _pam=seq[20:]
            gcgc ='%.2f' % ((N20.count('C')+N20.count('G'))/20)
            if seqdic.get(N20):
                off_list=seqdic[N20]
                total=sum(off_list)
                out.write('%s,%s,%s,%s,%s,%s,%s,%s,'%(name,seq,N20,_pam,'+',start,end,gcgc))
                for i in off_list:
                    out.write(str(i)+',')
                out.write(str(total)+',')
                if 'TTTT' in N20:
                    out.write('TTTT'+'\n')
                else:
                    out.write('-'+'\n')
            else:
                out.write('%s,%s,%s,%s,%s,%s,%s,%s,'%(name,seq,N20,_pam,'+',start,end,gcgc))
                for i in range(int(args.mismatch)+2):
                    out.write('0,')
                if 'TTTT' in N20:
                    out.write('TTTT'+'\n')
                else:
                    out.write('-'+'\n')
        count=0
        for sgRNA in sgRNA_F:
            count+=1
            name=key+'_a_'+str(count)
            start=sgRNA.start()
            end=sgRNA.end()
            seq=Fasta_reverse(sgRNA.group())
            N20=seq[:20]
            _pam=seq[20:]
            gcgc ='%.2f' % ((N20.count('C')+N20.count('G'))/20)
            if seqdic.get(N20):
                off_list=seqdic[N20]
                total=sum(off_list)
                out.write('%s,%s,%s,%s,%s,%s,%s,%s,'%(name,seq,N20,_pam,'-',start,end,gcgc))
                for i in off_list:
                    out.write(str(i)+',')
                out.write(str(total)+',')
                if 'TTTT' in N20:
                    out.write('TTTT'+'\n')
                else:
                    out.write('-'+'\n')
            else:
                out.write('%s,%s,%s,%s,%s,%s,%s,%s,'%(name,seq,N20,_pam,'-',start,end,gcgc))
                for i in range(int(args.mismatch)+2):
                    out.write('0,')
                if 'TTTT' in N20:
                    out.write('TTTT'+'\n')
                else:
                    out.write('-'+'\n')
    out.close()
    if  args.active !='':
        model=load_model('./model/%s.h5'%(args.active))
        data=pd.read_csv(args.output)
        seq=encode(data['sgR_seq+PAM'])
        prediction=model.predict(seq)
        data['sgR_efficiency']=prediction
        data.to_csv(args.output,index=False)
    info('Done')

############################################################################

def main():
    print(get_header())
    parser = argparse.ArgumentParser(description='sgRNAcas9-AI Parameters')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    parser.add_argument('-i', '--input', type=str,  help='gene fast file', default='')
    parser.add_argument('-g', '--genome', type=str,  help='genome fast file', default='')
    parser.add_argument('-a', '--active', type=str,  help='Cas9 Type (default: SpCas9)', default='SpCas9')
    parser.add_argument('-m', '--mismatch', type=str,  help='mismatch of sgRNA off-target prediction(default: 3)', default='3')
    parser.add_argument('-p', '--pam', type=str, help='Specify the pam to use for off-target analysis(default: NGG).', default='NGG')
    parser.add_argument('-t', '--thread', type=str, help='Specify the number of processes to use for analysis.(default: 1)', default='1')
    parser.add_argument('-o', '--output',  help='Output to use for the analysis (default: result.csv)', default='result.csv')
    args = parser.parse_args()
    sys.stdout.flush()

    #check parse_args
    if not args.genome or not args.input:
        parser.print_help()
        exit(1)

    #check model
    if args.active !='' and args.active not in cas_name:
        caninfo=''
        for i in cas_name:
            caninfo +='   - '+i+'\n'
        error('Model name :%s unsupport\n Support:\n %s'%(args.active,caninfo))
        exit(1)

    #check outfile
    if os.path.exists(args.output):
        warn('%s exists,removing'%args.output)
        os.remove(args.output)

    # check dedependencies
    info('Checking dependencies...')
    if check_src() and check_model():
        info('All the required dependencies are present!')

    #All check is ok
    #cal off-target
    cal_off(args.input,args.genome,args.pam,args.mismatch,args.thread)
    # #run
    run(args)


if __name__ == '__main__':
    main()
