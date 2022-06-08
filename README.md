## INTRODUCTION

The sgRNACas9-AI command-line software tool can be used to optimize the design of sgRNAs in CRISPR/Cas9 systems, including deep learning-based strategies for sgRNA activity prediction and assessment of potential off-target sites. In a nutshell, the features of sgRNACas9-AI:

  - sgRNA activity prediction based on deep learning strategy

  - Evaluation of sgRNA potential off-target sites based on sequence alignment strategy

## INSTALLATION

This tool can only be run on Linux systems, the steps to install the necessary dependencies for the program are as follows:

1.python3

  Linux[https](https://www.python.org/downloads/source/) ://www.python.org/downloads/source/

2.Install pip (pip should come with Python, but this should ensure it is present; see https://pip.pypa.io/en/stable/installation/)

3.Install package

```python
pip install regex
pip install numpy
pip install pandas
pip install tensorflow
```

4.Verify that sgRNACas9-AI is installed using the command:

```python
python sgRNACas9-AI.py --help
```

## USAGE INFORMATION

sgRNACas9-AI requires main two parameters:

  - input sequences in the form of fasta files (given by the -i or --input)

  - reference genome in the form of fasta files (given by the -g or --genome).

For example:

```python
python sgRNACas9-AI.py -i gene.fa -g genome.fa -a -m 3 -p NGG -t 8 -o result.txt
```

Parameter list：

```
-h or --help: show a help message and exit.

--version: show program version

-i or --input: input sequences in fasta files

-g or --genome: genome sequences in fasta files

-a or --active: predict sgRNA activity

-m or --mismatch: mismatch number of sgRNA off-target prediciton

-p or --pam: specify th pam to use for off-target analysis

-t or --thread: the number of processes to use for program

-o or --output: output file for result
```

Anyone can use the source codes, documents or the excutable file of sgRNACas9-AI free of charge 
for non-commercial use. For commercial use, please contact the author.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use this program in your research, please cite:

  1. Nie et al., sgRNAcas9-AI: a program for prediction of CRISPR/Cas9 and its variant sgRNA activity using deep learning, In prepare, 2022.

  2. Xie S, Shen B, Zhang C, Huang X, Zhang Y. sgRNAcas9: a software package for designing CRISPR sgRNA and evaluating potential off-target cleavage sites. PLoS One. 2014 Jun 23;9(6):e100448. doi: 10.1371/journal.pone.0100448.

Contact Us:
  College of Animal Sciences & Technology / College of Veterinary Medicine, Huazhong Agricultural University, Wuhan 430070, P. R. China.

  Please send bug reports to: ssxie@sibcb.ac.cn

  Other software resources: http://www.biootools.com

  Prof. Shengsong Xie

