# -*- coding: utf-8 -*- 
'''
Created on 28 May 2018

@author: Nicolas Villeriot
'''

import BioSeqTools
import re
import MincedLauncher
import gzip
from decimal import Decimal

#-------------------------------Function library------------------------------------------
'''
Function allowing to create a new fasta file with only contig with potential CRISPR
@param fastaFilePath: Absolute path of the fasta or fasta.gz file
@param outputPath: Absolute path where will be registered all results
@param mincedPath: Absolute path of minced executable
@param fasta: Input file format (True if fasta, False if fasta.gz)
@param minRL : Minimal repeat length allowed  (CRISPRCasFinder parameter)
@param maxRL : Maximal repeat length allowed (CRISPRCasFinder parameter)
@type fastaFilePath, outputPath, mincedPath, minRL, maxRL: String
@type fasta: Boolean
'''
def MetaSummary (fastaFilePath,outputPath,mincedPath,fasta,minRL,maxRL):
    
    ##Search potential CRISPRs, using CRT Minced
    if not fasta : 
        fastaFilePath = UnzipGZ(fastaFilePath, outputPath) ##gz format treatment
        
    minSS = int(round(Decimal(minRL)*Decimal(1.09))) #minced parameters
    maxSS = int(round(Decimal(maxRL)*Decimal(1.09))) #minced parameters
    error = MincedLauncher.CrisprSearch(mincedPath, fastaFilePath, outputPath, minRL, maxRL, minSS, maxSS)
    
    if error == "" :
        print("\nCrispr report creation...")
        (crisprs, contigsCr) = MincedLauncher.MincedReport(outputPath + "/MincedFiles/mincedRes", outputPath + "CrisprResult.csv")
        print("\nCrispr search completed !")
    else :
        print(error)
    ##
    
    ##List Id of the contigs with potential CRISPRs
    contigs = []
    for contig in contigsCr :
        contig = "_".join(contig.split("_")[:-2])
        if not (contig in contigs):
            contigs.append(contig)
    ##
    
    ##Create a fasta file of contigs with potential CRISPRs
    with open(outputPath + "Summary.fasta", "w") as fi :
        fi.write("")
    
    with open(outputPath + "Summary.fasta", "a") as fi :
        
        x = 0 
        contig = contigs[x] 
        for name, seq in BioSeqTools.ReadFasta(fastaFilePath):
            name = name[1:]
            
            if re.search(contig, name) != None : 
                ##avoid problematic characters
                name = name.replace("|","_") 
                ##
                #print(name)
                fi.write(">" + name + "\n" + seq + "\n")
                x += 1 
                if x < len(contigs) :
                    contig = contigs[x]
                else :
                    break
    ##  

'''
Function allowing to unzip a gz file 11/07/2018
@param fastaGzFile: Name and directory of the fasta gz file used as reference
@param out: Absolute path of the output directory
@type fastaGzFile, out: String
@return: String, Absolute path of the decompressed fasta file
'''
def UnzipGZ (fastaGzFile,out):
    with open(out + "Complete.fasta","w") as fastF :
        fastF.write("")    
    
    with gzip.open(fastaGzFile, 'rb') as gzF:
        with open(out + "Complete.fasta","a") as fastF :
            for line in gzF :
                fastF.write(line)
    
    return(out + "Complete.fasta")

'''
Function allowing to generate a formatted database of repeat, uses nameForm parameter to format sequence name
@param repeatLib : Absolute path of a repeat library
@param nameForm : Reformat fasta repeat sequence header to recognize which sequence is from the database
@param out: Absolute path of the output directory
@param fasta: Input file format (True if fasta, False if fasta.gz)
@type repeatLib, nameForm, out : String
@type fasta : Boolean
'''
def RepeatDatabase (repeatLib,nameForm, out, fasta):
    if not fasta : 
        repeatLib = UnzipGZ(repeatLib, out) ##gz format treatment
        
    if repeatLib != "" :
        for name, seq in BioSeqTools.ReadFasta(repeatLib):
            name = nameForm + name[1:]
            with open(out + "RepeatLib.fasta", "a") as repLib:
                repLib.write(name + "\n" + seq + "\n")
            
        
#-------------------------------Main------------------------------------------
'''
start = time.time()
MetaSummary("/home/Nicolas/Test/SequenceTest/NODE01_31386.fasta", "Output/", "/home/Nicolas/Test/minced-master/minced")
end = time.time()
print("Execution time : " + str(end-start) + " s")
'''