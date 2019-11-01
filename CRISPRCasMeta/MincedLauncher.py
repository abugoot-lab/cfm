'''
Created on 29 mars 2018

@author: Nicolas Villeriot
'''

import subprocess
import re

#-------------------------------Function library------------------------------------------
'''
Function allowing search CRISPR using Minced on a multi fasta file
@param mincedPath: Absolute path of minced executable
@param fastaFilePath: Absolute path of the fasta or fasta.gz file
@param outputPath: Absolute path where will be registered all results
@param minRL : Minimal repeat length allowed  (CRISPRCasFinder parameter)
@param maxRL : Maximal repeat length allowed (CRISPRCasFinder parameter)
@param minSS : Minimal spacer length allowed
@param maxSS : Maximal spacer length allowed 
@type mincedPath, fastaFilePath, outputPath, minRL, maxRL : String
@type minSS, maxSS : Integer
@return: String, error message
'''
def CrisprSearch (mincedPath,fastaFilePath,outputPath,minRL,maxRL,minSS,maxSS) :
    
    error = ""
    
    print("\nMinced...")
    
    try :
        with open(outputPath + "Minced.out", 'w') as fi :
            fi.write("")
        with open(outputPath + "Minced.out", 'a') as fi :
            fi.write(subprocess.check_output([mincedPath, "-minNR", "2", "-minRL", str(minRL), "-maxRL", str(maxRL), "-minSL", str(minSS), "-maxSL", str(maxSS), "-spacers", fastaFilePath, outputPath + "MincedFiles/mincedRes"]))
    except :
        error = "Error 101 : Command bash to run Minced failed !\n Check file Output/Minced.out if it exists for more details."
    
    return error

'''
Function allowing to read a Minced result file and create a csv file with readable results
@param reportPath: Absolute path where will be registered the CRT(Minced) results
@param csvFile: Absolute path of the csv file which will be generated
@type reportPath, csvFile: String
@return : Hash, CRT(Minced) results of each contig in a String ; String Array, contigs by alphabetical order
'''
def MincedReport (reportPath,csvFile) :
    with open(csvFile, "w") as fi :
        fi.write("Contig;Spacers;DRLength;CRISPR start;CRISPR end;Length;Array\n") ##CRISPR array 27/06/2018
    
    contigs = {}
    currentCont = []
    nbSeq = 0
    nbCR = 0
    nbRe = 0
    CRISPRarray = "" 
    nameOrder = [] 
    
    with open(reportPath, "r") as rep :
        for line in rep :
            listArg = line.split(" ")
            
            ##Build the CRISPR array 
            listSeq = line.split("\t")
            if re.search('^[0-9]+$',listSeq[0]) != None :
                CRISPRarray += "R_" + listSeq[2].rstrip("\n") + " S_" + listSeq[3].rstrip("\n") + " "
            ##
            
            
            if listArg[0] == "Sequence" :
                ##New sequence, last results recording in contigs hash
                if currentCont != [] :
                    crispr = 0
                    for i in range(1,len(currentCont)):
                        if (i % 5) == 0 : 
                            crispr += 1
                            nameOrder.append(currentCont[0] + "_crispr_" + str(crispr)) #conserves order of the input file
                            
                            contigs[currentCont[0] + "_crispr_" + str(crispr)] = currentCont[i-2] + ";" + currentCont[i-1] + ";" + currentCont[i-4] + ";" + currentCont[i-3] + ";" + str(int(currentCont[i-3]) - int(currentCont[i-4]) + 1) + ";" + currentCont[i] #CRISPR array 
                            
                            #print(currentCont[0] + "_crispr_" + str(crispr))
                
                currentCont = []
                currentCont.append(listArg[1].split("'")[1])
                nbSeq += 1
                ##
                
            elif listArg[0] == "CRISPR" :
                currentCont.append(listArg[5])
                currentCont.append(listArg[7].rstrip())
                nbCR += 1
            
            elif listArg[0] == "Repeats:" :
                ##End of this CRISPR array
                currentCont.append(str(int(listArg[1].split("\t")[0]) - 1))
                currentCont.append(listArg[3].split("\t")[0])
                currentCont.append(CRISPRarray) #CRISPR array
                CRISPRarray = "" #CRISPR array
                nbRe += 1
                ##
        
        ##Last sequence result recording
        if currentCont != [] :
            crispr = 0
            for i in range(1,len(currentCont)):
                if (i % 5) == 0 : #CRISPR array
                    crispr += 1
                    nameOrder.append(currentCont[0] + "_crispr_" + str(crispr)) #conserves order of the input file
                    contigs[currentCont[0] + "_crispr_" + str(crispr)] = currentCont[i-2] + ";" + currentCont[i-1] + ";" + currentCont[i-4] + ";" + currentCont[i-3] + ";" + str(int(currentCont[i-3]) - int(currentCont[i-4]) + 1) + ";" + currentCont[i] #CRISPR array
                    #print(currentCont[0] + "_crispr_" + str(crispr))
        ##
           
    print("Sequences : " + str(nbSeq) + "\nCRISPRs : " + str(nbCR) + "\nRepeats : "  + str(nbRe))
    
    ##Write results in a csv file 
    with open(csvFile, "a") as fi :
        for key in nameOrder:
            fi.write(key + ";")
            for crispr in contigs[key] :
                fi.write(crispr)
            fi.write("\n")
    ##
    return (contigs,nameOrder)
            
#-------------------------------Main------------------------------------------
'''
error = CrisprSearch("/home/Nicolas/Test/minced-master/minced", "/home/Nicolas/Test/SequenceTest/NODE01_31386.fasta")

if error == "" :
    print("\nCrispr report creation...")
    MincedReport("Output/MincedFiles/mincedAll.out", "Output/CrisprResult.csv")
    print("\nCrispr search completed !")
'''