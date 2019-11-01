# -*- coding: utf-8 -*- 
'''
Created on 9 August 2018

@author: Nicolas Villeriot
'''

import BioSeqTools
import json
import re
import subprocess

#-------------------------------Function library------------------------------------------

'''
Function allowing to read a .json file and creating a formatted .csv file
@param jsonFileName: Absolute path to the result.json file obtained with MetaFinder
@param outPath : Absolute path where will be registered all results
@type jsonFileName, outPath: String
@return : Hash, Summary of the result.json file
'''
def JsonSummary(jsonFileName, outPath):
    
    results = {} 
    ##Load json file into a tree
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    ##reinitialize the output file if it exists and write the header
    with open(outPath + "CRISPRCasSummary.csv", 'w') as result :
        result.write("READ\tSpacers\tEvidence Level\tDR Length\tDR Consensus\tDirection\tPotential Orientation\tCas\tCRISPR Start\tCRISPR End\tLength\tDr Group\n")
    ##
    
    ##Create a csv file which summarize the result.json file important data
    for contig in datas["Sequences"] :

            contigResult = []

            if contig["Crisprs"] != [] :
                for crispr in contig["Crisprs"]:
                    
                    if contig["Cas"] != []: 
                        contigResult.append(contig["Id"] + "\t" + str(crispr["Spacers"])+ "\t" + str(crispr["Evidence_Level"]) + "\t" + str(crispr["DR_Length"])+ "\t" + crispr["DR_Consensus"] + "\t" + crispr["CRISPRDirection"]+ "\t" + crispr["Potential_Orientation"] + "\tNa\t" + str(crispr["Start"]) + "\t" + str(crispr["End"])+ "\t" + str(crispr["End"] - crispr["Start"] + 1))
                    
                    crisprSeq = ""
                    for region in crispr["Regions"] :
                        if region["Sequence"] != "UNKNOWN" :
                            crisprSeq += region["Sequence"] + "\n"
                            
                    
            if contig["Cas"] != []:
                genes = ""
                
                for cas in contig["Cas"]:
                    
                    for subT in cas["Genes"] :
                        genes += subT["Sub_type"] + ","
                    genes = genes.rstrip(",")
                    genes += "|"
                
                genes = genes.rstrip("|")
                     
                if contigResult != [] :
                    listArg = []
                    for res in range(0,len(contigResult)) :
                        listArg = contigResult[res].split("\t")
                        listArg[7] = genes
                        contigResult[res] = "\t".join(listArg)
                else :
                    contigResult.append(contig["Id"] + "\t0\t0\t0\tNo DR\tNo Dir\tNo Ori\t" + genes + "\t0\t0\t0")
            
            if contigResult != [] :
                
                results[contig["Id"]] = []
                
                for res in range(0,len(contigResult)) :
                    results[contig["Id"]].append(contigResult[res])
    ##               
    return(results)
                
'''
Function allowing to create csv file from results of CRISPRCasMeta.
@param results : Data hash tree which summarizes the result.json file (output of JsonSummary function)
@param compRes : Data hash tree which contains dr groups (output of DRComparisonFast function)
@param outPath : Absolute path where will be registered all results
@type results, outPath : String
@type compRes : Hash
'''
def CsvParser(results,compRes,outPath): 
    nameOrder = []
    
    for name in results.keys() :
        nameOrder.append(name) 
    
    ##Write the csv file
    with open(outPath + "CRISPRCasSummary.tsv", 'a') as resFile :
            for name in nameOrder :
                gr = "dr_group_Na"
                if compRes != {} :
                    for drGrp in compRes.keys() :
                        if re.search(name,drGrp)!= None : 
                            #with open(outPath + "RepeatGroup.out", 'a') as grFile :
                            gr = "" 
                            #grFile.write(drGrp + "\n")
                            for dr in compRes[drGrp].keys() :
                                gr += dr + " " 
                                    #grFile.write(dr + "\n" + compRes[drGrp][dr] + "\n")
                                #grFile.write("\n")
                for crispr in results[name] :
                    resFile.write(crispr + "\t" + gr +"\n")
    ##
        
    print("Work finished ! Thank you !")


'''
Function allowing to create an human readable summary result file from result.json
@param jsonFileName : Path to the result.json file obtained with MetaFinder
@param outPath : Absolute path where will be registered all results
@param compRes : Data hash tree which contains dr groups (output of DRComparisonFast function)
@type jsonFileName, outPath : String
@type compRes : Hash
'''
def ResultSummary (jsonFileName, outPath, compRes): 
    
    jSum = JsonSummary(jsonFileName,outPath)
    CsvParser(jSum,compRes,outPath)

#-------------------------------      Main      ------------------------------------------
#compRes = ByDrComparisonLight("reLib", 1)
#compRes2 = 
##Test
#DbComparison ("", "", 1, "")
'''
compDr = DrComparisonFast("Compare/RepeatLib_B2F11_Ev3.fasta", 0)
with open("Compare/RepeatLib_B2F11_Ev3.out", "w") as resFile :
    for drGroup in compDr.keys():
        resFile.write(drGroup + "\n")
        for seq in compDr[drGroup].keys() :
            resFile.write(seq + "\n" + compDr[drGroup][seq] + "\n")
        
        resFile.write("\n")

CasClustersComparison("","","")
'''       

##
'''
jRes = glob.glob("Input/*.json")

jSum = {}

compDr = DrComparisonFast("Input/RepeatLib.fasta", 1)

i = 0
for jR in jRes :
    jSum = JsonSummary(jR, jSum,"Input/")
    
CsvParser(jSum,compDr,"Input/")'''
#CasRepeatComb(jSum,compDr,"Input/")
