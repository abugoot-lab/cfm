# -*- coding: utf-8 -*- 
'''
Created on 4 May 2018

This script can be useful when you want to do further analysis on CAS results of CRISPRCasFinder.
It allows to create some formatted files with the different CAS clusters of your sample.

@author: Nicolas Villeriot
'''

import BioSeqTools
import json
import re
from string import rstrip

#----------------------------Function library------------------------------#

'''
Function allowing to create a fasta file of clusters found in the rawCas.fna file from CRISPRCasFinder 4.2.11
@param rawCasFile: Absolute path of the rawCas.fna result
@param outDir : Absolute path where will be registered all results
@type rawCasFile, outDir: String
@return : Hash, Cas Type and the related strains or contigs
'''
def ClusterMaker(rawCasFile, outDir):
    
    casByGene = {} #Separated Cas
    
    with open(outDir + "CasClusters.fasta", "w") as clusterFile :
        clusterFile.write("")
    
    with open(outDir + "CasClusters.fasta", "a") as clusterFile :
        
        clustersCas = {}
        
        for name, seq in BioSeqTools.ReadFasta(rawCasFile):
            print(name)
            name = name.split("|")
            idHead = "_".join(name[1].split("_")[:-1])
            
            ##Separated Cas
            if not (name[3].split(" ")[0] in casByGene.keys()) :
                casByGene[name[3].split(" ")[0]] = []
            
            #casByGene[name[3].split(" ")[0]].append(name[0][1:] + "|" + seq) #genome version
            casByGene[name[3].split(" ")[0]].append(idHead + "|" + seq) #CRISPRCasMeta version 13/08/2018
            ##
            
            ##Build oriented Cas clusters with genes found in the rawCas.fna file
            if clustersCas != {} :
                if (idHead + name[2]) in clustersCas.keys() :
                    clustersCas[idHead + name[2]]["Header"] += " " + name[3].split(" ")[0]
                    clustersCas[idHead + name[2]]["Pos"].append(int(name[3].split(" ")[1].split(",")[0]))
                    clustersCas[idHead + name[2]]["Pos"].append(int(name[3].split(" ")[1].split(",")[1]))
                    clustersCas[idHead + name[2]]["Seq"] += seq
                    
                    
                else :
                    clustersCas[idHead + name[2]] = {}
                    clustersCas[idHead + name[2]]["Header"] = (">" + idHead + "|" + name[2] + "|" + name[3].split(" ")[0])
                    clustersCas[idHead + name[2]]["Pos"] = []
                    clustersCas[idHead + name[2]]["Pos"].append(int(name[3].split(" ")[1].split(",")[0]))
                    clustersCas[idHead + name[2]]["Pos"].append(int(name[3].split(" ")[1].split(",")[1]))
                    
                    if seq[:3] == "TTA" or seq[:3] == "CTA" or seq[:3] == "TCA" :
                        seq = BioSeqTools.ReverseComp(seq)
                        clustersCas[idHead + name[2]]["Dir"] = "rev comp"
                    else :
                        clustersCas[idHead + name[2]]["Dir"] = "forward"
                    
                    clustersCas[idHead + name[2]]["Seq"] = seq
                    
            else :
                clustersCas[idHead + name[2]] = {}
                clustersCas[idHead + name[2]]["Header"] = (">" + idHead + "|" + name[2] + "|" + name[3].split(" ")[0])
                clustersCas[idHead + name[2]]["Pos"] = []
                clustersCas[idHead + name[2]]["Pos"].append(int(name[3].split(" ")[1].split(",")[0]))
                clustersCas[idHead + name[2]]["Pos"].append(int(name[3].split(" ")[1].split(",")[1]))
                
                if seq[:3] == "TTA" or seq[:3] == "CTA" or seq[:3] == "TCA" :
                    seq = BioSeqTools.ReverseComp(seq)
                    clustersCas[idHead + name[2]]["Dir"] = "rev comp"
                else :
                    clustersCas[idHead + name[2]]["Dir"] = "forward"
                    
                clustersCas[idHead + name[2]]["Seq"] = seq
            ##
        
        ##Write oriented Cas cluster sequences in a fasta file
        for idH in clustersCas.keys() :
            newSeq = ""
            it = 0
            for nuc in clustersCas[idH]["Seq"] :
                it += 1
                newSeq += nuc
                if (it % 100) == 0 :
                    newSeq += "\n"
                    
            clusterFile.write(clustersCas[idH]["Header"] + "|" + str(min(clustersCas[idH]["Pos"])) + "_" + str(max(clustersCas[idH]["Pos"])) + "|" + clustersCas[idH]["Dir"] + "\n" + newSeq + "\n")            
        ##           
        
    return (casByGene)
        
'''
Function allowing to create specific fasta file for each type of cluster Cas found
@param outDir : Absolute path where will be registered all results
@type outDir: String
'''
def FileCasTyper(outDir):
    for name, seq in BioSeqTools.ReadFasta(outDir + "CasClusters.fasta"):
        name = name.split("|")
        
        newSeq = ""
        it = 0
        for nuc in seq :
            it += 1
            newSeq += nuc
            if (it % 100) == 0 :
                newSeq += "\n"
        
        with open(outDir + name[1] + ".fasta","a") as fi :
            fi.write("|".join(name) + "\n" + newSeq + "\n")


'''
Function allowing to create a tsv file from the CasClusters.fasta
@param outDir : Absolute path where will be registered all results
@param form : File format of the results (csv or tsv)
@param casTree : Cas Type and the related strains or contigs
@type outDir, form : String
@type casTree : Hash
'''
def CsvSummary(outDir,form,casTree):
    
    spacer = "\t"
    end = ".tsv"
    if form == "csv" :
        spacer = ";"
        end = ".csv"
    
    ##Separated Cas
    listCas = casTree.keys()
    ##
    
    ##Summary initialization
    with open(outDir + "Cas_Clusters_Summary"+ end,"w") as results :
        results.write("Id" + spacer +"Strain" + spacer +"Cas_Type" + spacer +"Location" + spacer + "Orientation" + spacer +"Cas_Cluster")
        for cas in listCas :
            results.write(spacer + cas)
        results.write("\n")
    ##
    
    with open(outDir + "Cas_Clusters_Summary"+ end,"a") as results :
        for name, seq in BioSeqTools.ReadFasta(outDir + "CasClusters.fasta"):
            name = name[1:].split("|")
            
            ##Strain Name
            stra = " "
            with open(outDir + "Metadata.tsv", "r") as meta :
                for line in meta :
                    if re.search(line.split("\t")[0],name[0]) != None :
                        stra = rstrip(line.split("\t")[1])
                        break       
            ##
            
            ##Separated Cas
            results.write(name[0] + spacer + stra + spacer + name[1] + spacer + name[3] + spacer + name[4] + spacer + seq)
            for cas in listCas :
                present = False
                
                if cas in name[2].split(" ") :
                    for strain in casTree[cas] :
                        if strain.split("|")[0] == name[0] :
                            results.write(spacer + strain.split("|")[1])
                            present = True
                            break
                
                if not present :
                    results.write(spacer + "None")
            
            results.write("\n")
            ##
            
    state = "Completed !"      
    return state

'''
Function allowing to take strain description from json file
@param outDir : Absolute path where will be registered all results
@param jsonFile : Absolute path to the result.json file obtained with MetaFinder
@type outDir, jsonFile : String
'''
def Lexique(outDir,jsonFile):
    with open(jsonFile, 'r') as fi :
        datas = json.load(fi)
    
    with open(outDir + "Metadata.tsv", "w") as lex :
        lex.write("Id\tStrain\n")
        
    with open(outDir + "Metadata.tsv", "a") as lex :
        for contig in datas["Sequences"] :
            acc = contig["Id"]
            desc = contig["Description"]
            
            if acc != "" and acc != None and desc != "" and desc != None and desc != "Unknown":
                lex.write(acc + "\t")
                
                if re.search("str\.",desc) != None :
                    desc = desc.split("str.")[1]
                elif re.search("strain",desc) != None :
                    desc = desc.split("strain")[1]
                elif len(desc.split(" ")) > 2 :
                    desc = " ".join(desc.split(" ")[2:])
                    
                if re.search(",",desc) != None :
                    lex.write(desc.split(",")[0] + "\n")
                else :
                    lex.write(desc + "\n")
            
            elif desc == "Unknown":
                lex.write(acc + "\t" + acc + "\n")
                
                    
#----------------------------      Main      ------------------------------#
'''
if len(sys.argv) == 4 :
    inpFile = sys.argv[1]
    outDir = sys.argv[2]
    jsonFile = sys.argv[3]
    
    ##Workflow
    Lexique(outDir, jsonFile)
    casTree = Cluster_Maker(inpFile, outDir)
    File_Cas_Typer(outDir)
    state = CSV_Summary(outDir,"tsv", casTree)
    print(state)
    ##
        
elif len(sys.argv) == 5 :
    inpFile = sys.argv[1]
    outDir = sys.argv[2]
    jsonFile = sys.argv[3]
    form = sys.argv[4]
    
    ##Workflow
    Lexique(outDir, jsonFile)
    casTree = Cluster_Maker(inpFile, outDir)
    File_Cas_Typer(outDir)
    state = CSV_Summary(outDir,form, casTree)
    print(state)
    ##
        
else :
    print("Clusters_CAS_Summary.py needs 3 or 4 arguments to work !\nClusters_Cas_Summary.py is a script which treats CAS results of CRISPRCasFinder for further analysis\n")
    print("Arguments : \n\t input/Directory/rawCas.fna : Path of the rawCas.fna result file obtained from CRISPRCasFinder\n\t output/Directory/ : Directory where will be stored all result files\n\t Directory/jsonFile.json : Path to the json result obtained from CRISPRCasFinder\n\t format : csv (or tsv as default value) (Optional argument)\n")
    print("Clusters_CAS_Summary.py example of command :\n Clusters_CAS_Summary.py CRISPRCasFinder_4.2.11/ResultDirectory/rawCas.fna Output/ CRISPRCasFinder_4.2.11/ResultDirectory/result.json csv")
'''