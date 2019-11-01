# -*- coding: utf-8 -*- 
'''
Created on 31 May 2018

@author: Nicolas Villeriot
'''

import BioSeqTools
import json
import re
import subprocess

#-------------------------------Function library------------------------------------------
'''
This function allows a subsequent treatment of MetaFinder results.
Its main purpose is to associate contigs (sequences) by there DR consensus.
@param jsonFileName : Path to the result.json file obtained with MetaFinder
@param misMax : Maximum number of mismatch allowed for DR consensus association
@param outPath : Absolute path where will be registered all results
@param vers : current version of CRISPRCasMeta
@param commandLine : complete command line sent to CRISPRCasMeta
@param nameForm : Reformat fasta repeat sequence header to recognize which sequence is from the database
@type jsonFileName, outPath, vers, commandLine, nameForm : String
@type misMax : Integer
@return : Hash, results of the dr comparison
'''
def JsonStat (jsonFileName,misMax,outPath,vers,commandLine,nameForm): 
    
    sta = {}
    #with open(outPath + "RepeatLib.fasta","w") as iniFi :
    #    iniFi.write("")
    
    ##Load json file into a tree
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    for contig in datas["Sequences"] :
        evLev = 0
        nbrCasType = 0
        nbrCas = 0
        cr = []
        casPos = []
        
        ##Repeat library
        repList = []
        if contig["Crisprs"] != []:
            for crispr in contig["Crisprs"] :
                if int(crispr["Evidence_Level"]) > 2 or contig["Cas"] != []: #only reliable results
                    repList.append(crispr["DR_Consensus"])
                
            AddRepeatLibrary(">" + contig["Id"],repList,outPath + "RepeatLib.fasta")
        ##
        
        if contig["Crisprs"] != [] and contig["Cas"] != [] :
            ##Number of evidence level > 2, CRISPR position
            for crispr in contig["Crisprs"] :
                if int(crispr["Evidence_Level"]) > 2 :
                    evLev += 1
                    cr.append(str(crispr["Start"]) + "_" + str(crispr["End"]))
            ##
            ##Number of cluster Cas and genes, cluster Cas position     
            for casSys in contig["Cas"] :
                nbrCasType += 1
                casPos.append(str(casSys["Start"]) + "_" + str(casSys["End"]))
                for cas in casSys["Genes"] :
                    nbrCas += 1
            ##
            ##Number of CRISPR/Cas systems with a CRISPR - Cas distance < 5 Kb
            prox = 0
            for locCR in cr:
                locCR = [int(locCR.split("_")[0]),int(locCR.split("_")[1])]
                for locCas in casPos :
                    locCas = [int(locCas.split("_")[0]),int(locCas.split("_")[1])]
                    if (locCas[0] > locCR[1] and locCas[0] < locCR[1] + 5000) or (locCas[1] > locCR[0] and locCas[1] < locCR[0] + 5000) :
                        prox += 1
            ##
                    
            sta[contig["Id"]] = [evLev,nbrCasType,nbrCas,nbrCas/nbrCasType,prox]
    

    #DrComparisonFast(outPath + "RepeatLib.fasta", misMax,jsonFileName,outPath,vers) #algorithm for DR comparison
    compRes = ByDrComparison(outPath + "RepeatLib.fasta", misMax,jsonFileName,outPath,vers,commandLine,nameForm) #algorithm for repeat comparison
    EvidenceLevelUpgrader(compRes, outPath + "result.json",nameForm)
        
    with open(outPath + "stat.tsv","w") as statFi :
        statFi.write("Id\tEvidence Level\tNbr CAS Type\tNbr CAS\tCAS Per Type\tCompleteSystem\n")
    
    with open(outPath + "stat.tsv","a") as statFi :
        for key in sta.keys() :
            statFi.write(key + "\t" + str(sta[key][0])+ "\t" + str(sta[key][1]) + "\t" + str(sta[key][2]) + "\t" + str(sta[key][3]) + "\t" + str(sta[key][4]) + "\n")
    return(compRes)

'''
Function allowing to add a new DR consensus to a fasta format library of repeat.
@param seqId : Name/Login of the current contig (sequence)
@param repeatList : List of DR consensus found in the current contig 
@param reLib : Path of the fasta format library of repeat.
@type seqId, reLib : String
@type repeatList : Array of String
'''       
def AddRepeatLibrary (seqId,repeatList,reLib):
    drNum = 0
    with open(reLib,"a") as reLibFile :
        if repeatList != [] :
            for repSeq in repeatList :
                drNum += 1
                reLibFile.write(seqId + "|" + str(drNum) + "\n" + repSeq + "\n")
                
        

'''
Function allowing to compare each DR consensus of a result.json file against a fasta format library of repeat.
@param jsonFileName : Path of the summary_Result.json file obtained with MetaFinder.
@param reLib : Path of the fasta format library of repeat.
@param misMax : Maximum number of mismatch allowed for DR consensus association
@param coverMin : Minimum alignment coverage allowed for DR consensus association
@type jsonFileName, reLib : String
@type misMax, coverMin : Integer
'''            
def ContigAssociation (jsonFileName,reLib,misMax,coverMin):
    listPotContig = []
    with open("Output/ContigsBlast.out","w") as bRes:
        bRes.write("")
    ##Load json file into a tree
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    for contig in datas["Sequences"] :
        if contig["Crisprs"] != [] :
            for crispr in contig["Crisprs"] :
                with open("Output/Temp.fasta","w") as temp :
                    temp.write("")
                
                with open("Output/Temp.fasta","w") as temp :
                    temp.write(">" + crispr["Name"] + "\n" + crispr["DR_Consensus"] + "\n")
                
                blastR = BlastAgainstFasta("Output/Temp.fasta", reLib)
                
                #print("\nCRISPR : " + crispr["Name"])
                with open("Output/ContigsBlast.out","a") as bRes:
                    bRes.write("\nCRISPR :\t" + crispr["Name"] + "\n")
                if blastR != None :
                        for res in blastR.split("\n") :
                            if res != "" :
                                #print(res)
                                with open("Output/ContigsBlast.out","a") as bRes:
                                    if (int(res.split("\t")[-1]) <= misMax) and (int(res.split("\t")[-4]) >= coverMin) and (re.search("_".join(crispr["Name"].split("_")[:-1]),res) == None) :
                                        bRes.write(res + "\n")
                                        #print("Potential : " + res.split("\t")[0].split("|")[0] + "\n")
                                        if not (res.split("\t")[0].split("|")[0] in listPotContig):
                                            listPotContig.append(res.split("\t")[0].split("|")[0])
    
    return (listPotContig)
 
'''
Function allowing to create groups of DR consensus with a fasta format library of repeat.
@param reLib : Path of the fasta format library of repeat.
@param misMax : Maximum number of mismatch allowed for DR consensus association
@param jsonFileName : Path to the result.json file obtained with MetaFinder
@param outPath : Absolute path where will be registered all results
@param vers : current version of CRISPRCasMeta
@type reLib, jsonFileName, outPath, vers : String
@type misMax : Integer
'''
def DrComparisonFast (reLib, misMax,jsonFileName, outPath,vers):
    compRes = {}
    drGroup = 0
    ##Repeat comparison, group repeats depending of the number of mismatch allowed
    if misMax != -1 :
        for name, seq in BioSeqTools.ReadFasta(reLib) :
            name = "_".join(name[1:].split("|"))
            lenSeq = len(seq)
            exist = False
            if compRes != {} :
                for drG in compRes.keys() :
                    for rep in compRes[drG].keys() :
                        if lenSeq >= len(compRes[drG][rep]) :
                            resE = BioSeqTools.CompareSeq(compRes[drG][rep], seq, misMax)
                        
                        else :
                            resE = BioSeqTools.CompareSeq(seq, compRes[drG][rep], misMax)
                    
                        if resE :
                            exist = True
                            compRes[drG][name] = seq
                            break
            
                if not exist :
                    drGroup += 1
                    compRes["Repeat_" + str(drGroup)] = {}
                    compRes["Repeat_" + str(drGroup)][name] = seq
                
            else :
                drGroup += 1
                compRes["Repeat_" + str(drGroup)] = {}
                compRes["Repeat_" + str(drGroup)][name] = seq
    ##
    
    ##Load json file into a tree ##10/07/2018
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    ##Add repeat comparison results to the result.json file
    datas["Version"] = vers
    for contig in datas["Sequences"] :
        if misMax != -1 :
            contig["RepeatGroups"] = [] 
            print("Dr comparison for : " + contig["Id"])
            
            for drG in compRes.keys() :
                groupRes = []
                find = False
                currentRes = {} #current CRISPR inclusion
                for rep in compRes[drG].keys() :
                    if re.search(contig["Id"],rep)!=None :
                        find = True
                        currentRes = {"Name" : rep,"Consensus" : compRes[drG][rep]} #current CRISPR inclusion
                    else :
                        groupRes.append({"Name" : rep,"Consensus" : compRes[drG][rep]})
                        
                if find and groupRes != [] :
                    groupRes.append(currentRes) #current CRISPR inclusion
                    contig["RepeatGroups"].append({"Name" : drG, "Crisprs" : groupRes})
                    
    with open(outPath + "result.json","w") as resFi :
        resFi.write("")
        
    with open(outPath + "result.json", "a") as resFi :
        resFi.write(json.dumps(datas, sort_keys=True,indent=1, separators=(',', ': ')))
    ##
    
'''
Function allowing to compare each DR consensus of a fasta format library of repeat.
@param reLib : Path of the fasta format library of repeat.
@param misMax : Maximum number of mismatch allowed for DR consensus association
@param jsonFileName : Path to the result.json file obtained with MetaFinder
@param outPath : Absolute path where will be registered all results
@param vers : current version of CRISPRCasMeta
@param commandLine : complete command line sent to CRISPRCasMeta
@param nameForm : Reformat fasta repeat sequence header to recognize which sequence is from the database
@type reLib, jsonFileName, outPath, vers, commandLine, nameForm : String
@type misMax : integer
@return : Hash, results of the dr comparison
'''
def ByDrComparison (reLib, misMax,jsonFileName, outPath,vers,commandLine,nameForm): 
    ##Repeat Comparison
    compRes = {}
    if misMax != -1 :
        print(" - Repeat comparison.")
        for name, seq in BioSeqTools.ReadFasta(reLib) :
            ##Only repeat of the studied metagenome are compared to the others
            follow = False
            
            if nameForm == "" :
                follow = True
            elif re.search(nameForm,name) == None :
                follow = True
            ##
            
            if follow :
                name = name[1:]
                #lenSeq = len(seq)
                for name1, seq1 in BioSeqTools.ReadFasta(reLib) :
                    name1 = name1[1:]
                    exist = False
                    if len(seq) >= len(seq1) :
                        exist = BioSeqTools.CompareSeq(seq1, seq, misMax)
                    else :
                        exist = BioSeqTools.CompareSeq(seq, seq1, misMax)
            
                    if exist :
                        if not (name in compRes.keys()) :
                            compRes[name] = {}
                            compRes[name][name] = seq
                
                        if name != name1 :
                            compRes[name][name1] = seq1
    ##
    
    ##Load json file into a tree ##10/07/2018
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    ##Add repeat comparison results to the result.json file
    datas["Version"] = vers
    datas["Command"] = commandLine
    drGr = 0
    for contig in datas["Sequences"] :
        contig["RepeatGroups"] = [] 
        #print("Dr comparison for : " + contig["Id"])
        if misMax != -1 :
            
            if compRes != {} :
                for crispr in compRes.keys() :
                    crisprSim = []
                    if re.search(contig["Id"],crispr)!= None :
                        for rep in compRes[crispr].keys() :
                            if rep != crispr :
                                crisprSim.append({"Name" : "_".join(rep.split("|")), "Consensus" : compRes[crispr][rep]})
                
                    if crisprSim != [] :
                        drGr += 1
                        contig["RepeatGroups"].append({"Name" : "Repeat_" + str(drGr), "Crisprs" : crisprSim})
    
    with open(outPath + "result.json","w") as resFi :
        resFi.write("")
        
    with open(outPath + "result.json", "a") as resFi :
        resFi.write(json.dumps(datas, sort_keys=True,indent=1, separators=(',', ': ')))    
    
    return (compRes)

'''
Function allowing to run a local blast (local database) on fasta file
@param querySeq : Absolute path of the query fasta file
@param dbFilePath: Absolute path of the local database used for the Blast
@type querySeq, dbFilePath: String
@return: String, error message
'''
def BlastAgainstFasta(querySeq, dbFilePath):
    try :
        res = ""
        #print("blastn -query " + querySeq + " -subject " + dbFilePath + " -outfmt 6 sacc stitle qcovs evalue pident")
        subprocess.check_output(["makeblastdb", "-dbtype", "nucl", "-in", dbFilePath])
        #res = subprocess.check_output(["blastn",  "-query", querySeq, "-subject", dbFilePath,"-word_size", "11" ,"-outfmt", "6 sacc stitle qseq sseq qcovs evalue pident sstart send"])
        res = subprocess.check_output(["blastn",  "-query", querySeq, "-subject", dbFilePath,"-word_size", "11" ,"-outfmt", "6 sacc stitle qseq sseq qcovs evalue pident mismatch"])
        return (res)
    except :
        print("Error: Command bash to run Blast failed ! ")
        
'''
Function allowing to complete MetaFinder CRISPR results adding those of Minced.
@param jsonFileName: Absolute path of the result.json file obtained with MetaFinder
@param out: Absolute path where will be registered all results
@type jsonFileName, out: String
'''            
def JsonComplete (jsonFileName,out):
    print(" - Complete MetaFinder results with CRT (Minced) results.")
    with open(out + "Completion.tsv","w") as resFile :
        resFile.write("CONTIG\tSpacers\tCRISPR Start\tCRISPR End\tCAS Cluster\tCAS Start\tCAS End\n")
    
    ##Load json file into a tree
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    for contig in datas["Sequences"] :
        if contig["Crisprs"] == [] and contig["Cas"] != [] :
            for casCluster in contig["Cas"] :
                listCas = ""
                for cas in casCluster["Genes"] :
                    listCas += cas["Sub_type"] + " "
                    
                with open(out + "CrisprResult.csv") as mcsv :
                    for line in mcsv :
                        line = line.split(";")
                        ##avoid problematic characters
                        line[0] = line[0].replace("|","_")
                        ##
                        if re.search(contig["Id"],line[0]) != None :
                            
                            with open(out + "Completion.tsv","a") as resFile :
                                resFile.write(contig["Id"] + "\t" + line[1] + "\t" + line[3] +"\t" + line[4] +"\t" + listCas +"\t" + str(casCluster["Start"]) + "\t" + str(casCluster["End"]) + "\n")
                                ##Complete the result.json file with the potential CRISPRs found by CRT (Minced)
                                crisprRes = {}
                                notACrispr = True
                                crisprRes["Name"] = line[0]
                                crisprRes["Start"] = int(line[3])
                                crisprRes["End"] = int(line[4])
                                crisprRes["DR_Consensus"] = "CRT results"
                                crisprRes["Repeat_ID"] = "CRT results"
                                crisprRes["DR_Length"] = int(line[2]) 
                                crisprRes["Spacers"] = int(line[1])
                                crisprRes["Potential_Orientation"] = "ND"
                                crisprRes["CRISPRDirection"] = "ND"
                                crisprRes["Evidence_Level"] = 1
                                crisprRes["Conservation_DRs"] = 0
                                crisprRes["Conservation_Spacers"] = 0
                                    
                                crisprRes["Regions"] = []
                                leftFlank = {}
                                leftFlank["Type"] = "LeftFLANK"
                                leftFlank["Start"] = int(line[3])
                                leftFlank["End"] = int(line[3])
                                leftFlank["Sequence"] = "UNKNOWN"
                                leftFlank["Leader"] = 0
                                crisprRes["Regions"].append(leftFlank)
                                    
                                currentPos = int(line[3])
                                for reg in line[6].split(" ") :
                                    regInfo = {}
                                    notACrispr = False
                                    if len(reg.split("_")) == 2 :
                                        regType = reg.split("_")[0]
                                        regSeq = reg.split("_")[1]
                                        regSize = len(regSeq)
                                        
                                        if regSeq != "" :
                                            if regType == "R" :
                                                regInfo["Type"] = "DR"
                                                regInfo["Start"] = currentPos
                                                currentPos += regSize
                                                regInfo["End"] = currentPos
                                                regInfo["Sequence"] = regSeq
                                                crisprRes["Regions"].append(regInfo)
                                            elif regType == "S" :
                                                regInfo["Type"] = "Spacer"
                                                regInfo["Start"] = currentPos
                                                currentPos += regSize
                                                regInfo["End"] = currentPos
                                                regInfo["Sequence"] = regSeq
                                                crisprRes["Regions"].append(regInfo)
                                            
                                        
                                    
                                rightFlank = {}
                                rightFlank["Type"] = "RightFLANK"
                                rightFlank["Start"] = int(line[4])
                                rightFlank["End"] = int(line[4])
                                rightFlank["Sequence"] = "UNKNOWN"
                                rightFlank["Leader"] = 0
                                crisprRes["Regions"].append(rightFlank)
                                ##
                                    
                                ##Avoid to repeat the same CRISPR for each CAS cluster near
                                isHere = False
                                if contig["Crisprs"] != [] :
                                    for crispr in contig["Crisprs"] :
                                        if re.search(crispr["Name"],crisprRes["Name"]) != None :
                                            isHere = True
                                        
                                    
                                if (not notACrispr) and (not isHere):
                                    contig["Crisprs"].append(crisprRes)
                                ##
                            
                            
        ##minimal distance between a crispr and a cluster Cas
        if contig["Crisprs"] != [] and contig["Cas"] != [] :
            posCR = []
            for crispr in contig["Crisprs"] :
                posCR.append(int(crispr["Start"]))
                posCR.append(int(crispr["End"]))
            
            minDist = -1
            for cas in contig["Cas"] :
                for gene in cas["Genes"] :
                    for pos in posCR :
                        startDist = abs(pos - int(gene["Start"]))
                        endDist = abs(pos - int(gene["End"]))
                        if startDist < minDist and startDist < endDist :
                            minDist = startDist
                        elif endDist < minDist and endDist < startDist :
                            minDist = endDist
                        elif minDist == -1 :
                            if startDist < endDist :
                                minDist = startDist
                            elif endDist < startDist :
                                minDist = endDist
            
            contig["CrisprCasMinimalDistance"] = minDist
        else :
            contig["CrisprCasMinimalDistance"] = -1   
        ##
        
    ##New result.json file                              
    with open(out + "result.json", "w") as compFi :
        compFi.write("")
    
    with open(out + "result.json", "a") as compFi :
        compFi.write(json.dumps(datas, sort_keys=True,indent=1, separators=(',', ': ')))
    ##
    
'''
Function allowing to create a fasta file with only CRISPR/Cas systems sequence + 1 Kb flanking regions
@param summaryFasta: Absolute path of the fasta file containing only sequences with CRISPRs found with CRT(Minced)
@param jsonFileName: Absolute path of the result.json file obtained with MetaFinder
@param out: Absolute path where will be registered all results
@type summaryFasta, jsonFileName, out: String
'''
def CrtFiltered(summaryFasta, jsonFileName, out):
    print(" - Create a fasta file of CRISPR/Cas systems.")
    ##Load json file into a tree
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    with open(out + "CRTFilteredContigs.fasta","w") as filtFi :
        filtFi.write("")
    
    with open(out + "CRTFilteredContigs.fasta","a") as filtFi :
          
        for name, seq in BioSeqTools.ReadFasta(summaryFasta):
            name = name[1:]
            for contig in datas["Sequences"] :
                if name == contig["Id"] :
                    completeName = ">" + name
                    zoomSeq = ""
                    listPos = []
                    spacers = 0
                    casGenes = 0
                    sumEvLev = 0
                
                    ##Resize the sequence
                    if contig["Crisprs"] != [] :
                        for crispr in contig["Crisprs"] :
                            if int(crispr["Evidence_Level"]) > 2 or contig["Cas"] != [] : ## only reliable results 17/08/2018
                                listPos.append(crispr["Start"])
                                listPos.append(crispr["End"])
                                spacers += crispr["Spacers"]
                                sumEvLev += crispr["Evidence_Level"]
                
                    if contig["Cas"] != [] :
                        for cas in contig["Cas"] :
                            listPos.append(cas["Start"])
                            listPos.append(cas["End"])
                            casGenes += len(cas["Genes"])
                    
                    if listPos != [] :
                        maxPos = max(listPos)
                        minPos = min(listPos)
                        if (minPos - 1001) >= 0 and (maxPos + 1000) <= len(seq) :
                            zoomSeq = seq[minPos-1001:maxPos+1000]
                        elif (minPos - 1001) >= 0  :
                            zoomSeq = seq[minPos-1001:]
                        elif (maxPos + 1000) <= len(seq):
                            zoomSeq = seq[:maxPos+1000]
                        else :
                            zoomSeq = seq
                        ##Resize the sequence
                    
                        completeName += "|Crisprs_" + str(len(contig["Crisprs"])) + "|Spacers_" + str(spacers) + "|SumEvidenceLevel_" + str(sumEvLev) + "|Cas_" + str(len(contig["Cas"])) + "|Genes_" + str(casGenes) + "\n"
                        filtFi.write(completeName + zoomSeq + "\n")
                    break
            
'''
Function allowing to upgrade evidence level (1,2 and 3) of a CRISPR when a similar repeat has an evidence level 4
@param compRes : Data hash tree which contains dr groups (output of DRComparisonFast function)
@param jsonFileName: Absolute path of the result.json file obtained with MetaFinder
@param nameForm : Reformat fasta repeat sequence header to recognize which sequence is from the database
@type compRes: Hash
@type jsonFileName, nameForm: String
'''
def EvidenceLevelUpgrader(compRes,jsonFileName, nameForm):
    ##Load json file into a tree
    with open(jsonFileName, 'r') as fi :
        datas = json.load(fi)
    ##
    
    if compRes != {} :
        print(" - Evidence level upgrade.")
        for repe1 in compRes.keys() :
            cont1 = repe1.split("|")[0]
            crisprName = "_".join(repe1.split("|"))
            for contig in datas["Sequences"] :
                if contig["Id"] == cont1 and contig["Crisprs"] != []:
                    for crispr in contig["Crisprs"] :
                        
                        ##Case where current crispr can be upgraded
                        if crispr["Evidence_Level"] < 4 :
                            for repe2 in compRes[repe1].keys():
                                cont2 = repe2.split("|")[0]
                                repe2 = "_".join(repe2.split("|"))
                                if re.search(nameForm,repe2) != None and nameForm != "" : ##Repeat database for evidence level upgrade
                                    crispr["Evidence_Level"] = 4
                                    break
                                else :
                                    for contig2 in datas["Sequences"] :
                                        if contig2["Id"] == cont2 and contig2["Crisprs"] != []:
                                            for crispr2 in contig2["Crisprs"] :
                                                if crispr2["Evidence_Level"] == 4 :
                                                    crispr["Evidence_Level"] = 4
                                                    break
                                            break
                        ##
                    break
              
    
    ##New result.json file                              
    with open(jsonFileName, "w") as compFi :
        compFi.write("")
    
    with open(jsonFileName, "a") as compFi :
        compFi.write(json.dumps(datas, sort_keys=True,indent=1, separators=(',', ': ')))
    ##
               
#-------------------------------Main------------------------------------------
'''
start = time.time()
JsonStat("/home/villeriot/Test/MetaFinder/B2F05_Summary_1000_25_6_2018_11_2_55/result.json","/home/villeriot/Test/MetaFinder/B2F05_Summary_1000_25_6_2018_11_2_55/summary_Result.json",1,50)
JsonComplete("/home/villeriot/Test/MetaFinder/B2F05_Summary_1000_25_6_2018_11_2_55/summary_Result.json", "Output/CrisprResult.csv", "Output/Completion.tsv")
end = time.time()
print("Execution time : " + str(end-start) + " s")
'''