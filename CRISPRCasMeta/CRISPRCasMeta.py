# -*- coding: utf-8 -*- 
'''
Created on 26 June 2018

@author: Nicolas Villeriot
@version: 1.4.2
'''

import FastaSummary
import StatResult
import ContigAnalyser 
import ClustersCasSummary 
import time
import subprocess
import re
import os
import shutil
from optparse import OptionParser #command line argument parser
import BioSeqTools

#-------------------------------Function library------------------------------------------
'''
Function allowing to clean workspace before running.
@param out : Absolute path where will be registered all results
@type out : String
@return: String, error message
'''
def Cleaner(out):
    out = out.rstrip("/")
    cleanerFailed = ""
    print("Step 1/4 : Cleaner ...")
    start = time.time()
    try :
        #shutil.rmtree("/" + out, ignore_errors=True)
        outDir = subprocess.check_output(["mkdir",out])
        outDir += subprocess.check_output(["mkdir",out + "/MincedFiles"])
        print("Workspace preparation :\n" + outDir)
        end = time.time()
        print("Execution time (1/4) : " + str(end-start) + " s\n")
    except :
        cleanerFailed = "Error 001 : Workspace issues. \nMake sure the absolute path of your \"Output/\" folder exist and is accessible."
        print("Error 001 : Workspace issues. \nMake sure the absolute path of your \"Output/\" folder exist and is accessible.")
        end = time.time()
        print("Execution time (1/4) : " + str(end-start) + " s\n")
        
    return(cleanerFailed)
        
'''
Function allowing to create a summary fasta with only sequences for which CRT(Minced) found a CRISPR array.
@param fastaFile : Absolute path of the fasta or fasta.gz file
@param mincedPath : Absolute path of minced executable
@param out : Absolute path where will be registered all results
@param minRL : Minimal repeat length allowed  (CRISPRCasFinder parameter)
@param maxRL : Maximal repeat length allowed (CRISPRCasFinder parameter)
@param fasta : Input file format (True if fasta, False if fasta.gz)
@param repeatLib : Absolute path of a repeat library
@param nameForm : Reformat fasta repeat sequence header to recognize which sequence is from the database
@type fastaFile, mincedPath, out, minRL, maxRL, repeatLib, nameForm : String
@type fasta : boolean
@return: String, error message
'''
def StartMinced (fastaFile,mincedPath,out,fasta,minRL,maxRL, repeatLib,nameForm): 
    mincedOut = ""
    print("Step 2/4 : Minced ...")
    start = time.time()
    ##check input format
    if fasta :
        try :
            for name, seq in BioSeqTools.ReadFasta(fastaFile):
                if name[0] != ">" :
                    print("Error 003 : Fasta file, sequence name format not supported.")
                    mincedOut = "Error 003 : Fasta file, sequence name format not supported."
                if re.search("[^a-zA-Z]",seq) != None :
                    print("Error 004 : Fasta file, sequence format not supported.")
                    mincedOut = "Error 004 : Fasta file, sequence format not supported."
                    
                break
        except:
            print("Error 002 : Unreadable fasta file.")
            mincedOut = "Error 002 : Unreadable fasta file."
    
    elif re.search(".gz$",fastaFile) == None :
        print("Error 005 : Unreadable fasta gz file.")
        mincedOut = "Error 005 : Unreadable fasta gz file."
    ##
    if mincedOut == "" :
        try :
            FastaSummary.MetaSummary(fastaFile,out,mincedPath,fasta,minRL,maxRL)
            FastaSummary.RepeatDatabase(repeatLib, nameForm, out, fasta)
        except :
            mincedOut = "Failed"
    end = time.time()
    print("Execution time (2/4) : " + str(end-start) + " s\n")
    return(mincedOut)

'''
Function allowing to search CRISPR/Cas systems using CRISPRCasFinder metagenome version.
@param ccfPath : Absolute path of the CRISPRCasFinders metagenomic version
@param minRL : Minimal repeat length allowed  (CRISPRCasFinder parameter)
@param maxRL : Maximal repeat length allowed (CRISPRCasFinder parameter)
@param arMis : Repeat mismatch sensitivity (CRISPRCasFinder parameter)
@param minSS : Minimal spacer length allowed  (CRISPRCasFinder parameter)
@param maxSS : Maximal spacer length allowed (CRISPRCasFinder parameter)
@param maxPS : Maximal allowed percentage of similarity between Spacers (CRISPRCasFinder parameter)
@param maxPMR : Percentage mismatchs allowed between Repeats (CRISPRCasFinder parameter)
@param maxPMTR : Percentage mismatchs allowed for truncated Repeat (CRISPRCasFinder parameter)
@param betDetectDR : Option allowing to better detect the truncated DR (default: False) (CRISPRCasFinder parameter)
@param flankSize : Size of Flanking regions in base pairs (bp) for each analyzed CRISPR array (CRISPRCasFinder parameter)
@param cpuM : Number of cpu alocated for the script execution (CRISPRCasFinder parameter)
@param out : Absolute path where will be registered all results
@type ccfPath, minRL, maxRL, arMis, minSS, maxSS, maxPS, maxPMR, maxPMTR, betDetectDR, flankSize, cpuM, out : String
@return: String, error message
''' 
def CRISPRCasSearch (ccfPath,minRL,maxRL,arMis,minSS,maxSS,maxPS,maxPMR,maxPMTR,betDetectDR,flankSize,cpuM,out):
    bdt = ""
    if betDetectDR :
        bdt = "-bDT"
    
    n =""
    if arMis :
        n = "-n"
        
    ccfFailed = ""
    print("Step 3/4 : MetaFinder ...")
    start = time.time()
    try :
        print(" ".join(["perl",ccfPath + "MetaFinder.pl","-cf",ccfPath + "CasFinder-2.0.2","-def","General","-cas","-i",out + "Summary.fasta","-cpuM",cpuM,"-meta","-keep","-mr",minRL,"-xr",maxRL, n,"-pm", minSS,"-px",maxSS,"-s",maxPS,"-md", maxPMR,"-PMT",maxPMTR,bdt,"-fl",flankSize,"-so",ccfPath + "sel392v2.so"]))
        outCCF = subprocess.Popen(["perl",ccfPath + "MetaFinder.pl","-cf",ccfPath + "CasFinder-2.0.2","-def","General","-cas","-i",out + "Summary.fasta","-cpuM",cpuM,"-meta","-keep" ,"-mr",minRL,"-xr",maxRL, n,"-pm", minSS,"-px",maxSS,"-s",maxPS,"-md", maxPMR,"-PMT",maxPMTR,bdt,"-fl",flankSize,"-so",ccfPath + "sel392v2.so"], cwd=out)
        
        outCCF.wait()
        end = time.time()
        print("Execution time (3/4) : " + str(end-start) + " s\n")
    except :
        ccfFailed = "Error 102 : Command for MetaFinder failed !"
        print("Error 102 : Command for MetaFinder failed !")
        end = time.time()
        print("Execution time (3/4) : " + str(end-start) + " s\n")
        
        ##Clean the workspace to avoid issues and memory space consumption 
        print("Emergency workspace cleaning ...")
        resultPath = ""
        for root, dirnames, filenames in os.walk(ccfPath):
            for dirname in dirnames:
                if re.search("Summary",dirname) != None :
                    resultPath = ccfPath + dirname + "/"
        shutil.rmtree(resultPath, ignore_errors=True)
        print("A manual cleaning of the MetaFinder directory can be necessary.")
        ##
    
    return (ccfFailed)
        
'''
Function allowing to do further analysis on CRISPRCasFinder metagenome results.
@param misMax : Maximum number of mismatch allowed for DR comparison
@param out : Absolute path where will be registered all results
@param vers : current version of CRISPRCasMeta
@param commandLine : complete command line sent to CRISPRCasMeta
@param inputDir : input directory which contains the fasta or fasta.gz file
@param nameForm : Reformat fasta repeat sequence header to recognize which sequence is from the database
@type out, vers, commandLine, inputDir, nameForm : String
@type misMax : Integer
''' 
def MetaFinalisation (misMax,out,vers,commandLine,inputDir,nameForm): 
    print("Step 4/4 : Analysis ...")
    start = time.time()
    resultPath = ""
    for root, dirnames, filenames in os.walk(out):
        for dirname in dirnames:
            if re.search("Summary",dirname) != None :
                resultPath = out + dirname + "/"
            
    compRes = StatResult.JsonStat(resultPath + "result.json", misMax, out, vers, commandLine, nameForm)
    StatResult.JsonComplete(out + "result.json", out)
    print(" > Create a fasta file with for CRISPR/Cas systems.")
    StatResult.CrtFiltered(out + "Summary.fasta", out + "result.json", out) 
    print(" > Create an human readable Summary file (.csv format) from result.json.")
    ContigAnalyser.ResultSummary(out + "result.json",out,compRes)
    try :
        print(" > Cleaning the workspace.")
        shutil.move(resultPath + "rawCas.fna",out)
        shutil.move(resultPath + "rawCRISPRs.fna",out)
        shutil.rmtree(resultPath, ignore_errors=True)
        shutil.rmtree(out + "MincedFiles/", ignore_errors=True) 
        os.remove(out + "Minced.out") 
        os.remove(out + "Summary.fasta") 
        
    except :
        print("Error 006 : CRISPRCasMeta workspace cleaner failed !\nPossible issue : environment constraint (directories, files).\n")
    print(" > Creating a Cas clusters fasta file and its summary (.csv format).")
    ClustersCasSummary.Lexique(out, out + "result.json") 
    casTree = ClustersCasSummary.ClusterMaker(out + "rawCas.fna", out) 
    ClustersCasSummary.CsvSummary(out, "tsv", casTree) 
    try :
        os.remove(out + "Metadata.tsv")
        shutil.move(out,inputDir) #move output directory in the input directory
    except :
        print("Error 006 : CRISPRCasMeta workspace cleaner failed !\nPossible issue : environment constraint (directories, files).\n")
    end = time.time()
    print("Execution time (4/4) : " + str(end-start) + " s\n")

#-------------------------------      Main      ------------------------------------------
vers = "1.4.2"

parser = OptionParser() #command line parser
cwd = os.path.dirname(os.path.realpath(__file__)) + "/" #current working directory

##command line arguments 
parser.add_option("-i", "--input", action="store", type = "string", dest="fastaFile", default="", help="Absolute path of the fasta or fasta.gz input file (metagenome)")
parser.add_option("-o", "--output", action="store", type = "string", dest="out", default=(cwd + "Output/"), help="Absolute path of the output folder where the results will be registered")
parser.add_option("-m", "--minced", action="store", type = "string", dest="mincedPath", default=(cwd + "minced-0.3.0/minced"), help="Absolute path of CRT(minced) program")
parser.add_option("--ccm", action="store", type = "string", dest="ccfPath", default=(cwd + "MetaFinder/"), help="Absolute path of MetaFinder metagenomic version folder")
parser.add_option("--mis", action="store", type = "string", dest="misMax", default="-1", help="Maximum number of mismatch allowed for DR comparison. By default, CRISPRCasMeta do not perform DR comparison. (default: -1, recommended: 1, limits: 0 to 5)") #06/08/2018
parser.add_option("--mr","--mindr", action="store", type = "string", dest="minRL", default="23", help="Minimal repeat length (default: 23)")
parser.add_option("--xr","--maxdr", action="store", type = "string", dest="maxRL", default="55", help="Maximal repeat length (default: 55)")
parser.add_option("-n","--nomism", action="store_true", dest="arMis", help="Option used to do not allow mismatches (default value is $mismOne when this option is not called. i.e. mismatches are allowed by default)")
parser.add_option("--pm","--minsp", action="store", type = "string", dest="minSS", default="0.6", help="Minimal Spacers size in function of Repeat size (default: 0.6)")
parser.add_option("--px","--maxsp", action="store", type = "string", dest="maxSS", default="2.5", help="Maximal Spacers size in function of Repeat size (default: 2.5)")
parser.add_option("-s","--spsim", action="store", type = "string", dest="maxPS", default="60", help="Maximal allowed percentage of similarity between Spacers (default: 60)")
parser.add_option("--md","--mismdrs", action="store", type = "string", dest="maxPMR", default="20", help="Percentage mismatchs allowed between Repeats (default: 20)")
parser.add_option("-t", "--truncdr", action="store", type = "string", dest="maxPMTR", default="33", help="Percentage mismatchs allowed for truncated Repeat (default: 33)")
parser.add_option("--bdt", action="store_true", dest="betDetectDR", help="Option allowing to better detect the truncated DR (default: False)")
parser.add_option("--fl","--flank", action="store", type = "string", dest="flankSize", default="100", help="Size of Flanking regions in base pairs (bp) for each analyzed CRISPR array (default: 100)")
parser.add_option("--cpum", action="store", type = "string", dest="cpuM", default="1", help="Option allowing to set number of CPUs to use for MacSyFinder (default: 1)")
parser.add_option("--fa", action="store_true", dest="fasta", help="Option allowing to set the input format to fasta (by default input is a fasta.gz format)")
parser.add_option("-l", "--replib", action="store", type = "string", dest="repeatLib", default="", help="Absolute path of a repeat library (--mis option mandatory, fasta or fasta.gz format depending of --fa option)")
parser.add_option("--forlib", action="store", type = "string", dest="libFormat", default="dbHeader", help="Reformat fasta repeat sequence header to recognize which sequence is from the database(--mis and --reLib options mandatory)")


options, args = parser.parse_args()

fastaFile = options.fastaFile
out = cwd + options.out + "/"
mincedPath = options.mincedPath
ccfPath = options.ccfPath
misMax = options.misMax
minRL = options.minRL
maxRL = options.maxRL
arMis = options.arMis
minSS = options.minSS
maxSS = options.maxSS
maxPS = options.maxPS
maxPMR = options.maxPMR
maxPMTR = options.maxPMTR
betDetectDR = options.betDetectDR
flankSize = options.flankSize
cpuM = options.cpuM
fasta = options.fasta
##
##Use repeat database for repeat comparison
repeatLib = options.repeatLib
libFormat = ">" + options.libFormat + "_"
if misMax == "-1" and repeatLib != "" :
    print("Error 201 : use repeat library option (--repLib) whereas --mis option not used!\nCannot perform repeat comparison, process end!\n")
    fastaFile = ""
    
##

## Command line
optAr = ""
optBdt = ""
optFa = ""
if arMis :
    optAr = " -n "
if betDetectDR :
    optBdt = " --bdt "
if fasta :
    optFa = " --fa "
commandLine = "python CRISPRCasMeta.py -i " + fastaFile + " -o "  + out + " -m " + mincedPath + " --ccm " + ccfPath + " --mis " + misMax + " --mr " + minRL + " --xr " + maxRL + optAr + " --pm " + minSS + " --px " + maxSS + " -s " + maxPS + " --md " + maxPMR + " -t " + maxPMTR + optBdt + " --fl " + flankSize + " --cpum " + cpuM + optFa + " --repLib " + repeatLib + " --forLib " + libFormat
##

print(betDetectDR)
print(arMis)
##

stopProcess = ""

startTot = time.time()
stopProcess = Cleaner(out)
if stopProcess == "" and fastaFile != "" :
    stopProcess = StartMinced(fastaFile,mincedPath,out,fasta,minRL,maxRL,repeatLib,libFormat)
stopProcess = ""    
if stopProcess == "" and fastaFile != "" :
    stopProcess = CRISPRCasSearch(ccfPath,minRL,maxRL,arMis,minSS,maxSS,maxPS,maxPMR,maxPMTR,betDetectDR,flankSize,cpuM,out)
if stopProcess == "" and fastaFile != "" :
    inputDir = "/".join(fastaFile.split("/")[:-1]) + "/"
    MetaFinalisation(int(misMax),out,vers,commandLine,inputDir,libFormat)
#shutil.rmtree(out, ignore_errors=True)
    
endTot = time.time()

print("Execution time (Total) : " + str(endTot-startTot) + " s\n")
