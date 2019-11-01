# -*- coding: utf-8 -*- 
'''
Created on 24 August 2018

@author: Nicolas Villeriot
'''

import regex

#----------------------------Function library------------------------------#
'''
Function allowing to read a fasta file and returning 
for each iteration a sequence name and its complete sequence
@param fastaFileName: Name and directory of the fasta file used as reference
@type fastaFileName: String
@return: a sequence and its name for each iteration 
'''
def ReadFasta(fastaFileName):
    with open(fastaFileName) as fastaFile:
        #Initialization
        name = None
        seq = []
        
        # the file line by line
        for line in fastaFile:
            #Discard "\n" at the line end
            line = line.rstrip()
            
            #Takes the chromosome name if the line begin by '>' and return it with the corresponding sequence
            if line.startswith(">"):
                if name: 
                    yield (name, ''.join(seq))
                name = line
                seq = []
            else:
                seq.append(line)
                
        #Return the last chromosome name and sequence
        if name: 
            yield (name, ''.join(seq))


'''
Function allowing to reverse complement DNA locus
@param dnaSeq : DNA sequence to reverse complement
@type dnaSeq : String
@return : String, reverse DNA sequence
'''      
def ReverseComp (dnaSeq):
    comp = ""
    for nuc in reversed(dnaSeq) :
        if nuc == "T" :
            comp += "A"
        elif nuc == "A" :
            comp += "T"
        elif nuc == "C" :
            comp += "G"
        elif nuc == "G" :
            comp += "C"
        else :
            comp += nuc
    
    return (comp)

'''
Function allowing to determine if two sequences are similar with a maximum number of match
@param seqMin: Smallest sequence
@param seqMax: Largest sequence
@param misMax: Maximum number of mismatch allowed for DR consensus association
@type seqMin, seqMax: String
@type misMax: Integer
@return: boolean, True if sequences are similar
'''
def CompareSeq (seqMin,seqMax,misMax):
    similar = False
    if (regex.search("("+seqMin+"){e<="+str(misMax)+"}",seqMax) != None) or (regex.search("("+seqMin+"){e<="+str(misMax)+"}",ReverseComp(seqMax)) != None) : #Repeat comparison depending to number of mismatch allowed
        similar = True
    return(similar)
#----------------------------      Main      ------------------------------#
