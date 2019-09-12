from Scripts.BEseqType import *
from Bio.Seq import Seq
import csv
from Scripts.baseEditors import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter


def Main(snp):
    # sort DNA. First, determine which bases we wish to replace. 4 cases:
    # 1. C to T: use CBE list       2. A to G: use ABE list
    # 3. T to C: switch to reverse complement and use ABE
    # 4. G to A: switch to RC and use CBE
    for snp in SNPS:
        if snp.mutation == "C" and snp.wt == "T":
            BElist = CBElist
            MinorBElist=CBElistMinor
        elif snp.mutation == "A" and snp.wt == "G":
            BElist = ABElist
            MinorBElist = ABElistMinor
        elif snp.mutation == "T" and snp.wt == "C":
            snp=getRevComp(snp)
            snp.mutation = "A"
            snp.wt="G"
            BElist = ABElist
            MinorBElist = ABElistMinor
        elif snp.mutation == "G" and snp.wt == "A":
            snp = getRevComp(snp)
            snp.mutation = "C"
            snp.wt="T"
            BElist = CBElist
            MinorBElist = CBElistMinor
        else:
            BElist=None
            MinorBElist=None


        try:
          # check for matches in major window
            check_match = matchBE(snp, BElist)
            matches.update(check_match)

            #check for matches in minor window
            check_match_minor = matchBE(snp, MinorBElist)
            matchesMinor.update(check_match_minor)

            #check for clean match
            clean, quiet= cleanMatch(snp,check_match,BElist)
            cleanMatchdic.update(clean)
            quietMatchdic.update(quiet)




        except:
            print ("Error: %s is not an appropriate snp" %snp.snpID)