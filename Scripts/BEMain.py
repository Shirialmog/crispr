from BEseqType import *
from Bio.Seq import Seq
import csv
from baseEditors import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter


def importSNPS(SNPfile):
    rslts=[]
    with open(SNPfile,"r") as csv_file:
        try:
            csv_read = csv.reader(csv_file, delimiter=",")
        except MemoryError:
            print("MemoryError, try splitting your file into smaller files")
            exit()

        for row in csv_read:
            snpID=row[0]
            find_mut_len=len(row[9])
            """
            the sequence is given as a string, with 25 bases before mutation and 25 after. 
            need to find length of mutation to splice properly
            """
            mut_len=(find_mut_len-1)/2
            mutation=row[9][2]
            wt=row[9][0]
            sequence5 = row[7][0:25]
            sequence3 = row[7][26:]
            readingFrame=row[14]
            aaPosition=row[15]
            geneName=row[12]
            geneID=row[13]
            newSNP=BEseqType(snpID,sequence5, sequence3, wt, mutation, readingFrame, aaPosition, geneName, geneID)
            rslts.append(newSNP)
        rsltsDic={}
        num=0
        for snp in rslts:
            rsltsDic[snp.snpID]=snp.geneName
            num=num+1

        return rslts, rsltsDic

def matchBE(snp, BElist):
    matches_list = []
    matches = {}
    for BE in BElist:
        if BElist[BE][5] == "U":
            sequence=snp.seq3
        else:
            sequence=snp.seq5
            sequence=sequence[::-1] # get reverse


        PAM = BElist[BE][0]
        # find if PAM matches in correct place
        if len(sequence) > BElist[BE][2]:
            match = False
            start = BElist[BE][1] -1
            end = BElist[BE][2]-1
            window = (end-start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                j=0
                while (j<=lenPAM):
                    if sequence[temp] in BEletter[PAM[j]]:
                        j=j+1
                        temp=temp+1
                    else:
                        break
                    if j == len(PAM):
                        match = True
                        break
            if match is True:
                matches_list.append(BE)

        matches[snp.snpID] = matches_list
    return matches

def cleanMatch(snp, Matches, BElist):
    fnum=0
    cleanMatches={}
    clean_list=[]
    quietDic={}
    quiet_list=[]
    for match in Matches:
        for BE in Matches[match]:
            if BElist[BE][5] == "U":
                seq3 = snp.seq3
                seq5=snp.seq5
            else:
                seq3 = snp.seq5
                seq5 = snp.seq3
                seq3 = seq3[::-1]  # get reverse
                seq5 = seq5[::-1]  # get reverse

            PAM=BElist[BE][0]
            start=BElist[BE][1]-1
            end=BElist[BE][2]-1
            #direction=BElist[BE][5]
            locations=[]
            #find locations of PAM. for each location, check if clean
            window = (end - start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                j=0
                while (j<=lenPAM):
                    if seq3[temp] in BEletter[PAM[j]]:
                        j=j+1
                        temp=temp+1
                    else:
                        break
                    if j == len(PAM):
                        locations.append(i + start) #maybe j is not neccesary?
                        break

            totalSeq1=seq5+snp.wt+seq3
            diff=0
            if snp.readingFrame=="1":
                if len(snp.seq5)%3==1:
                    totalSeq1=totalSeq1[1:]
                    diff = 1
                elif len(snp.seq5)%3==2:
                    totalSeq1=totalSeq1[2:]
                    diff = 2
            elif snp.readingFrame=="2":
                if len(snp.seq5) % 3 == 0:
                    totalSeq1 = totalSeq1[2:]
                    diff = 2
                elif len(snp.seq5) % 3 == 2:
                    totalSeq1 = totalSeq1[1:]
                    diff = 1
            elif snp.readingFrame=="3":
                if len(snp.seq5) % 3 == 0:
                    totalSeq1 = totalSeq1[1:]
                    diff = 1
                elif len(snp.seq5) % 3 == 1:
                    totalSeq1 = totalSeq1[2:]
                    diff = 2
            if len(totalSeq1)%3==1:
                totalSeq1=totalSeq1[:-1]
            if len(totalSeq1)%3==2:
                totalSeq1=totalSeq1[:-2]
            totalSeq2=Seq(totalSeq1)
            protein_seq=totalSeq2.translate()
            max_num = 0
            protein_match=False
            for loc in locations:
                temp_seq=totalSeq1
                loc=loc+len(snp.seq5)
                activation_window = totalSeq1[loc-end-diff:loc-start]
                num = 0  # number of times the variant appears within the activation window
                new_AW = []
                for i in range(len(activation_window)):
                    if activation_window[i]==snp.mutation:
                        num=num+1
                        new_AW.append(snp.wt)
                    else:
                        new_AW.append(activation_window[i])
                new_AW=''.join(new_AW)
                if num>max_num:
                    max_num=num
                finalSeq=totalSeq1[0:loc-end-diff]+new_AW+totalSeq1[loc-start:]
                protein_seq_new=Seq(finalSeq).translate()

                if protein_seq==protein_seq_new:
                    protein_match=True

            if max_num==0:
                    clean_list.append(BE)
            if protein_match==True:
                    quiet_list.append(BE)
        cleanMatches[snp.snpID]=clean_list
        quietDic[snp.snpID]=quiet_list

    return cleanMatches, quietDic

def isStopCodon (seq5, x, seq3, readingFrame):
    #take one reading frame containing the mutation
    totalSequence=seq5+x+seq3
    start=len(seq5)
    if readingFrame=="1":
        codon=totalSequence[start:start+3]
    elif readingFrame=="2":
        codon=totalSequence[start-1:start+2]
    else:
        codon=totalSequence[start-2:start+1]
    if codon== "TAG" or codon== "TAA" or codon== "TGA":
        return True


def getRevComp(snp):
    len5=len(snp.seq5)
    totalSeq= snp.seq5+snp.mutation+snp.seq3
    total_seq=Seq(totalSeq)
    rc_seq=total_seq.reverse_complement()
    seq5=str(rc_seq[0:len5])
    seq3=str(rc_seq[len5+1:])
    new_snp=BEseqType(snp.snpID,seq5 ,seq3,snp.mutation,snp.wt, snp.readingFrame, snp.aaPosition,snp.geneName,snp.geneID)
    return new_snp


def Main(DB):
    totnum=0
    SNPS,rsltsDic=importSNPS(DB) #parsing cvs file
    matches={}
    matchesMinor={}
    cleanMatchdic = {}
    quietMatchdic={}
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

        #if isStopCodon(snp.seq5, snp.mutation, snp.seq3, snp.readingFrame)==True:
            #print ("mutation in snp %s resulted in a stop codon " % snp.snpID)
        #if isStopCodon(snp.seq5, snp.wt, snp.seq3, snp.readingFrame)==True:
            #print ("wt of snp %s results in a stop codon " % snp.snpID)
    empty_keys = [k for k, v in matches.items() if v == []]
    print("matches:", 27257 - (len(empty_keys)))
    empty_keys = [k for k, v in cleanMatchdic.items() if v == []]
    print ("clean:", 27257-(len(empty_keys)))
    empty_keys = [k for k, v in quietMatchdic.items() if v == []]
    print("quiet:", 27257 - (len(empty_keys)))
    return matches, matchesMinor, cleanMatchdic, quietMatchdic, rsltsDic


matches, matchesMinor, cleanMatchdic, quietMatchdic ,rsltsDic=Main("snps_with_clinvar.csv")
with open('matches_file.csv', mode='w') as f:
    f.write('snpID, matches\n')
    for key in matches.keys():
        f.write("%s,%s\n" % (key, matches[key]))
with open('cleanMatches_file.csv', mode='w') as f:
    f.write('snpID, matches\n')
    for key in cleanMatchdic.keys():
        f.write("%s,%s\n" % (key, cleanMatchdic[key]))
with open('quietMatches_file.csv', mode='w') as f:
    f.write('snpID, matches\n')
    for key in quietMatchdic.keys():
        f.write("%s,%s\n" % (key, quietMatchdic[key]))

with open('full_results.csv', mode='w') as f:
    f.write('snpID, Gene Name, BE1, BE2, BE3,HF-BE3,BE4(max),BE4-Gam,YE1-BE3,YEE-BE3, VQR-BE3,VRER-BE3,SaBE3, SaBE4,SaBE4-Gam, Sa(KKH)-BE3,Cas12a-BE,Target-AID,Target-AID-NG,xBE3,eA3A-BE3,BE-PLUS,ABE 7.9, ABE 7.10,ABE 7.10*,xABE,ABESa,VQR-ABE,VRER-ABE, Sa(KKH)-ABE\n')
    keyList=matches.keys()
    hh=0
    beList=["BE1", "BE2", "BE3","HF-BE3","BE4(MAX)","BE4-Gam","YE1-BE3","YEE-BE3", "VQR-BE3","VRER-BE3","SaBE3", "SaBE4","SaBE4-Gam", "Sa(KKH)-BE3","Cas12a-BE","Target-AID","Target-AID-NG","xBE3","eA3A-BE3","BE-PLUS", "ABE 7.9","ABE 7.10","ABE 7.10*","xABE" ,"ABESa","VQR-ABE","VRER-ABE","Sa(KKH)-ABE"]
    for key in keyList:
        gg=0
        mlist=[]
        for BE in beList:
            temp=[]
            if BE in matches[key]:
                temp.append("Match")
            else:
                temp.append("0")
            if BE in cleanMatchdic[key]:
                temp.append("Clean Match")
            else:
                temp.append("0")
            if BE in quietMatchdic[key]:
                temp.append("Quiet Match")
            else:
                temp.append("0")
            temp=' '.join(temp)
            if temp=='0 0 0' or temp=='Match 0 0':
                gg=gg+1
            mlist.append(temp)
        if gg==28:
            print (key, BE )
            hh=hh+1

        f.write("%s,%s,%s\n" %(key,rsltsDic[key], mlist))
    print(hh)