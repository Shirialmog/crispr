from Scripts.BEseqType import *
from Bio.Seq import Seq
from Scripts.baseEditors import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter

def matchBE(snp, BElist):
    matches_list = []
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
    return matches_list

def cleanMatch(snp, Matches, BElist):
    clean_list=[]
    quiet_list=[]
    for BE in Matches:
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


    return clean_list, quiet_list

def getRevComp(snp):
    len5=len(snp.seq5)
    totalSeq= snp.seq5+snp.mutation+snp.seq3
    total_seq=Seq(totalSeq)
    rc_seq=total_seq.reverse_complement()
    seq5=str(rc_seq[0:len5])
    seq3=str(rc_seq[len5+1:])
    new_snp=BEseqType(snp.snpID,seq5 ,seq3,snp.mutation,snp.wt, snp.readingFrame, snp.aaPosition,snp.geneName,snp.geneID)
    return new_snp

def MainBE(upSeq, downSeq, mutation, wt, readingFrame=2):
    snp=BEseqType(1234,upSeq, downSeq, wt, mutation, readingFrame, 1, 12, 12)
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
        return ("BAD")


    # check for matches in major window
    check_match = matchBE(snp, BElist)

    # check for matches in minor window
    check_match_minor = matchBE(snp, MinorBElist)

    # check for clean match
    clean, quiet = cleanMatch(snp, check_match, BElist)

    return ', '.join(check_match), ', '.join(clean), ', '.join(quiet)


