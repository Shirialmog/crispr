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

def cleanMatch(snp, Matches, BElist,rev):
    printRevSeq=None
    printRevCorSeq=None
    clean_list=[]
    clean_dic={}
    quiet_list=[]
    quiet_dic={}
    locations_dic = {}
    locFromEndList = []
    locFromEndDic = {}
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
        locations_list=[]

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
                    locations_list.append(i + start)
                    break
        locations_dic[BE]=locations_list

        totalSeq1=seq5+snp.mutation+seq3
        printSeq5=seq5 #for results
        printSeq3 = seq3  # for results
        diff=0
        if snp.readingFrame=="1":
            if len(snp.seq5)%3==1:
                totalSeq1=totalSeq1[1:]
                printSeq5=printSeq5[1:]
                diff = 1
            elif len(snp.seq5)%3==2:
                totalSeq1=totalSeq1[2:]
                printSeq5 = printSeq5[2:]
                diff = 2
        elif snp.readingFrame=="2":
            if len(snp.seq5) % 3 == 0:
                totalSeq1 = totalSeq1[2:]
                printSeq5 = printSeq5[2:]
                diff = 2
            elif len(snp.seq5) % 3 == 2:
                totalSeq1 = totalSeq1[1:]
                printSeq5=printSeq5[1:]
                diff = 1
        elif snp.readingFrame=="3":
            if len(snp.seq5) % 3 == 0:
                totalSeq1 = totalSeq1[1:]
                printSeq5=printSeq5[1:]
                diff = 1
            elif len(snp.seq5) % 3 == 1:
                totalSeq1 = totalSeq1[2:]
                printSeq5=printSeq5[2:]
                diff = 2
        if len(totalSeq1)%3==1:
            totalSeq1=totalSeq1[:-1]
            printSeq3 = printSeq3[:-1]
        if len(totalSeq1)%3==2:
            totalSeq1=totalSeq1[:-2]
            printSeq3 = printSeq3[:-2]
        totalSeq2=Seq(totalSeq1)

        protein_seq=totalSeq2.translate()
        max_num = 0
        protein_match=False

        for loc in locations_dic[BE]:
            locFromEndList.append(len(printSeq3-loc))
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


            if rev==False:
                oldShowSeq = totalSeq1[0:loc - end - diff] + "[" + activation_window + "]" + totalSeq1[loc - start:]
                finalShowSeq=totalSeq1[0:loc-end-diff]+"["+new_AW+"]"+totalSeq1[loc-start:]
                origMutSeq = printSeq5 + "["+snp.mutation+"]" + printSeq3
                origRevSeq=''
                printPam = printPamSeq(printSeq5, snp.mutation, printSeq3, locations_dic, BE, PAM)
            else:
                beginningP=Seq(totalSeq1[0:loc - end - diff]).reverse_complement()
                old_AW=Seq(activation_window).reverse_complement()
                new_AW=Seq(new_AW).reverse_complement()
                endP=Seq(totalSeq1[loc - start:]).reverse_complement()
                oldShowSeq = endP + "<b>" + old_AW + "</b>" + beginningP
                finalShowSeq = endP + "<b>" + new_AW + "</b>" + beginningP
                origMutSeq = Seq(printSeq3).reverse_complement() + "<b>" + Seq(snp.mutation).reverse_complement() + "</b>" + Seq(printSeq5).reverse_complement()
                origRevSeq=printSeq5+"<b>"+snp.mutation+"</b>"+printSeq3
                printRevCorSeq= printPamSeq(beginningP.reverse_complement(), new_AW.reverse_complement(), endP.reverse_complement(), locations_dic, BE, PAM)
                printPam = printPamSeq(beginningP.reverse_complement(), old_AW.reverse_complement(), endP.reverse_complement(), locations_dic, BE, PAM)
            protein_seq_new=Seq(finalSeq).translate()

            if protein_seq==protein_seq_new:
                protein_match=True

        if max_num==1:
                clean_list.append(BE)
                clean_dic[BE]=[finalShowSeq,PAM]

        if protein_match==True:
                quiet_list.append(BE)
                quiet_dic[BE]=[finalShowSeq,PAM]





    return clean_dic,quiet_dic,origMutSeq,origRevSeq,printPam,printRevCorSeq

def printPamSeq(seq5,activationWindow,seq3, locations_dic,BE,PAM):
    for loc in locations_dic[BE]:
        locEnd=len(seq3)-loc
        printPAM=seq5+"<b>"+activationWindow+"</b>"+seq3[0:locEnd]+"<b><class style='color:#DD96F0'>"+seq3[locEnd:locEnd+len(PAM)]+"</b></class>"+seq3[locEnd+len(PAM):]
    return printPAM

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
    refSeq=snp.seq5+"<b>"+snp.wt+"</b>"+snp.seq3
    mutSeq = snp.seq5 + "<b>" + snp.mutation + "</b>" + snp.seq3
    rev=False
    if snp.mutation == "C" and snp.wt == "T":
        BElist = CBElist
        MinorBElist=CBElistMinor
    elif snp.mutation == "A" and snp.wt == "G":
        BElist = ABElist
        MinorBElist = ABElistMinor
    elif snp.mutation == "T" and snp.wt == "C":
        snp=getRevComp(snp)
        rev=True
        snp.mutation = "A"
        snp.wt="G"
        BElist = ABElist
        MinorBElist = ABElistMinor
    elif snp.mutation == "G" and snp.wt == "A":
        snp = getRevComp(snp)
        rev=True
        snp.mutation = "C"
        snp.wt="T"
        BElist = CBElist
        MinorBElist = CBElistMinor
    else:
        return ("Cannot be fixed with this tool")


    # check for matches in major window
    check_match = matchBE(snp, BElist)

    # check for clean match
    clean_dic,quiet_dic,origMutSeq, origRevSeq, printPam, printRevCorSeq= cleanMatch(snp, check_match, BElist,rev)

    return clean_dic, quiet_dic,refSeq,mutSeq, origMutSeq, origRevSeq, printPam, printRevCorSeq


