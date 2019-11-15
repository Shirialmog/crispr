from Scripts.BEseqType import *
from Bio.Seq import Seq
from Scripts.baseEditors import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter
from Scripts.transverse import SpecialCleanMatch,origPro

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
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        match = True
                        break
            if match is True:
                matches_list.append(BE)
    return matches_list

def cleanMatch(snp,Matches, BElist,rev):
    printRevSeq=None
    printRevCorSeq=None
    origMutSeq = ''
    origRevSeq = ''
    clean_dic={}
    quiet_dic={}
    locations_dic = {}
    locFromEndList = []
    locFromEndDic = {}
    originalProtein=origPro(snp,rev)
    print (originalProtein)
    for BE in Matches:
        if BElist[BE][5] == "U":
            seq3 = snp.seq3
            seq5=snp.seq5
        else:
            seq3 = snp.seq5[::-1] # get reverse
            seq5 = snp.seq3[::-1]  # get reverse

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
                    j+=1
                    temp+=1
                else:
                    break
                if j == len(PAM):
                    locations_list.append(i + start)
                    break
        locations_dic[BE]=locations_list

        totalSeq1 = seq5 + snp.mutation + seq3
        totalSeq2 = seq5 + snp.wt + seq3
        printSeq5=seq5 #for results
        printSeq3=seq3  # for results
        diff=0
        if snp.readingFrame=="1":
            print ("lenseq5",len(seq5))
            if len(snp.seq5)%3==1:
                print ("ff")
                totalSeq1=totalSeq1[1:]
                totalSeq2 = totalSeq2[1:]
                printSeq5=printSeq5[1:]
                diff = 1
                print (totalSeq2)
            elif len(snp.seq5)%3==2:
                totalSeq1=totalSeq1[2:]
                totalSeq2 = totalSeq2[2:]
                printSeq5 = printSeq5[2:]
                diff = 2
        elif snp.readingFrame=="2":
            if len(snp.seq5) % 3 == 0:
                totalSeq1 = totalSeq1[2:]
                totalSeq2 = totalSeq2[2:]
                printSeq5 = printSeq5[2:]
                diff = 2
            elif len(snp.seq5) % 3 == 2:
                totalSeq1 = totalSeq1[1:]
                totalSeq2 = totalSeq1[2:]
                printSeq5=printSeq5[1:]
                diff = 1
        elif snp.readingFrame=="3":
            if len(snp.seq5) % 3 == 0:
                totalSeq1 = totalSeq1[1:]
                totalSeq2 = totalSeq2[1:]
                printSeq5=printSeq5[1:]
                diff = 1
            elif len(snp.seq5) % 3 == 1:
                totalSeq1 = totalSeq1[2:]
                totalSeq2 = totalSeq2[2:]
                printSeq5=printSeq5[2:]
                diff = 2
        if len(totalSeq1)%3==1:
            totalSeq1=totalSeq1[:-1]
            totalSeq2 = totalSeq2[:-1]
            printSeq3 = printSeq3[:-1]
        if len(totalSeq1)%3==2:
            totalSeq1=totalSeq1[:-2]
            totalSeq2 = totalSeq2[:-2]
            printSeq3 = printSeq3[:-2]
        print ("totalseq", totalSeq2)
        totalSeq2=Seq(totalSeq2)
        if rev==False:
            protein_seq=totalSeq2.translate()
        else:
            protein_seq=totalSeq2.reverse_complement().translate()

        clean_list = []
        quiet_list=[]
        origMutSeq = {}
        printPam = {}
        printRevCorSeq = {}
        finalShowSeq = {}
        for loc in locations_dic[BE]:
            protein_match = False
            locFromEnd=len(printSeq3)-loc
            loc=loc+len(printSeq5)
            activation_window=totalSeq1[len(printSeq3)-locFromEnd+len(printSeq5)-end:len(printSeq3)-locFromEnd+len(printSeq5)-start+1]
            num = 0  # number of times the variant appears within the activation window
            max_num=0
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


            if rev==False:
                new = new_AW
                finalSeq = Seq(totalSeq1[0:loc - end] + str(new) + totalSeq1[loc - start+1:])
                beginningP = Seq(
                    totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end])
                endP = Seq(totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:])
                oldShowSeq = totalSeq1[0:loc - end - diff] + "<b>" + activation_window + "</b>" + totalSeq1[loc - start:]
                finalShowSeq[loc]=totalSeq1[0:loc-end-diff]+"<b>"+new_AW+"</b>"+totalSeq1[loc-start:]
                origMutSeq[loc] = printSeq5 + "<b>"+snp.mutation+"</b>" + printSeq3
                origRevSeq=''
                #printPam[loc] = printPamSeq(printSeq5, new_AW, printSeq3, locFromEnd, PAM)
                printPam[loc] = printPamSeq(beginningP, new_AW, endP, locFromEnd, PAM)
            else:
                finalSeq2 = totalSeq1[0:loc - end] +new_AW + totalSeq1[loc - start+1:]
                #beginningP=Seq(totalSeq1[0:loc - end - diff]).reverse_complement()
                beginningP=Seq(totalSeq1[0:len(printSeq3)-locFromEnd+len(printSeq5)-end]).reverse_complement()
                old_AW=Seq(activation_window).reverse_complement()
                new_AW=Seq(new_AW).reverse_complement()
                #endP=Seq(totalSeq1[loc - start:]).reverse_complement()
                endP=Seq(totalSeq1[len(printSeq3)-locFromEnd+len(printSeq5)-start+1:]).reverse_complement()
                finalSeq = endP + new_AW + beginningP
                oldShowSeq = endP + "<b>" + old_AW + "</b>" + beginningP
                finalShowSeq[loc] = endP + "<class style='color:blue'><b>" + new_AW + "</b></class>" + beginningP
                origMutSeq[loc] = Seq(printSeq3).reverse_complement() + "<b>" + Seq(snp.mutation).reverse_complement() + "</b>" + Seq(printSeq5).reverse_complement()
                origRevSeq=printSeq5+"<b>"+snp.mutation+"</b>"+printSeq3
                printRevCorSeq[loc]= printPamSeq(beginningP.reverse_complement(), new_AW.reverse_complement(), endP.reverse_complement(), locFromEnd, PAM)
                printPam[loc]= printPamSeq(beginningP.reverse_complement(), old_AW.reverse_complement(), endP.reverse_complement(), locFromEnd, PAM)
            protein_seq_new=finalSeq.translate()
            if protein_seq==protein_seq_new:
                protein_match=True
            print (protein_seq_new,protein_seq)
            print (finalSeq)
            if max_num==1:
                clean_list.append(loc)
            if protein_match == True:
                quiet_list.append(loc)
        print ("list:", quiet_list,BE)
        clean_dic[BE]=[origMutSeq,printPam,printRevCorSeq,finalShowSeq,PAM,clean_list,rev]

        #quiet_list=quietCheck(snp,seq5,seq3,locations_dic,rev,BE,start,end)
        quiet_dic[BE]=[origMutSeq,printPam,printRevCorSeq,finalShowSeq,PAM, quiet_list,rev]

    return clean_dic,quiet_dic,origMutSeq,origRevSeq,locations_dic,originalProtein



def printPamSeq(seq5,activationWindow,seq3, locfromEnd,PAM):
    len3=len(seq3)
    printPAM=seq5+"<class style='color:blue'><b>"+activationWindow+"</b></class>"+seq3[0:len3-locfromEnd]+"<b><class style='color:#DD96F0'>"+seq3[len3-locfromEnd:len3-locfromEnd+len(PAM)]+"</b></class>"+seq3[len3-locfromEnd+len(PAM):]
    return printPAM

def getRevComp(snp):
    len5=len(snp.seq5)
    len3 = len(snp.seq3)
    totalSeq= snp.seq5+snp.mutation+snp.seq3
    total_seq=Seq(totalSeq)
    rc_seq=total_seq.reverse_complement()
    seq5=str(rc_seq[0:len3])
    start=len(rc_seq)-len5
    seq3=str(rc_seq[start:])
    new_snp=BEseqType(snp.snpID,seq5 ,seq3,snp.mutation,snp.wt, snp.readingFrame, snp.aaPosition,snp.geneName,snp.geneID)
    return new_snp


def checkRF(snp):
    #this function will check the other 2 bases in the reading frame to see whether fixing them may result in a quiet result
    if snp.readingFrame==1:
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(12, zero_seq5, zero_seq3, zero_wt, zero_mutation, 1, 1, 1, 1)
        first_seq5=snp.seq5+snp.mutation
        first_mutation=snp.seq3[0]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq3[1:]
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,2,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation+snp.seq3[0]
        second_mutation = snp.seq3[1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[2:]
        snp2 = BEseqType(12, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp0, snp1,snp2

    elif snp.readingFrame==2:
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(12, zero_seq5, zero_seq3, zero_wt, zero_mutation, 2, 1, 1, 1)
        first_seq5=snp.seq5[:-1]
        first_mutation=snp.seq5[-1]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.mutation+snp.seq3
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation
        second_mutation = snp.seq3[0]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[1:]
        snp2 = BEseqType(12, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp0,snp1,snp2

    # elif snp.readingFrame==3:
    else:
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(12, zero_seq5, zero_seq3, zero_wt, zero_mutation, 3, 1, 1, 1)
        first_seq5=snp.seq5[:-2]
        first_mutation=snp.seq5[-2]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq5[-1]+snp.mutation+snp.seq3
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5[:-1]
        second_mutation = snp.seq5[-1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.mutation+snp.seq3
        snp2 = BEseqType(12, second_seq5, second_seq3,second_wt, second_mutation,2,1,1,1)
        return snp0, snp1,snp2

def find_cor(base):
    if base=="C":
        return "T"
    if base=="T":
        return "C"
    if base=="A":
        return "G"
    if base=="G":
        return "A"

def MainBE(upSeq, downSeq, mutation, wt, readingFrame=2):
    upSeq=upSeq.upper()
    downSeq=downSeq.upper()
    mutation=mutation.upper()
    wt=wt.upper()
    snp=BEseqType(1234,upSeq, downSeq, wt, mutation, readingFrame, 1, 12, 12)
    refSeq=snp.seq5+"<b>"+snp.wt+"</b>"+snp.seq3
    mutSeq = snp.seq5 + "<b>" + snp.mutation + "</b>" + snp.seq3
    rev,rev0,rev2,rev3=False,False,False,False
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
        BElist=[]

    snp0,snp2,snp3=checkRF(snp)
    if snp0.mutation == "C":
        BElist0 = CBElist
    elif snp0.mutation == "A":
        BElist0 = ABElist
    elif snp0.mutation == "T":
        snp0 = getRevComp(snp0)
        rev0 = True
        snp0.mutation = "A"
        snp0.wt = "G"
        BElist0 = ABElist
    elif snp0.mutation == "G":
        snp0 = getRevComp(snp0)
        rev0 = True
        snp0.mutation = "C"
        snp0.wt = "T"
        BElist0 = CBElist
    if snp2.mutation == "C":
        BElist2 = CBElist
    elif snp2.mutation == "A":
        BElist2 = ABElist
    elif snp2.mutation == "T":
        snp2 = getRevComp(snp2)
        rev2 = True
        snp2.mutation = "A"
        snp2.wt = "G"
        BElist2 = ABElist
    elif snp2.mutation == "G":
        snp2 = getRevComp(snp2)
        rev2 = True
        snp2.mutation = "C"
        snp2.wt = "T"
        BElist2 = CBElist

    if snp3.mutation == "C":
        BElist3 = CBElist
    elif snp3.mutation == "A":
        BElist3 = ABElist
    elif snp3.mutation =="T":
        snp3 = getRevComp(snp3)
        rev3 = True
        snp3.mutation = "A"
        snp3.wt = "G"
        BElist3 = ABElist
    elif snp3.mutation == "G":
        snp3 = getRevComp(snp3)
        rev3 = True
        snp3.mutation = "C"
        snp3.wt = "T"
        BElist3 = CBElist

    check_match = matchBE(snp, BElist)
    print ("match:", check_match)
    try:
        check_match0 = matchBE(snp0, BElist0)
        print ("check_match0", check_match0)
    except:
        check_match0=[]
        BElist0=[]
    try:
        check_match2=matchBE(snp2,BElist2)
        print("check_match2", check_match2)
    except:
        BElist2=[]
        check_match2=[]
    try:
        check_match3=matchBE(snp3,BElist3)
        print ("check_match3", check_match3)
    except:
        BElist3=[]
        check_match3=[]
    #check for clean match
    clean_dic,quiet_dic,origMutSeq, origRevSeq,locations_dic,originalProtein= cleanMatch(snp, check_match, BElist,rev)
    # clean_dic0, quiet_dic0, origMutSeq0, origRevSeq0, locations_dic0 = SpecialCleanMatch(snp0, check_match0, BElist0, rev0,originalProtein)
    # clean_dic2, quiet_dic2, origMutSeq2, origRevSeq2, locations_dic2 = SpecialCleanMatch(snp2, check_match2, BElist2, rev2,originalProtein)
    # clean_dic3, quiet_dic3, origMutSeq3, origRevSeq3, locations_dic3 = SpecialCleanMatch(snp3, check_match3, BElist3, rev3,originalProtein)
    # for key in quiet_dic0:
    #     if key not in quiet_dic:
    #         quiet_dic[key]=quiet_dic0[key]
    # for key in quiet_dic2:
    #     if key not in quiet_dic:
    #         quiet_dic[key]=quiet_dic2[key]
    # for key in quiet_dic3:
    #     if key not in quiet_dic:
    #         quiet_dic[key]=quiet_dic3[key]
    # except:
    #     return [],[],[],[],[],[],[]
    origMutSeq='1'
    origRevSeq='1'
    try:
        return clean_dic, quiet_dic,refSeq,mutSeq, origMutSeq, origRevSeq,locations_dic
    except:
        return [],[],[],[],[],[],[]


