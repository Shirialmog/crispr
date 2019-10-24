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
    clean_dic={}
    quiet_dic={}
    locations_dic = {}
    locFromEndList = []
    locFromEndDic = {}
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
            if len(snp.seq5)%3==1:
                totalSeq1=totalSeq1[1:]
                totalSeq2 = totalSeq2[1:]
                printSeq5=printSeq5[1:]
                diff = 1
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
                print (finalSeq2)
                #beginningP=Seq(totalSeq1[0:loc - end - diff]).reverse_complement()
                beginningP=Seq(totalSeq1[0:len(printSeq3)-locFromEnd+len(printSeq5)-end]).reverse_complement()
                old_AW=Seq(activation_window).reverse_complement()
                new_AW=Seq(new_AW).reverse_complement()
                #endP=Seq(totalSeq1[loc - start:]).reverse_complement()
                endP=Seq(totalSeq1[len(printSeq3)-locFromEnd+len(printSeq5)-start+1:]).reverse_complement()
                finalSeq = endP + new_AW + beginningP
                print (finalSeq)
                oldShowSeq = endP + "<b>" + old_AW + "</b>" + beginningP
                finalShowSeq[loc] = endP + "<class style='color:blue'><b>" + new_AW + "</b></class>" + beginningP
                origMutSeq[loc] = Seq(printSeq3).reverse_complement() + "<b>" + Seq(snp.mutation).reverse_complement() + "</b>" + Seq(printSeq5).reverse_complement()
                origRevSeq=printSeq5+"<b>"+snp.mutation+"</b>"+printSeq3
                printRevCorSeq[loc]= printPamSeq(beginningP.reverse_complement(), new_AW.reverse_complement(), endP.reverse_complement(), locFromEnd, PAM)
                printPam[loc]= printPamSeq(beginningP.reverse_complement(), old_AW.reverse_complement(), endP.reverse_complement(), locFromEnd, PAM)
            print ("f:",finalSeq)
            protein_seq_new=finalSeq.translate()
            if protein_seq==protein_seq_new:
                protein_match=True

            if max_num==1:
                clean_list.append(loc)
            if protein_match == True:
                quiet_list.append(loc)
        clean_dic[BE]=[origMutSeq,printPam,printRevCorSeq,finalShowSeq,PAM,clean_list,rev]

        quiet_list=quietCheck(snp,seq5,seq3,locations_dic,rev,BE,start,end)
        quiet_dic[BE]=[origMutSeq,printPam,printRevCorSeq,finalShowSeq,PAM, quiet_list,rev]

    return clean_dic,quiet_dic,origMutSeq,origRevSeq,locations_dic

def quietCheck(snp,seq5,seq3,locations_dic,rev,BE,start,end):
    quiet_list = []
    totalSeq1 = seq5 + snp.mutation + seq3
    totalSeq2 = seq5 + snp.wt + seq3
    printSeq5 = seq5  # for results
    printSeq3 = seq3  # for results
    diff = 0
    if snp.readingFrame == "1":
        if len(snp.seq5) % 3 == 1:
            totalSeq1 = totalSeq1[1:]
            totalSeq2 = totalSeq2[1:]
            printSeq5 = printSeq5[1:]
            diff = 1
        elif len(snp.seq5) % 3 == 2:
            totalSeq1 = totalSeq1[2:]
            totalSeq2 = totalSeq2[2:]
            printSeq5 = printSeq5[2:]
            diff = 2
    elif snp.readingFrame == "2":
        if len(snp.seq5) % 3 == 0:
            totalSeq1 = totalSeq1[2:]
            totalSeq2 = totalSeq2[2:]
            printSeq5 = printSeq5[2:]
            diff = 2
        elif len(snp.seq5) % 3 == 2:
            totalSeq1 = totalSeq1[1:]
            totalSeq2 = totalSeq1[2:]
            printSeq5 = printSeq5[1:]
            diff = 1
    elif snp.readingFrame == "3":
        if len(snp.seq5) % 3 == 0:
            totalSeq1 = totalSeq1[1:]
            totalSeq2 = totalSeq2[1:]
            printSeq5 = printSeq5[1:]
            diff = 1
        elif len(snp.seq5) % 3 == 1:
            totalSeq1 = totalSeq1[2:]
            totalSeq2 = totalSeq2[2:]
            printSeq5 = printSeq5[2:]
            diff = 2
    if len(totalSeq1) % 3 == 1:
        totalSeq1 = totalSeq1[:-1]
        totalSeq2 = totalSeq2[:-1]
        printSeq3 = printSeq3[:-1]
    if len(totalSeq1) % 3 == 2:
        totalSeq1 = totalSeq1[:-2]
        totalSeq2 = totalSeq2[:-2]
        printSeq3 = printSeq3[:-2]
    totalSeq2 = Seq(totalSeq2)


    for loc in locations_dic[BE]:
        protein_match = False
        locFromEnd = len(printSeq3) - loc
        loc = loc + len(printSeq5)
        activation_window = totalSeq1[
                            len(printSeq3) - locFromEnd + len(printSeq5) - end:len(printSeq3) - locFromEnd + len(
                                printSeq5) - start + 1]

        new_AW = []
        for i in range(len(activation_window)):
            if activation_window[i] == snp.mutation:
                new_AW.append(snp.wt)
            else:
                new_AW.append(activation_window[i])
        new_AW = ''.join(new_AW)


        if rev == False:
            new = new_AW
            finalSeq = Seq(totalSeq1[0:loc - end] + str(new) + totalSeq1[loc - start + 1:])
            protein_seq = totalSeq2.translate()
        else:
            beginningP = Seq(totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end]).reverse_complement()
            old_AW = Seq(activation_window).reverse_complement()
            new_AW = Seq(new_AW).reverse_complement()
            # endP=Seq(totalSeq1[loc - start:]).reverse_complement()
            endP = Seq(totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:]).reverse_complement()
            finalSeq = endP + new_AW + beginningP
            protein_seq = totalSeq2.reverse_complement().translate()

        protein_seq_new = finalSeq.translate()
        if protein_seq == protein_seq_new:
            protein_match = True
        if protein_match == True:
            quiet_list.append(loc)

    return quiet_list

def printPamSeq(seq5,activationWindow,seq3, locfromEnd,PAM):
    len3=len(seq3)
    printPAM=seq5+"<class style='color:blue'><b>"+activationWindow+"</b></class>"+seq3[0:len3-locfromEnd]+"<b><class style='color:#DD96F0'>"+seq3[len3-locfromEnd:len3-locfromEnd+len(PAM)]+"</b></class>"+seq3[len3-locfromEnd+len(PAM):]
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


def checkRF(snp):
    #this function will check the other 2 bases in the reading frame to see whether fixing them may result in a quiet result
    if snp.readingFrame==1:
        first_seq5=snp.seq5+snp.mutation
        first_mutation=snp.seq3[0]
        first_wt=str(Seq(snp.mutation).reverse_complement())
        first_seq3=snp.seq3[1:]
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,2,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation+snp.seq3[0]
        second_mutation = snp.seq3[1]
        second_wt = str(Seq(snp.mutation).reverse_complement())
        second_seq3 = snp.seq3[2:]
        snp2 = BEseqType(12, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp1,snp2

    elif snp.readingFrame==2:
        first_seq5=snp.seq5[:-1]
        first_mutation=snp.seq5[-1]
        first_wt=str(Seq(snp.mutation).reverse_complement())
        first_seq3=snp.mutation+snp.seq3
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation
        second_mutation = snp.seq3[0]
        second_wt = str(Seq(snp.mutation).reverse_complement())
        second_seq3 = snp.seq3[1:]
        snp2 = BEseqType(12, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp1,snp2

    # elif snp.readingFrame==3:
    else:
        first_seq5=snp.seq5[:-2]
        first_mutation=snp.seq5[-2]
        first_wt=str(Seq(snp.mutation).reverse_complement())
        first_seq3=snp.seq5[-1]+snp.mutation+snp.seq3
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5[:-1]
        second_mutation = snp.seq3[-1]
        second_wt = str(Seq(snp.mutation).reverse_complement())
        second_seq3 = snp.mutation+snp.seq3
        snp2 = BEseqType(12, second_seq5, second_seq3,second_wt, second_mutation,2,1,1,1)
        return snp1,snp2


def MainBE(upSeq, downSeq, mutation, wt, readingFrame=2):
    upSeq=upSeq.upper()
    downSeq=downSeq.upper()
    mutation=mutation.upper()
    wt=wt.upper()
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

    snp2,snp3=checkRF(snp)
    print (2,3)
    if snp2.mutation == "C":
        BElist2 = CBElist
    elif snp2.mutation == "A":
        BElist2 = ABElist
    elif snp2.mutation == "T":
        snp2 = getRevComp(snp)
        rev = True
        snp2.mutation = "A"
        snp2.wt = "G"
        BElist2 = ABElist
    elif snp2.mutation == "G":
        snp2 = getRevComp(snp)
        rev = True
        snp2.mutation = "C"
        snp2.wt = "T"
        BElist2 = CBElist

    if snp3.mutation == "C":
        BElist3 = CBElist
    elif snp3.mutation == "A":
        BElist3 = ABElist
    elif snp3.mutation == "T":
        snp3 = getRevComp(snp)
        rev = True
        snp3.mutation = "A"
        snp3.wt = "G"
        BElist3 = ABElist
    elif snp3.mutation == "G":
        snp3 = getRevComp(snp)
        rev = True
        snp3.mutation = "C"
        snp3.wt = "T"
        BElist3 = CBElist
    try:
        #check for matches in major window
        check_match = matchBE(snp, BElist)
        check_match2=matchBE(snp2,BElist2)
        check_match3=matchBE(snp3,BElist3)
        #check for clean match
        clean_dic,quiet_dic,origMutSeq, origRevSeq,locations_dic= cleanMatch(snp, check_match, BElist,rev)
        clean_dic2, quiet_dic2, origMutSeq2, origRevSeq2, locations_dic2 = cleanMatch(snp2, check_match2, BElist, rev)
        clean_dic3, quiet_dic3, origMutSeq3, origRevSeq3, locations_dic3 = cleanMatch(snp3, check_match3, BElist, rev)
        for key in quiet_dic2:
            if key not in quiet_dic:
                quiet_dic[key]=quiet_dic2[key]
        for key in quiet_dic3:
            if key not in quiet_dic:
                quiet_dic[key]=quiet_dic3[key]
    except:
        return [],[],[],[],[],[],[]
    return clean_dic, quiet_dic,refSeq,mutSeq, origMutSeq, origRevSeq,locations_dic


