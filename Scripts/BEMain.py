from Scripts.BEseqType import *
from Bio.Seq import Seq
import csv
from Scripts.baseEditorsTable import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter

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
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        match = True
                        break
            if match is True:
                matches_list.append(BE)

        matches[snp.snpID] = matches_list
    return matches

def cleanMatch(snp, Matches, BElist,rev):
    cleanMatches={}
    clean_list=[]
    quietDic={}
    quiet_list=[]
    locations_dic={}
    for match in Matches:
        for BE in Matches[match]:
            if BElist[BE][5] == "U":
                seq3 = snp.seq3
                seq5=snp.seq5
            else:
                seq3 = snp.seq5[::-1]  # get reverse
                seq5 = snp.seq3[::-1]  # get reverse

            PAM=BElist[BE][0]
            start=BElist[BE][1]-1
            end=BElist[BE][2]-1
            locations=[]  #find locations of PAM. for each location, check if clean
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
                        locations.append(i + start) #maybe j is not neccesary?
                        break
            locations_dic[BE] = locations

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
            if rev==False:
                protein_seq = Seq(totalSeq2).translate()
            else:
                protein_seq = Seq(totalSeq2).reverse_complement().translate()

            #clean_list = []
            #quiet_list = []
            max_num = 0
            protein_match = False
            for loc in locations_dic[BE]:
                protein_match = False
                temp_seq=totalSeq1
                # loc=loc+len(snp.seq5)
                # activation_window = totalSeq1[loc-end-diff:loc-start]
                locFromEnd = len(printSeq3) - loc
                loc = loc + len(printSeq5)
                activation_window = totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - end:len(
                    printSeq3) - locFromEnd + len(printSeq5) - start + 1]
                num = 0  # number of times the variant appears within the activation window
                max_num = 0
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

                # finalSeq=totalSeq1[0:loc-end-diff]+new_AW+totalSeq1[loc-start:]
                # protein_seq_new=Seq(finalSeq).translate()

                if rev==False:
                    finalSeq = Seq(totalSeq1[0:loc - end] + str(new_AW) + totalSeq1[loc - start + 1:])
                else:
                    beginningP = Seq(
                        totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end]).reverse_complement()
                    new_AW = Seq(new_AW).reverse_complement()
                    # endP=Seq(totalSeq1[loc - start:]).reverse_complement()
                    endP = Seq(
                        totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:]).reverse_complement()
                    finalSeq = endP + new_AW + beginningP

                protein_seq_new = finalSeq.translate()
                if max_num==1:
                    clean_list.append(BE)
                if protein_seq==protein_seq_new:
                    protein_match=True

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

def Main(DB):
    SNPS,rsltsDic=importSNPS(DB) #parsing cvs file
    matches={}
    matches2 = {}
    matches3 = {}
    matchesMinor={}
    cleanMatchdic = {}
    quietMatchdic={}
    rev = False
    # sort DNA. First, determine which bases we wish to replace. 4 cases:
    # 1. C to T: use CBE list       2. A to G: use ABE list
    # 3. T to C: switch to reverse complement and use ABE
    # 4. G to A: switch to RC and use CBE
    for snp in SNPS:
        rev = False
        if snp.mutation == "C" and snp.wt == "T":
            BElist = CBElist
            MinorBElist=CBElistMinor
        elif snp.mutation == "A" and snp.wt == "G":
            BElist = ABElist
            MinorBElist = ABElistMinor
        elif snp.mutation == "T" and snp.wt == "C":
            snp=getRevComp(snp)
            rev = True
            snp.mutation = "A"
            snp.wt="G"
            BElist = ABElist
            MinorBElist = ABElistMinor
        elif snp.mutation == "G" and snp.wt == "A":
            snp = getRevComp(snp)
            rev = True
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
            try:
                snp2, snp3 = checkRF(snp)
                check_match2 = matchBE(snp2, BElist)
                matches2.update(check_match2)
                check_match3 = matchBE(snp3, BElist)
                matches3.update(check_match3)
            except:
                g=1


            #check for matches in minor window
            check_match_minor = matchBE(snp, MinorBElist)
            matchesMinor.update(check_match_minor)

            #check for clean match
            clean, quiet= cleanMatch(snp,check_match,BElist,rev)
            try:
                clean2, quiet2 = cleanMatch(snp2, check_match2, BElist)
                clean3, quiet3 = cleanMatch(snp3, check_match3, BElist)
                for key in quiet2:
                    if key not in quiet:
                        quiet[key] = quiet2[key]
                for key in quiet3:
                    if key not in quiet:
                        quiet[key] = quiet3[key]
            except:
                g=1
            cleanMatchdic.update(clean)
            quietMatchdic.update(quiet)

        except:
            g=1
            #print ("Error: %s is not an appropriate snp" %snp.snpID)

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
    return matches, cleanMatchdic, quietMatchdic, rsltsDic


matches, cleanMatchdic, quietMatchdic ,rsltsDic=Main("snps_with_clinvar.csv")
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

with open('snps_with_clinvar.csv',mode='r') as file:
    file_read = csv.reader(file, delimiter=",")
    total_table={}
    for row in file_read:
        table={}
        gene_rsID=row[0]
        table['gene_rsID']=row[0]
        table['genome_assembly']=row[2]
        table['chromosome']=row[3]
        table['orientation']=row[4]
        table['start_position']=row[5]
        table['end_position']=row[6]
        table['xmer']=row[7]
        table['clinical_sig']=row[8]
        table['snp']=row[9]
        table['residue']=row[10]
        table['var_type']=row[11]
        table['gene_name']=row[12]
        table['gene_ID']=row[13]
        table['reading_frame']=row[14]
        table['aaPosition']=row[15]
        table['condition']=str(row[16])
        total_table[gene_rsID]=table




with open('full_results.csv', mode='w') as f:
    f.write('snpID,Genome Assembly,Chromosome,Orientation, Start Position, End Position,5Xmer,SNP,Residue,Var Type, Gene Name, Gene ID, Reading Frame, aaPosition, BE1, BE2, BE3,HF-BE3,BE4(max),BE4-Gam,YE1-BE3,YEE-BE3, VQR-BE3,VRER-BE3,SaBE3, SaBE4,SaBE4-Gam, Sa(KKH)-BE3,Cas12a-BE,Target-AID,Target-AID-NG,xBE3,eA3A-BE3,BE-PLUS,CP-CBEmax variants,evoAPOBEC1-BE4max, evoFERNY-BE4max,evoCDA1-BE4max,ABE 7.9, ABE 7.10,ABE 7.10*,xABE,NG-ABEmax,ABESa,VQR-ABE,VRER-ABE, Sa(KKH)-ABE,CP-ABEmax variants,Clinical Significance,condition\n')
    keyList=matches.keys()
    beList=["BE1", "BE2", "BE3", "HF-BE3", "BE4(max)", "BE4-Gam","YE1-BE3","YEE-BE3", "VQR-BE3","VRER-BE3","SaBE3", "SaBE4", "SaBE4-Gam", "Sa(KKH)-BE3","Cas12a-BE","Target-AID","Target-AID-NG","xBE3","eA3A-BE3","BE-PLUS","CP-CBEmax variants","evoAPOBEC1-BE4max", "evoFERNY-BE4max","evoCDA1-BE4max", "ABE 7.9","ABE 7.10","ABE 7.10*","xABE","NG-ABEmax" ,"ABESa","VQR-ABE","VRER-ABE","Sa(KKH)-ABE","CP-ABEmax variants"]
    for key in keyList:
        print (key)
        mlist=[]
        for BE in beList:
            temp=0
            if BE in quietMatchdic[key]:
                temp="Quiet"
            if BE in cleanMatchdic[key]:
                temp="Clean"

            mlist.append(temp)
        if mlist!=[0]*34:
            f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %(key,total_table[key]['genome_assembly'],total_table[key]['chromosome'],
                               total_table[key]['orientation'],total_table[key]['start_position'],
                               total_table[key]['end_position'],total_table[key]['xmer'],
                                total_table[key]['snp'],
                                total_table[key]['residue'],total_table[key]['var_type'],
                               total_table[key]['gene_name'],total_table[key]['gene_ID'],
                                total_table[key]['reading_frame'],total_table[key]['aaPosition'],mlist,total_table[key]['clinical_sig'],total_table[key]['condition']))
