#Modules:
from Scripts.casDBFile_new import casDB, CASletter, CASnames,casDic #Scripts.
from Scripts.seqType_main import *
from Scripts.seqTools import *
#Libraries
import xml.etree.ElementTree as ET
import csv
import time

DBLOCATION = "snpDB.xml" #SNP DB location

#statistics global variables
badParse = 0
AVERAGERSLTS = 0
totalSequances = 0


def importSNPS(SNPfile):
    global badParse
    global totalSequances
    rslts=[]
    with open(SNPfile,"r") as csv_file:
        try:
            csv_read = csv.reader(csv_file, delimiter=",")
        except MemoryError:
            print("MemoryError, try splitting your file into smaller files")
            exit()
        except Exception as e:
            print (e)
        for row in csv_read:
            try:
                if len(row[9])>3:
                    continue
                snpID=row[0]
                mutation=row[9][2]
                wt=row[9][0]
                sequence5 = row[7][0:25]
                sequence3 = row[7][26:]
                readingFrame=row[14]
                aaPosition=row[15]
                geneName=row[12]
                geneID=row[13]
                genechromo=row[2]
                repl=""
                orientation=row[4]
                mutatype=''
                fromto="From "+wt+" to "+ mutation
                clinical=(row[8]).split(',')
                clinical=('.').join(clinical)
                condition=(row[16]).split(',')
                condition=('.').join(condition)
                # newSNP=seqType(snpID,geneName,geneID,sequence5, sequence3, wt, mutation, genechromo,fromto,orientation,
                #                mutatype,clinical,condition,readingFrame=readingFrame)
                newSNP = seqType(snpID, geneName, geneID, sequence5, sequence3, wt, mutation, genechromo, fromto,
                                 orientation,
                                 mutatype, clinical, condition,aaPosition=aaPosition, readingFrame=readingFrame)
                rslts.append(newSNP)
                totalSequances+=1
            except:
                badParse+=1

        return rslts
def importDB():
    """
    opens DB to extract XML data using DBLOCATION constant
    :return: XML data from file, type: xml.etree.ElementTree
    """
    with open(DBLOCATION,"r",encoding='utf-8') as f:
        dbRawData = ""
        print("loading DB")
        try:
            dbRawData = f.read()
        except MemoryError:
            print("Oh no! we got a MemoryError! try splitting your file into smaller files")
            exit()
        print("Finished Successfully")
        print("Parsing XML")
        dbXmlRoot = ET.fromstring(dbRawData)
        print("Finished XML Parsing")
        return dbXmlRoot


def parse(DB):
    """
    Gets XML db, returns seqType DB
    :param DB: XML elementTree root
    :return: rslt dict with SeqType variables
    """
    global badParse
    global totalSequances
    rsltDB = []
    print("Parsing DB")
    for seq in DB:
        totalSequances += 1
        try:
            orientation,mutatype,geneId,geneName,clinical,chromosome,aaPosition,readingFrame,residue = "NA","NA","NA","NA","NA","NA","NA","NA","NA"
            rsId= seq.get("rsId")
            rsId= seq.get("rsId")
            for element in seq:
                if element.tag == "Sequence":
                    seq5 = element.find("Seq3").text.upper() #sequence
                    seq3 = element.find("Seq5").text.upper()
                    variants = element.find("Observed").text.upper().split("/") #variants
                if element.tag == "Assembly":
                    refallele = element.find("Component").find("MapLoc").get("refAllele").upper()
                    try:
                        geneId = element.find("Component").find("MapLoc").find("FxnSet").get("geneId")  # geneId
                        geneName = element.find("Component").find("MapLoc").find("FxnSet").get("symbol")  # geneName
                        clinical = seq.find("Phenotype").find("ClinicalSignificance").text  # <<<<<>>>>> add to field
                        chromosome = element.find("Component").get("chromosome")
                        orientation = element.find("Component").get("orientation")
                        aaPosition = element.find("Component").find("MapLoc").find("FxnSet").get("aaPosition")
                        readingFrame = element.find("Component").find("MapLoc").find("FxnSet").get("readingFrame")
                        residue = element.find("Component").find("MapLoc").find("FxnSet").get("residue")
                        for set in element.find("Component").findall("MapLoc"):
                            if set.find("FxnSet").get("fxnClass") != "reference":
                                mutatype = set.find("FxnSet").get("fxnClass")  #mutation type
                    except Exception:
                        True #skipping, non critical data
            mutations = variants.copy()
            try:
                mutations.remove(refallele)
            except ValueError: #refallele not in variants
                True
            try:
                aaPosition = str(int(aaPosition) - 1)
            except Exception:
                True
            fromto = ""
            for mutation in mutations:
                fromto += "from "+refallele+" to "+mutation+";"
            seq3 = seq3[-25:] #cutting sequences to minimum to avoid as many illegal nucs as possible
            seq5 = seq5[:25]
            if legalSeq(seq3+refallele+seq5+"".join(variants),"-AGTC") and (60 > len(seq3+refallele+seq5) > 40): #sequence size limitations
                a = seqType(rsId, geneId, geneName, seq3, seq5, refallele, mutations, chromosome, fromto, orientation, mutatype, clinical,aaPosition=aaPosition,readingFrame=readingFrame,residue=residue)
                rsltDB.append(a)
            else:
                badParse += 1
            #     print("illegal sequence:",seq3+"<"+refallele+">"+seq5)
        except AttributeError as e:
            badParse += 1
            print(e, rsId)
    print("Parsing finished successfully")
    print("Initial parsing, ",len(rsltDB)," sequences were loaded")
    return rsltDB

def match(CAS,sequence,seq3len,rtrnLoc=False,specific=False):
    """
    gets full Sequence and CAS's PAM, returns Match if the PAM matches the Variant somewhere
    if rtrnLoc = True, returns location of the match instead of True, or False if none found
    :param CAS: String
    :param Variant: String
    :param specific: Boolean - check only the specific location given
    :return: if rtrnLoc -> list with int locations; else Boolean
    """
    global snps_with_pam
    window = len(CAS)
    movement = window
    rslt = []
    start = seq3len - window + 1
    if specific == True: #don't move
        movement = 1
        start = seq3len
    for i in range(movement): #window size options
        match = True #flag
        for j in range(window): #compare single chars in window
            try:
                if sequence[start+i+j] not in CASletter[CAS[j]]:
                    match = False
                    break
            except Exception as e:
                match = False
                print("Corrupt Data, skipping variant",e)
        if match:
            if rtrnLoc:
                rslt.append(start + i)
            else:
                return True
    if rtrnLoc:
        return rslt
    else:
        return False

def rank(DB):
    """
    <<<<UNAVAILABLE FUNCTION>>>>>
    Loops on rslt DB and adds rank using external ranking module
    :param DB: seqType database
    :return: None
    """
    return DB

def parseRslts(rslts):
    """
    gets seq.rslts returns string ready for csv
    :param rslts:
    :return: String
    """
    parsed = ["" for x in range(len(casDB))] #return variable for function template
    for column in range(len(casDB)):
        for rslt in rslts: #strings
            if casDB[column] in rslt[1].split(" ") != -1:
                if rslt[1].find("comp") != -1:
                    parsed[column] += rslt[0] + " comp; "
                else:
                    parsed[column] += rslt[0] + "; "
    return parsed

def exportSeq(DB):
    """
    outputs a seqTypeDB to a CSV file Using CSV module, e.g. used to show which sequences were tested
    :param DB: {[key - string]:(seqType,[CAS1,CAS2...CASn],...)}
    :return: None, OUTput is a CSV file in the project dir
    """
    print("exporting Imported sequences to file")
    try:
        with open('AnalysedList.csv', 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            # Write column titles
            spamwriter.writerow(["srID","geneID","geneName","Clinical Significance","WT sequance","Variants","chromosome"])
            #Build data rows
            rows = []
            print("Building rows")
            for seq in DB: #data
                if seq.geneName == None:
                    seq.geneName = "NA"
                #',' in csv seperates columns
                stringified = str(seq.snpID+','+seq.geneID+','+seq.geneName+','+seq.clinical+','
                                  +seq.seq3[-25:] + "<" + seq.wt + ">" + seq.seq5[:25]+','+"|".join(seq.varSeqList)+','+seq.chromo)
                rows.append([stringified])
            print("Creating CSV File")
            spamwriter.writerows(rows)
        print("File created successfully")
    except PermissionError as e:
        print("Writing CSV Failed, Please close file and run again")

def export(DB):
    """
    outputs a seqTypeDB to a CSV file Using CSV module
    :param DB: {[key - string]:(seqType,[CAS1,CAS2...CASn],...)}
    :return: None, OUTput is a CSV file in the project dir
    """
    print("exporting to file")
    try:
        with open('outputs.csv', 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            spamwriter.writerow([",".join(["srID","geneID","geneName","Clinical Significance","Condition",
                                           "WT sequance","Variants","chromosome","base replacements",
                                            "readingFrame",
                                           ",".join(casDB)]).replace('"','')]) #column titles
            rows = []
            print("Building rows")
            DBkeys = DB.keys()
            for key in DBkeys: #data
                seq = DB[key]
                try:
                    rslts = parseRslts(seq.rslts)
                except AttributeError as e:
                    print("Error:",seq)
                    rslts = seq
                # if seq.geneName == None:
                #     seq.geneName = ""
                # if (seq.residue != None and seq.aaPosition != None and seq.wt != None):
                #     base_replacement = str(seq.residue)+str(seq.aaPosition)+str(seq.wt)
                stringified = str(seq.snpID+','+seq.geneID+','+seq.geneName+','+seq.clinical+','+seq.condition+','
                                  +seq.seq3[-25:] + "<" + seq.wt + ">" + seq.seq5[:25]+','+"|".join(seq.varSeqList)+','+seq.chromo+','+seq.baseRepl+','
                                  +str(seq.readingFrame)+','+",".join(rslts))
                rows.append([stringified])
            print("Creating CSV File")
            spamwriter.writerows(rows)
        print("File created successfully")
    except PermissionError as e:
        print("Writing CSV Failed, Please close file and run again")

def matchOnDB(DB):
    """

    :param DB: gets seqType DB
    :return: rsltsDB - {key-String wtSequence:[seqType,[CAS1,CAS2...CASn],...]}
    """
    snps_with_pam=0
    count=0
    print("Matching and building results DB")
    global AVERAGERSLTS
    rsltDB = {}
    for seq in DB: #Go over all objects in originPyDB
        currentWT = seq.seq3[-25:]+seq.wt+seq.seq5[:25] #DB keys are strings, for good hash function in dictionary instead of seqtype hash function
        tempSeq = duplicate(seq)

        for variant in seq.varSeqList: #Go over all .patho
            for CAS in casDB: #Go over all PAMseqs
                try:
                   currentVar = seq.seq3+variant+seq.seq5
                   compwt = reverseSeq(currentWT)[::-1] #on complementary we work with the reversed antisense, so the PAM remains the same, but the window starting point changes to seq5 length
                   compVar = reverseSeq(currentVar)[::-1]
                   matchPositions = match(CAS, currentVar, len(seq.seq3), rtrnLoc=True)

                   if matchPositions != []: #if there are any locations that match on variant
                       count+=1
                       for matchIndex in matchPositions: #go over each of these locations
                           casname=getCas(CAS)
                           for name in casname:
                                valid = True
                                for cass in casDic[name]:
                                   if (match(cass,currentWT,matchIndex,specific=True)):
                                       valid=False
                                       break
                                   if (match(CAS, currentWT, matchIndex, specific=True)):
                                       valid = False
                                       break
                                   if name=="SpCas9" and (match("NAG",currentWT,matchIndex,specific=True)):
                                       valid=False
                                       break
                                   if valid==True:
                                      tempSeq.rslts.append((seq.seq3[-25:] + "<" + variant + ">" + seq.seq5[:25], CAS+" on "+str(matchIndex))) #rslts:[("VAR1 DNA",[CAS1,...,CASn]),...,("VARn DNA",[CAS1,...,CASn])]
                                      if rsltDB.get(currentWT) == None: #update DB - new entry
                                          rsltDB[currentWT] = tempSeq
                                      break

                   matchPositions = match(CAS, compVar, len(seq.seq5), rtrnLoc=True) #again for complementary
                   if matchPositions != []:
                       count+=1
                       for matchIndex in matchPositions:
                           casname = getCas(CAS)
                           for name in casname:
                               valid = True
                               for cass in casDic[name]:
                                   if (match(cass, compwt, matchIndex,specific=True)):
                                       valid=False
                                       break
                                   if (match(CAS, compwt, matchIndex, specific=True)):
                                       valid = False
                                       break
                                   if name=="SpCas9" and (match("NAG", compwt, matchIndex,specific=True)):
                                       valid = False
                                       break

                                   if valid==True:
                                        tempSeq.rslts.append((seq.seq3[-25:] + "<" + variant + ">" + seq.seq5[:25], "comp. " + CAS+" on "+str(matchIndex)))
                                        if rsltDB.get(currentWT) == None:
                                            rsltDB[currentWT] = tempSeq
                                        break
                except Exception as e:
                    print("failed to match due to an unexpected error:",e,"\n Statistics unreliable")
        if count>0:
            snps_with_pam+=1
    return rsltDB,snps_with_pam

def gRNA(seq,PAM,loc=True):
    """
    Gets a sequance, PAM and Variant location, returns calculated guideRNA
    :param loc: Variant location, if none passed, assume seq is seqType
    :param seq: seqType, or String and loc provided
    :param PAM: String
    :return:
    """
    if loc: #seqType
        loc = match(PAM, seq, len(seq.seq3), loc=True)
        return reverseSeq(seq[loc:loc + len(PAM)])
    loc = match(PAM,seq,loc,loc=True)
    return reverseSeq(seq[loc:loc+len(PAM)])

def getCasName(Cas):
    if CASnames[Cas] != None:
        return CASnames[Cas]
    else:
        return "Na"
def getCas(Cas):
    li=[]
    for cass in casDic:
        if Cas in casDic[cass]:
            li.append(cass)
    return li


def main():
    t0 = time.time() #for statistics
    rawDB = importSNPS("snps_with_clinvar2.csv") #*import whole XML DB and *parse to originPyDB <<<<>>>>maybe XML isnt persistant in the way it's written
    #DB = parse(rawDB) #pyDB's will contain seqType objects with .seqID .original-string .patho-list of strings .workingCas list of strings)
    #exportSeq(DB) #creates a file of analysed sequences (in addition to results file which contains only hits)
    rsltDB,snps_with_pam = matchOnDB(rawDB) #results DB is a dict, by definition duplicates are not possible\allowed in a set
    export(rsltDB)#export to CSV (for every ID: the seq, which CASes work, future imp: what is the ranking between them, future imp: what is the affect of the pathogen)
    t1 = time.time()

    print("---------statistics:---------")
    print ('Number of SNPS checked: ',  totalSequances)
    print ('Bad parse: ',badParse)
    print ('Number of SNPS with PAM match: ' ,snps_with_pam)
    print ('Number of results: ', len(rsltDB))
    #print("Debug, Hits found:",int(len(rsltDB))," ",int((len(rsltDB)/(totalSequances-badParse))*100),"%")
    # if len(rsltDB) == 0:
    #     print("Average Hits per sequence: 0")
    # else:
    #     print("Average Hits per sequence: ",AVERAGERSLTS/len(rsltDB))
    # print("Debug: failed RS's:",badParse," out of: ",totalSequances," ",int((badParse/totalSequances)*100),"%")
    # print("Debug, script ran in: ",t1-t0," seconds over ",totalSequances," sequences")

if __name__ == "__main__":
    main()