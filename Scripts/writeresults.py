def write_results(matches, cleanMatchdic, quietMatchdic,rsltsDic,filepath):
    with open(filepath, mode='w') as f:
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


            f.write("%s,%s,%s\n" %(key,rsltsDic[key], mlist))