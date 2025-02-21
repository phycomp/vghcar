'''整併期別 國民健康署癌症存活率AJCC整併期別:
1. 若病理分期PSD字根/字首為4或6者 (AJCC y prefix), 亦即有手術且手術前有接受放射治療/化學治療/標靶治療者 -> 選臨床期別(2nd)TCDB_CLG
2. 有手術且術前沒有接受放射治療/化學治療/標靶治療 -> 選病理TCDB_PSG
手術之術式PRIST為20-90 (肝癌包含術式13 攝護腺癌不包含術式21-23、25)
3. 其他特殊情境
  選病理 但病理不詳且pn in ('X' '0' '99') & opln_h in ('0' '1') 改選臨床TCDB_CLG
  選病理 但病理為888/BBB者改選臨床
  選臨床 但臨床為888/BBB者改選病理
  選臨床 但未手術、cM0且pM1者改選病理必須寫入 daily batch ETL for SMART
肝癌術式=['13']+[str(x) for x in range(20, 91)]
攝護腺癌術式=[str(x) for x in range(20, 91)]
for v in ['21', '22', '23', '25']: 攝護腺癌術式.remove(v)
#Alogorithm for HPA workStage group assigment
'''
from streamlit import session_state
from numpy import nan
from re import search

def rdcCLG(dSeries):
  valueCLG={"0":'1', "0  ":'1', "0A":'2', "0A ":'2', "0IS":'3', "1":'4', "1  ":'4', "1A":'5', "1A ":'6', "1A1":'6', "1A2":'7', "1B":'8', "1B ":'9', "1B1":'10', "1B2":'11', "1C":'12', "1C ":'13', "1S":'14', "2":'15', "2  ":'16', "2A":'17', "2A ":'18', "2A1":'19', "2A2":'20', "2B":'21', "2B ":'21', "2C":'22', "3":'23', "3  ":'23', "3A":'24', "3A ":'24', "3B":'25', "3B ":'25', "3C":'26', "3C ":'26', "3C1":'27', "3C2":'28', "4":'29', "4  ":'29', "4A":'30', "4A ":'30', "4B":'31', "4B ":'31', "4C":'32', "4C ":'32'}
  return dSeries.map(valueCLG)#.astype(str)

'''
PM(病理M分類) 與 CM(臨床M分類)的選擇:
1. CM=0, PM=BBB -> working M stage (STAGE_M)=CM=0, working STAGE
2. CM=1, PM=BBB -> working M stage (STAGE_M)=CM=1
3. CM=0, PM=1 ->-> working M stage (STAGE_M)=PM=1
4. CM=0, PM=anything else ->  working M stage (STAGE_M)=1
'''
def rtrvFullStage(row):
    tcdbPRIST, tcdbCLM, tcdbCLG, tcdbPAN, tcdbPAM, tcdbPSG, tcdbPSD, tcdbDOMSUPS, tcdbSPPSTF=row[113], row[128], row[129], row[133], row[134], row[135], row[136], row[147], row[149]
    if tcdbPSD in ["4", "6"]:
      if tcdbCLG in ['888','BBB']: stageG=tcdbPSG
      elif tcdbSPPSTF=='00' and tcdbCLM[0]=='0' and tcdbPAM[0]=='1' and tcdbPSG[0]=='4': stageG=tcdbPSG
      else: stageG=tcdbCLG
    elif (tcdbDOMSUPS not in ['','00000000','99999999']) and (tcdbSPPSTF>='20' and tcdbSPPSTF<='99' and  not (tcdbSPPSTF in ['21','22','23','25'] and tcdbPRIST=='C619') or (tcdbPRIST=='C220' and tcdbSPPSTF=='13')):
      if tcdbPAN in ['X','0','99'] and tcdbPSD=='999': stageG=tcdbCLG
      elif tcdbPSD in ['888','BBB']:stageG=tcdbCLG
      else: stageG=tcdbPSG
    else: stageG=tcdbCLG
    return stageG

def 整併期別(row, ndx癌症, ndx臨床, ndx病理, ndx字根, ndx原發部位, ndx手術類別, ndx手術日期, ndxCM, ndxPM):
  rowNdx, 癌症, 臨床期別, 病理期別, 字根=row.name, row[ndx癌症], row[ndx臨床], row[ndx病理], row[ndx字根]
  原發部位, 手術類別, 手術日期=row[ndx原發部位], row[ndx手術類別], row[ndx手術日期]
  CM, PM=row[ndxCM], row[ndxPM]
  if 字根 in ['4', '6']:
    if 臨床期別 in ['888','BBB']: HPA期別=病理期別
    elif 手術類別=='00' and CM[0]=='0' and PM[0]=='1' and 病理期別[0]=='4': HPA期別=病理期別
    else: HPA期別=臨床期別
  elif 手術日期 not in ['','00000000','99999999'] and 手術類別 in range(20, 91) and not (手術類別 in ['21','22','23','25'] and 原發部位=='C619') or (原發部位=='C220' and 手術類別=='13'):
      if 病理期別=='999' and 病理N in ['X','0','99']: HPA期別=臨床期別
      elif 病理期別 in ['888','BBB']: HPA期別=臨床期別
      else: HPA期別=病理期別
  else:
    HPA期別=臨床期別
  allCncr=session_state['AllCancer']
  allCncr['HPA期別'][rowNdx]=HPA期別

def 合併期別(row, ndx癌症, ndx臨床, ndx病理, ndx字根):#, ndxHPA):
  #去除期別=['888', '000', 'BBB'] #3_13 in 888, 000, BBB
  #43 PBC_PGROUP1, 129 TCDB_CLG, 135 TCDB_PSG, 136 TCDB_PSD
  去除期別=['888', '999', 'BBB']
  rowNdx, 癌症, 臨床期別, 病理期別, 病理字根=row.name, row[ndx癌症], row[ndx臨床], row[ndx病理], row[ndx字根]
  #rowNdx, 癌症, 臨床期別, 病理期別, 病理字根=row.name, row[43], row[129], row[135], row[136]
  if 癌症 is not nan: 癌症=癌症.strip()
  if 病理字根 in ['4', '6']:  #病理字根==>PSD 3_14=4 or 3_14=6
    HPA期別=臨床期別
  else:
    HPA期別=臨床期別 if 病理期別 in 去除期別 else 病理期別
    if 病理期別 in 去除期別: HPA期別=臨床期別
    else: HPA期別=病理期別
  df=session_state['AllCancer']
  df['HPA期別'][rowNdx]=HPA期別

def isValid(value):
  return False if value is nan else True
rmTAGs=['BBB','999','888','OC','-9', nan, '']

def AJCC期別(row, ajcc7, gndrAJCC, sexNdx, pristNdx, O3Mndx):
    rowNdx, tcdbGndr, tcdbPRIST, icdO3M=row.name, row[sexNdx], row[pristNdx], row[O3Mndx]    #TCDB_HISTGY
    #rowNdx, tcdbGndr, tcdbPRIST, icdO3M=row.name, row[104], row[113], row[115]    #TCDB_HISTGY
    allCncr=session_state['AllCancer']
    if not tcdbPRIST:# or not icdO3M:
      allCncr['ajccID'][rowNdx]=None
      allCncr['ajccTitle'][rowNdx]=None
      allCncr['ajccSchmID'][rowNdx]=None
      allCncr['ajccSchmTitle'][rowNdx]=None
    else:
      mainO3T, trailO3T=tcdbPRIST[:2], tcdbPRIST[3]
      icdO3T=f'{mainO3T}.{trailO3T}'
      fltrAJCC=ajcc7[ajcc7.O3T.apply(lambda vl:(not search('C7[012]', mainO3T) and vl =='LymphomaNotBrain') or vl!='LymphomaNotBrain')] #"C70", "C71", "C72"
      rsltAJCC=fltrAJCC.query(f'O3MBegin<="{icdO3M}" & O3MEnd>="{icdO3M}"')
      if len(rsltAJCC)==1:
          ajccID=rsltAJCC.ChapterId.values[0]
          ajccTitle=rsltAJCC.ChapterTitle.values[0]
          ajccSchmID=rsltAJCC.SchemaId.values[0]
          ajccSchmTitle=rsltAJCC.SchemaName.values[0]
          allCncr['ajccID'][rowNdx]=f"{ajccID}"#Series(data=f"{ajccID}", index=[rowNdx])
          allCncr['ajccTitle'][rowNdx]=f"{ajccTitle}"#Series(data=f"{ajccTitle}", index=[rowNdx])
          allCncr['ajccSchmID'][rowNdx]=f"{ajccSchmID}"#Series(data=f"{ajccSchmID}", index=[rowNdx])
          allCncr['ajccSchmTitle'][rowNdx]=f"{ajccSchmTitle}"#Series(data=f"{ajccSchmTitle}", index=[rowNdx])
          #print('rowNdx',allCncr['ajccTitle'][rowNdx])
          #print('rsltAJCC existed', rowNdx)
          #print('rsltAJCC existed', rowNdx, ajccID, ajccTitle, ajccSchmID, ajccSchmTitle)
          '''
          return ajccID, ajccTitle, ajccSchmID, ajccSchmTitle
          try:
              #allCncr.insert(lastClmn, 'ajccID', Series(ajccID, index=[rowNdx]))
              allCncr['ajccID']=Series(data=ajccID, index=[rowNdx])
          except:
              allCncr.insert(lastClmn, 'ajccID', '')
              allCncr.ajccID[0]=ajccID
          #print('AJCC', ajccID)
          #print('rsltAJCC', '|'.join(map(lambda x:str(len(x)), [ajccID, ajccTitle, ajccSchmID, ajccSchmTitle])))
          '''
      else:
          if tcdbGndr=='2':
              rsltGndrAJCC=gndrAJCC.query(f'ChapterId=="37" and O3T=="{icdO3T}" and O3MBegin<="{icdO3M}" and O3MEnd>="{icdO3M}"')
          else:
              rsltGndrAJCC=gndrAJCC.query(f'ChapterId=="28" and O3T=="{icdO3T}" and O3MBegin<="{icdO3M}" and O3MEnd>="{icdO3M}"')
              #rsltGndrAJCC=gndrAJCC.apply(lambda vl:True if vl.ChapterId=='28' and vl.O3T==icdO3M and vl.O3MBegin<=icdO3M<=vl.O3MEnd')
          #allCncr.ajccID=rsltGndrAJCC.ChapterId
          #allCncr.ajccTitle=rsltGndrAJCC.ChapterTitle
          #allCncr.ajccSchmID=rsltGndrAJCC.SchemaId
          #allCncr.ajccSchmTitle=rsltGndrAJCC.SchemaName
          if len(rsltGndrAJCC)==1:
              ajccID=rsltGndrAJCC.ChapterId.values[0]
              ajccTitle=rsltGndrAJCC.ChapterTitle.values[0]
              ajccSchmID=rsltGndrAJCC.SchemaId.values[0]
              ajccSchmTitle=rsltGndrAJCC.SchemaName.values[0]
              allCncr['ajccID'][rowNdx]=f"{ajccID}"   #Series(data=f"{ajccID}", index=[rowNdx])
              allCncr['ajccTitle'][rowNdx]=f"{ajccTitle}"
              allCncr['ajccSchmID'][rowNdx]=f"{ajccSchmID}"
              allCncr['ajccSchmTitle'][rowNdx]=f"{ajccSchmTitle}"
              print('rsltGndrAJCC existed', rowNdx, ajccID, ajccTitle, ajccSchmID, ajccSchmTitle)
              '''
              try:
                  #allCncr.insert(lastClmn, 'ajccID', Series(ajccID, index=[rowNdx]))
                  allCncr['ajccID']=Series(data=ajccID, index=[rowNdx])
              except:
                  #allCncr.insert(lastClmn, 'ajccID', '')
                  allCncr['ajccID']=Series(data=ajccID, index=[rowNdx])
                  allCncr.ajccID[0]=ajccID
              return ajccID, ajccTitle, ajccSchmID, ajccSchmTitle
              #print('gndrAJCC', ajccID)
              #print('genderAJCC', ajccID, ajccTitle, ajccSchmID, ajccSchmTitle)
              #print('genderAJCC', '|'.join(map(lambda x:str(len(x)), [ajccID, ajccTitle, ajccSchmID, ajccSchmTitle])))
              '''
          else:
              #allCncr.insert(0, 'ajccID', None)
              allCncr['ajccID'][rowNdx]=None  #Series(data=False, index=[rowNdx])
              allCncr['ajccTitle'][rowNdx]=None   #Series(data=False, index=[rowNdx])
              allCncr['ajccSchmID'][rowNdx]=None  #Series(data=False, index=[rowNdx])
              allCncr['ajccSchmTitle'][rowNdx]=None   #Series(data=False, index=[rowNdx])
              #print('NA existed')
              '''
              try:
                  #allCncr.insert(lastClmn, 'ajccID', Series('', index=[rowNdx]))
                  allCncr['ajccID']=Series(data=ajccID, index=[rowNdx])
              except: allCncr.insert(lastClmn, 'ajccID', '')
              return None, None, None, None
              print(len(allCncr.ajccID))
    #allCncr=concat(allCncr, df)
#del allCncr.ajccID
#df=DataFrame()
anthrCncr.apply(AJCC期別, axis=1)#, args=(df, allCncr))
allCncr
            '''

def rtrvFullStage(row, ndxCLG, ndxPSG):
    rmvTags=['8','B','9','N','-']
    rowNdx, TCDB_CLT, TCDB_CLN, TCDB_CLM, TCDB_CLG, TCDB_PSG, TCDB_PSD, TCDB_COTHSTG, ajccID, CLG_VALUE, PSG_VALUE=row.name, row[126], row[127], row[128], row[129], row[135], row[136], row[141], row[330], row[ndxCLG], row[ndxPSG]
    #print('CLG_VALUE, PSG_VALUE', CLG_VALUE, PSG_VALUE)
    if TCDB_PSD in ['4', '6']:
      if CLG_VALUE is nan: fullStage=TCDB_PSG
      elif PSG_VALUE is nan: fullStage=TCDB_CLG
      else: fullStage=TCDB_CLG if CLG_VALUE > PSG_VALUE else TCDB_PSG
      #else: fullStage=TCDB_CLG if isValid(CLG_VALUE) and isValid(PSG_VALUE) and CLG_VALUE > PSG_VALUE else TCDB_PSG
      #fullStage=TCDB_CLG if isValid(CLG_VALUE) and isValid(PSG_VALUE) and CLG_VALUE > PSG_VALUE else TCDB_PSG
    else:
      fullStage=TCDB_PSG if TCDB_CLG in rmTAGs else TCDB_CLG
    fullStage='' if fullStage in rmTAGs else fullStage
    fullStage=TCDB_COTHSTG if '57' == ajccID and not fullStage else fullStage
    fullStage='' if fullStage in ['0000', '8888', '9999'] else fullStage
    allCncr=session_state['AllCancer']
    allCncr['fullStage'][rowNdx]=fullStage
'''
    tStage=TCDB_CLT[0] if TCDB_CLT is not nan else TCDB_CLT
    nStage=TCDB_CLN[0] if TCDB_CLN is not nan else TCDB_CLN
    mStage=TCDB_CLM[0] if TCDB_CLM is not nan else TCDB_CLM
    #print('fullStage', fullStage)
    if fullStage and fullStage is not nan:
        firstStage=fullStage[0]
        simpleStage='' if firstStage in rmvTags else firstStage
        stage4C='C' if firstStage in ["4C ", "4C"] else firstStage
        stage4C='' if stage4C in rmvTags else stage4C
        allCncr['fullStage'][rowNdx]=fullStage
        allCncr['simpleStage'][rowNdx]=simpleStage
        allCncr['stage4C'][rowNdx]=stage4C
    if tStage and tStage is not nan:
        tStage='' if tStage in rmvTags else tStage
        allCncr['tStage'][rowNdx]=tStage
    if nStage and nStage is not nan:
        nStage='' if nStage in rmvTags else nStage
        allCncr['nStage'][rowNdx]=nStage
    if mStage and mStage is not nan:
        mStage='' if mStage in rmvTags else mStage
        allCncr['mStage'][rowNdx]=mStage
allCncr['fullStage']=None
allCncr['simpleStage']=None
allCncr['stage4C']=None
allCncr['tStage']=None
allCncr['nStage']=None
allCncr['mStage']=None
allCncr.apply(rtrvStage, axis=1, args=(allCncr,))

'''

def nvoStage(row):
  去除期別=['888', '000', 'BBB'] #3_13 in 888, 000, BBB
  癌症, 原發部位, 臨床期別, 病理N, 病理期別, 病理字根, 淋巴結範圍=row[43].strip(), row[113], row[129], row[133], row[135], row[136], row[151]
  if 病理字根 in ['4', '6']:    #病理字根==>PSD 3_14=4 or 3_14=6
    workStage=臨床期別   #3_7 CLG
    if workStage in 去除期別: workStage=病理期別 # else 臨床期別 ['888', '000', 'BBB']
  else:
    workStage=病理期別  #3_7 CLG 臨床期別 3_13 PSG 病理期別
    if workStage in 去除期別: workStage=臨床期別 #else 病理期別
  if 臨床期別 in ['000', '888'] and 病理期別 in ['BBB', '888']: workStage = 'X'
  if 臨床期別 in 去除期別 and 病理期別 in 去除期別: workStage = 'X'

def mergeStage(row):
  癌症, 原發部位, 臨床期別, 病理N, 病理期別, 病理字根, 淋巴結範圍=row[43].strip(), row[113], row[129], row[133], row[135], row[136], row[151]
  術式=原發部位[1:3]
  if 癌症=='HCC_18' and 術式 in 肝癌術式: 癌症期別=病理期別    #肝癌
  elif 癌症=='PROSTATE_41' and 術式 in 攝護腺癌術式: #攝護腺癌
    癌症期別=病理期別
  #PBC_PGROUP1
  if 病理字根 in ['4', '6']:
    if 有手術:
      fullStage=臨床期別 if 治療 else 病理期別  #沒有治療
      #病理N in ['X' '0' '99'] and 淋巴結範圍 in ['0' '1']:#淋巴結手術範圍 TCDB_SRLNSOF
    elif 有手術: fullStage=病理期別
    fullStage=臨床期別 if CLG_VALUE > PSG_VALUE else 病理期別
  rmvTags=['8','B','9','N','-']
  rowNdx, TCDB_CLT, TCDB_CLN, TCDB_CLM, TCDB_CLG, TCDB_PSG, TCDB_PSD, TCDB_COTHSTG, ajccID, CLG_VALUE, PSG_VALUE=row.name, row[126], row[127], row[128], row[129], row[135], row[136], row[141], row[330], row[331], row[332]
