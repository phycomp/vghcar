#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) stop("At least one argument must be supplied (input file).n", call.=FALSE)
fullpath=args[1]
fname=basename(fullpath)
abspath=dirname(fullpath)
library(stringr, grep)	
DATA_DATE<-str_match(fname, 'AllCancer20(\\d+-\\d+-\\d+)')[2]	#AllCancer2018-12-01.csv
DATA_DATE<-gsub('-', '', DATA_DATE)
print(DATA_DATE)
#DATA_DATE <- match	#"180820" # 資料日期
HPA_DATE <- "2019-12-31" # 國健署回壓日期
##### Path Setting Area #######################################################
if ("Windows" == Sys.info()['sysname']) {
    DATASET_DIR <- "D:/SMART/VGHCAR_DATASET/"
    CODE_DIR <- "D:/SMART/VGHCAR/"
} else if ("Linux" == Sys.info()['sysname']) {
    ## AllCancer.csv
    DATASET_DIR <- "UUID"	#../VGHCAR_DATASET/"
    CODE_DIR <- ""	#../shiny-server/"
} else {
    print("Please check your OS environment!!!")
    print("Please check your OS environment!!!")
    print("Please check your OS environment!!!")
    return;
}
###############################################################################


library(dplyr)

# 外部函式
source(paste0(CODE_DIR, "libs/time/dateConversion.r"), encoding = "UTF-8")

# 因為執行順序的關係，把主程式包成 main
..main <- function() {

    # 測試使用 TRUE
    passExtendAJCC7 <- FALSE
    # passExtendAJCC7 <- TRUE

    if (passExtendAJCC7) {
        print(paste("Loading extend AJCC7 for testing.", Sys.time()))
        rawData <- readRDS(file = paste0(DATASET_DIR, "extendAJCC7", DATA_DATE, ".rds"))
    } else {
        print(paste("Read all data.", Sys.time()))
        rawData <- .readRawData()
        uuidList <- .readUuidList()

        print(paste("Join UUID List", Sys.time()))
        rawData <- .joinUuidList(rawData, uuidList)

        print(paste("Make and save new UUID List", Sys.time()))
        uuidList <- .makeUuidList(rawData)
        .saveUuidList(uuidList)

        print("Start extend raw data.")
        print(paste("Start extend raw data.", Sys.time()))
        print(paste(".extendAJCC7", Sys.time()))
        rawData <- .extendAJCC7(rawData)
        saveRDS(rawData, file = paste0(DATASET_DIR, "extendAJCC7", DATA_DATE, ".rds"))
    }

    print(paste(".extendStage", Sys.time()))
    rawData <- .extendStage(rawData)
    print(paste(".extendDates", Sys.time()))
    rawData <- .extendDates(rawData)
    print(paste(".extendStatus", Sys.time()))
    rawData <- .extendStatus(rawData)
    print(paste(".extendClassify", Sys.time()))
    rawData <- .extendClassify(rawData)
    print(paste(".extendRename", Sys.time()))
    rawData <- .extendRename(rawData)

    # print(paste("", Sys.time()))
    # rawData <- (rawData)

    print('Wait for saving data!')
    return(rawData)
}

.readRawData <- function() {
    # baseTablePath <- paste0(DATASET_DIR, "AllCancer.csv")
    baseTablePath <- file.path(fullpath)#abspath, "AllCancer20", substr(DATA_DATE, 1, 2), "-", substr(DATA_DATE, 3, 4), "-", substr(DATA_DATE, 5, 6), ".csv")
    rawData <- read.csv(baseTablePath, fileEncoding = "big5", encoding="utf8", stringsAsFactors=FALSE)

    return(rawData)
}

.readUuidList <- function() {
    baseTablePath <- paste0(CODE_DIR, "UUID/uuidList.csv")
    uuidList <- read.csv(baseTablePath, fileEncoding = "big5", encoding="utf8", stringsAsFactors=FALSE)

    return(uuidList)
}

.firstMakeUuidList <- function(rawData) {
    chooseColumn <- c("LV_UUID", "TCDB_CISTCSER")
    uuidList <- rawData[,chooseColumn]
    colnames(uuidList) <- c("UNIQ_UUID", "TCDB_CISTCSER")

    return(uuidList)
}

.makeUuidList <- function(rawData) {
    chooseColumn <- c("UNIQ_UUID", "TCDB_CISTCSER")
    uuidList <- rawData[,chooseColumn]

    return(uuidList)
}

.joinUuidList <- function(rawData, uuidList) {
    extendedData <- rawData %>%
                    left_join(uuidList, by = c("TCDB_CISTCSER" = "TCDB_CISTCSER")) %>%
                    mutate(UNIQ_UUID = if_else(is.na(UNIQ_UUID), LV_UUID, UNIQ_UUID))
    return(extendedData)
}

.saveUuidList <- function(uuidList) {
    baseTablePath <- paste0(DATASET_DIR, "uuidList.csv")
    write.csv(uuidList, file = baseTablePath, row.names = FALSE)
}

.extendAJCC7 <- function(rawData) {
    # O3T 的大項目
    rawData <- rawData %>% mutate(MainO3T = substr(TCDB_PRIST, 1, 3))
    # O3T
    rawData <- rawData %>% mutate(ICD_O3T = paste0(substr(TCDB_PRIST, 1, 3), ".", substr(TCDB_PRIST, 4, 4)))
    # O3M
    rawData <- rawData %>% mutate(ICD_O3M = TCDB_HISTGY)

    # 增加 AJCC7 的編碼
    AJCC7 <- read.csv(paste0(CODE_DIR, "AJCC7/AJCC7.csv"), encoding = "utf8", strip.white = TRUE, stringsAsFactors = FALSE)
    AJCC7_Gender <- read.csv(paste0(CODE_DIR, "AJCC7/AJCC7_Gender.csv"), encoding = "utf8", strip.white = TRUE, stringsAsFactors = FALSE)
print("AJCC7 for loop")
    totalCount = dim(rawData)[1]
    for (row in 1:totalCount) {
        if (0 == row %% 100) {
            print(paste("In loop:", row, row / 100, "/", totalCount, ", now time:", Sys.time()))
        }
        specificAJCC7 <- AJCC7 %>%
            dplyr::filter("LymphomaNotBrain" == O3T | "Lymphoma" == O3T | O3T == rawData$ICD_O3T[row]) %>%
            dplyr::filter("LymphomaNotBrain" != O3T | ("LymphomaNotBrain" == O3T & !(substr(rawData$ICD_O3T[row], 1, 3) %in% c("C70", "C71", "C72")))) %>%
            dplyr::filter(O3MBegin <= rawData$ICD_O3M[row]) %>%
            dplyr::filter(rawData$ICD_O3M[row] <= O3MEnd)

        if (1 == dim(specificAJCC7)[1]) {
            rawData[row, "AJCC7_ID"] <- specificAJCC7["ChapterId"]
            rawData[row, "AJCC7_NAME"] <- specificAJCC7["ChapterTitle"]
            rawData[row, "AJCC7_SUBID"] <- specificAJCC7["SchemaId"]
            rawData[row, "AJCC7_SUBNAME"] <- specificAJCC7["SchemaName"]
        } else {
            if (is.na(rawData[row, "TCDB_SEX"])) {
                next
            }
            # Female
            if ("2" == rawData[row, "TCDB_SEX"]) {
                genderSpecificAJCC7 <- AJCC7_Gender %>%
                    dplyr::filter("37" == ChapterId) %>%
                    dplyr::filter(O3T == rawData$ICD_O3T[row]) %>%
                    dplyr::filter(O3MBegin <= rawData$ICD_O3M[row]) %>%
                    dplyr::filter(rawData$ICD_O3M[row] <= O3MEnd)
            } else {
                genderSpecificAJCC7 <- AJCC7_Gender %>%
                    dplyr::filter("28" == ChapterId) %>%
                    dplyr::filter(O3T == rawData$ICD_O3T[row]) %>%
                    dplyr::filter(O3MBegin <= rawData$ICD_O3M[row]) %>%
                    dplyr::filter(rawData$ICD_O3M[row] <= O3MEnd)
            }
            if (1 == dim(genderSpecificAJCC7)[1]) {
                rawData[row, "AJCC7_ID"] <- genderSpecificAJCC7["ChapterId"]
                rawData[row, "AJCC7_NAME"] <- genderSpecificAJCC7["ChapterTitle"]
                rawData[row, "AJCC7_SUBID"] <- genderSpecificAJCC7["SchemaId"]
                rawData[row, "AJCC7_SUBNAME"] <- genderSpecificAJCC7["SchemaName"]
            }
        }
    }

    return(rawData)
}

.extendStage <- function(rawData) {
    # 定義癌症期別
    rawData <- rawData %>%
        # 病理判斷 (1st) 臨床期別組合 PSG Clinical Stage Group
        mutate(Severity = substr(TCDB_PSG, 1, 1)) %>%
        # 醫師判斷 (2nd) 臨床期別組合 CLG Clinical Stage Group
        mutate(Severity = if_else(Severity %in% c('8','B','0','9','','N','-'), substr(TCDB_CLG, 1, 1), Severity)) %>%
        # 不屬於期別的，清空
        mutate(Severity = if_else(Severity %in% c('8','B','0','9','','N','-'), '', Severity))


    # 為了判斷期別所增加的數字(病理期別)
    rawData <- rawData %>%
        mutate(CLG_VALUE = 0) %>%
        mutate(CLG_VALUE = if_else("0"    == TCDB_CLG, 1, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("0  "  == TCDB_CLG, 1, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("0A"   == TCDB_CLG, 2, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("0A "  == TCDB_CLG, 2, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("0IS"  == TCDB_CLG, 3, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1"    == TCDB_CLG, 4, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1  "  == TCDB_CLG, 4, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1A"   == TCDB_CLG, 5, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1A "  == TCDB_CLG, 6, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1A1"  == TCDB_CLG, 6, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1A2"  == TCDB_CLG, 7, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1B"   == TCDB_CLG, 8, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1B "  == TCDB_CLG, 9, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1B1"  == TCDB_CLG, 10, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1B2"  == TCDB_CLG, 11, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1C"   == TCDB_CLG, 12, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1C "  == TCDB_CLG, 13, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("1S"   == TCDB_CLG, 14, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2"    == TCDB_CLG, 15, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2  "  == TCDB_CLG, 16, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2A"   == TCDB_CLG, 17, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2A "  == TCDB_CLG, 18, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2A1"  == TCDB_CLG, 19, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2A2"  == TCDB_CLG, 20, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2B"   == TCDB_CLG, 21, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2B "  == TCDB_CLG, 21, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("2C"   == TCDB_CLG, 22, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3"    == TCDB_CLG, 23, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3  "  == TCDB_CLG, 23, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3A"   == TCDB_CLG, 24, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3A "  == TCDB_CLG, 24, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3B"   == TCDB_CLG, 25, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3B "  == TCDB_CLG, 25, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3C"   == TCDB_CLG, 26, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3C "  == TCDB_CLG, 26, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3C1"  == TCDB_CLG, 27, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("3C2"  == TCDB_CLG, 28, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4"    == TCDB_CLG, 29, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4  "  == TCDB_CLG, 29, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4A"   == TCDB_CLG, 30, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4A "  == TCDB_CLG, 30, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4B"   == TCDB_CLG, 31, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4B "  == TCDB_CLG, 31, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4C"   == TCDB_CLG, 32, CLG_VALUE)) %>%
        mutate(CLG_VALUE = if_else("4C "  == TCDB_CLG, 32, CLG_VALUE))


    # 為了判斷期別所增加的數字(臨床期別)
    rawData <- rawData %>%
        mutate(PSG_VALUE = 0) %>%
        mutate(PSG_VALUE = if_else("0"    == TCDB_PSG, 1, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("0  "  == TCDB_PSG, 1, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("0A"   == TCDB_PSG, 2, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("0A "  == TCDB_PSG, 2, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("0IS"  == TCDB_PSG, 3, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1"    == TCDB_PSG, 4, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1  "  == TCDB_PSG, 4, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1A"   == TCDB_PSG, 5, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1A "  == TCDB_PSG, 6, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1A1"  == TCDB_PSG, 6, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1A2"  == TCDB_PSG, 7, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1B"   == TCDB_PSG, 8, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1B "  == TCDB_PSG, 9, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1B1"  == TCDB_PSG, 10, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1B2"  == TCDB_PSG, 11, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1C"   == TCDB_PSG, 12, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1C "  == TCDB_PSG, 13, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("1S"   == TCDB_PSG, 14, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2"    == TCDB_PSG, 15, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2  "  == TCDB_PSG, 16, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2A"   == TCDB_PSG, 17, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2A "  == TCDB_PSG, 18, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2A1"  == TCDB_PSG, 19, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2A2"  == TCDB_PSG, 20, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2B"   == TCDB_PSG, 21, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2B "  == TCDB_PSG, 21, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("2C"   == TCDB_PSG, 22, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3"    == TCDB_PSG, 23, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3  "  == TCDB_PSG, 23, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3A"   == TCDB_PSG, 24, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3A "  == TCDB_PSG, 24, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3B"   == TCDB_PSG, 25, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3B "  == TCDB_PSG, 25, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3C"   == TCDB_PSG, 26, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3C "  == TCDB_PSG, 26, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3C1"  == TCDB_PSG, 27, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("3C2"  == TCDB_PSG, 28, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4"    == TCDB_PSG, 29, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4  "  == TCDB_PSG, 29, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4A"   == TCDB_PSG, 30, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4A "  == TCDB_PSG, 30, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4B"   == TCDB_PSG, 31, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4B "  == TCDB_PSG, 31, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4C"   == TCDB_PSG, 32, PSG_VALUE)) %>%
        mutate(PSG_VALUE = if_else("4C "  == TCDB_PSG, 32, PSG_VALUE))

    rawData$TCDB_VOTHSTG <- as.character(rawData$TCDB_VOTHSTG)

    # 定義完整癌症期別
    rawData <- rawData %>%
        # 3-14 病理分期字根/字首* (TCDB_PSD)
        # 4_Y-首次治療期間或治療後進行的病理分期
        # 6_M&Y-多原發腫瘤併首次治療期間進行之病理分期
        mutate(FullStage = if_else((TCDB_PSD %in% c('4', '6')), if_else((CLG_VALUE > PSG_VALUE), TCDB_CLG, TCDB_PSG), if_else((TCDB_CLG %in% c('BBB','999','888','OC','-9','NA','')), TCDB_PSG, TCDB_CLG))) %>%
        # 不屬於期別的，清空
        mutate(FullStage = if_else(FullStage %in% c('BBB','999','888','OC','-9','NA',''), '', FullStage)) %>%
        # 3-19 其他分期系統期別(臨床分期)* (Lymphoma 才要)
        mutate(FullStage = if_else(('57' == AJCC7_ID & '' == FullStage), TCDB_COTHSTG, FullStage)) %>%
        # 不屬於期別的，清空
        mutate(FullStage = if_else(FullStage %in% c('0000', '8888', '9999'), '', FullStage))

    # 定義癌症期別 Stage I, II, III, IV
    rawData <- rawData %>%
        # 只取第一碼
        mutate(SimpleStage = substr(FullStage, 1, 1)) %>%
        # 不屬於期別的，清空
        mutate(SimpleStage = if_else(SimpleStage %in% c('8','B','9','','N','-'), '', SimpleStage))

    # 定義癌症期別 Stage I, II, III, IV/IVA/IVB, IVC
    rawData <- rawData %>%
        # 4C 獨立，其餘只取第一碼
        mutate(Stage4C = if_else(("4C" == FullStage | "4C " == FullStage), 'C', substr(FullStage,1,1))) %>%
        # 不屬於期別的，清空
        mutate(Stage4C = if_else(Stage4C %in% c('8','B','9','','N','-'), '', Stage4C))

    rawData <- rawData %>%
        # 只取第一碼
        mutate(tStage = substr(TCDB_CLT, 1, 1)) %>%
        # 不屬於期別的，清空
        mutate(tStage = if_else(tStage %in% c('8','B','9','','N','-'), '', tStage))

    rawData <- rawData %>%
        # 只取第一碼
        mutate(nStage = substr(TCDB_CLN, 1, 1)) %>%
        # 不屬於期別的，清空
        mutate(nStage = if_else(nStage %in% c('8','B','9','','N','-'), '', nStage))

    rawData <- rawData %>%
        # 只取第一碼
        mutate(mStage = substr(TCDB_CLM, 1, 1)) %>%
        # 不屬於期別的，清空
        mutate(mStage = if_else(mStage %in% c('8','B','9','','N','-'), '', mStage))


    return(rawData)
}



.extendDates <- function(rawData) {
    # 診斷日期採用 2-5 最初診斷日
    rawData <- rawData %>% mutate(StartDate = dateConversion(TCDB_DOID))

    # 2018.7 修改 startDate 抓法，先抓最初病理診斷日期 pbcancer.PB3MDATE → 最初臨床診斷日期 pbcancer.PB3TDATE → 最初診斷日期 CISTCRM.DOID → 首次就診日期 CISTCRM.DOFC
    #rawData <- rawData %>%
    #mutate(startDate_1 = if_else(!is.na(TCDB_DOID) & dateConversion(TCDB_DOID)>dateConversion(TCDB_DOFC),dateConversion(TCDB_DOID),dateConversion(TCDB_DOFC))) %>%
    #mutate(startDate_2 = if_else(!is.na(PBC_PB3TDATE) & dateConversion(PBC_PB3TDATE)>dateConversion(TCDB_DOID),dateConversion(PBC_PB3TDATE),dateConversion(TCDB_DOID))) %>%
    #mutate(startDate = if_else(!is.na(dateConversion(PBC_PB3MDATE)) & dateConversion(PBC_PB3MDATE)>dateConversion(PBC_PB3TDATE),dateConversion(PBC_PB3MDATE),dateConversion(PBC_PB3TDATE)))


    # 因為撈資料的問題，院外死亡日期改在這邊做 replace
    rawData <- rawData %>% mutate(LV_OUTDIE = gsub("院外死亡日期[：:]", "", PBAB_PDESC))

    # 增加結束時間  
#message('LV_OUTDIE',rawData$LV_OUTDIE)
    rawData <- rawData %>%
        # 最後聯絡或死亡日期 (1st) DATATYPE IN ('1','4') AND VITSS IN ('0','2') --5-4 死亡狀態 0死亡 2死亡但日期不詳
        mutate(EndDate = if_else(TCDB_DATATYPE %in% c('1', '4') & TCDB_VITSS %in% c('0','2'), dateConversion(TCDB_DLCOD), NULL)) %>%
        # 醫師註記院外死亡   (2nd)
        mutate(EndDate = if_else(is.na(EndDate), dateConversion(LV_OUTDIE), EndDate)) %>%
        # 院內註記死亡日期   (3rd) (PMEDSTAT = '3')
        mutate(EndDate = if_else(is.na(EndDate) & 3 == PBAS_PMEDSTAT, dateConversion(PBAS_PLSTVDT), EndDate)) %>%
        # 如果沒有值，放癌登的最後日期，當做 Lost Follow-Up
        mutate(EndDate = if_else(is.na(EndDate), dateConversion(TCDB_DLCOD), EndDate)) %>%
        # 如果還是沒有值，就去看看他有沒有掛其他科
        mutate(EndDate = if_else(is.na(EndDate) & 3 != PBAS_PMEDSTAT, dateConversion(PBAS_PLSTVDT), EndDate))

    # 結束時間 in VGH
    rawData <- rawData %>%
        # 院內註記日期
        mutate(EndDateInVgh = if_else(dateConversion(TCDB_DLCOD) > dateConversion(PBAS_PLSTVDT), dateConversion(TCDB_DLCOD), dateConversion(PBAS_PLSTVDT)))

    # 結束時間 in VGH then 國建署
    rawData <- rawData %>%
        # 院內註記日期
        mutate(EndDateInVghThenOutdie = EndDateInVgh) %>%
        # 醫師註記院外死亡
        mutate(EndDateInVghThenOutdie = if_else('' != LV_OUTDIE & EndDateInVghThenOutdie < dateConversion(LV_OUTDIE), dateConversion(LV_OUTDIE), EndDateInVghThenOutdie))

        # # 先.....不要節外生枝好了
        # mutate(EndDateInVghThenOutdie = if_else('' != LV_OUTDIE & EndDateInVghThenOutdie < dateConversion(LV_OUTDIE), dateConversion(LV_OUTDIE), EndDateInVghThenOutdie)) %>%
        # # 如果沒有院外死亡日期，然後現有的存活日期比國健署回壓日期還要早，就要改成國健署的日期
        # mutate(EndDateInVghThenOutdie = if_else('' == LV_OUTDIE & EndDateInVghThenOutdie < dateConversion(HPA_DATE), dateConversion(HPA_DATE), EndDateInVghThenOutdie))

    # 結束時間 in VGH then 國建署 then VGH
    rawData <- rawData %>%
        # 結束時間 in VGH then 國建署
        mutate(EndDateInVghThenOutdieThenVgh = EndDateInVghThenOutdie) %>%
        # 院內註記日期
        mutate(EndDateInVghThenOutdieThenVgh = if_else(EndDateInVghThenOutdieThenVgh < dateConversion(PBAS_PLSTVDT), dateConversion(PBAS_PLSTVDT), EndDateInVghThenOutdieThenVgh)) %>%
        mutate(EndDateInVghThenOutdieThenVgh = if_else(is.na(EndDateInVghThenOutdieThenVgh), dateConversion(TCDB_DLCOD), EndDateInVghThenOutdieThenVgh))

    # 以開始時間排序，判斷是否為此病人的第一顆癌症
        rawData <- rawData %>% arrange(StartDate)
        rawData$IsFirstCancer <- !duplicated(rawData$TCDB_PHISTNUM)

    # 計算 EndDate
    # 診斷日期採用 2-5 最初診斷日
    rawData <- rawData %>% mutate(startDate = dateConversion(TCDB_DOID))

    # 2018.7 修改 startDate 抓法，先抓最初病理診斷日期 pbcancer.PB3MDATE → 最初臨床診斷日期 pbcancer.PB3TDATE → 最初診斷日期 CISTCRM.DOID → 首次就診日期 CISTCRM.DOFC
    #rawData <- rawData %>%
    #mutate(startDate_1 = if_else(!is.na(TCDB_DOID) & dateConversion(TCDB_DOID)>dateConversion(TCDB_DOFC),dateConversion(TCDB_DOID),dateConversion(TCDB_DOFC))) %>%
    #mutate(startDate_2 = if_else(!is.na(PBC_PB3TDATE) & dateConversion(PBC_PB3TDATE)>dateConversion(TCDB_DOID),dateConversion(PBC_PB3TDATE),dateConversion(TCDB_DOID))) %>%
    #mutate(startDate = if_else(!is.na(dateConversion(PBC_PB3MDATE)) & dateConversion(PBC_PB3MDATE)>dateConversion(PBC_PB3TDATE),dateConversion(PBC_PB3MDATE),dateConversion(PBC_PB3TDATE)))



    # 復發日期
    rawData <- rawData %>% mutate(recurrentDate  = dateConversion(PBC_DFSDATE))
    rawData <- rawData %>% mutate(tRecurrentDate = dateConversion(PBC_TDATE))
    rawData <- rawData %>% mutate(nRecurrentDate = dateConversion(PBC_NDATE))
    rawData <- rawData %>% mutate(mRecurrentDate = dateConversion(PBC_MDATE))


    # 最大可能存活日期，拿來算 Lost Follow-Up Date
    rawData <- rawData %>%
        # 國健署回報日期 2016-12-31 -> 改成設定檔
        mutate(maxPossibleLiveDate = dateConversion(HPA_DATE)) %>%
        # 如果癌登的最後日期比較大，則放癌登的最後日期
        mutate(maxPossibleLiveDate = if_else(!is.na(TCDB_DLCOD) & dateConversion(TCDB_DLCOD) > maxPossibleLiveDate, dateConversion(TCDB_DLCOD), maxPossibleLiveDate)) %>%
        # 如果最後就診日比較大，則放最後就診日
        mutate(maxPossibleLiveDate = if_else(!is.na(PBAS_PLSTVDT) & dateConversion(PBAS_PLSTVDT) > maxPossibleLiveDate, dateConversion(PBAS_PLSTVDT), maxPossibleLiveDate))
        # 如果門診記錄的日期比較大，則放門診記錄日期 >>>還沒撈<<<
        # mutate(maxPossibleLiveDate = if_else(!is.na(DTA_DTPDATE) & dateConversion(DTA_DTPDATE) > maxPossibleLiveDate, dateConversion(DTA_DTPDATE), maxPossibleLiveDate))

    # 計算死亡日期
    rawData <- rawData %>%
        mutate(deathDate = EndDateInVghThenOutdie)

    # OLD 計算死亡日期
    # rawData <- rawData %>%
    #     # 最後聯絡或死亡日期 (1st) DATATYPE IN ('1','4') AND VITSS IN ('0','2') --5-4 死亡狀態 0死亡 2死亡但日期不詳
    #     mutate(deathDate = if_else(TCDB_DATATYPE %in% c('1', '4') & TCDB_VITSS %in% c('0','2'), dateConversion(TCDB_DLCOD), NULL)) %>%
    #     # 醫師註記院外死亡   (2nd)
    #     mutate(deathDate = if_else(is.na(deathDate), dateConversion(LV_OUTDIE), deathDate)) %>%
    #     # 院內註記死亡日期   (3rd) (PMEDSTAT = '3')
    #     mutate(deathDate = if_else(is.na(deathDate) & 3 == PBAS_PMEDSTAT, dateConversion(PBAS_PLSTVDT), deathDate)) %>%
    #     # 如果仍然沒有值，則放最大可能存活日期當做 Lost Follow-Up
    #     mutate(deathDate = if_else(is.na(deathDate), maxPossibleLiveDate, deathDate))


    # 計算個別的 duration
    rawData <- rawData %>%
        mutate(recurrentDuration = if_else(!is.na(recurrentDate), recurrentDate - startDate, NA_real_)) %>% 
        mutate(recurrentDuration = as.numeric(recurrentDuration)) %>%
        mutate(recurrentDuration = if_else(recurrentDuration < 0, 0, recurrentDuration))
    rawData <- rawData %>%
        mutate(tDuration = if_else(!is.na(tRecurrentDate), tRecurrentDate - startDate, NA_real_)) %>% 
        mutate(tDuration = as.numeric(tDuration)) %>%
        mutate(tDuration = if_else(tDuration < 0, 0, tDuration))
    rawData <- rawData %>%
        mutate(nDuration = if_else(!is.na(nRecurrentDate), nRecurrentDate - startDate, NA_real_)) %>% 
        mutate(nDuration = as.numeric(nDuration)) %>%
        mutate(nDuration = if_else(nDuration < 0, 0, nDuration))
    rawData <- rawData %>%
        mutate(lrControlDuration = if_else(tDuration < nDuration ,tDuration ,nDuration)) %>% 
        mutate(lrControlDuration = as.numeric(lrControlDuration)) %>%
        mutate(lrControlDuration = if_else(lrControlDuration < 0, 0, lrControlDuration))
    rawData <- rawData %>%
        mutate(mDuration = if_else(!is.na(mRecurrentDate), mRecurrentDate - startDate, NA_real_)) %>% 
        mutate(mDuration = as.numeric(mDuration)) %>%
        mutate(mDuration = if_else(mDuration < 0, 0, mDuration))
    rawData <- rawData %>% 
        mutate(deathDuration = deathDate - startDate) %>% 
        mutate(deathDuration = as.numeric(deathDuration)) %>%
        mutate(deathDuration = if_else(deathDuration < 0, 0, deathDuration))

    # 真正的 Duration 計算
    # Local Control rate
    rawData <- rawData %>%
        mutate(LocalControlRateDuration = if_else(!is.na(tDuration), tDuration, deathDuration))
    
    # Locoregional Control rate
    rawData <- rawData %>%
        mutate(LocoregionalControlRateDuration = if_else(!is.na(lrControlDuration), lrControlDuration, deathDuration)) %>%
        mutate(LocoregionalControlRateDuration = if_else((LocoregionalControlRateDuration > LocalControlRateDuration), LocalControlRateDuration, LocoregionalControlRateDuration)) %>%
        mutate(LocoregionalControlRateDuration = if_else(is.na(LocoregionalControlRateDuration), deathDuration, LocoregionalControlRateDuration))

    # DistantMetasis Control rate
    rawData <- rawData %>%
        mutate(DistantMetastaseControlRateDuration = if_else(!is.na(mDuration), mDuration, deathDuration))

    # nodal Control Rate
    rawData <- rawData %>%
        mutate(NodalControlRateDuration = if_else(!is.na(nDuration), nDuration, deathDuration))


    # Overall survival
    rawData <- rawData %>% 
        mutate(OverallSurvivalDuration = deathDuration)
    # rawData <- rawData %>% 
    #     mutate(OverallSurvivalDuration = deathDuration) %>% 
    #     mutate(OverallSurvivalDuration = as.numeric(OverallSurvivalDuration)) %>%
    #     mutate(OverallSurvivalDuration = if_else(OverallSurvivalDuration < 0, 0, OverallSurvivalDuration))

    # Disease-Specific Survival
    rawData <- rawData %>% mutate(DiseaseSpecificSurvivalDuration = deathDuration)

    # Progression Free Survival
    rawData <- rawData %>%
        mutate(ProgressionFreeSurvivalDuration = if_else(!is.na(recurrentDuration), recurrentDuration, deathDuration))

    # Local Progression Free Survival
    rawData <- rawData %>%
        mutate(LocalProgressionFreeSurvivalDuration = LocalControlRateDuration)

    # Locoregional Progression Free Survival
    rawData <- rawData %>%
        mutate(LocoregionalProgressionFreeSurvivalDuration = LocoregionalControlRateDuration)

    # Distant Metastase Free Survival
    rawData <- rawData %>%
        mutate(DistantMetastaseFreeSurvivalDuration = DistantMetastaseControlRateDuration)


    return(rawData)
}


.extendStatus <- function(rawData) {
    # 增加存活狀態
    rawData <- rawData %>%
        # 先看癌登
        mutate(VitalStatus = TCDB_VITSS) %>%
        # 院外死亡日，註記死亡
        mutate(VitalStatus = if_else('' != LV_OUTDIE, as.integer(0), VitalStatus)) %>%
        # 院內死亡日，註記死亡
        mutate(VitalStatus = if_else(!is.na(PBAS_PLSTVDT) & 3 == PBAS_PMEDSTAT, as.integer(0), VitalStatus)) %>%
        # 沒有存活的資料，假設他還活著
        mutate(VitalStatus = if_else(is.na(VitalStatus), as.integer(1), VitalStatus))

    # 增加存活狀態 in VGH
    rawData <- rawData %>%
        # 院內存活狀態
        mutate(VitalStatusInVgh = if_else(!is.na(PBAS_PLSTVDT) & 3 == PBAS_PMEDSTAT, as.integer(0), as.integer(1)))

    # 結束時間 in VGH then 國建署
    rawData <- rawData %>%
        # 院內註記日期
        mutate(VitalStatusInVghThenOutdie = if_else(!is.na(PBAS_PLSTVDT) & 3 == PBAS_PMEDSTAT, as.integer(0), as.integer(1))) %>%
        # 醫師註記院外死亡
        mutate(VitalStatusInVghThenOutdie = if_else('' != LV_OUTDIE, as.integer(0), VitalStatusInVghThenOutdie))

    # 結束時間 in VGH then 國建署 then VGH
    rawData <- rawData %>%
        # 院內註記日期
        mutate(VitalStatusInVghThenOutdieThenVgh = if_else(!is.na(PBAS_PLSTVDT) & 3 == PBAS_PMEDSTAT, as.integer(0), as.integer(1))) %>%
        # 醫師註記院外死亡
        mutate(VitalStatusInVghThenOutdieThenVgh = if_else('' != LV_OUTDIE, as.integer(0), VitalStatusInVghThenOutdieThenVgh)) %>%
        # 院內註記日期
        mutate(VitalStatusInVghThenOutdieThenVgh = if_else(dateConversion(LV_OUTDIE) < dateConversion(PBAS_PLSTVDT), if_else(!is.na(PBAS_PLSTVDT) & 3 == PBAS_PMEDSTAT, as.integer(0), as.integer(1)), VitalStatusInVghThenOutdieThenVgh))


    rawData <- rawData %>%
        mutate(deathByCancer = if_else(TCDB_COAOD == TCDB_PRIST, '1', '0'))

    rawData <- rawData %>%
        mutate(deathByOthers = if_else((TCDB_COAOD != TCDB_PRIST), '1', '0')) %>%
        mutate(deathByOthers = if_else(("0000" == TCDB_COAOD), '0', deathByOthers)) %>%
        mutate(deathByOthers = if_else(("0" == TCDB_COAOD), '0', deathByOthers))


    rawData <- rawData %>%
        mutate(tRecurrence = if_else(PBC_TFLAG %in% c('1', '2'), '1', '0'))

    rawData <- rawData %>%
        mutate(nRecurrence = if_else(PBC_NFLAG %in% c('1', '2'), '1', '0'))

    rawData <- rawData %>%
        mutate(mRecurrence = if_else(PBC_MFLAG %in% c('1', '2'), '1', '0'))


    # 計算 event (event = 1, censor = 0)
    rawData <- rawData %>%
        mutate(LocalControlRateEvent = if_else(('1' == tRecurrence), '1', '0'))

    rawData <- rawData %>%
        mutate(LocoregionalControlRateEvent = if_else(('1' == tRecurrence | '1' == nRecurrence), '1', '0'))

    rawData <- rawData %>%
        mutate(DistantMetastaseControlRateEvent = if_else(('1' == mRecurrence), '1', '0'))

    rawData <- rawData %>%
        mutate(NodalControlRateEvent = if_else(('1' == nRecurrence), '1', '0'))

    # 理論上要用這個方法
    # rawData <- rawData %>%
    #     mutate(OverallSurvivalEvent = if_else(('1' == deathByCancer | '1' == deathByOthers), '1', '0'))

    # 用舊的年報的方法，反向就正常了
    rawData <- rawData %>%
        mutate(OverallSurvivalEvent = if_else('0' == VitalStatusInVghThenOutdie, '1', '0'))

    rawData <- rawData %>%
        mutate(DiseaseSpecificSurvivalEvent = if_else(('1' == deathByCancer), '1', '0'))

    rawData <- rawData %>%
        mutate(ProgressionFreeSurvivalEvent = if_else(('1' == deathByCancer | '1' == deathByOthers | '1' == tRecurrence | '1' == nRecurrence | '1' == mRecurrence), '1', '0'))

    rawData <- rawData %>%
        mutate(LocalProgressionFreeSurvivalEvent = if_else(('1' == deathByCancer | '1' == deathByOthers | '1' == tRecurrence), '1', '0'))

    rawData <- rawData %>%
        mutate(LocoregionalProgressionFreeSurvivalEvent = if_else(('1' == deathByCancer | '1' == deathByOthers | '1' == tRecurrence | '1' == nRecurrence), '1', '0'))

    rawData <- rawData %>%
        mutate(DistantMetastaseFreeSurvivalEvent = if_else(('1' == deathByCancer | '1' == deathByOthers | '1' == mRecurrence), '1', '0'))


    return(rawData)
}


.extendClassify <- function(rawData) {
    # 原本的分類
    rawData <- rawData %>% mutate(VGHTPE_GROUP = trimws(PBC_PGROUP1))

    # 年齡區間
    rawData$AgeRange <- floor(rawData$TCDB_DISAGE/10)
    rawData$AgeRange[rawData$AgeRange < 4] <- 3
    rawData$AgeRange[rawData$AgeRange >= 8] <- 8
    factor(rawData$AgeRange)

    # 因為我們有兒癌團隊好些年了，他們一直想看存活率
    # 分Child (< 19歲) / adult (>= 19歲) / all ages
    rawData <- rawData %>% mutate(IsChild = if_else(TCDB_DISAGE >= 19, 0, 1))

    return(rawData)
}

.extendRename <- function(rawData) {
    # 為了方便辨識而改名字
    ## for Dr. Shiao
    rawData <- rawData %>% mutate(histno = TCDB_PHISTNUM)
    rawData <- rawData %>% mutate(sex = TCDB_SEX) # 1.5
    rawData <- rawData %>% mutate(age = TCDB_DISAGE) # 2.1

    rawData <- rawData %>% mutate(AJCC_7th_Stage = SimpleStage)
    rawData <- rawData %>% mutate(AJCC_7th_Stage_4C = Stage4C)


    ## for SMART 年報
    rawData <- rawData %>% mutate(START_DATE = StartDate)
    rawData <- rawData %>% mutate(SEX = TCDB_SEX) # 1.5
    rawData <- rawData %>% mutate(AGE = TCDB_DISAGE) # 2.1
    rawData <- rawData %>% mutate(AGE_RANGE = AgeRange)

    rawData <- rawData %>% mutate(MAIN_O3T = MainO3T)
    rawData <- rawData %>% mutate(IS_CHILD = IsChild)
    rawData <- rawData %>% mutate(IS_FIRST = if_else(IsFirstCancer, 1, 0))

    rawData <- rawData %>% mutate(STAGE = SimpleStage)
    rawData <- rawData %>% mutate(STAGE_4C = Stage4C)
    rawData <- rawData %>% mutate(STAGE_T = tStage)
    rawData <- rawData %>% mutate(STAGE_N = nStage)
    rawData <- rawData %>% mutate(STAGE_M = mStage)

    # rawData <- rawData %>% mutate( = )

    return(rawData)
}

.transferFactorVariable <- function(rawData) {
    # base table 的 類別變量 對應表
    baseTableFactorMapping <- read.csv(paste0(CODE_DIR, "libs/mapping/baseTableFactorMapping.csv"), stringsAsFactors=FALSE, fileEncoding="utf8")

    # 對全部的變數做考慮
    rawDataColumnNameList <- colnames(rawData)
    for (columnName in rawDataColumnNameList) {
        # 有在 類別變量 對應表 的變數，轉 factor
        if (1 == sum(columnName == baseTableFactorMapping$COLNAME)) {
            keyString <- baseTableFactorMapping[columnName == baseTableFactorMapping$COLNAME, "LISTKEY"]
            levelsList <- unlist(strsplit(keyString, "[$]"))

            valueString <- baseTableFactorMapping[columnName == baseTableFactorMapping$COLNAME, "LISTVAL"]
            labelsList <- unlist(strsplit(valueString, "[$]"))

            rawData[,columnName] <- factor(rawData[,columnName], levels = levelsList, labels = labelsList)
        }
    }

    return(rawData)
}

..makeMccCarData <- function(extendedRawData) {
    # 2016 年報版本的欄位
    chooseColumn <- c(
        "TCDB_SEX",
        "AgeRange",
        "Severity",
        "MainO3T",
        "AJCC7_ID",
        "AJCC7_SUBID",
        "VGHTPE_GROUP",
        "IsChild",
        "IsFirstCancer",
        "FullStage",
        "SimpleStage",
        "Stage4C",
        "Severity",
        "TCDB_DLCOD",
        "StartDate",
        "EndDate",
        "EndDateInVgh",
        "EndDateInVghThenOutdie",
        "EndDateInVghThenOutdieThenVgh",
        "VitalStatus",
        "VitalStatusInVgh",
        "VitalStatusInVghThenOutdie",
        "VitalStatusInVghThenOutdieThenVgh"
    )

    allCancerRawDataLocal <- extendedRawData[,chooseColumn]
    allCancerRawDataLocal <- .transferFactorVariable(allCancerRawDataLocal)

    return(allCancerRawDataLocal)
}


..makeSmartCarData <- function(extendedRawData) {
    # SMART 版本的年報欄位
    chooseColumn <- c(
        "UNIQ_UUID",
        "START_DATE",
        "SEX",
        "AGE",
        "AGE_RANGE",
        "AJCC7_ID",
        "AJCC7_SUBID",
        "VGHTPE_GROUP",
        "MAIN_O3T",
        "IS_CHILD",
        "IS_FIRST",
        "STAGE",
        "STAGE_4C",
        "STAGE_T",
        "STAGE_N",
        "STAGE_M",
        "LocalControlRateDuration",
        "LocalControlRateEvent",
        "LocoregionalControlRateDuration",
        "LocoregionalControlRateEvent",
        "DistantMetastaseControlRateDuration",
        "DistantMetastaseControlRateEvent",
        "NodalControlRateDuration",
        "NodalControlRateEvent",
        "OverallSurvivalDuration",
        "OverallSurvivalEvent",
        "DiseaseSpecificSurvivalDuration",
        "DiseaseSpecificSurvivalEvent",
        "ProgressionFreeSurvivalDuration",
        "ProgressionFreeSurvivalEvent",
        "LocalProgressionFreeSurvivalDuration",
        "LocalProgressionFreeSurvivalEvent",
        "LocoregionalProgressionFreeSurvivalDuration",
        "LocoregionalProgressionFreeSurvivalEvent",
        "DistantMetastaseFreeSurvivalDuration",
        "DistantMetastaseFreeSurvivalEvent"
    )

    smartCarData <- extendedRawData[,chooseColumn]

    return(smartCarData)
}



# 執行主程式
extendedRawData <- ..main()
saveRDS(extendedRawData, file = paste0(DATASET_DIR, "extendedRawData", DATA_DATE, ".rds"))
#print(rawData[2,'PHISTNUM', 'PBCSEQNO'])

# 製作 MCC 版本的年報資料
allCancerRawDataLocal <- ..makeMccCarData(extendedRawData)
saveRDS(allCancerRawDataLocal, file = paste0(DATASET_DIR, "allCancerRawDataLocal", DATA_DATE, ".rds"))
write.csv(allCancerRawDataLocal, file = paste0(DATASET_DIR, "allCancerRawDataLocal", DATA_DATE, ".csv"), row.names = FALSE)

# 製作 SMART 版本的年報資料
smartCarData <- ..makeSmartCarData(extendedRawData)
saveRDS(smartCarData, file = paste0(DATASET_DIR, "smartCarData", DATA_DATE, ".rds"))
write.csv(smartCarData, file = paste0(DATASET_DIR, "smartCarData", DATA_DATE, ".csv"), row.names = FALSE)


print("Please RUN THE COMMAND in MobaXTerm")
#print("  sh D:/SMART/VGHCAR/sh/copyUuidList.sh        ")
#print(paste0("  sh D:/SMART/VGHCAR/sh/uploadDataToServer.sh bhlee ", DATA_DATE, "        "))



### 舊的注解，仍具參考價值 ###
# 延伸欄位並沒有因為檔案變小，執行速度變快。
# 所以還是要全部欄位一起轉，然後產生年報用的小檔案。

# 目前在 Windows 上面轉資料的方法：
## Step 1:
# 一開始要用 SQL 撈出來的資料 (AllCancerUseful.sql)
# 目前已經設定在 Trinity 做排程，
# 除非 Trinity 當天沒做出來，不然不需要自己運作，
# 目前設定 Trinity 每天 23:45 開始撈資料，
# 並且產生 AllCancer_日期.csv 到 10.97.249.120 的 /tmp 資料夾
# 產檔路徑已更改 IP 為 10.97.235.41 
# 請把 AllCancer_日期.csv 複製到 DATASET_DIR
# 並且設定這隻檔案的 DATA_DATE 變數

## Step 2:
# 打開 RStudio
# 打開這隻 transferData.r
# 按下 Source 的按鈕讓他跑

## Step 3:
# 做完複製特定檔案到特定資料夾：
#   uuidList.csv => CODE_DIR/UUID/uuidList.csv

## Step 4:
# 把 smartCarData日期.csv 上傳到 10.97.249.120
# 資料夾路徑是 /docker/vghcar/dataset

## 附註：
# 這隻程式尚待完成的事項 (非必要，只是做完比較方便而已)
# [ ] 將 DATA_DATE 改成參數，執行 source 的時候從外面傳進來。
# [ ] 放到 Linux Server 上面，讓他可以用一行指令執行。
# [ ] 用 Trinity 做排程，接在撈資料的後面放 Step，讓他自動執行。
