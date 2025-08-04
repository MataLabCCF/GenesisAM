import os
import gzip
import argparse

def isToAdd(ID, covarData):
    if ID in covarData:
        return True
    else:
        return False


def execute(command):
    print("========================================================")
    print(command)
    os.system(command)
    #input()

def createVCF(begin, end, msp, outputFolder, outputName, covarName):
    if covarName:
        covarData = open(covarName)
        header = True
        dictIDs = {}
        for line in covarData:
            if header:
                header = False
            else:
                IID = line.split()[0]
                dictIDs[IID] = 0



    for chrom in range(begin, end+1):
        mspWithChrom = msp.replace("*", str(chrom))
        print(f"Building the VCF for chromosome {chrom} ({mspWithChrom})")
        fileMSP = open(mspWithChrom)
        lineCount = 0
        dictAnc = {}
        vcfOut = {}

        for line in fileMSP:
            lineCount = lineCount+1

            if lineCount == 1:
                ancs=line.strip().replace("#Subpopulation order/codes: ", "").split()
                for anc in ancs:
                    POP, ID = anc.split("=")
                    dictAnc[int(ID)] = POP
                for ID in dictAnc:
                    POP = dictAnc[ID]
                    vcfOut[ID] = open(f"{outputFolder}/LA_Results/{outputName}_{POP}_chrom{chrom}.vcf", "w")

            elif lineCount == 2:
                headerVCF = ("##fileformat=VCFv4.2\n##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">\n"
                          "##INFO=<ID=FBP,Number=1,Type=String,Description=\"Coordinate (in bp) of the first variant used in this window\">\n"
                          "##INFO=<ID=LBP,Number=1,Type=String,Description=\"Coordinate (in bp) of the last variant used in this window\">\n"
                          "##INFO=<ID=NV,Number=1,Type=String,Description=\"Number of variants used in this window\">\n"
                          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Local ancestry genotype\">\n"
                          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

                headerMSP = line.strip().split("\t")
                for i in range(6, len(headerMSP), 2):
                    ID = headerMSP[i].replace(".0", "")
                    if covarName != "":
                        if isToAdd(ID, dictIDs):
                            headerVCF = f"{headerVCF}\t{ID}"
                    else:
                        headerVCF = f"{headerVCF}\t{ID}"

                for ID in vcfOut:
                    vcfOut[ID].write(f"{headerVCF}\n")

            else:
                split = line.strip().split()
                CHROM = f"chr{chrom}"
                FBP = split[1]
                LBP = split[2]
                NV= split[5]
                POS = int((int(FBP)+int(LBP))/2)
                REF = "T" #oTher ancestry
                ALT = "A" #Ancestry
                QUAL = FILTER = "."
                IDVCF = f"chrom{chrom}:Window{lineCount-2}"
                INFO = f"FBP={FBP};LBP={LBP};NV={NV}"
                FORMAT = "GT"

                for ID in vcfOut:
                    vcfOut[ID].write(f"{CHROM}\t{POS}\t{IDVCF}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}")

                for i in range(6, len(split), 2):
                    sampleID = headerMSP[i].replace(".0", "")
                    if covarName != "":
                        if isToAdd(sampleID, dictIDs):
                            for ID in vcfOut:
                                GT = getGT(split[i], split[i+1], str(ID))
                                vcfOut[ID].write(f"\t{GT}")
                    else:
                        for ID in vcfOut:
                            GT = getGT(split[i], split[i+1], str(ID))
                            vcfOut[ID].write(f"\t{GT}")


                for ID in vcfOut:
                    vcfOut[ID].write(f"\n")

        for ID in vcfOut:
            vcfOut[ID].close()
    return dictAnc

def getGT(A1, A2, AncID):
    if A1 == AncID:
        GT = "1"
    else:
        GT = "0"
    if A2 == AncID:
        GT = f"{GT}|1"
    else:
        GT = f"{GT}|0"
    return GT

def createAMFile(dictAnc,outputFolder, covarName):
    print(f"Creating AM: {outputFolder}/AM.R")
    fileAM = open(f"{outputFolder}/AM.R", "w")

    #Headers
    fileAM.write("library(GENESIS)\n")
    fileAM.write("library(GWASTools)\n")
    fileAM.write("library(SNPRelate)\n")
    fileAM.write("library(SeqArray)\n")
    fileAM.write("options <- commandArgs(trailingOnly = TRUE)\n")

    #GetOption
    fileAM.write(f"\n\n")
    for i in dictAnc:
        POP = dictAnc[i]
        fileAM.write(f"{POP} = options[{i+1}]\n")
    fileAM.write(f"kingFile = options[{i+2}]\n")
    fileAM.write(f"outputName = options[{i+3}]\n")
    fileAM.write(f"covarName = options[{i+4}]\n")

    #GDS conversion
    fileAM.write(f"\n\n")
    for i in dictAnc:
        POP = dictAnc[i]
        fileAM.write(f"{POP}GDS = paste(outputName, \"_{POP}.gds\", sep = \"\")\n")
        fileAM.write(f"snpgdsVCF2GDS({POP}, {POP}GDS, verbose=TRUE)\n")

    #Read King
    fileAM.write(f"\n\n")
    fileAM.write(f"GRM <- kingToMatrix(kingFile, estimator=\"Kinship\")\n")

    # Read GDS
    fileAM.write(f"\n\n")
    for i in dictAnc:
        POP = dictAnc[i]
        if i == 0:
            fileAM.write(f"files <- list({POP.lower()}={POP}GDS")
        else:
            fileAM.write(f", {POP.lower()}={POP}GDS")
    fileAM.write(f")\n")
    fileAM.write(f"gdsList <- lapply(files, GdsGenotypeReader)\n")

    #Creating ScanAnnot
    fileAM.write(f"\n\n")
    fileAM.write(f"scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=getScanID(gdsList[[1]]), stringsAsFactors=FALSE))\n")

    #Add covar on scanAnnot
    fileAM.write(f"\n\n")
    fileAM.write(f"covar = read.table(covarName, sep = \"\t\", header = T)\n")
    fileAM.write(f"scanID = getScanID(scanAnnot)\n")
    fileAM.write(f"X = scanID[which(!covar$IID %in% scanID)]\n")
    fileAM.write(f"covarInCommon = covar[!covar$IID %in% X,]\n")
    fileAM.write(f"\n")


    covar = open(covarName)
    for line in covar:
        header = line.strip().split()
        break
    covar.close()

    covarParam = ""
    for field in header:
        if field.upper() == "STATUS" or field.upper() == "DISEASE":
            fileAM.write(f"scanAnnot$PHENO = covarInCommon${field}\n")
        elif field.upper() == "SEX" or field.upper() == "GENDER":
            if covarParam == "":
                covarParam=f"\"{field}\""
            else:
                covarParam=f"{covarParam},\"{field}\""
            fileAM.write(f"covar${field} = as.factor(covar${field})\n")
            fileAM.write(f"scanAnnot$SEX = covarInCommon${field}\n")
        else:
            if field != "IID" and field != "SampleID" and field != "ID":
                if covarParam == "":
                    covarParam=f"\"{field}\""
                else:
                    covarParam=f"{covarParam},\"{field}\""
                fileAM.write(f"scanAnnot${field} = covarInCommon${field}\n")

    # Fit Null
    fileAM.write(f"\n\n")
    fileAM.write(f"nullModel = fitNullModel(scanAnnot, outcome=\"PHENO\", cov.mat=GRM, family=\"binomial\",covars=c({covarParam}))\n")
    fileAM.write(f"genoDataList <- lapply(gdsList, GenotypeData, scanAnnot=scanAnnot)\n")
    fileAM.write(f"genoIterators <- lapply(genoDataList[1:{len(dictAnc)-1}], GenotypeBlockIterator)\n")

    # Association
    fileAM.write(f"myassoc <- admixMap(genoIterators, nullModel)\n")
    fileAM.write(f"outputDF = subset(myassoc, select=-c(variant.id))\n")
    fileAM.write(f"write.csv(outputDF, paste(outputName, \".csv\", sep=\"\"), row.names = FALSE, quote = FALSE)\n")

    for i in dictAnc:
        POP = dictAnc[i]
        fileAM.write(f"genoIterators <- lapply(genoDataList[{i+1}], GenotypeBlockIterator)\n")
        fileAM.write(f"myassoc <- admixMap(genoIterators, nullModel)\n")
        fileAM.write(f"outputDF = subset(myassoc, select=-c(variant.id))\n")
        fileAM.write(f"write.csv(outputDF, paste(outputName, \"_{POP}.csv\", sep=\"\"), row.names = FALSE, quote = FALSE)\n")
    fileAM.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Admixture Mapping with GENESIS')

    requiredGeneral = parser.add_argument_group("Required arguments for all steps")
    requiredGeneral.add_argument('-m', '--msp', help='MSP file with chromosome replaced by *',required=True)
    requiredGeneral.add_argument('-o', '--outputName', help='Name of output name', required=True)
    requiredGeneral.add_argument('-O', '--outputFolder', help='Name of output folder', required=True)
    requiredGeneral.add_argument('-b', '--begin', help='First chromosome (default = 1)', default=1, type=int)
    requiredGeneral.add_argument('-e', '--end', help='Last chromosome  (default = 22)', default=22, type=int)
    requiredGeneral.add_argument('-p', '--plinkFile', help='Plink prefix file to infer KING matrix', required= True)

    software = parser.add_argument_group("Required softwares")
    software.add_argument('--plink1', help='Plink 1 path (calculate number of independent tests)', required=True)
    #software.add_argument('--plink2', help='Plink 2 path (king matrix calculation)', required=True)
    software.add_argument('-R', '--rscript', help='Rscript path', required=True)
    software.add_argument('-K', '--king', help='King path', required=True)

    other = parser.add_argument_group("Other arguments")
    #other.add_argument('-v', '--vcfOnly', help='Just convert MSP to VCF', required=False, default=False, action="store_true")
    other.add_argument('-c', '--covar', help='Covar file to build a MSP with only covar data', default="", required = False)


    args = parser.parse_args()
    execute(f"mkdir {args.outputFolder}")
    execute(f"mkdir {args.outputFolder}/LA_Results")
    execute(f"mkdir {args.outputFolder}/PLINK_Results")
    execute(f"mkdir {args.outputFolder}/AM_Results")


    #Make VCS
    dictAnc = createVCF(args.begin, args.end, args.msp, args.outputFolder, args.outputName, args.covar)

    #Prepare PLINK to KING relationship
    execute(f"cp {args.plinkFile}.bed {args.outputFolder}/PLINK_Results/{args.outputName}_KING.bed")
    execute(f"cp {args.plinkFile}.bim {args.outputFolder}/PLINK_Results/{args.outputName}_KING.bim")
    newFamFile = open(f"{args.outputFolder}/PLINK_Results/{args.outputName}_KING.fam", "w")
    famFile = open(f"{args.plinkFile}.fam")
    for line in famFile:
        ID = line.split()[1]
        sex = line.split()[4]
        newFamFile.write(f"{ID}\t{ID}\t0\t0\t{sex}\t-9\n")
    newFamFile.close()
    famFile.close()

    #King inference and GRM setup (set cutoff < 0.03 (fourth degree) as 0
    execute(f"{args.king} -b {args.outputFolder}/PLINK_Results/{args.outputName}_KING.bed --kinship --prefix {args.outputFolder}/PLINK_Results/{args.outputName}_KING")
    kingFile = open(f"{args.outputFolder}/PLINK_Results/{args.outputName}_KING.kin0")
    kingFileFinal = open(f"{args.outputFolder}/PLINK_Results/{args.outputName}_KING_Fixed.kin0", "w")

    header = True
    for line in kingFile:
        if header:
            kingFileFinal.write(line)
            header = False
        else:
            FID1, ID1, FID2, ID2, N_SNP, HetHet, IBS0, Kinship= line.strip().split()
            if float(Kinship) < 0.03:
                Kinship = 0.00
            kingFileFinal.write(f"{FID1}	{ID1}	{FID2}	{ID2}	{N_SNP}	{HetHet}	{IBS0}	{Kinship}\n")
    kingFileFinal.close()
    kingFile.close()

    #GENESIS LA
    createAMFile(dictAnc, args.outputFolder, args.covar)
    VCFPath = f"{args.outputFolder}/LA_Results/"
    for chrom in range(args.begin, args.end+1):

        commandLine = f"{args.rscript} {args.outputFolder}/AM.R"
        for i in dictAnc:
            POP = dictAnc[i]
            commandLine = f"{commandLine} {VCFPath}{args.outputName}_{POP}_chrom{chrom}.vcf"

        commandLine = (f"{commandLine} {args.outputFolder}/PLINK_Results/{args.outputName}_KING_Fixed.kin0 "
                       f"{args.outputFolder}/AM_Results/{args.outputName}_CHROM{chrom} {args.covar}")
        execute(commandLine)


    # Number of independent tests
    # Calculation based on Gao X, Becker LC, Becker DM, Starmer J, Province MA (2009) Avoiding the high
    # Bonferroni penalty in genome-wide association studies. Genetic Epidemiology

    #PVT = open(f"{args.outputFolder}/PLINK_Results/{args.outputName}_PVT.txt", "w")
    for ID in dictAnc:
        POP = dictAnc[ID]
        fileOut = open(f"{args.outputFolder}/LA_Results/Dosage_{POP}_chromAll.tsv", "w")
        for chrom in range(args.begin, args.end+1):
            VCF = f"{args.outputFolder}/LA_Results/{args.outputName}_{POP}_chrom{chrom}.vcf"
            fileIn = open(VCF)

            header = True
            for line in fileIn:
                if header:
                    if line.startswith("#CHROM"):
                        header = False
                else:
                    split = line.strip().split()

                    for i in range(9, len(split)):
                        A1, A2 = split[i].split("/")
                        dosage = int(A1) + int(A2)
                        if i == 9:
                            fileOut.write(f"{dosage}")
                        else:
                            fileOut.write(f"\t{dosage}")
                    fileOut.write("\n")

        fileOut.close()
        print(f"Rscript simpleM.R {args.outputFolder}/LA_Results/Dosage_{POP}_chromAll.tsv > {args.outputFolder}/LA_Results/PVT_{POP}.txt")
        os.system(f"Rscript simpleM.R {args.outputFolder}/LA_Results/Dosage_{POP}_chromAll.tsv > {args.outputFolder}/LA_Results/PVT_{POP}.txt")