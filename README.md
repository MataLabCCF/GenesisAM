# GenesisAM

![alt text](GENESIS_AM.png)

## Example of command line

```
python GENESIS_AM.py --msp /my/folder/to/local/ancestry/result/chr\*/query_results.msp -b 1 -e 22 \
--plinkFile /my/folder/to/plink/file/LPD_Phase2 --plink1 plink -R Rscript \
--KING king -o fileName -O /my/output/folder/
```

## Parameters

```
GenesisAM
-h, --help: shows this menu
options:

Required arguments:
-m, --msp: MSP file with chromosome replaced by *
-o, --outputName: Name of output file
-O, --outputFolder: Name of the output folder
-b, --begin: First chromosome (default = 1)
-e, --end: Last chromosome  (default = 22)
-p, --plinkFile: Plink prefix file to infer KING matrix

Optional arguments:
-c, --covar: Covar file to build a MSP with only covar data

Programs:
--plink1: Plink 1 path (calculate number of independent tests)
-R, --rscript: Rscript path
-K, --king: King path

```

## Acknowledgements
This work is supported by NIH Grant R01 1R01NS112499-01A1, MJFF Grant ID: 18298, ASAP-GP2 and Parkinson's Foundation

# Contact
Created by Thiago Peixoto Leal. PhD ([PEIXOTT@ccf.org](PEIXOTT@ccf.org) or [thpeixotol@hotmail.com](thpeixotol@hotmail.com))
