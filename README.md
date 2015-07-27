# cufflinks-cibersort
Running a collection of cufflinks gene FPKM files through CIBERSORT


### Dependencies

- R 
- `Rserve` for R
- CIBERSORT jar

### Running example

```sh

python 
    cufflinks-cibersort.py 
        --cufflinks-files cufflinks-all-samples/*/gene*fpkm* 
        --mixture-file CIBERSORT_package/LM22.txt 
        --cibersort-jar CIBERSORT_package/CIBERSORT.jar 
        --output cibersort-output.tsv 
        --print-top-n
```


### Usage
```sh
python cufflinks-cibersort.py
```

```sh
usage: cufflinks-cibersort.py [-h] --cufflinks-files CUFFLINKS_FILES
                              [CUFFLINKS_FILES ...] [--gep-output GEP_OUTPUT]
                              --output OUTPUT
                              [--cibersort-output CIBERSORT_OUTPUT]
                              --mixture-file MIXTURE_FILE --cibersort-jar
                              CIBERSORT_JAR [--min-fpkm MIN_FPKM]
                              [--print-top-n] [--n N]
cufflinks-cibersort.py: error: the following arguments are required: --cufflinks-files, --output, --mixture-file, --cibersort-jar
```

