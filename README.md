# qc
this script compiles metrics that are important for CNVqc with platfrom and date range information

## DEPENDENCIES
dependency.yml

## auto_cnv_qc.py
usage: auto_cnv_qc.py [-h] --platform {Nextseq,Novaseq} --date DATE --panel {Exome,Neuropathy,LynchHRD} [--outdir OUTDIR]

## args
```
usage: auto_cnv_qc.py [-h] --platform {Nextseq,Novaseq} --date DATE --panel
                      {Exome,Neuropathy,LynchHRD} [--outdir OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  --platform {Nextseq,Novaseq}
                        [Input] Platform
  --date DATE           [Input] date prefix for runs: eg., (2205)
  --panel {Exome,Neuropathy,LynchHRD}
                        [Input] Panel
  --outdir OUTDIR       [Optional] Path to output directory
```