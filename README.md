# CNVqc Metrics Compilation

This script is designed to compile and analyze CNVqc metrics based on the specified platform and date range. It provides a comprehensive overview of CNV quality control metrics, coverage statistics, and other relevant data.

## Dependencies

Ensure you have the required dependencies installed. Refer to `dependency.yml` for a list of necessary packages and libraries.

## Usage

To run the script, use the following command:

```
auto_cnv_qc.py --platform {Nextseq,Novaseq} --date DATE --panel {Exome,Neuropathy,LynchHRD} [--outdir OUTDIR]
```

### Arguments

- `--platform {Nextseq,Novaseq}`: Specify the sequencing platform. Choose between 'Nextseq' and 'Novaseq'.
  
- `--date DATE`: Provide the date prefix for runs, e.g., `2205`.
  
- `--panel {Exome,Neuropathy,LynchHRD}`: Define the panel. Options include 'Exome', 'Neuropathy', and 'LynchHRD'.
  
- `--outdir OUTDIR`: (Optional) Set the path to the output directory. If not provided, the current working directory will be used.

## Features

- **Run ID Retrieval**: The script fetches run IDs based on the platform, prefix, and panel combination.
  
- **QC Status Extraction**: Extracts the QC status from the provided QC file.
  
- **CNV Metrics Compilation**: Gathers CNV metrics from the metrics file.
  
- **CNV Calls Count**: Counts the number of CNV calls in the cnv.bed file.
  
- **Coverage QC Status**: Determines the coverage QC status based on specific criteria.
  
- **Coverage Statistics**: Computes coverage statistics for the specified runs.
  
- **Data Merging**: Merges CNV QC and Coverage Statistics data for comprehensive analysis.
  
- **Visualization**: Generates regression plots for specific metrics for visual analysis.

## Output

The script will generate a series of output files, including:

- CNV QC dataframes
- Coverage statistics dataframes
- Merged dataframes with CNV QC and coverage statistics
- Regression plots for specific metrics
