#!/usr/bin/python

import sys
import argparse
import json
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


"""
    This script is used to generate CNVqc metrics with platfrom and date range information

"""


def get_args():
    """Parse Arguments"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--platform",
                        action="store",
                        dest="platform",
                        choices=['Nextseq', 'Novaseq'],
                        required=True,
                        help="[Input] Platform")
    parser.add_argument("--date",
                        action="store",
                        dest="date",
                        required=True,
                        help="[Input] date prefix for runs: eg., (2205)")
    parser.add_argument("--panel",
                        action="store",
                        dest="panel",
                        choices=['Exome', 'Neuropathy', 'LynchHRD'],
                        required=True,
                        help="[Input] Panel")
    parser.add_argument("--outdir",
                        action="store",
                        dest="outdir",
                        type=Path,
                        default=Path.cwd(),
                        help="[Optional] Path to output directory")

    return parser.parse_args()


def get_run_ids(platform, prefix, panel):
    """ Get Run IDs for platform, prefix and panel combination """

    run_ids = []
    path_map = {
        "Exome": f"{prefix}*/Aligned_Exomes_SGE/Project_98_Exome",
        "Neuropathy": f"{prefix}*/Aligned_Panel_9801_SGE/Project_9801_Neuropathy",
        "LynchHRD": f"{prefix}*/Aligned_Panel_283_SGE/Project_283_lynchHRD"
    }
    root_dir = Path(f'/NGS/{platform}')
    panel_path = path_map[panel]
    if 'Novaseq' in platform:
        panel_path = panel_path.replace('_98_', '_298_')
    panel_dirs = root_dir.glob(panel_path)
    for run_path in panel_dirs:
        run_ids.append(run_path.parent.parent.stem)
    return run_ids


def get_qc_status(qc_file):
    """ Get QC status from QC file """

    qc_status = ''
    with open(qc_file, 'r') as qc_fh:
        for line in qc_fh:
            if line.startswith('#QC'):
                qc_status = line.strip().split(" ")[1]
    return qc_status


def get_cnv_metrics(metrics_file):
    """ Get CNV metrics from the metrics file """

    with open(metrics_file, 'r') as met_fh:
        for line in met_fh:
            if not line.startswith('sample'):
                return list(map(float, line.strip().split("\t")[2:]))


def get_cnv_calls(bed_file):
    """ Get number of CNV calls in cnv.bed file """

    bdf = pd.read_csv(bed_file,
                      comment='#',
                      sep='\t',
                      names=[
                          '#Chromosome', 'start', 'end', 'type', 'genes',
                          'coverage ratio', 'copy number', 'method'
                      ])
    return [len(bdf)]


def get_cnv_qc(run_ids, platform, panel, outdir):
    """ Get CNV QC, returns QC dataframe """

    path_map = {
        "Exome": "Aligned_Exomes_SGE/Project_98_Exome",
        "Neuropathy": "Aligned_Panel_9801_SGE/Project_9801_Neuropathy",
        "LynchHRD": "Aligned_Panel_283_SGE/Project_283_lynchHRD"
    }

    cnv_qc_df = pd.DataFrame()
    for run_id in run_ids:
        run_path = Path(f'/NGS/{platform}/{run_id}/{path_map[panel]}')
        if 'Novaseq' in platform:
            run_path = Path(str(run_path).replace('_98_', '_298_'))
        summary_file = outdir.joinpath(run_id + '_sample_cnv_qc_summary.csv')
        df_header = [
            'sample', 'CNV QC', 'stdev', 'mad', 'iqr', 'bivar', 'cnv_calls'
        ]
        run_qc_data = []
        for run_dir in run_path.iterdir():
            if run_dir.stem.startswith('Sample_') and ('NEG' not in run_dir.stem) and (
                    'Neg' not in run_dir.stem):
                accession = run_dir.stem.replace('Sample_', '')
                qc_file = run_dir.joinpath(accession + '.cnv.qc')
                metrics_file = run_dir.joinpath(accession + '.metrics')
                bed_file = run_dir.joinpath(accession + '.cnv.bed')

                if qc_file.is_file():
                    qc_status = [get_qc_status(qc_file)]
                else:
                    qc_status = []
                if metrics_file.is_file():
                    cnv_metrics = get_cnv_metrics(metrics_file)
                else:
                    cnv_metrics = []
                if bed_file.is_file():
                    cnv_calls = get_cnv_calls(bed_file)
                else:
                    cnv_calls = [0]
                sample_qc_data = ([accession] + qc_status + cnv_metrics +
                                  cnv_calls)
                run_qc_data.append(sample_qc_data)
        run_qc_df = pd.DataFrame(run_qc_data, columns=df_header)
        run_qc_df.to_csv(summary_file, index=False)
        run_qc_df['run_id'] = run_id
        if cnv_qc_df.empty:
            cnv_qc_df = run_qc_df
        else:
            cnv_qc_df = pd.concat([cnv_qc_df, run_qc_df], ignore_index=True)
    return cnv_qc_df


def get_cov_status(row):
    """ Get Coverage QC status, returns PASS or FAIL """

    mean_qs = float(row['mean_QS'])
    q30 = float(row['q30'])
    x10 = float(row['X10'])
    x20 = float(row['X20'])
    mean_cov = float(row['mean_coverage'])
    status = "FAIL"

    if mean_qs >= 30 and q30 >= 75 and x10 >= 98 and x20 >= 97 and mean_cov >= 70:
        status = "PASS"

    return status


def get_cov_stats(run_ids, platform, panel, outdir):
    """ Get Coverage statistics, returns dataframe """

    path_map = {
        "Exome": f"Aligned_Exomes_SGE/Project_98_Exome",
        "Neuropathy": "Aligned_Panel_9801_SGE/Project_9801_Neuropathy",
        "LynchHRD": "Aligned_Panel_283_SGE/Project_283_lynchHRD"
    }

    col_filter = [
        'sample', 'specificity', 'mean_QS', 'perfect_index', 'q30', 'X10',
        'X20', 'X50'
    ]

    cov_stats_df = pd.DataFrame()
    for run_id in run_ids:
        run_path = Path(f'/NGS/{platform}/{run_id}/{path_map[panel]}')
        if 'Novaseq' in platform:
            run_path = Path(str(run_path).replace('_98_', '_298_'))
        summary_file = outdir.joinpath(run_id + '_sample_cnv_qc_summary.csv')
        run_cov_df = []
        coverage_file = outdir.joinpath(run_id + '_cov_stats.xlsx')
        for run_dir in run_path.iterdir():
            if run_dir.stem.startswith('Sample_') and ('NEG' not in run_dir.stem) and (
                    'Neg' not in run_dir.stem):
                accession = run_dir.stem.replace('Sample_', '')
                cov_file = run_path.joinpath(run_dir.stem,
                                             accession + '.coverage.json')
                if panel == "Neuropathy":
                    cov_file = Path(str(cov_file).replace('coverage.json','trim.coverage.json'))
                if cov_file.is_file():
                    with open(cov_file, 'r') as cov_fh:
                        sample_cov_df = pd.json_normalize(json.load(cov_fh))
                        if 'X1000' in sample_cov_df.columns:
                            filter_columns = col_filter + [
                                'X150', 'X500', 'X1000', 'yield',
                                'mean_coverage'
                            ]
                        else:
                            filter_columns = col_filter + [
                                'X100', 'pcr_dup_rate', 'mean_coverage'
                            ]
                        run_cov_df.append(sample_cov_df[filter_columns])
        run_cov_df = pd.concat(run_cov_df, axis=0)
        run_cov_df.fillna(0, inplace=True)
        run_cov_df['coverage_qc'] = run_cov_df.apply(get_cov_status, axis=1)
        run_cov_df.to_excel(coverage_file,
                            sheet_name='coverage_stats',
                            index=False,
                            header=True)
        run_cov_df.drop('coverage_qc', axis=1, inplace=True)
        run_cov_df['sample'] = run_cov_df['sample'].str.replace(
            '.trim', '', regex=True)  # Handle 9801
        if cov_stats_df.empty:
            cov_stats_df = run_cov_df
        else:
            cov_stats_df = pd.concat([cov_stats_df, run_cov_df],
                                     ignore_index=True)

    return cov_stats_df


def ngs_qc_old(row):
    """ Get NGS QC based on old criteria """

    status = "PASS"
    if float(row['X10']) < 0.95 or float(row['X20'] < 0.9
                                         or row['mean_coverage'] < 50):
        status = "FAIL"
    
    return status


def get_merged_df(cnv_qc_df, cov_stats_df):
    """ Merge CNV QC and Coverage Statistics dataframe """

    df_headers = [
        'run_id', 'sample', 'X10', 'X20', 'X50', 'X100', 'mean_coverage',
        'Exome_QC_old', 'Exome_QC_updated', 'CNV QC', 'cnv_calls', 'stdev', 'mad',
        'iqr', 'bivar', 'MadIQRBivarSum'
    ]
    merged_df = pd.merge(cnv_qc_df,
                         cov_stats_df,
                         left_on='sample',
                         right_on='sample')
    merged_df['Exome_QC_old'] = merged_df.apply(ngs_qc_old, axis=1)
    merged_df['Exome_QC_updated'] = merged_df.apply(get_cov_status, axis=1)
    merged_df['MadIQRBivarSum'] = merged_df['mad'] + merged_df[
        'iqr'] + merged_df['bivar']
    merged_df = merged_df[df_headers]
    return merged_df


def get_run_stats(df):
    """ Get Run statistics """

    passing_df = df[df['CNV QC'] == 'PASS']
    grouped_df = df.groupby('run_id').mean()
    grouped_df['num_samples'] = df.groupby('run_id').size().to_frame('counts')
    grouped_df['num_samples_passing_CNVqc'] = passing_df.groupby(
        'run_id').size().to_frame('counts')
    grouped_df['avg_cnv_calls_in_passing_samples'] = passing_df.groupby(
        'run_id')['cnv_calls'].mean().round(2)
    grouped_df.fillna(0, inplace=True)
    grouped_df['CNV_failure_rate'] = (
        1 -
        (grouped_df['num_samples_passing_CNVqc'] / grouped_df['num_samples']))
    grouped_df.reset_index(inplace=True)
    grouped_df.rename(columns={'cnv_calls': 'avg_cnv_calls'}, inplace=True)
    return grouped_df


def get_summary_df(df):
    """ Get Summary dataframe """

    passing_NGSqc_df = df[df['Exome_QC_updated'] == 'PASS'].copy()
    passing_NGSqc_df["group"] = 'Total'
    summary_df = passing_NGSqc_df.groupby('group').mean()
    summary_df['num_samples_passing_Exomeqc'] = passing_NGSqc_df.groupby(
        'group').size().to_frame('counts')
    summary_df.reset_index(inplace=True)
    return summary_df


def make_regplots(outdir, df, x_axis, y_axis):
    """ Make Plots """

    sns.regplot(x=x_axis, y=y_axis, data=df, line_kws={"color": "red"}, ci=95)
    plt.savefig(outdir.joinpath(f"{x_axis}VS{y_axis}"))
    plt.clf()


def compile_results(outdir, platform, prefix, panel, df):
    """ Compile results """

    writer = pd.ExcelWriter(outdir.joinpath(f'{platform}_{prefix}_{panel}_stats.xlsx'),
                            engine='xlsxwriter')

    # Make all_samples sheet
    df.to_excel(writer, sheet_name='all_samples', index=False)

    # Make clinical_samples sheet
    clinical_df = df[~df['sample'].str.contains('RD')]
    clinical_df.to_excel(writer, sheet_name='clinical_samples', index=False)

    # Make run_stats sheet
    run_stats_all_df = get_run_stats(df)
    run_stats_all_df.to_excel(writer, sheet_name='run_stats_all', index=False)

    run_stats_clinical_df = get_run_stats(clinical_df)
    run_stats_clinical_df.to_excel(writer,
                                   sheet_name='run_stats_clinical',
                                   index=False)

    # Make summary sheet
    summary_df = get_summary_df(clinical_df)
    summary_df.to_excel(writer, sheet_name='summary', index=False)

    writer.save()

    make_regplots(outdir, run_stats_clinical_df, "CNV_failure_rate",
                  "MadIQRBivarSum")
    make_regplots(outdir, run_stats_clinical_df,
                  "avg_cnv_calls_in_passing_samples", "MadIQRBivarSum")

def run():
    try:
        args = get_args()
        if not args.outdir.is_dir():
            try:
                args.outdir.mkdir(parents=True, exist_ok=True)
            except:
                print(
                    f"ERROR: Unable to create output directory {args.outdir}. Exiting!\n"
                )

    except KeyboardInterrupt:
        pass

    outdir = args.outdir.joinpath(args.platform + '_' + args.date + '_' +
                                  args.panel)
    outdir.mkdir(parents=True, exist_ok=True)

    run_ids = get_run_ids(args.platform, args.date, args.panel)

    if len(run_ids) > 0:
        cnv_qc_df = get_cnv_qc(run_ids, args.platform, args.panel, outdir)
        cov_stats_df = get_cov_stats(run_ids, args.platform, args.panel,
                                     outdir)
        merged_df = get_merged_df(cnv_qc_df, cov_stats_df)
        compile_results(outdir, args.platform, args.date, args.panel, merged_df)

    else:
        print("No runs to process. Exiting! \n")
        sys.exit(0)


if __name__ == "__main__":
    run()
    