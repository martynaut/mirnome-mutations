import click
import pandas as pd
import os
import re
import numpy as np
from post_extraction_analysis import post_analyse


def filter_and_combine(output_folder, include_merger, include_filtering):

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(output_folder)]
    files = [file for sublist in files for file in sublist]

    results_files = []

    for file in files:
        if re.match('.*results_[a-zA-Z0-9]*\.csv$', file):
            results_files.append(file)

    if not include_merger and len(results_files) > 1:
        click.echo("Multiple algorithms used and merger not set. Closing")
        return 1

    if not include_filtering:
        to_merge = results_files
    else:
        to_merge = [file.replace('.csv', '_eval.csv') for file in results_files]

    if include_filtering:
        for file in results_files:
            post_analyse(file,
                         file.replace('.csv', '_eval.csv'))

    if include_merger:

        results_df = pd.DataFrame()

        for file in to_merge:
            temp_df = pd.read_csv(file)
            alg = re.search('results_([a-zA-Z0-9]+)', file)
            temp_df['alg'] = alg.group(1)
            results_df = pd.concat([results_df, temp_df])
    else:
        results_df = pd.read_csv(results_files[0])

    columns = ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format',
               'tumor', 'normal', 'indiv_name', 'indiv_id', 'sample_id_tumor_name',
               'sample_id_tumor_aliQ', 'sample_id_normal_name',
               'sample_id_normal_aliQ', 'norm_ref_count', 'norm_alt_count',
               'tumor_ref_count', 'tumor_alt_count', 'BQ_ref_tum', 'BQ_alt_tum',
               'BQ_ref_norm', 'BQ_alt_norm', 'QSS_ref_tum', 'QSS_alt_tum',
               'QSS_ref_nor', 'QSS_alt_nor', 'SSC', 'SPV', 'eval', 'alg']

    if include_filtering:
        results_df = results_df[results_df['eval']]
        results_df = results_df[columns]
    else:
        columns.remove('eval')
        results_df = results_df[columns]

    results_df['control:mut/norm'] = results_df['norm_alt_count'].div(results_df['norm_ref_count'])
    results_df['control:mut/norm'] = results_df['control:mut/norm'].astype(float)
    results_df['tumor:mut/norm'] = results_df['tumor_alt_count'].div(results_df['tumor_ref_count'])
    results_df['tumor:mut/norm'] = results_df['tumor:mut/norm'].astype(float)
    results_df['ratio'] = results_df['tumor:mut/norm'].div(results_df['control:mut/norm'])

    results_df.replace(np.inf, 0, inplace=True)

    results_df.to_csv(output_folder + '/all_mutations_filtered.csv',
                      sep=',',
                      index=False)
    return 0


@click.command()
@click.argument('output_folder')
@click.option('--include_merger', '-m')
@click.option('--include_filtering', '-f')
def main(output_folder,
         include_merger='',
         include_filtering=''
         ):
    filter_and_combine(output_folder, include_merger, include_filtering)


if __name__ == "__main__":
    main()
