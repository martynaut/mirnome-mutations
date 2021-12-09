import click
import pandas as pd
import os
import re
import numpy as np
from post_extraction_analysis import post_analyse
from distinct_occur_helpers import concat_alg


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
        for file in results_files:
            post_analyse(file,
                         file.replace('.csv', '_eval.csv'))

    if include_merger:
        results_df = pd.DataFrame()

        for file in to_merge:
            temp_df = pd.read_csv(file)
            alg = re.search('results_([a-zA-Z0-9]+)', file)
            temp_df['alg'] = alg.group(1)
            if include_filtering:
                temp_df = temp_df[temp_df['eval']]
            results_df = pd.concat([results_df, temp_df])

            results_df = results_df.groupby(['chrom', 'pos',
                                             'indiv_name', 'ref', 'alt']).agg({'alg': concat_alg,
                                                                               'norm_ref_count': sum,
                                                                               'norm_alt_count': sum,
                                                                               'tumor_ref_count': sum,
                                                                               'tumor_alt_count': sum
                                                                               }).reset_index()

    else:
        results_df = pd.read_csv(to_merge[0])
        alg = re.search('results_([a-zA-Z0-9]+)', to_merge[0])
        try:
            results_df['alg'] = alg.group(1)
        except AttributeError:
            results_df['alg'] = 'unknown'

    columns = ['chrom', 'pos', 'ref', 'alt', 'indiv_name', 'norm_ref_count', 'norm_alt_count',
               'tumor_ref_count', 'tumor_alt_count', 'alg']

    if results_df.shape[0] == 0:
        print('no mutations found')
        return 1

    results_df = results_df[columns]
    results_df.replace(np.inf, 0, inplace=True)
    results_df.to_csv(output_folder + '/all_mutations.csv',
                      sep=',',
                      index=False)

    return 0


@click.command()
@click.argument('output_folder')
@click.option('--include_merger', '-m', '/m')
@click.option('--include_filtering', '-f', '/f')
def main(output_folder,
         include_merger='',
         include_filtering=''
         ):
    filter_and_combine(output_folder, include_merger, include_filtering)


if __name__ == "__main__":
    main()
