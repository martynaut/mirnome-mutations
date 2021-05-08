import click
import pandas as pd
import os
import re
import numpy as np
from post_extraction_analysis import post_analyse
from distinct_occure_helpers import type_of_mutation, concat_alg, \
    subst_type


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

    columns_to_drop = ['control:mut/norm', 'tumor:mut/norm', 'ratio', 'eval',
                       'qual',
                       'filter', 'info', 'format', 'normal', 'tumor',
                       'indiv_id', 'sample_id_tumor_name',
                       'sample_id_tumor_aliq', 'sample_id_normal_name',
                       'sample_id_normal_aliq'
                       ]

    if include_filtering:
        results_df = results_df[results_df['eval']]
        results_df = results_df[columns]
    else:
        columns.remove('eval')
        columns_to_drop.remove('eval')
        results_df = results_df[columns]

    results_df['control:mut/norm'] = results_df['norm_alt_count'].div(results_df['norm_ref_count'])
    results_df['control:mut/norm'] = results_df['control:mut/norm'].astype(float)
    results_df['tumor:mut/norm'] = results_df['tumor_alt_count'].div(results_df['tumor_ref_count'])
    results_df['tumor:mut/norm'] = results_df['tumor:mut/norm'].astype(float)
    results_df['ratio'] = results_df['tumor:mut/norm'].div(results_df['control:mut/norm'])

    results_df.replace(np.inf, 0, inplace=True)

    results_df.to_csv(output_folder + '/all_mutations.csv',
                      sep=',',
                      index=False)

    all_mutations_raw = results_df.copy()

    all_mutations_raw.columns = all_mutations_raw.columns.str.lower()

    all_mutations_raw.drop(columns_to_drop, axis=1, inplace=True)
    if all_mutations_raw.shape[0] == 0:
        print('no mutations found')
        return 1
    all_mutations_raw['mutation_type'] = all_mutations_raw.apply(lambda x: type_of_mutation(x), axis=1)

    all_mutations_raw.fillna(-1, inplace=True)

    all_mutations = all_mutations_raw.groupby(['chrom', 'pos',
                                               'indiv_name', 'ref', 'alt',
                                               'mutation_type']).agg({'alg': concat_alg,
                                                                      'norm_ref_count': sum,
                                                                      'norm_alt_count': sum,
                                                                      'tumor_ref_count': sum,
                                                                      'tumor_alt_count': sum
                                                                      }).reset_index()
    all_mutations['type_of_subst'] = all_mutations.apply(lambda x: subst_type(x), axis=1)

    all_mutations.to_csv(output_folder + '/all_mutations_algorithms_merged.csv',
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
