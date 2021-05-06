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
        if re.match('.*/results_[a-zA-Z0-9]*\.csv$', file):
            results_files.append(file)

    if include_filtering:

        post_analyse(output_folder + '/results_muse.csv',
                     output_folder + '/results_muse_eval.csv')

        post_analyse(output_folder + '/results_mutect2.csv',
                     output_folder + '/results_mutect2_eval.csv')
        post_analyse(output_folder + '/results_varscan2.csv',
                     output_folder + '/results_varscan2_eval.csv')
        post_analyse(output_folder + '/results_somaticsniper.csv',
                     output_folder + '/results_somaticsniper_eval.csv')

    if include_merger and include_filtering:
        pass

    if include_merger:

        df_muse = pd.read_csv(output_folder + '/results_muse_eval.csv')
        df_mutect2 = pd.read_csv(output_folder + '/results_mutect2_eval.csv')
        df_varscan2 = pd.read_csv(output_folder + '/results_varscan2_eval.csv')
        df_ss = pd.read_csv(output_folder + '/results_somaticsniper_eval.csv')

        df_muse = df_muse[df_muse['eval']]
        df_mutect2 = df_mutect2[df_mutect2['eval']]
        df_varscan2 = df_varscan2[df_varscan2['eval']]
        df_ss = df_ss[df_ss['eval']]

        df_muse['alg'] = 'muse'
        df_mutect2['alg'] = 'mutect2'
        df_varscan2['alg'] = 'varscan2'
        df_ss['alg'] = 'somaticsniper'

        df = pd.concat([df_muse, df_mutect2, df_varscan2, df_ss])

        df = df[['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format',
                 'tumor', 'normal', 'indiv_name', 'indiv_id', 'sample_id_tumor_name',
                 'sample_id_tumor_aliQ', 'sample_id_normal_name',
                 'sample_id_normal_aliQ', 'norm_ref_count', 'norm_alt_count',
                 'tumor_ref_count', 'tumor_alt_count', 'BQ_ref_tum', 'BQ_alt_tum',
                 'BQ_ref_norm', 'BQ_alt_norm', 'QSS_ref_tum', 'QSS_alt_tum',
                 'QSS_ref_nor', 'QSS_alt_nor', 'SSC', 'SPV', 'eval', 'alg']]

        df['control:mut/norm'] = df['norm_alt_count'] / df['norm_ref_count']
        df['tumor:mut/norm'] = df['tumor_alt_count'] / df['tumor_ref_count']

        df['ratio'] = df['tumor:mut/norm'] / df['control:mut/norm']

        df.replace(np.inf, 0, inplace=True)

        df.to_csv(output_folder + '/all_mutations_filtered.csv',
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
