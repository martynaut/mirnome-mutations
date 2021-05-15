import click
import os
import pandas as pd
import warnings
from distinct_occure_helpers import type_of_mutation, \
    subst_type, seq_type, from_end, from_start

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None


def add_info(output_folder, localization_file, input_file=''):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if input_file:
        all_mutations_raw = pd.read_csv(input_file)
    else:
        all_mutations_raw = pd.read_csv(output_folder + '/all_mutations.csv')

    all_mutations_raw.columns = all_mutations_raw.columns.str.lower()

    all_mutations_raw.dropna(axis=0, inplace=True, how='all')

    columns_to_drop = ['control:mut/norm', 'tumor:mut/norm', 'ratio', 'eval',
                       'qual',
                       'filter', 'info', 'format', 'normal', 'tumor',
                       'indiv_id', 'sample_id_tumor_name',
                       'sample_id_tumor_aliq', 'sample_id_normal_name',
                       'sample_id_normal_aliq'
                       ]

    for column in columns_to_drop:
        if column in all_mutations_raw.columns:
            all_mutations_raw.drop(column, axis=1, inplace=True)

    all_mutations_raw['mutation_type'] = all_mutations_raw.apply(lambda x: type_of_mutation(x), axis=1)

    all_mutations_raw['type_of_subst'] = all_mutations_raw.apply(lambda x: subst_type(x), axis=1)

    all_mutations_raw.fillna(-1, inplace=True)

    localizations = pd.read_csv(localization_file, sep=',')

    all_mutations = all_mutations_raw.join(localizations.set_index('chrom'), on='chrom', how='left')
    all_mutations = all_mutations[(all_mutations['pos'] >= all_mutations['start'])
                                  & (all_mutations['pos'] <= all_mutations['stop'])]
    column = all_mutations.columns[-1]
    double_check = all_mutations.groupby(['chrom', 'pos', 'indiv_name'])[[column]].count()
    double_check.columns = ['no_of_loc']
    all_mutations = all_mutations.join(double_check, on=['chrom', 'pos', 'indiv_name'], how='left')

    all_mutations['seq_type'] = all_mutations['name'].apply(lambda x: seq_type(x))
    if all_mutations.shape[0] > 0:
        all_mutations['from_start'] = all_mutations.apply(lambda x: from_start(x, 'start', 'stop'), axis=1)
        all_mutations['from end'] = all_mutations.apply(lambda x: from_end(x, 'stop', 'start'), axis=1)
    all_mutations.to_csv(output_folder + '/all_mutations_with_localization.csv',
                         sep=',',
                         index=False)
    return 0


@click.command()
@click.argument('output_folder')
@click.argument('localization_file')
@click.option('--input_file', '-f')
def main(output_folder,
         localization_file,
         input_file=''
         ):
    add_info(output_folder, localization_file, input_file)


if __name__ == "__main__":
    main()
