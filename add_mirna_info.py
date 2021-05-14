import click
import os
import pandas as pd
from distinct_occure_helpers import seq_type, from_end, from_start
import warnings

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None


def add_info(output_folder, localization_file, input_file=''):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if input_file:
        all_mutations = pd.read_csv(input_file)
    else:
        all_mutations = pd.read_csv(output_folder + '/all_mutations_algorithms_merged.csv')

    localizations = pd.read_csv(localization_file, sep=',')
    all_mutations = all_mutations.join(localizations.set_index('chrom'), on='chrom', how='left')
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
