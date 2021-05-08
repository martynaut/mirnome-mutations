import click
import pandas as pd
from distinct_occure_helpers import seq_type, from_end, from_start
import warnings

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None


def add_info(output_folder, localization_file):
    all_mutations = pd.read_csv(output_folder + '/all_mutations_algorithms_merged.csv')

    localizations = pd.read_csv(localization_file, sep=',')
    all_mutations = all_mutations.join(localizations.set_index('chrom'), on='chrom', how='left')
    all_mutations = all_mutations[(all_mutations['pos'] >= all_mutations['start'])
                                  & (all_mutations['pos'] <= all_mutations['stop'])]
    double_check = all_mutations.groupby(['chrom', 'pos', 'indiv_name'])[['ref']].count()
    double_check.columns = ['no_of_loc']
    all_mutations = all_mutations.join(double_check, on=['chrom', 'pos', 'indiv_name'], how='left')

    all_mutations['seq_type'] = all_mutations['name'].apply(lambda x: seq_type(x))

    all_mutations['from_start'] = all_mutations.apply(lambda x: from_start(x, 'start', 'stop'), axis=1)
    all_mutations['from end'] = all_mutations.apply(lambda x: from_end(x, 'stop', 'start'), axis=1)

    all_mutations.to_csv(output_folder + '/all_mutations_with_localization.csv',
                         sep=',',
                         index=False)
    return 0


@click.command()
@click.argument('output_folder')
@click.argument('localization_file')
def main(output_folder,
         localization_file
         ):
    add_info(output_folder, localization_file)


if __name__ == "__main__":
    main()
