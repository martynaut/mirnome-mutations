import click
import pandas as pd
from distinct_occure_helpers import concat_ints, if_complex
import warnings


warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None


def dist_occur(output_folder):

    all_mutations = pd.read_csv(output_folder + '/all_mutations_with_weights.csv')

    df_complex = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name']).agg({
        'pos': ['nunique', 'count']
    }).reset_index()
    df_complex.columns = ['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name',
                          'pos_nunique', 'pos_count']
    df_complex['complex'] = df_complex['pos_count'].apply(lambda x: 0 if x < 2 else 1)

    df_complex.to_csv(output_folder + '/complex_mutations.csv',
                      sep=',',
                      index=False)

    df_by_gene = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type']).agg(
        {'indiv_name': ['nunique',
                        'count'],
         'pos': 'nunique'
         }).reset_index()
    df_by_gene.columns = ['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'indiv_name_nunique',
                          'indiv_name_count', 'pos_nunique']
    df_by_gene['if_complex'] = df_by_gene.apply(lambda x: if_complex(x, df_complex), axis=1)
    df_by_gene.to_csv(output_folder + '/occur.csv',
                      sep=',',
                      index=False)

    if 'norm_ref_count' in all_mutations.columns:

        df_by_mutation = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'pos', 'ref', 'alt',
                                                'mutation_type', 'type_of_subst']).agg(
            {'indiv_name': 'nunique',
             'norm_ref_count': [sum, concat_ints],
             'norm_alt_count': [sum, concat_ints],
             'tumor_ref_count': [sum, concat_ints],
             'tumor_alt_count': [sum, concat_ints],
             'no_of_loc': concat_ints
             }).reset_index()

    else:
        df_by_mutation = all_mutations.groupby(['chrom', 'pre_name', 'id', 'start_pre', 'seq_type', 'pos', 'ref', 'alt',
                                                'mutation_type', 'type_of_subst']).agg(
            {'indiv_name': 'nunique',
             'no_of_loc': concat_ints
             }).reset_index()
    
    df_by_mutation.columns = [' '.join(col).strip() for col in df_by_mutation.columns.values]

    df_by_mutation.to_csv(output_folder + '/distinct_mutations.csv',
                          sep=',',
                          index=False)


@click.command()
@click.argument('output_folder')
def main(output_folder
         ):
    dist_occur(output_folder)


if __name__ == "__main__":
    main()
