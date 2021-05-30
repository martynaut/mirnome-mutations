import pandas as pd
import click
from hgvs_helpers import var_c_p_prep, type_of_mutation, mut_lenght


def hgvs_nomenclature(output_folder):

    table = pd.read_csv(output_folder + '/all_mutations_with_weights.csv')

    # create HGVS from a scratch
    # create 4 new tables for subst, ins, del _1nt, i del_longer, creates new column "to_hgvs_g" with e.g
    # g.423274G>T or g.514_515delAT; at the end joins it to new table 'table'

    table['mutation_type'] = table.apply(lambda x: type_of_mutation(x), axis=1)
    table['delins'] = None
    table['delins'] = table['alt'].apply(lambda x: mut_lenght(x))

    df_subst = table[((table['mutation_type'] == 'subst') & (table['delins'] != 'delins'))].copy()
    if df_subst.shape[0] >= 1:
        df_subst['to_hgvs_g'] = df_subst.apply(lambda x: 'g.' + str(x['pos']) + x['ref'] + '>' + x['alt'], axis=1)

    df_ins = table[(table['mutation_type'] == 'ins')].copy()
    if df_ins.shape[0] >= 1:
        df_ins['to_hgvs_g'] = df_ins.apply(
            lambda x: 'g.' + str(x['pos']) + '_' + str(x['pos'] + 1) + 'ins' + str(x['alt'][1:]), axis=1)

    df_del_1 = table[(table['mutation_type'] == 'del') & (table['ref'].str.len() == 2)].copy()
    if df_del_1.shape[0] >= 1:
        df_del_1['to_hgvs_g'] = df_del_1.apply(lambda x: 'g.' + str(x['pos'] + 1) + 'del' + str(x['ref'][1:]), axis=1)

    df_del_2 = table[(table['mutation_type'] == 'del') & (table['ref'].str.len() > 2)].copy()
    if df_del_2.shape[0] >= 1:
        df_del_2['to_hgvs_g'] = df_del_2.apply(
            lambda x: 'g.' + str(x['pos'] + 1) + '_' + str(x['pos'] + len(x['ref']) - 1) + 'del' + str(x['ref'][1:]),
            axis=1)

    df_delins = table[((table['mutation_type'] == 'subst') & (table['delins'] == 'delins'))].copy()
    if df_delins.shape[0] >= 1:
        df_delins['to_hgvs_g'] = df_delins.apply(
            lambda x: 'g.' + str(x['pos']) + '_' + str(x['pos'] + len(x['ref']) - 1) + 'delins' + x['alt'], axis=1
        )

    table = pd.concat([df_subst, df_ins, df_del_1, df_del_2, df_delins], sort=False)

    table['hgvs_c_p'] = table.apply(lambda x: var_c_p_prep(x), axis=1)

    table.to_csv(output_folder + '/all_mutations_with_hgvs.csv', sep=',', index=False)


@click.command()
@click.argument('output_folder')
def main(output_folder
         ):
    hgvs_nomenclature(output_folder)


if __name__ == "__main__":
    main()
