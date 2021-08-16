import pandas as pd
import click
from hgvs_helpers import var_c_p_prep, type_of_mutation, mut_lenght, rev_comp


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

    table['hgvs_g'] = table.apply(lambda x: var_c_p_prep(x), axis=1)
    table.drop(['to_hgvs_g'], axis=1, inplace=True)

    table.to_csv(output_folder + '/all_mutations_with_hgvs.csv', sep=',', index=False)

    # N. HGVS

    table['shuffle'] = table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1

    table.loc[table['orientation'] == '-', 'shuffle'] = table['stop_pre_build'].astype(int) - table['pos']. \
        astype(int) + 1

    table.loc[table['shuffle'] <= 0, 'shuffle'] = '1' + (table['shuffle'] - 1).astype(str)

    table.loc[(table['pos'].astype(int) > table['stop_pre_build'].astype(int)) &
              (table['orientation'] == '+'), 'shuffle'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int)).astype(str)
    table.loc[(table['pos'].astype(int) < table['start_pre_build'].astype(int)) &
              (table['orientation'] == '-'), 'shuffle'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int)).astype(str)

    table['hgvs'] = ''

    table.loc[((table['mutation_type'] == 'subst') & (table['orientation'] == '+')), 'hgvs'] = \
        'n.' + table['shuffle'].astype(str) + table['ref'] + '>' + table['alt']

    table.loc[((table['mutation_type'] == 'subst') & (table['orientation'] == '-')), 'hgvs'] = \
        'n.' + table['shuffle'].astype(str) + table['ref'].apply(lambda x: rev_comp(x)) + '>' + \
        table['alt'].apply(lambda x: rev_comp(x))

    table['shuffle_ins'] = (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '_' + \
                           (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1 + 1).astype(str)

    table.loc[table['orientation'] == '-', 'shuffle_ins'] = (table['stop_pre_build'].astype(int) - table['pos'].
                                                             astype(int) + 1).astype(str) + '_' + \
                                                            (table['stop_pre_build'].astype(int) - table['pos'].
                                                             astype(int) + 1 + 1).astype(str)

    table.loc[(table['pos'] <= table['start_pre_build']) & (table['orientation'] == '+'), 'shuffle_ins'] = \
        '1' + (table['pos'].astype(int) - table['start_pre_build'].astype(int)).astype(str) + \
        '_' + '1' + (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str)

    table.loc[(table['pos'] <= table['start_pre_build']) & (table['orientation'] == '-'), 'shuffle_ins'] = \
        '1' + (table['stop_pre_build'].astype(int) - table['pos'].astype(int)).astype(str) + \
        '_' + '1' + (table['stop_pre_build'].astype(int) - table['pos'].astype(int) + 1).astype(str)

    table.loc[(table['pos'].astype(int) > table['stop_pre_build']) &
              (table['orientation'] == '+'), 'shuffle_ins'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int)).astype(str) + '_' + \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int) + 1).astype(str)
    table.loc[(table['pos'].astype(int) < table['start_pre_build']) &
              (table['orientation'] == '-'), 'shuffle_ins'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int)).astype(str) + '_' + \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) + 1).astype(str)

    table.loc[((table['mutation_type'] == 'ins') & (table['orientation'] == '+')), 'hgvs'] = \
        'n.' + table['shuffle_ins'] + 'ins' + \
        table['alt'].apply(lambda x: x[1:])

    table.loc[((table['mutation_type'] == 'ins') & (table['orientation'] == '-')), 'hgvs'] = \
        'n.' + table['shuffle_ins'] + 'ins' + \
        table['alt'].apply(lambda x: rev_comp(x[:-1])[::-1])

    # table['shuffle_del1'] = ''
    # table['shuffle_del2'] = ''
    # table['shuffle_del'] = ''

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '+'), 'shuffle_del1'] = \
        (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1 + 1)

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '+'), 'shuffle_del2'] = \
        (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1 + 1 +
         table['alt'].apply(lambda x: len(x) - 1))

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '-'), 'shuffle_del1'] = \
        (table['stop_pre_build'].astype(int) - table['pos'].
         astype(int) + 1 + 1)

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '-'), 'shuffle_del2'] = \
        (table['stop_pre_build'].astype(int) - table['pos'].
         astype(int) + 1 + 1 +
         table['alt'].apply(lambda x: len(x) - 1))

    table.loc[(table['shuffle_del1'] <= 0) & (table['mutation_type'] == 'del'), 'shuffle_del'] = \
        '1' + (table['shuffle_del1'] - 1).astype(str)

    table.loc[(table['shuffle_del2'] <= 0) & (table['mutation_type'] == 'del'), 'shuffle_del'] = \
        '1' + (table['shuffle_del2'] - 1).astype(str)

    table.loc[((table['pos'].astype(int) + 1) > table['stop_pre_build'].astype(int)) &
              (table['orientation'] == '+') & (table['mutation_type'] == 'del'), 'shuffle_del1'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
              table['pos'].astype(int) - table['stop_pre_build'].astype(int) + 1).astype(str)

    table.loc[((table['pos'].astype(int) + table['alt'].apply(lambda x: len(x) - 1) + 1) > table['stop_pre_build'].
               astype(int)) &
              (table['orientation'] == '+') & (table['mutation_type'] == 'del'), 'shuffle_del2'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int) + 1 + table['alt'].
                apply(lambda x: len(x) - 1)).astype(str)

    table.loc[((table['pos'].astype(int) + 1) < table['start_pre_build'].astype(int)) &
              (table['orientation'] == '-') & (table['mutation_type'] == 'del'), 'shuffle_del1'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) + 1).astype(str)

    table.loc[((table['pos'].astype(int) + table['alt'].apply(lambda x: len(x) - 1) + 1) < table['start_pre_build'].
               astype(int)) &
              (table['orientation'] == '-') & (table['mutation_type'] == 'del'), 'shuffle_del2'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) + 1 + table['alt'].
                apply(lambda x: len(x) - 1)).astype(str)

    table.loc[(table['mutation_type'] == 'del') & (table['alt'].apply(lambda x: len(x)) == 2), 'shuffle_del'] = \
        table['shuffle_del1'].astype(str)

    table.loc[(table['mutation_type'] == 'del') & (table['alt'].apply(lambda x: len(x)) > 2), 'shuffle_del'] = \
        table['shuffle_del1'].astype(str) + '_' + table['shuffle_del2'].astype(str)

    table.loc[(table['mutation_type'] == 'del') &
              (table['orientation'] == '+'), 'hgvs'] = \
        table['shuffle_del'].astype(str) + 'del' + table['alt'].apply(lambda x: x[1:])

    table.loc[(table['mutation_type'] == 'del') &
              (table['orientation'] == '-'), 'hgvs'] = \
        table['shuffle_del'].astype(str) + 'del' + table['alt'].apply(lambda x: rev_comp(x[:-1])[::-1])

    table.drop(['shuffle', 'shuffle_ins', 'shuffle_del',
                'norm_ref_count', 'norm_alt_count', 'tumor_ref_count', 'tumor_alt_count',
                'seq_type'],
               axis=1, inplace=True)
    table.to_csv(output_folder + '/all_mutations_with_n_hgvs.csv', sep=',', index=False)


@click.command()
@click.argument('output_folder')
def main(output_folder
         ):
    hgvs_nomenclature(output_folder)


if __name__ == "__main__":
    main()
