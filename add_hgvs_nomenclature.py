import pandas as pd
import click
from hgvs_helpers import var_c_p_prep, rev_comp, tryconvert


def hgvs_nomenclature(output_folder, weight_filter):

    table = pd.read_csv(output_folder + '/all_mutations_with_weights.csv')

    table = table[table['weight'] >= weight_filter]

    # create HGVS from a scratch
    # create 4 new tables for subst, ins, del _1nt, i del_longer, creates new column "to_hgvs_g" with e.g
    # g.423274G>T or g.514_515delAT; at the end joins it to new table 'table'

    df_subst = table[table['mutation_type'] == 'subst'].copy()
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

    df_delins = table[table['mutation_type'] == 'indel'].copy()
    if df_delins.shape[0] >= 1:
        df_delins['to_hgvs_g'] = df_delins.apply(
            lambda x: 'g.' + str(x['pos']) + '_' + str(x['pos'] + len(x['ref']) - 1) + 'delins' + x['alt'], axis=1
        )

    df_multi = table[table['mutation_type'] == 'multimutation'].copy()
    if df_multi.shape[0] >= 1:
        df_multi['to_hgvs_g'] = ''

    table = pd.concat([df_subst, df_ins, df_del_1, df_del_2, df_delins, df_multi], sort=False)

    table['hgvs_g'] = table.apply(lambda x: var_c_p_prep(x), axis=1)
    table.drop(['to_hgvs_g'], axis=1, inplace=True)

    table.to_csv(output_folder + '/all_mutations_with_hgvs.csv', sep=',', index=False)


def hgvs_n_nomenclature(output_folder):

    table = pd.read_csv(output_folder + '/all_mutations_with_hgvs.csv')
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

    table['shuffle_ins'] = (table['pos'].astype(int) - table['start_pre_build'].astype(int)).astype(str) + '_' + \
                           (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str)

    table.loc[table['orientation'] == '-', 'shuffle_ins'] = (table['stop_pre_build'].astype(int) - table['pos'].
                                                             astype(int)).astype(str) + '_' + \
                                                            (table['stop_pre_build'].astype(int) - table['pos'].
                                                             astype(int) + 1).astype(str)

    table.loc[(table['pos'] <= table['start_pre_build']) & (table['orientation'] == '+'), 'shuffle_ins'] = \
        '1' + (table['pos'].astype(int) - table['start_pre_build'].astype(int)).astype(str) + \
        '_' + '1' + (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str)

    table.loc[(table['pos'] <= table['start_pre_build']) & (table['orientation'] == '-'), 'shuffle_ins'] = \
        '1' + (table['stop_pre_build'].astype(int) - table['pos'].astype(int)).astype(str) + \
        '_' + '1' + (table['stop_pre_build'].astype(int) - table['pos'].astype(int) + 1).astype(str)

    table.loc[(table['pos'].astype(int) > table['stop_pre_build']) &
              (table['orientation'] == '+'), 'shuffle_ins'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int) - 1).astype(str) + '_' + \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int)).astype(str)
    table.loc[(table['pos'].astype(int) < table['start_pre_build']) &
              (table['orientation'] == '-'), 'shuffle_ins'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) - 1).astype(str) + '_' + \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(str) + '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int)).astype(str)

    table.loc[((table['mutation_type'] == 'ins') & (table['orientation'] == '+')), 'hgvs'] = \
        'n.' + table['shuffle_ins'] + 'ins' + \
        table['alt'].apply(lambda x: x[1:])

    table.loc[((table['mutation_type'] == 'ins') & (table['orientation'] == '-')), 'hgvs'] = \
        'n.' + table['shuffle_ins'] + 'ins' + \
        table['alt'].apply(lambda x: rev_comp(x[1:])[::-1])

    # del

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '+'), 'shuffle_del1'] = \
        (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1 + 1).astype(int)

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '+'), 'shuffle_del2'] = \
        (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1 +
         table['ref'].apply(lambda x: len(x) - 1)).astype(int)

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '-'), 'shuffle_del1'] = \
        (table['stop_pre_build'].astype(int) - table['pos'].
         astype(int) + 1 + 1).astype(int)

    table.loc[(table['mutation_type'] == 'del') & (table['orientation'] == '-'), 'shuffle_del2'] = \
        (table['stop_pre_build'].astype(int) - table['pos'].
         astype(int) + 1 +
         table['ref'].apply(lambda x: len(x) - 1)).astype(int)

    table['shuffle_del1'].apply(lambda x: f'{x:.0f}')
    table['shuffle_del2'].apply(lambda x: f'{x:.0f}')

    table['shuffle_del1'].fillna(10000, inplace=True)
    table['shuffle_del2'].fillna(10000, inplace=True)

    table.loc[(table['shuffle_del1'] <= 0) & (table['mutation_type'] == 'del'), 'shuffle_del1'] = \
        '1' + (table['shuffle_del1'] - 1).astype(int).astype(str)

    table.loc[(table['shuffle_del2'] <= 0) & (table['mutation_type'] == 'del'), 'shuffle_del2'] = \
        '1' + (table['shuffle_del2'] - 1).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + 1) > table['stop_pre_build'].astype(int)) &
              (table['orientation'] == '+') & (table['mutation_type'] == 'del'), 'shuffle_del1'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
              table['pos'].astype(int) - table['stop_pre_build'].astype(int) + 1).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + table['ref'].apply(lambda x: len(x) - 1) + 1) > table['stop_pre_build'].
               astype(int)) &
              (table['orientation'] == '+') & (table['mutation_type'] == 'del'), 'shuffle_del2'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int) + 1 + table['ref'].
                apply(lambda x: len(x) - 1)).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + 1) < table['start_pre_build'].astype(int)) &
              (table['orientation'] == '-') & (table['mutation_type'] == 'del'), 'shuffle_del1'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) + 1).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + table['alt'].apply(lambda x: len(x) - 1) + 1) < table['start_pre_build'].
               astype(int)) &
              (table['orientation'] == '-') & (table['mutation_type'] == 'del'), 'shuffle_del2'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) + 1 + table['ref'].
                apply(lambda x: len(x) - 1)).astype(int).astype(str)

    table.loc[(table['mutation_type'] == 'del') & (table['ref'].apply(lambda x: len(x)) == 2), 'shuffle_del'] = \
        table['shuffle_del1'].apply(lambda x: tryconvert(x))

    table.loc[(table['mutation_type'] == 'del') & (table['ref'].apply(lambda x: len(x)) > 2), 'shuffle_del'] = \
        table['shuffle_del1'].apply(lambda x: tryconvert(x)) + '_' \
        + table['shuffle_del2'].apply(lambda x: tryconvert(x))

    table.loc[(table['mutation_type'] == 'del') &
              (table['orientation'] == '+'), 'hgvs'] = \
        'n.' + table['shuffle_del'].astype(str) + 'del' + table['ref'].apply(lambda x: x[1:])

    table.loc[(table['mutation_type'] == 'del') &
              (table['orientation'] == '-'), 'hgvs'] = \
        'n.' + table['shuffle_del'].astype(str) + 'del' + table['ref'].apply(lambda x: rev_comp(x[1:])[::-1])

    # indel

    table.loc[(table['mutation_type'] == 'indel') & (table['orientation'] == '+'), 'shuffle_indel1'] = \
        (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int)

    table.loc[(table['mutation_type'] == 'indel') & (table['orientation'] == '+'), 'shuffle_indel2'] = \
        (table['pos'].astype(int) - table['start_pre_build'].astype(int) + 1 +
         table['ref'].apply(lambda x: len(x) - 1)).astype(int)

    table.loc[(table['mutation_type'] == 'indel') & (table['orientation'] == '-'), 'shuffle_indel1'] = \
        (table['stop_pre_build'].astype(int) - table['pos'].
         astype(int) + 1).astype(int)

    table.loc[(table['mutation_type'] == 'indel') & (table['orientation'] == '-'), 'shuffle_indel2'] = \
        (table['stop_pre_build'].astype(int) - table['pos'].
         astype(int) + 1 +
         table['ref'].apply(lambda x: len(x) - 1)).astype(int)

    table['shuffle_indel1'].apply(lambda x: f'{x:.0f}')
    table['shuffle_indel2'].apply(lambda x: f'{x:.0f}')

    table['shuffle_indel1'].fillna(10000, inplace=True)
    table['shuffle_indel2'].fillna(10000, inplace=True)

    table.loc[(table['shuffle_indel1'] <= 0) & (table['mutation_type'] == 'indel'), 'shuffle_indel1'] = \
        '1' + (table['shuffle_indel1'] - 1).astype(int).astype(str)

    table.loc[(table['shuffle_indel2'] <= 0) & (table['mutation_type'] == 'indel'), 'shuffle_indel2'] = \
        '1' + (table['shuffle_indel2'] - 1).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + 1) > table['stop_pre_build'].astype(int)) &
              (table['orientation'] == '+') & (table['mutation_type'] == 'indel'), 'shuffle_indel1'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
              table['pos'].astype(int) - table['stop_pre_build'].astype(int)).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + table['ref'].apply(lambda x: len(x) - 1) + 1) > table['stop_pre_build'].
               astype(int)) &
              (table['orientation'] == '+') & (table['mutation_type'] == 'indel'), 'shuffle_indel2'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
                table['pos'].astype(int) - table['stop_pre_build'].astype(int) + 1 + table['ref'].
                apply(lambda x: len(x) - 1)).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + 1) < table['start_pre_build'].astype(int)) &
              (table['orientation'] == '-') & (table['mutation_type'] == 'indel'), 'shuffle_indel1'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int)).astype(int).astype(str)

    table.loc[((table['pos'].astype(int) + table['alt'].apply(lambda x: len(x) - 1) + 1) < table['start_pre_build'].
               astype(int)) &
              (table['orientation'] == '-') & (table['mutation_type'] == 'indel'), 'shuffle_indel2'] = \
        (table['stop_pre_build'].astype(int) - table['start_pre_build'].astype(int) + 1).astype(int).astype(str) + \
        '+' + (
                table['start_pre_build'].astype(int) - table['pos'].astype(int) + 1 + table['ref'].
                apply(lambda x: len(x) - 1)).astype(int).astype(str)

    table.loc[(table['mutation_type'] == 'indel'), 'shuffle_indel'] = \
        table['shuffle_indel1'].apply(lambda x: tryconvert(x)) + '_' \
        + table['shuffle_indel2'].apply(lambda x: tryconvert(x))

    table.loc[(table['mutation_type'] == 'indel') &
              (table['orientation'] == '+'), 'hgvs'] = \
        'n.' + table['shuffle_indel'].astype(str) + 'delins' + table['alt']

    table.loc[(table['mutation_type'] == 'indel') &
              (table['orientation'] == '-'), 'hgvs'] = \
        'n.' + table['shuffle_indel'].astype(str) + 'delins' + table['alt'].apply(lambda x: rev_comp(x)[::-1])

    table.drop(['shuffle', 'shuffle_ins', 'shuffle_del', 'shuffle_del1', 'shuffle_del2',
                'shuffle_indel', 'shuffle_indel1', 'shuffle_indel2'],
 #               'seq_type'],
               axis=1, inplace=True)
    try:
        table.drop(['norm_ref_count', 'norm_alt_count', 'tumor_ref_count', 'tumor_alt_count'],
                   axis=1, inplace=True)
    except KeyError:
        pass
    try:
        table.drop(['start_pre'],
                   axis=1, inplace=True)
    except KeyError:
        pass
    table.sort_values(by=['chrom', 'pos'], axis=0, inplace=True)
    table.to_csv(output_folder + '/all_mutations_with_n_hgvs.csv', sep=',', index=False)


@click.command()
@click.argument('output_folder')
def main(output_folder
         ):
    hgvs_nomenclature(output_folder)
    hgvs_n_nomenclature(output_folder)


if __name__ == "__main__":
    main()
