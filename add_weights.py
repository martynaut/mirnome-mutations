import pandas as pd
import click
from add_weights_helpers import rev_comp, change_in_motifs, iupac_dict, motifs, \
    motif_check


def add_mutation_weights(output_folder, coordinates_with_seq):

    df_all_mut = pd.read_csv(
        output_folder + '/all_mutations_with_localization.csv'
    )

    df_all_mut['cut_region'] = df_all_mut.apply(
        lambda x: 1 if
        ((x['type'] == 'flanking-5' and x['from end'] == -1)
            or (x['type'] == 'flanking-3' and x['from_start'] == 1)
            or (x['type'] == 'loop' and x['from_start'] == 1)
            or (x['type'] == 'loop' and x['from end'] == -1)
            or (x['type'] == 'pre-seed' and x['from_start'] == 1)
            or (x['type'] == 'post-seed' and x['from end'] == -1)
         ) else 0,
        axis=1
    )

    df_all_mut['duplex'] = df_all_mut['type'].apply(
        lambda x: 1 if x in ['pre-seed', 'seed', 'post-seed'] else 0
    )

    df_all_mut['seed'] = df_all_mut.apply(
        lambda x: 1 if x['type'] == 'seed' and ((x['arm'] == x['balance'])
                                                or x['balance'] == 'both') else 0,
        axis=1
    )

    df_coordinates = pd.read_csv(
        coordinates_with_seq,
        sep='\t',
        header=None
    )

    df_all_mut = df_all_mut.join(df_coordinates.set_index([0, 3]),
                                 on=['chrom', 'pre_name'],
                                 how='left')

    df_all_mut = df_all_mut[(df_all_mut['pos'] >= df_all_mut[1])
                            & (df_all_mut['pos'] <= df_all_mut[2])]

    df_all_mut['whole_seq'] = df_all_mut[4].str.upper()

    df_all_mut['whole_seq'] = df_all_mut.apply(
        lambda x: x['whole_seq'].replace('T', 'U') if x['orientation'] == '+'
        else ''.join([rev_comp[c] for c in list(x['whole_seq'])])[::-1],
        axis=1
    )

    for key in motifs.keys():
        motif = key
        for letter, reg in iupac_dict.items():
            motif = motif.replace(letter, reg)
        df_all_mut['motifs_{}'.format(key)] = df_all_mut.apply(lambda x: motif_check(x, motif, key), axis=1)

    motifs_cols = [col for col in df_all_mut.columns if 'motifs_' in str(col)]

    df_all_mut.head()

    df_all_mut['motifs'] = df_all_mut.apply(
        lambda x: change_in_motifs(x, motifs_cols),
        axis=1
    )

    df_all_mut['weight'] = df_all_mut.apply(
        lambda x: 2 if x['seed'] == 1 else (
            1.5 if (x['duplex'] == 1)
            or (x['motifs'] == 1) or (x['cut_region'] == 1) else 1
        ),
        axis=1
    )

    df_all_mut.to_csv(
        output_folder + '/all_mutations_with_weights.csv'
    )


@click.command()
@click.argument('output_folder')
@click.argument('coordinates_with_seq')
def main(output_folder,
         coordinates_with_seq
         ):
    add_mutation_weights(output_folder, coordinates_with_seq)


if __name__ == "__main__":
    main()
