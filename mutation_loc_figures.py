import pandas as pd
import matplotlib
import click
matplotlib.use('TkAgg')
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
import numpy as np
from matplotlib import lines

plt.rcParams['svg.fonttype'] = 'none'

image_path = 'input_files/primirna_background.tiff'

SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

rc('font', size=MEDIUM_SIZE)          # controls default text sizes
rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def create_plot(data_df, output_name, mutations=0, genes=0, mirna_type='both'):

    df_plot_5, loop_value = prepare_data_5p(data_df)

    df_plot_3 = prepare_data_3p(data_df)

    max_value = max(df_plot_5['pos'].max(), df_plot_3['pos'].max())

    rc("pdf", fonttype=42)
    sns.set_style(style='white')

    palette = {'flanking-5': 'grey',
               'flanking-3': 'grey',
               'pre-seed': '#5481A6',
               'seed': 'darkblue',
               'post-seed': '#5481A6',
               'loop': 'grey',
               'silent-pre': '#5481A6',
               'silent-post': '#5481A6',
               'silent-seed': '#5481A6'
               }

    fig = plt.figure(figsize=(25, 10))
    ax = fig.add_axes([0, 0, 1, 1])

    ax.axis('off')
    im = plt.imread(image_path)
    plt.imshow(im)
    plt.xticks([])
    plt.yticks([])

    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()

    plt.text(x_max * 0.92, y_min * 1, str(mutations), horizontalalignment='center',
             verticalalignment='center', fontdict={'size': '42'})
    plt.text(x_max * 0.92, y_min * 1.09, str(genes), horizontalalignment='center',
             verticalalignment='center', fontdict={'size': '42'})
    line = lines.Line2D([x_max * 0.89, x_max * 0.95], [y_min * 1.04, y_min * 1.04], lw=0.5, color='black')
    ax.add_line(line)
    line.set_clip_on(False)

    plt.text(x_max * 0.24, y_min * -0.09, 'flanking region', horizontalalignment='center', fontdict={'size': '42'})
    line = lines.Line2D([x_max * 0.02, x_max * 0.46], [y_min * -0.07, y_min * -0.07], lw=0.5, color='grey', alpha=0.8)
    ax.add_line(line)
    line.set_clip_on(False)

    plt.text(x_max * 0.65, y_min * -0.09, 'miRNA', horizontalalignment='center', fontdict={'size': '42'})
    line = lines.Line2D([x_max * 0.47, x_max * 0.83], [y_min * -0.07, y_min * -0.07], lw=0.5, color='#5481A6',
                        alpha=0.8)
    ax.add_line(line)
    line.set_clip_on(False)

    plt.text(x_max * 0.91, y_min * -0.09, 'loop', horizontalalignment='center', fontdict={'size': '42'})
    line = lines.Line2D([x_max * 0.84, x_max * 0.98], [y_min * -0.07, y_min * -0.07], lw=0.5, color='grey', alpha=0.8)
    ax.add_line(line)
    line.set_clip_on(False)

    if mirna_type != '3p':
        plt.text(x_max * 0.54, y_min * 0, 'seed', horizontalalignment='center', fontdict={'size': '42'})
        line = lines.Line2D([x_max * 0.49, x_max * 0.59], [y_min * 0.03, y_min * 0.03], lw=0.5, color='darkblue',
                            alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)

    if mirna_type != '5p':
        line = lines.Line2D([x_max * 0.675, x_max * 0.775], [y_min * 1.125, y_min * 1.125], lw=0.5, color='darkblue',
                            alpha=0.8)
        ax.add_line(line)
        line.set_clip_on(False)
        plt.text(x_max * 0.725, y_min * 1.195, 'seed', horizontalalignment='center', fontdict={'size': '42'})

    if loop_value > 0:
        plt.text(
            x_max * 0.92, y_min * 0.6, '+ {} loop\nmutations'.format(loop_value),
            horizontalalignment='center',
            verticalalignment='center',
            fontdict={'size': '30'}
        )

    a = plt.axes([.058, .52, .82, .35])

    hue_order = ['flanking-5', 'pre-seed', 'seed', 'post-seed', 'loop']

    if mirna_type == '3p':
        hue_order = ['flanking-5', 'silent-pre', 'silent-seed', 'silent-post', 'loop']

        df_plot_5['type'] = df_plot_5['type'].apply(
            lambda x: 'silent-seed' if x == 'seed' else (
                'silent-pre' if x == 'pre-seed' else (
                    'silent-post' if x == 'post-seed' else x
                )
            )
        )
    labels = [str(x) for x in range(-25, 0)] + \
             [str(x) for x in range(1, 23)] + ['+' + str(x) for x in range(1, 4)]

    labels = ['L' if x == '+4' else x for x in labels]

    ax = sns.barplot(x="from_start", y="pos", hue="type",
                     data=df_plot_5, dodge=False,
                     hue_order=hue_order,
                     palette=palette,
                     ax=a)

    for loc in ['right', 'top', 'left', 'bottom']:
        ax.spines[loc].set_visible(False)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend(loc='upper right', ncol=2)

    plt.setp(ax.get_legend().get_texts(), fontsize='18')
    plt.setp(ax.get_legend().get_title(), fontsize='22')
    plt.grid(b=True, which='major', axis='y', color='lightgrey',
             linestyle='-', linewidth=0.75, zorder=2, alpha=0.5)

    ax.set_xticklabels(labels)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(30)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(30)

    if max_value > 5:
        plt.yticks(np.arange(0, max_value + 1, np.floor(max_value/3)))
        plot_limit = max_value + 2
    else:
        plt.yticks(np.arange(0, max_value + 1, 1))
        plot_limit = max_value + 1

    ax.set_ylim([0, plot_limit])
    ax.xaxis.tick_bottom()
    for label in ax.get_xticklabels():
        if label.get_text() not in ['-25', '-20', '-15', '-10', '-5', '1', '5', '10', '15', '20',
                                    '+1']:
            label.set_visible(False)

    new_width = 0.35
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_width

        # we change the bar width
        patch.set_width(new_width)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

    ax.tick_params(axis='both', which='major', pad=8, width=0.5)

    plt.setp(ax.patches, linewidth=0)

    ax.get_legend().remove()

    b = plt.axes([.025, .02, .82, .35], facecolor='w')
    hue_order = ['flanking-3', 'pre-seed', 'seed', 'post-seed', 'loop']
    labels = ['+' + str(x) for x in range(1, 26)][::-1] + [str(x) for x in range(1, 23)][::-1] + \
             ['-' + str(x) for x in range(1, 4)]

    if mirna_type == '5p':
        hue_order = ['flanking-3', 'silent-pre', 'silent-seed', 'silent-post', 'loop']

        df_plot_3['type'] = df_plot_3['type'].apply(
            lambda x: 'silent-seed' if x == 'seed' else (
                'silent-pre' if x == 'pre-seed' else(
                    'silent-post' if x == 'post-seed' else x
                )
            )
        )

    ax = sns.barplot(x="from_start", y="pos", hue="type",
                     data=df_plot_3, dodge=False,
                     hue_order=hue_order,
                     palette=palette,
                     ax=b)

    for loc in ['right', 'top', 'left', 'bottom']:
        ax.spines[loc].set_visible(False)

    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.grid(b=True, which='major', axis='y', color='lightgrey',
             linestyle='-', linewidth=0.75, zorder=1, alpha=0.5)

    ax.set_xticklabels(labels)

    if max_value > 5:
        plt.yticks(np.arange(0, max_value + 1, np.floor(max_value/3)))
        plot_limit = max_value + 2
    else:
        plt.yticks(np.arange(0, max_value + 1, 1))
        plot_limit = max_value + 1

    ax.get_legend().remove()
    ax.set_ylim([plot_limit, 0])

    ax.xaxis.tick_top()

    for tick in ax.xaxis.get_major_ticks():
        tick.label2.set_fontsize(30)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(30)

    for label in ax.get_xticklabels():
        if label.get_text() not in ['+25', '+20', '+15', '+10', '+5', '+1', '1', '5', '10', '15', '20',
                                    '-5']:
            label.set_visible(False)

    new_width = 0.35
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_width

        # we change the bar width
        patch.set_width(new_width)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

    plt.setp(ax.patches, linewidth=0)

    plt.savefig(output_name, format='svg', dpi=300, transparent=True, compression="lzw", bbox_inches='tight')
    plt.close()


def prepare_data_5p(df_temp):

    add_loop = df_temp[(df_temp['arm'] == 'loop') & ((df_temp['from_start'] >= 4) & (df_temp['from end'] <= -4))]

    add_loop_value = add_loop.shape[0]

    # add_loop_df = pd.DataFrame([['loop', 'loop', 51, add_loop_value]], columns=['arm', 'type',
    #                                                                             'from_start',
    #                                                                             'pos'])

    dataframe = df_temp[(df_temp['arm'] == '5p') |
                        ((df_temp['arm'] == 'loop') & (df_temp['from_start'] < 4))].groupby(['arm', 'type',
                                                                                             'from_start'],
                                                                                            as_index=False
                                                                                            )[['pos']]. \
        count()[[
            'arm', 'type',
            'from_start',
            'pos'
        ]]

    dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '5p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '5p'), 'from_start'] + 25

    dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '5p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '5p'), 'from_start'] + 26

    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '5p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '5p'), 'from_start'] + 33

    dataframe.loc[(dataframe['type'] == 'loop') & (dataframe['arm'] == 'loop'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'loop') & (dataframe['arm'] == 'loop'), 'from_start'] + 47

    for x in range(1, 26):
        if dataframe[(dataframe['type'] == 'flanking-5') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'flanking-5', x, 0]], columns=['arm', 'type',
                                                                          'from_start',
                                                                          'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(26, 27):
        if dataframe[(dataframe['type'] == 'pre-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'pre-seed', x, 0]], columns=['arm', 'type',
                                                                        'from_start',
                                                                        'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(27, 34):
        if dataframe[(dataframe['type'] == 'seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'seed', x, 0]], columns=['arm', 'type',
                                                                    'from_start',
                                                                    'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(34, 48):
        if dataframe[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['5p', 'post-seed', x, 0]], columns=['arm', 'type',
                                                                         'from_start',
                                                                         'pos'])
            dataframe = pd.concat([dataframe, new_row])
    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] > 47), 'from_start'] = 47

    for x in range(48, 53 - 2):
        if dataframe[(dataframe['type'] == 'loop') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['loop', 'loop', x, 0]], columns=['arm', 'type',
                                                                      'from_start',
                                                                      'pos'])
            dataframe = pd.concat([dataframe, new_row])

    # dataframe = pd.concat([dataframe, add_loop_df])

    dataframe = dataframe.groupby(['arm', 'type',
                                   'from_start'],
                                  as_index=False)[['pos']].sum()[['arm', 'type',
                                                                  'from_start',
                                                                  'pos']]

    return dataframe, add_loop_value


def prepare_data_3p(df_temp):
    dataframe = df_temp[(df_temp['arm'] == '3p')].groupby(['arm', 'type', 'from_start'],
                                                          as_index=False)[['pos']].count()[['arm', 'type',
                                                                                            'from_start',
                                                                                            'pos']]

    try:

        dataframe_loop = df_temp[(df_temp['arm'] == 'loop') &
                                 (df_temp['from end'] > -4)].groupby(['arm', 'type', 'from end'],
                                                                     as_index=False)[['pos']].count()[['arm', 'type',
                                                                                                       'from end',
                                                                                                       'pos']]
        dataframe_loop['from_start'] = dataframe_loop['from end'].apply(
            lambda start: start + 4
        )
        dataframe_loop.drop('from end', inplace=True, axis=1)
    except KeyError:

        dataframe_loop = df_temp[(df_temp['arm'] == 'loop') &
                                 (df_temp['from_end'] > -4)].groupby(['arm', 'type', 'from_end'],
                                                                     as_index=False)[['pos']].count()[['arm', 'type',
                                                                                                       'from_end',
                                                                                                       'pos']]
        dataframe_loop['from_start'] = dataframe_loop['from_end'].apply(
            lambda start: start + 4
        )
        dataframe_loop.drop('from_end', inplace=True, axis=1)

    dataframe = pd.concat([dataframe, dataframe_loop], sort=False)

    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['arm'] == '3p'), 'from_start'] + 13 - 2

    dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'seed') & (dataframe['arm'] == '3p'), 'from_start'] + 6 - 2

    dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'pre-seed') & (dataframe['arm'] == '3p'), 'from_start'] + 5 - 2

    dataframe.loc[(dataframe['type'] == 'flanking-3') & (dataframe['arm'] == '3p'), 'from_start'] = \
        dataframe.loc[(dataframe['type'] == 'flanking-3') & (dataframe['arm'] == '3p'), 'from_start'] + 27 - 2

    for x in range(28 - 2, 53 - 2):
        if dataframe[(dataframe['type'] == 'flanking-3') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'flanking-3', x, 0]], columns=['arm', 'type',
                                                                          'from_start',
                                                                          'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(14 - 2, 28 - 2):
        if dataframe[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'post-seed', x, 0]], columns=['arm', 'type',
                                                                         'from_start',
                                                                         'pos'])
            dataframe = pd.concat([dataframe, new_row])

    dataframe.loc[(dataframe['type'] == 'post-seed') & (dataframe['from_start'] > 27 - 2), 'from_start'] = 27 - 2

    for x in range(7 - 2, 14 - 2):
        if dataframe[(dataframe['type'] == 'seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'seed', x, 0]], columns=['arm', 'type',
                                                                    'from_start',
                                                                    'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(6 - 2, 7 - 2):
        if dataframe[(dataframe['type'] == 'pre-seed') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['3p', 'pre-seed', x, 0]], columns=['arm', 'type',
                                                                        'from_start',
                                                                        'pos'])
            dataframe = pd.concat([dataframe, new_row])

    for x in range(1, 6 - 2):
        if dataframe[(dataframe['type'] == 'loop') & (dataframe['from_start'] == x)].shape[0] == 0:
            new_row = pd.DataFrame([['loop', 'loop', x, 0]], columns=['arm', 'type',
                                                                      'from_start',
                                                                      'pos'])
            dataframe = pd.concat([dataframe, new_row])

    dataframe['from_start'] = dataframe['from_start'].apply(lambda start: start * -1)

    dataframe = dataframe.groupby(['arm', 'type',
                                   'from_start'],
                                  as_index=False)[['pos']].sum()[['arm', 'type',
                                                                  'from_start',
                                                                  'pos']]

    return dataframe


def prepare_figure(output_folder):

    if not os.path.exists(output_folder + '/plots'):
        os.makedirs(output_folder + '/plots')

    df_temp = pd.read_csv(output_folder + '/all_mutations_with_localization.csv')
    # df_temp = df_temp[df_temp['mutation_type'] == 'subst']

    mutations = df_temp.shape[0]
    genes = df_temp['pre_name'].nunique()

    create_plot(df_temp, output_folder + '/plots/plot_miRNA.svg',
                mutations, genes, mirna_type='both')

    # 5' dominant miRNAs

    df_5_dominant = df_temp[df_temp['balance'] == '5p']

    mutations = df_5_dominant.shape[0]
    genes = df_5_dominant['pre_name'].nunique()

    create_plot(df_5_dominant, output_folder + '/plots/plot_5p_balance_miRNA.svg',
                mutations, genes, mirna_type='5p')

    # 3p dominant miRNAs

    df_3_dominant = df_temp[df_temp['balance'] == '3p']

    mutations = df_3_dominant.shape[0]
    genes = df_3_dominant['pre_name'].nunique()

    create_plot(df_3_dominant, output_folder + '/plots/plot_3p_balance_miRNA.svg',
                mutations, genes, mirna_type='3p')

    # balanced miRNAs

    df_no_dominant = df_temp[df_temp['balance'] == 'both']

    mutations = df_no_dominant.shape[0]
    genes = df_no_dominant['pre_name'].nunique()

    create_plot(df_no_dominant, output_folder + '/plots/plot_balanced_miRNA.svg',
                mutations, genes, mirna_type='both')


@click.command()
@click.argument('output_folder')
def main(output_folder):
    prepare_figure(output_folder)


if __name__ == "__main__":
    main()
