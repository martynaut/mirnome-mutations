import click
import gzip
import pandas as pd
import os
import shutil
import numpy as np
from extract_results_for_mirnaome_helpers import retract_counts, retract_info
from prepare_vcf_files_helpers import update_dict_with_file


pd.options.mode.chained_assignment = None


def each_file_processing(filename, coordinates, dict_with_files):

    print(filename)

    dataframe_records = pd.DataFrame()
    if filename.endswith('vcf.gz'):
        f = gzip.open(filename)
    else:
        f = open(filename)

    columns = []
    tumor_name = 'tumor'
    normal_name = 'normal'
    for line in f.readlines():
        try:
            line = line.decode('ascii')
        except AttributeError:
            pass
        if line[:2] == '##':
            if line.startswith('##tumor_sample='):
                tumor_name = line.replace("##tumor_sample=", '').strip().lower()
            elif line.startswith("##normal_sample="):
                normal_name = line.replace("##normal_sample=", '').strip().lower()
        elif line[:1] == '#':
            columns = line.replace('#', '').strip().lower().split('\t')
        else:
            position = line.split('\t')[:5]
            if position[0].startswith('chr'):
                chromosome = position[0]
            else:
                chromosome = 'chr' + position[0]
            if coordinates[(coordinates['chr'] == chromosome) &
                           (coordinates['start_ref'] < int(position[1])) &
                           (coordinates['stop_ref'] > int(position[1]))].shape[0] > 0:
                new_record = pd.DataFrame([line.replace('\n',
                                                        '').replace(';',
                                                                    ':').replace('"',
                                                                                 '').split('\t')],
                                          columns=columns)
                new_record['chrom'] = chromosome
                new_record['indiv_name'] = dict_with_files[filename]['indiv_name']
                new_record['indiv_id'] = dict_with_files[filename]['indiv_id']
                new_record['sample_id_tumor_name'] = dict_with_files[filename]['sample_id_tumor_name']
                new_record['sample_id_tumor_aliQ'] = dict_with_files[filename]['sample_id_tumor_aliQ']
                new_record['sample_id_normal_name'] = dict_with_files[filename]['sample_id_normal_name']
                new_record['sample_id_normal_aliQ'] = dict_with_files[filename]['sample_id_normal_aliQ']
                try:
                    (new_record['norm_ref_count'], new_record['norm_alt_count'],
                        new_record['tumor_ref_count'], new_record['tumor_alt_count'],
                        new_record['BQ_ref_tum'],
                        new_record['BQ_alt_tum'],
                        new_record['BQ_ref_norm'],
                        new_record['BQ_alt_norm'],
                        new_record['QSS_ref_tum'],
                        new_record['QSS_alt_tum'],
                        new_record['QSS_ref_nor'],
                        new_record['QSS_alt_nor'],
                        new_record['SSC']) = retract_counts(
                            new_record[normal_name], new_record[tumor_name], new_record['format'],
                            new_record['ref'], new_record['alt'])
                except KeyError:
                    (new_record['norm_ref_count'], new_record['norm_alt_count'],
                     new_record['tumor_ref_count'], new_record['tumor_alt_count'],
                     new_record['BQ_ref_tum'],
                     new_record['BQ_alt_tum'],
                     new_record['BQ_ref_norm'],
                     new_record['BQ_alt_norm'],
                     new_record['QSS_ref_tum'],
                     new_record['QSS_alt_tum'],
                     new_record['QSS_ref_nor'],
                     new_record['QSS_alt_nor'],
                     new_record['SSC']) = (np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan,
                                           np.nan)

                if np.isnan(float(new_record['SSC'].values[0])):
                    new_record['SSC'], new_record['SPV'] = retract_info(
                        new_record['info']
                    )
                else:
                    new_record['SPV'] = np.nan
                dataframe_records = pd.concat([dataframe_records, new_record])
    f.close()
    return dataframe_records


def all_files_processing(input_folder, output_folder, coordinates_file, pass_arg=False):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    else:
        click.echo("Cleaning output folder")
        for filename in os.listdir(output_folder):
            file_path = os.path.join(output_folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete {}. Reason: {}'.format(file_path, e))

    coordinates = pd.read_table(coordinates_file,
                                names=['chr', 'start', 'stop', 'miRNA_gene_ID'])

    coordinates['start_ref'] = coordinates['start']
    coordinates['stop_ref'] = coordinates['stop']

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
    files = [file for sublist in files for file in sublist]

    vcf_files = []
    for file in files:
        if file.endswith('vcf.gz') or file.endswith('vcf'):
            vcf_files.append(file)

    # files summary
    dict_with_files = {}
    for gz_file in vcf_files:
        dict_with_files = update_dict_with_file(gz_file, dict_with_files)

    df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')

    df_dict_with_files.index.name = 'filename'
    df_dict_with_files.to_csv(output_folder + '/files_summary.csv', sep=',')

    counter = 0
    all_files = df_dict_with_files.shape[0]
    if all_files == 0:
        print('No files found in the input directory')
        return 1
    for file_type in list(df_dict_with_files['type_of_file'].unique()):

        results_df = pd.DataFrame(columns=['chrom',
                                           'pos',
                                           'id',
                                           'ref',
                                           'alt',
                                           'qual',
                                           'filter',
                                           'info',
                                           'format',
                                           'normal',
                                           'tumor',
                                           'indiv_name',
                                           'indiv_id',
                                           'sample_id_tumor_name',
                                           'sample_id_tumor_aliQ',
                                           'sample_id_normal_name',
                                           'sample_id_normal_aliQ',
                                           'norm_ref_count',
                                           'norm_alt_count',
                                           'tumor_ref_count',
                                           'tumor_alt_count',
                                           'BQ_ref_tum',
                                           'BQ_alt_tum',
                                           'BQ_ref_norm',
                                           'BQ_alt_norm',
                                           'QSS_ref_tum',
                                           'QSS_alt_tum',
                                           'QSS_ref_nor',
                                           'QSS_alt_nor',
                                           'SSC',
                                           'SPV'
                                           ])

        for file in list(df_dict_with_files.loc[df_dict_with_files['type_of_file'] == file_type, :].index):
            counter += 1
            print(str(counter) + ' / ' + str(all_files) + '\n')
            dataframe = each_file_processing(file,
                                             coordinates,
                                             dict_with_files
                                             )
            results_df = pd.concat([results_df, dataframe])
        if pass_arg:
            try:
                results_df = results_df[results_df['filter'].str.contains('PASS')]
            except KeyError:
                results_df = results_df
        results_df.to_csv(output_folder + '/results_{}.csv'.format(file_type),
                          sep=',',
                          index=False)


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.option('--pass_arg', '-p', '/p')
def main(input_folder,  output_folder, coordinates_file, pass_arg=False
         ):
    all_files_processing(input_folder, output_folder, coordinates_file, pass_arg)


if __name__ == "__main__":
    main()
