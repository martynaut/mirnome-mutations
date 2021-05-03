import click
import gzip
import pandas as pd
import os
import numpy as np
from extract_results_for_mirnaome_helpers import retract_counts, retract_info
from prepare_vcf_files_helpers import update_dict_with_file


pd.options.mode.chained_assignment = None


def each_file_processing(filename, coordinates, dict_with_files):

    print(filename)

    dataframe_records = pd.DataFrame()

    with gzip.open(filename) as f:
        columns = []
        for line in f.readlines():
            line = line.decode('ascii')
            if line[:1] == '##':
                pass
            elif line[:1] == '#':
                columns = line.replace('#', '').strip().lower().split('\t')
            else:
                position = line.split('\t')[:5]
                if coordinates[(coordinates['chr'] == position[0]) &
                               (coordinates['start_ref'] < int(position[1])) &
                               (coordinates['stop_ref'] > int(position[1]))].shape[0] > 0:

                    new_record = pd.DataFrame([line.replace('\n',
                                                            '').replace(';',
                                                                        ':').replace('"',
                                                                                     '').split('\t')],
                                              columns=columns)
                    new_record['indiv_name'] = dict_with_files[filename]['indiv_name']
                    new_record['indiv_id'] = dict_with_files[filename]['indiv_id']
                    new_record['sample_id_tumor_name'] = dict_with_files[filename]['sample_id_tumor_name']
                    new_record['sample_id_tumor_aliQ'] = dict_with_files[filename]['sample_id_tumor_aliQ']
                    new_record['sample_id_normal_name'] = dict_with_files[filename]['sample_id_normal_name']
                    new_record['sample_id_normal_aliQ'] = dict_with_files[filename]['sample_id_normal_aliQ']
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
                            new_record['normal'], new_record['tumor'], new_record['format'],
                            new_record['ref'], new_record['alt'])

                    if np.isnan(float(new_record['SSC'].values[0])):
                        new_record['SSC'], new_record['SPV'] = retract_info(
                            new_record['info']
                        )
                    else:
                        new_record['SPV'] = np.nan
                    dataframe_records = pd.concat([dataframe_records, new_record])

    return dataframe_records


def all_files_processing(input_folder, output_folder, coordinates_file):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    coordinates = pd.read_table(coordinates_file,
                                names=['chr', 'start', 'stop', 'gene'])

    coordinates['start_ref'] = coordinates['start']
    coordinates['stop_ref'] = coordinates['stop']

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
    files = [file for sublist in files for file in sublist]

    files_temp = [[x[0] + '/' + y for y in x[2]] for x in os.walk(output_folder + '/temp')]
    files_temp = [file for sublist in files_temp for file in sublist]

    files.extend(files_temp)

    gz_files = []
    for file in files:
        if file.endswith('vcf.gz') or file.endswith('vcf'):
            gz_files.append(file)
        else:
            pass

    # files summary
    dict_with_files = {}
    for gz_file in gz_files:
        dict_with_files = update_dict_with_file(gz_file, dict_with_files)

    df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')

    df_dict_with_files.index.name = 'filename'
    df_dict_with_files.to_csv(output_folder + '/files_summary.csv', sep=',')

    counter = 0
    all_files = df_dict_with_files.shape[0]
    for file_type in list(df_dict_with_files['type_of_file'].unique()):

        results_df = pd.DataFrame()

        for file in list(df_dict_with_files.loc[df_dict_with_files['type_of_file'] == file_type, :].index):
            counter += 1
            print(str(counter) + ' / ' + str(all_files) + '\n')
            dataframe = each_file_processing(file,
                                             coordinates,
                                             dict_with_files
                                             )
            results_df = pd.concat([results_df, dataframe])
        results_df = results_df[results_df['filter'].str.contains('PASS')]
        results_df.to_csv(output_folder + '/results_{}.csv'.format(file_type),
                          sep=',',
                          index=False)


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
def main(input_folder,  output_folder, coordinates_file,
         ):
    all_files_processing(input_folder, output_folder, coordinates_file)


if __name__ == "__main__":
    main()
