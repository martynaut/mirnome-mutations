import click
import gzip
import pandas as pd
import os
import shutil
from prepare_vcf_files_helpers import update_dict_with_file, change_format, change_info


pd.options.mode.chained_assignment = None


def make_unique_files(input_folder, output_folder, copy_input):

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

    if not os.path.exists(output_folder + '/temp'):
        os.makedirs(output_folder + '/temp')

    files = [[x[0] + '/' + y for y in x[2]] for x in os.walk(input_folder)]
    flat_files = [file for sublist in files for file in sublist]
    gz_files = [file for file in flat_files if file.endswith('vcf.gz')]
    
    if not gz_files:
        print("No vcf.gz files in the directory")
        return 1

    dict_with_files = {}
    for gz_file in gz_files:
        dict_with_files = update_dict_with_file(gz_file, dict_with_files)

    df_dict_with_files = pd.DataFrame.from_dict(dict_with_files, orient='index')
    df_dict_with_files.index.name = 'filename'
    df_dict_with_files.to_csv(output_folder + '/files_summary_before_merging.csv', sep=',')
    df_dict_with_files_grouped = df_dict_with_files.reset_index().groupby(['indiv_name',
                                                                           'type_of_file']).agg('nunique')
    df_dict_with_files_grouped.to_csv(output_folder + '/files_summary_count_per_patient_before_merging.csv', sep=',')

    df_not_unique_patients = df_dict_with_files_grouped.loc[df_dict_with_files_grouped['filename'] != 1, :]
    df_not_unique_patients.to_csv(output_folder + '/not_unique_patients.csv', sep=',')
    with open(output_folder+'/do_not_use.txt', 'w+') as do_not_use_file:
        for patient in list(df_not_unique_patients.unstack().index.unique()):
            this_patient = df_not_unique_patients.xs(patient, level=0)
            for file_type in list(this_patient.index.unique()):
                first = True
                with open(output_folder + '/temp/' + patient+'_'+file_type+'.vcf', 'wb') as combined:
                    temp_df = df_dict_with_files.loc[(df_dict_with_files['indiv_name'] == patient) &
                                                     (df_dict_with_files['type_of_file'] == file_type), :]

                    lines_df = pd.DataFrame()
                    columns = []
                    for filename in list(temp_df.index.unique()):
                        print(filename)
                        do_not_use_file.write(filename+'\n')
                        with gzip.open(filename) as f:
                            for line in f.readlines():
                                dline = line.decode('ascii')
                                if dline.startswith('##') and first:
                                    combined.write(line)
                                elif dline.startswith('##'):
                                    pass
                                elif dline.startswith('#') and first:
                                    combined.write(line)
                                    columns = dline.replace('#', '').strip().split('\t')
                                elif dline.startswith('#'):
                                    columns = dline.replace('#', '').strip().split('\t')
                                else:
                                    new_record = \
                                        pd.DataFrame([dline.replace('\n',
                                                                    '').replace(';',
                                                                                ':').replace('"',
                                                                                             '').split('\t')],
                                                     columns=columns)

                                    new_columns_normal = new_record['NORMAL'].str.split(":", expand=True)

                                    normal_columns = list(map(lambda x: x + '_normal',
                                                              new_record['FORMAT'].str.strip().str.split(":").
                                                              values[0]))
                                    try:
                                        new_columns_normal.columns = normal_columns
                                    except ValueError:
                                        normal_columns.remove('SS_normal')
                                        new_columns_normal.columns = normal_columns

                                    new_record = pd.concat([new_record, new_columns_normal], axis=1)

                                    new_columns_tumor = new_record['TUMOR'].str.split(":", expand=True)
                                    tumor_columns = list(map(lambda x: x + '_tumor',
                                                             new_record['FORMAT'].str.strip().str.split(":").
                                                             values[0]))
                                    new_columns_tumor.columns = tumor_columns
                                    new_record = pd.concat([new_record, new_columns_tumor], axis=1)
                                    lines_df = pd.concat([lines_df, new_record])
                        first = False
                    lines_df = lines_df[lines_df['FILTER'].str.contains('PASS')]
                    filter_columns = list(lines_df.columns)
                    filter_columns.remove('CHROM')
                    filter_columns.remove('POS')
                    filter_columns.remove('ID')
                    filter_columns.remove('REF')
                    filter_columns.remove('ALT')
                    filter_columns_names = list(set(list(map(lambda x: x.replace('_tumor', '').replace('_normal', ''),
                                                    filter_columns))))
                    filter_columns_names.remove('FORMAT')
                    filter_columns_names.remove('FILTER')
                    filter_columns_names.remove('TUMOR')
                    filter_columns_names.remove('NORMAL')
                    filter_columns_names.remove('QUAL')
                    filter_columns_names.remove('INFO')
                    if 'SS' in filter_columns_names:
                        filter_columns_names.remove('SS')
                        filter_columns_names.append('SS')
                    format_all = ':'.join(filter_columns_names)

                    lines_df['FORMAT'] = format_all
                    group_dict = {
                            'INFO': lambda x: change_info(x, file_type),
                            'AD_normal': lambda x: change_format(x, file_type, 'AD_normal'),
                            'AD_tumor': lambda x: change_format(x, file_type, 'AD_tumor'),
                            'BQ_normal': lambda x: change_format(x, file_type, 'BQ_normal'),
                            'BQ_tumor': lambda x: change_format(x, file_type, 'BQ_tumor'),
                            'QSS_normal': lambda x: change_format(x, file_type, 'QSS_normal'),
                            'QSS_tumor': lambda x: change_format(x, file_type, 'QSS_tumor'),
                            'BCOUNT_normal': lambda x: change_format(x, file_type, 'BCOUNT_normal'),
                            'BCOUNT_tumor': lambda x: change_format(x, file_type, 'BCOUNT_tumor'),
                            'SSC_normal': lambda x: change_format(x, file_type, 'SSC_normal'),
                            'SSC_tumor': lambda x: change_format(x, file_type, 'SSC_tumor'),
                            'RD_normal': lambda x: change_format(x, file_type, 'RD_normal'),
                            'RD_tumor': lambda x: change_format(x, file_type, 'RD_tumor')
                        }
                    for col in filter_columns:
                        if col in ['AD_normal', 'AD_tumor',
                                   'BQ_normal', 'BQ_tumor',
                                   'QSS_normal', 'QSS_tumor',
                                   'BCOUNT_normal', 'BCOUNT_tumor',
                                   'SSC_normal', 'SSC_tumor',
                                   'RD_normal', 'RD_tumor']:
                            continue
                        else:
                            group_dict[col] = 'first'
                    for col in ['AD_normal', 'AD_tumor',
                                'BQ_normal', 'BQ_tumor',
                                'QSS_normal', 'QSS_tumor',
                                'BCOUNT_normal', 'BCOUNT_tumor',
                                'SSC_normal', 'SSC_tumor',
                                'RD_normal', 'RD_tumor']:
                        if col not in filter_columns:
                            del group_dict[col]
                    lines_df_grouped = lines_df.groupby(['CHROM', 'POS', 'ID',
                                                         'REF',
                                                         'ALT']).agg(group_dict).reset_index()

                    normal_filter_column_names = filter_columns_names.copy()
                    try:
                        normal_filter_column_names.remove('SS')
                    except ValueError:
                        pass
                    lines_df_grouped['NORMAL'] = lines_df_grouped.apply(lambda x: ':'.join(
                        map(lambda y: '' if str(y) == 'nan' else str(y),
                            [x[col_name] for col_name in list(map(
                                lambda y: y + '_normal', normal_filter_column_names))])
                    ), axis=1)
                    lines_df_grouped['TUMOR'] = lines_df_grouped.apply(lambda x: ':'.join(
                        map(lambda y: '' if str(y) == 'nan' else str(y),
                            [x[col_name] for col_name in list(map(lambda y: y + '_tumor', filter_columns_names))])
                    ), axis=1)
                    lines_df_grouped = lines_df_grouped[columns]

                    for index, row in lines_df_grouped.iterrows():
                        combined.write('\t'.join(row.tolist()).encode('utf-8')+b'\n')
                with open(output_folder + '/temp/' + patient+'_'+file_type+'.vcf', 'rb') as combined, \
                        gzip.open(output_folder + '/' + patient+'_'+file_type+'.vcf.gz', 'wb') as f_out:
                    shutil.copyfileobj(combined, f_out)
    with open(output_folder + '/do_not_use.txt', 'r') as file_dont:
        do_not_use = []
        for line in file_dont.readlines():
            do_not_use.append(line.strip())
        gz_files_single = []
        for file in gz_files:
            if file[-6:] == 'vcf.gz' and file not in do_not_use:
                gz_files_single.append(file)

    for file in gz_files_single:
        if copy_input == 1:
            shutil.copyfile(file, output_folder + '/' + file.split('/')[-1])
        else:
            shutil.move(file, output_folder + '/' + file.split('/')[-1])

    click.echo("Cleaning temp folder")
    for filename in os.listdir(output_folder + '/temp'):
        file_path = os.path.join(output_folder + '/temp', filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete {}. Reason: {}'.format(file_path, e))


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.option('--copy_input', '-c')
def main(input_folder,  output_folder, copy_input=''
         ):
    if not copy_input:
        copy_input = 0
    copy_input = int(copy_input)
    make_unique_files(input_folder, output_folder, copy_input)


if __name__ == "__main__":
    main()
