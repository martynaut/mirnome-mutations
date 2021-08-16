import click
from extract_results_for_mirnaome import all_files_processing
from merge_algorithms import filter_and_combine
from distinct_occure import dist_occur
from add_mirna_info import add_info
from mutation_loc_figures import prepare_figure, prepare_figures_per_mirna
from add_weights import add_mutation_weights
from add_hgvs_nomenclature import hgvs_nomenclature


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.argument('localization_file')
@click.argument('coordinates_with_seq')
@click.option('--from_step', '-s', '/s')
@click.option('--end_step', '-es', '/es')
@click.option('--include_merger', '-m', '/m')
@click.option('--include_filtering', '-f', '/f')
@click.option('--csv_file', '-c', '/c')
@click.option('--pass_arg', '-p', '/p')
def main(input_folder,  output_folder, coordinates_file,
         localization_file, coordinates_with_seq,
         from_step='',
         end_step='',
         include_merger='',
         include_filtering='',
         csv_file='',
         pass_arg=''):

    if not from_step:
        from_step = 0
    if not end_step:
        end_step = 6
    if csv_file:
        from_step = 3
    if pass_arg == 'True' or pass_arg == '1':
        pass_arg = True
    else:
        pass_arg = False
    if not include_merger:
        include_merger = 0
    if not include_filtering:
        include_filtering = 0
    from_step = int(from_step)
    end_step = int(end_step)
    include_merger = int(include_merger)
    include_filtering = int(include_filtering)

    if from_step <= 1 <= end_step:
        click.echo("Step 1: Extract results for mirnaome")
        all_files_processing(input_folder, output_folder, coordinates_file, pass_arg)
    else:
        click.echo("Skipping step 1")

    if from_step <= 2 <= end_step:
        click.echo("Step 2 (optional): Merge all algorithms and filter mutations")
        if include_merger:
            click.echo("Merger of algorithms enabled")
        if include_filtering:
            click.echo("Filtering of mutations included")
        result = filter_and_combine(output_folder, include_merger, include_filtering)
        if result == 1:
            return 1
    else:
        click.echo("Skipping step 2")
    if from_step <= 3 <= end_step:
        click.echo("Step 3: Add miRNA information")
        if csv_file:
            add_info(output_folder, localization_file, csv_file)
        else:
            add_info(output_folder, localization_file)
    else:
        click.echo("Skipping step 3")
    if from_step <= 4 <= end_step:
        click.echo("Step 4: Add mutation weights and HGVS nomenclature")
        add_mutation_weights(output_folder, coordinates_with_seq)
        hgvs_nomenclature(output_folder)
    else:
        click.echo("Skipping step 3")
    if from_step <= 5 <= end_step:
        click.echo("Step 5: Make complex, distinct and occure files")
        dist_occur(output_folder)
    else:
        click.echo("Skipping step 4")
    if from_step <= 6 <= end_step:
        click.echo("Step 6: All mutations visualizations")
        prepare_figure(output_folder)
        prepare_figures_per_mirna(output_folder)
    else:
        click.echo("Skipping step 5")

    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
