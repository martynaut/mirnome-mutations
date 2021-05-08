import click
from extract_results_for_mirnaome import all_files_processing
from merge_algorithms import filter_and_combine
from distinct_occure import dist_occur
from add_mirna_info import add_info


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.argument('localization_file')
@click.option('--from_step', '-s')
@click.option('--end_step', '-es')
@click.option('--include_merger', '-m')
@click.option('--include_filtering', '-f')
def main(input_folder,  output_folder, coordinates_file,
         localization_file,
         from_step='',
         end_step='',
         include_merger='',
         include_filtering=''):

    if not from_step:
        from_step = 0
    if not end_step:
        end_step = 3
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
        all_files_processing(input_folder, output_folder, coordinates_file)
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
        add_info(output_folder, localization_file)
    else:
        click.echo("Skipping step 3")
    if from_step <= 4 <= end_step:
        click.echo("Step 4: Make complex, distinct and occure files")
        dist_occur(output_folder)
    else:
        click.echo("Skipping step 3")

    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
