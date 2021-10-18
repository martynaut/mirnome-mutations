# miRMut - Pipeline for annotation of mutations in miRNA genes

The miRMut pipeline assigns mutation characteristics including miRNA gene IDs (miRBase and MirGeneDB),
mutation localization within miRNA precursor structure, potential RNA-binding motif disruption,
ascription of mutation according to HGVS nomenclature, as well as miRNA gene characteristics such as miRNA gene confidence
and miRNA arm balance. The pipeline includes creating a tabular and graphical summary of mutations.

Primarily, the scripts were used to analyse mutations in cancer patients using somatic mutation data from TCGA.
The scripts may be reused for other data sources. 

For conditions to reuse of these scripts please refer to `LICENSE` file.

## Using miRMut

Details on how to prepare environment and use the script can be found on [GitHub wiki](https://github.com/martynaut/mirnome-mutations/wiki) pages for this repository.

Example run on sample input files is prepared in `run_test.sh` bash script.

## Authors

Martyna O. Urbanek-Trzeciak, Paulina Galka-Marciniak, Piotr Kozlowski

Institute of Bioorganic Chemistry, Polish Academy of Sciences, Noskowskiego 12/14, 61-704, 
Poznan, Poland

## Citation

Code and execution:

[ref TBD]

Results of the analysis of TCGA data with the usage of the miRMut scripts:

**Pan-cancer analysis of somatic mutations in miRNA genes**

Martyna Olga Urbanek-Trzeciak, Paulina Galka-Marciniak, Paulina Maria Nawrocka, Ewelina Kowal, Sylwia Szwec, Maciej Giefing,  Piotr Kozlowski
doi: https://doi.org/10.1016/j.ebiom.2020.103051

https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(20)30427-8/fulltext

## Contact

For any issues, please create a GitHub Issue.

## Funding

This work was supported by research grants from the Polish National Science Centre [2016/22/A/NZ2/00184 and 2020/39/B/NZ5/01970 (to P.K.), 2015/17/N/NZ3/03629 (to M.O.U-T.) and 2020/39/D/NZ2/03106 (to P.G-M.)]
