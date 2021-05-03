# miRNome somatic mutations analysis Python scripts

Project containing Python scripts used to analyse somatic mutations in miRNA genes.
Primarily, the scripts were used to analyse mutations in cancer patients using somatic mutation data from TCGA.
The scripts may be reused for other data sources with input data prepared as somatic mutation data 
in TCGA (https://cancergenome.nih.gov/). 

For conditions to reuse of these scripts please refer to `LICENSE` file.

## Pre-run preparation

All needed Python libraries are gathered in `requirements.txt` (Python 3 is needed).

### Preparing ViennaRNA

Download ViennaRNA distribution from https://www.tbi.univie.ac.at/RNA/#download 
to chosen folder
and install using instructions:
```bash
tar -zxvf ViennaRNA-2.4.11.tar.gz
cd ViennaRNA-2.4.11
./configure
make
sudo make install
```

### Before runing mirnome analysis input files (localization and confidence scores) need to be prepared.
#### Examplary input files are included in the repository.

The files are prepared based on downloaded files from miRBase and mirgeneDB from chosen release :
- hairpin.fa (from miRBase)
- hsa.gff3.txt (from miRBase)
- confidence.txt (from miRBase)
- confidence_score.txt (from miRBase)
- aliases.txt (from miRBase)
- mirna_chromosome_build.txt (from miRBase)
- hsa.gff (from mirgeneDB)

ViennaRNA is needed (for installation see above). Important: add **absolute** path to Vienna package 
    
To prepare the files run script `prepare_localization_file.py` similar to (you can also find `run_localization.sh` file in the repository):

    python3 prepare_localization_file.py /path/to/ViennaRNA-2.4.11 \
    /path/to/hairpin.fa \
    /path/to/hsa.gff3.txt \
    /path/to/confidence.txt \
    /path/to/confidence_score.txt \
    /path/to/aliases.txt \
    /path/to/mirna_chromosome_build.txt \
    /path/to/hsa.gff \
    /path/to/calculated/input/files/

### How to use it

Example run is prepared in `run_lusc_test.sh` bash script.

To run analysis run

```
python3 run_mirnaome_analysis.py /path/to/input/LUSC \
/path/to/output/DATA_LUSC \
/path/to/calculated/input/files/coordinates.bed \
/path/to/calculated/input/files/localizations.csv \
 -s 1
```

See additional features running  

```bash
python3 run_mirnaome_analysis.py --help
```

To skip steps of the analysis (if first steps were already completed) use `-s` argument adding step
from which script should start.

## Authors

Martyna O. Urbanek-Trzeciak, Paulina Galka-Marciniak, Piotr Kozlowski

Institute of Bioorganic Chemistry, Polish Academy of Sciences, Noskowskiego 12/14, 61-704, 
Poznan, Poland

## Citation

**Pan-cancer analysis of somatic mutations in miRNA genes**

Martyna Olga Urbanek-Trzeciak, Paulina Galka-Marciniak, Paulina Maria Nawrocka, Ewelina Kowal, Sylwia Szwec, Maciej Giefing,  Piotr Kozlowski
doi: https://doi.org/10.1016/j.ebiom.2020.103051

https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(20)30427-8/fulltext

## Contact

For any issues, please create a GitHub Issue.

## Funding

This work was supported by research grants from the Polish National Science Centre [2016/22/A/NZ2/00184 (to P.K.) and 2015/17/N/NZ3/03629 (to M.O. U.-T.)]
