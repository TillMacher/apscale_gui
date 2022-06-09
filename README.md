# APSCALE graphical user interface
Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

[![Downloads](https://pepy.tech/badge/apscale)](https://pepy.tech/project/apscale)  - Apscale

[![Downloads](https://pepy.tech/badge/apscale_gui)](https://pepy.tech/project/apscale_gui)  - Apscale GUI

## Introduction
The APSCALE Graphical User Interface is a metabarcoding pipeline that handles the most common tasks in metabarcoding
pipelines like paired-end merging, primer trimming, quality filtering, otu clustering and
denoising. It uses a Graphical interface and is configured via a single configuration file.
It automatically uses the available ressources on the machine it runs on while still providing the option
to use less if desired.

For more information on the pipeline running in the background visit [APSCALE](https://github.com/DominikBuchner/apscale).

## Installation

APSCALE can be installed on all common operating systems (Windows, Linux, MacOS).
APSCALE requires Python 3.7 or higher and can be easily installed via pip in any command line:

`pip install apscale_gui`

To update apscale_gui run:

`pip install --upgrade apscale_gui`

### Further dependencies - vsearch

APSCALE calls vsearch for multiple modules. It should be installed and be in PATH to be executed
from anywhere on the system.

Check the vsearch Github page for further info:

https://github.com/torognes/vsearch

Support for compressed files with zlib is necessary. For Unix based systems this is shipped with
vsearch, for Windows the zlib.dll can be downloaded via:

[zlib for Windows](https://sourceforge.net/projects/mingw-w64/files/External%20binary%20packages%20%28Win64%20hosted%29/Binaries%20%2864-bit%29/zlib-1.2.5-bin-x64.zip/download)

The dll has to be in the same folder as the vsearch executable. If you need help with adding a folder to PATH in windows
please take a look at the first answer on this stackoverflow issue:

[How to add a folder to PATH Windows](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho)

To check if everything is correctly set up please type this into your command line:

`vsearch --version`

It should return a message similar to this:

```
vsearch v2.19.0_win_x86_64, 31.9GB RAM, 24 cores
https://github.com/torognes/vsearch

Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016)
VSEARCH: a versatile open source tool for metagenomics
PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584

Compiled with support for gzip-compressed files, and the library is loaded.
zlib version 1.2.5, compile flags 65
Compiled with support for bzip2-compressed files, but the library was not found.
```

### Further dependencies - cutadapt

APSCALE also calls cutadapt with some modules. Cutadapt should be downloaded and installed
automatically with the APSCALE installation. To check this, type:

`cutadapt --version`

and it should return the version number, for example:

`3.5`

### Further dependencies - blastn

APSCALE also calls blastn for the local blast modules. It should be installed and be in PATH to be executed
from anywhere on the system.

Check the BLAST Software home page:

https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

you can download it from here:

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

To check this, type:

`blastn -version`

and it should return the version number, for example:

`blastn: 2.12.0+ Package:
  blast 2.12.0, build Jun  4 2021 04:06:33`

## Tutorial
### Creating a new project
Create a new folder (e.g. on your desktop) and name it for example: 'APSCALE_projects'.

Now run the APSCALE GUI with:

`python -m apscale_gui` or simply `apscale_gui`

You will be asked to select an output directory.

Select the folder you just created ('APSCALE_projects').

Now create a new project using the GUI by typing your desired name of the project (e.g. My_new_project'). A new folder in your output directory will be created.

Already existing project folders can be loaded from here in the future.

In this case a new, blank project folder was created.

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_data/apscale_start.png" width="40%" height="40%">

### Data structure

APSCALE is organized in projects with the following structure:

<pre>
/YOUR_PROJECT_PATH/My_new_project/
├───1_raw data
│   └───data
├───2_demultiplexing
│   └───data
├───3_PE_merging
│   └───data
├───4_primer_trimming
│   └───data
├───5_quality_filtering
│   └───data
├───6_dereplication_pooling
│   └───data
│       ├───dereplication
│       └───pooling
├───7_otu_clustering
│   └───data
├───8_denoising
│   └───data
└───9_lulu_filtering
    ├───denoising
    │   └───data
    └───otu_filtering
        └───data

</pre>

### Input data
APSCALE expects *demultiplexed .fastq.gz files* in the 2_demultiplexing/data folder.

If you prefer to have your data all in one place you can copy the raw data into 1_raw_data/data. However, demultiplexing won't be handled by APSCALE directly, but the GUI version has a demultiplexing tool implemented (see https://github.com/DominikBuchner/demultiplexer).


### The interace
When loading a project you will be greeted by the APSCALE home window.

From here a multitude of DNA metabarcoding related tools can be started.

<img width="810" alt="image" src="https://user-images.githubusercontent.com/48299746/172805134-9a9166e7-8b88-472f-a742-f4b2e7d9526e.png">

### Running apscale: All-in-One Analysis
The APSCALE pipeline can easily be started via the All-in-One window.

First, the settings need to be adjusted.
Therefore, one can either adjust the settings from within the GUI and apply them via the green button.
Or one can open the settings file (either from within the GUI or from the project folder) and adjust all settings according to the data set.

Most settings can be left on default. However, following settings need to be adjusted:
- Forward primer sequence (in 5'-3' orientation)
- Reverse primer sequence (in 5'-3' orientation)
- Length of the target fragment (after primer trimming)

To run APSCALE, simply select the steps to perform, click on 'Run analysis' and sit back and enjoy!

<img width="958" alt="image" src="https://user-images.githubusercontent.com/48299746/172804927-b24730f0-530c-4653-9191-a0970ec95cde.png">

## Output

APSCALE will output an OTU table and an ESV table, as well as two .fasta files, which can be used for taxnomic assignment. For example, for COI sequences,
BOLDigger (https://github.com/DominikBuchner/BOLDigger) can be used directly with the output of APSCALE to assign taxomoy to the OTUs / ESVs using the Barcode of Life Data system (BOLD) database. Furthermore, the ESV and OTU tables are compatible with TaxonTableTools (https://github.com/TillMacher/TaxonTableTools), which can be used for DNA metabarcoding specific analyses.

## Local BLAST

The local BLAST tool is really simple to use.

1. Select your sequences (.fasta) and OTU table (.xlsx).
2. Build a new database from a source file (see available dabases below). This only needs to be done once.
3. Select your database to perform the BLAST against.
4. Run the BLAST (blastn is recommended)
5. Filter the BLAST results. The hits per OTU will be filtered as follows:
- By e-value (the e-value is the number of expected hits of similar quality which could be found just by chance):
- The hit(s) with the lowest e-value are kept (the lower the e-value the better).
- By taxonomy:
- Hits with the same taxonomy are dereplicated.
- Hits are adjusted according to thresholds (default: species >=98%, genus >=95%, family >=90%, order >=85%) and dereplicated.
- Hits with still conflicting taxonomy are set back to the most recent common taxonomy
- OTU without matches are collected from the OTU table

<img width="700" alt="image" src="https://user-images.githubusercontent.com/48299746/172805346-f35884da-8727-4efd-a2a2-e5465450939e.png">

## Available databases for local BLAST


### Diat.barcode database
Available from here: https://www6.inrae.fr/carrtel-collection_eng/Barcoding-database/Database-download

Please download the latest .xlsx file!

### Midori2 database
Available from here: http://www.reference-midori.info/download.php#

Please download the latest .fasta file!

For example: Databases/GenBank249/BLAST_AA_sp/fasta/MIDORI_LONGEST_AA_GB249_CO1_BLAST.fasta.zip

Unzip it to recieve the .fasta file!

### Custom NCBI database

Visit the Genbank homepage (https://www.ncbi.nlm.nih.gov/) and search for sequences to add to your database.

Then select
  * Send to:
  * Complete record
  * File
  * GenBank (full)

Then download the .gb file!

<details><summary>Click here to see an example</summary>

<img width="1175" alt="image" src="https://user-images.githubusercontent.com/48299746/172684354-a9645207-6703-4b5d-a23f-45b79efd722c.png">

</details>

Alternatively (for large datasets) one can use the Entrez-Direct tool: https://www.ncbi.nlm.nih.gov/books/NBK179288/

The following command would download all 12S reference sequences for vertebrates:

`esearch -db nuccore -query '12S[All Fields] AND ("Vertebrata"[Organism] OR "Vertebrata"[Organism] OR Vertebrata[All Fields])' | efetch -format gbc > Desktop/vertebrate_sequences.gb`

### My database is missing!

Just let us know if there is need for further databases and we will try to add them.

