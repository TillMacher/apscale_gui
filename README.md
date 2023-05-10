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
APSCALE expects *demultiplexed .fastq.gz files* in the 2_demultiplexing/data folder (see above).

APSCALE expects the paired-end reads to end on e.g. *_R1.fastq.gz* and *_R2.fastq.gz*! If APSCALE crashes and you need to rename your files you can simply use the rename tool integrated in APSCALE.

If you prefer to have your data all in one place you can copy the raw data into 1_raw_data/data. However, demultiplexing won't be handled by APSCALE directly, but the GUI version has a demultiplexing tool implemented (see https://github.com/DominikBuchner/demultiplexer).


### The interface
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

### Output

APSCALE will create following output files (that are relevant for downstream analyses):
- Lulu-filtered OTU table (.xlsx and .snappy)
- Lulu-filtered OTU sequences (.fasta)
- Lulu-filtered ESV table (.xlsx and .snappy)
- Lulu-filtered ESV sequences (.fasta)

These files can be used for taxonomic assignment. For example, for COI sequences, BOLDigger (https://github.com/DominikBuchner/BOLDigger) can be used directly with the output of APSCALE to assign taxomoy to the OTUs / ESVs using the Barcode of Life Data system (BOLD) database. Furthermore, the ESV and OTU tables are compatible with TaxonTableTools (https://github.com/TillMacher/TaxonTableTools), which can be used for DNA metabarcoding specific analyses.

<details><summary> Click here to an exemplary APSCALE project </summary>

<pre>
/YOUR_PROJECT_PATH/My_new_project/
├───1_raw data
│   └───data
│       ├───raw_data_R1.fastq.gz
│       └───raw_data_R2.fastq.gz
├───2_demultiplexing
│   └───data
│       ├───SAMPLE_1_a_R1.fastq.gz
│       ├───SAMPLE_1_a_R2.fastq.gz
│       ├───SAMPLE_1_b_R1.fastq.gz
│       ├───SAMPLE_1_b_R2.fastq.gz
│       ├───SAMPLE_2_a_R1.fastq.gz
│       ├───SAMPLE_2_a_R2.fastq.gz
│       ├───SAMPLE_2_b_R1.fastq.gz
│       ├───SAMPLE_2_b_R2.fastq.gz
│       └───...
├───3_PE_merging
│   └───data
│       ├───SAMPLE_1_a_PE.fastq.gz
│       ├───SAMPLE_1_b_PE.fastq.gz
│       ├───SAMPLE_2_a_PE.fastq.gz
│       ├───SAMPLE_2_b_PE.fastq.gz
│       └───...
├───4_primer_trimming
│   └───data
│       ├───SAMPLE_1_a_PE_trimmed.fastq.gz
│       ├───SAMPLE_1_b_PE_trimmed.fastq.gz
│       ├───SAMPLE_2_a_PE_trimmed.fastq.gz
│       ├───SAMPLE_2_b_PE_trimmed.fastq.gz
│       └───...
├───5_quality_filtering
│   └───data
│       ├───SAMPLE_1_a_PE_trimmed_filtered.fastq.gz
│       ├───SAMPLE_1_b_PE_trimmed_filtered.fastq.gz
│       ├───SAMPLE_2_a_PE_trimmed_filtered.fastq.gz
│       ├───SAMPLE_2_b_PE_trimmed_filtered.fastq.gz
│       └───...
├───6_dereplication_pooling
│   └───data
│       ├───dereplication
│       │   ├───SAMPLE_1_a_PE_trimmed_filtered_dereplicated.fastq.gz
│       │   ├───SAMPLE_1_b_PE_trimmed_filtered_dereplicated.fastq.gz
│       │   ├───SAMPLE_2_a_PE_trimmed_filtered_dereplicated.fastq.gz
│       │   ├───SAMPLE_2_b_PE_trimmed_filtered_dereplicated.fastq.gz
│       │   └───...
│       └───pooling
│           ├───pooled_sequences_dereplicated.fasta.gz
│           └───pooled_sequences.fasta.gz
├───7_otu_clustering
│   └───data
│   ├───tutorial_apscale_OTU_table.parquet.snappy
│   ├───tutorial_apscale_OTU_table.xlsx
│   └───tutorial_apscale_OTUs.fasta
├───8_denoising
│   └───data
│   ├───tutorial_apscale_ESV_table.parquet.snappy
│   ├───tutorial_apscale_ESV_table.xlsx
│   └───tutorial_apscale_ESVs.fasta
└───9_lulu_filtering
    ├───denoising
    │   └───data
    │   ├───tutorial_apscale_ESV_table_filtered.parquet.snappy
    │   ├───tutorial_apscale_ESV_table_filtered.xlsx
    │   └───tutorial_apscale_ESVs_filtered.fasta
    └───otu_clustering    
        └───data
        ├───tutorial_apscale_OTU_table_filtered.parquet.snappy
        ├───tutorial_apscale_OTU_table_filtered.xlsx
        └───tutorial_apscale_OTUs_filtered.fasta

</pre>

</details>

## APSCALE modules

### Demultiplexing

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173058747-e1589d38-7a7c-493d-8f87-427b18475378.png)

Raw reads are demultiplexed into individual files, based on indiced and/or tags (see [Bohmann et al., 2022](https://doi.org/10.1111/1755-0998.13512) for an overview)
  
</details>

### Paired-end merging

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173062975-18f7f8d6-b2fb-4ec7-9ab1-72597b620770.png)

Paired-end reads are merged into a single read.
  
</details>

### Primer trimming

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173063025-72ed4de0-d3f4-4a99-b7a0-22f3f10c4d2c.png)
  
Adapter or primer sequences are removed from each read.
  
</details>
  
### Quality & length filtering

<details><summary> Learn more </summary>
 
![image](https://user-images.githubusercontent.com/48299746/173065000-86700cf8-8f79-45f3-8d68-2e3e15b113dd.png)

Reads are filtered according to the expected length of the target fragment. Usually a certain threshold around the expected length is applied (e.g., +-10 of the target fragment length).
  
![image](https://user-images.githubusercontent.com/48299746/173064244-f70439ee-3cd2-4646-b290-d648c231ccb4.png)

Additionally reads are filtered by quality. APSCALE uses the 'maximum expected error' value for quality filtering, which is calculated based on Phred quality score. You can learn more about quality filtering in the [usearch documentation](https://www.drive5.com/usearch/manual/exp_errs.html).
  
</details>

### Dereplication & pooling

<details><summary> Learn more </summary>
  
![image](https://user-images.githubusercontent.com/48299746/173064538-293c9661-f239-4620-8735-9035fe091682.png)

Initially, reads are dereplicated per sample. Only reads with an abundance of at least 4 (default value) are kept.
  
![image](https://user-images.githubusercontent.com/48299746/173064752-082721fb-5afb-423c-af48-6ab9578e684d.png)

Then, reads are pooled into a single file and globally dereplicated. The pooled and dereplicated reads are used for clustering and denoising.
  
</details>

### OTU clustering

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173063360-7e624d90-bc6e-4ca8-aa62-60bc931e83ef.png)
  
Reads are clustered into Operational Taxonomic Units (OTUs), based on a similarity threshold (e.g., 97% similarity).
  
</details>

### Denoising (ESVs)

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173072263-847817e4-86c6-4722-99c2-e15ab5124b84.png)

Reads are denoised into Exact Sequence Variants (ESVs). Here, neighbours with small numbers of differences and small abundance compared to X are predicted to be bad reads of X (see [Edgar 2016](https://doi.org/10.1101/081257) for more details). Denoising is an error removal step.
  
</details>

### Chimera removal (both for OTUs and ESVs)

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173065405-45ad9244-6cbb-4906-9f96-a33626b42cfd.png)

Chimeras are artificial products derived from two biological sequences. They can occur through incomplete extension during PCR. You can learn more about chimeras in the [usearch documentation](https://drive5.com/usearch/manual/chimeras.html). Chimeras are removed from the OTUs and ESVs.
  
</details>

### LULU filtering

<details><summary> Learn more </summary>

The LULU filtering algorithm is used to reduce the number of erroneous OTUs/ESVs to achieve more realistic biodiversity metrics. More details can be found in [Frøslev et al., 2017](https://www.nature.com/articles/s41467-017-01312-x).
  
</details>

### Re-mapping

<details><summary> Learn more </summary>

![image](https://user-images.githubusercontent.com/48299746/173072350-7e661c26-0bc0-405e-bc93-d2c7fdd83baf.png)
 
Lastly, OTUs and ESVs are re-mapped to the sequences of each sample and read tables are created.
  
</details>


## Summary statistics

APSCALE will write all relevant statistics for each module to a project report file. In the ASPCALE-GUI version one can additionally calculate many relevant statistics for the processed dataset. All plots are stored as .pdf and interactive .html charts.

You can check out some examples below:

<details><summary> Boxplot of reads per sample for each module</summary>

![newplot (5)](https://user-images.githubusercontent.com/48299746/173040665-b15f9d71-e10a-4615-a9cc-7a4a03661f16.png)

</details>

<details><summary> Summary of reads per sample for each module (excel table)</summary>

![image](https://user-images.githubusercontent.com/48299746/173041399-57aba00c-dd0f-433e-a0e2-9916e4aeeadc.png)

</details>

<details><summary> OTU summary (all samples)</summary>

![newplot (2)](https://user-images.githubusercontent.com/48299746/173040814-48724104-8448-4877-962f-bfc5a5cf525f.png)

</details>

<details><summary> OTU summary (negative controls)</summary>

![newplot (3)](https://user-images.githubusercontent.com/48299746/173040993-7e3668c4-74de-4828-85c3-437d91841fed.png)

</details>

<details><summary> OTU summary (sample 1 consisting of 4 extraktion replicates with each 2 PCR replicates)</summary>

![newplot (4)](https://user-images.githubusercontent.com/48299746/173040942-36a3ae48-3c3c-4dbb-a7dd-5c1076c0df26.png)

</details>

<details><summary> OTU heatmap (all samples; log of reads) </summary>

![newplot (6)](https://user-images.githubusercontent.com/48299746/173051725-b3756ed0-34b8-4756-84ac-679dda9a1c66.png)

</details>

<details><summary> LULU filtering </summary>

![image](https://user-images.githubusercontent.com/48299746/173078145-e75d5401-343d-4484-a4b7-da52fc1410d3.png)

</details>


## Local BLAST

The local BLAST tool is really simple to use.

<img width="700" alt="image" src="https://user-images.githubusercontent.com/48299746/172805346-f35884da-8727-4efd-a2a2-e5465450939e.png">

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

The following exemplary BLAST results...

| ID  | Hit  | Phylum | Class  | Order | Family  | Genus | Species  | Similarity (%) | E-Value |
| -----  | -----  | ----- | -----  | ----- | -----  | ----- | -----  | ----- | ----- |
| OTU_1  | Hit_1  |  Chordata | Actinopteri  | Esociformes   | Esocidae     | Esox | Esox lucius  | 100 | 3.33e-68 |
| OTU_1  | Hit_2  |  Chordata | Actinopteri  | Esociformes   | Esocidae     | Esox | Esox lucius  | 100 | 3.33e-68 |
| OTU_2  | Hit_1  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Leuciscus | Leuciscus aspius  | 100 | 3.43e-59 |
| OTU_2  | Hit_2  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Squalius | Squalius cephalus  | 100 | 3.43e-59 |
| OTU_3  | Hit_1  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Rutilus | Rutilus rutilus | 95 | 4.77e-35 |
| OTU_3  | Hit_2  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Rutilus | Rutilus rutilus | 95 | 4.77e-35 |
| OTU_4  | Hit_1  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Leuciscus | Leuciscus aspius  | 100 | 1.05e-46 |
| OTU_4  | Hit_2  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Squalius | Squalius cephalus  | 99 | 9.27e-16 |
| OTU_4  | Hit_3  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Barbus | Barbus barbus  | 98 | 1.68e-12 |

... would be filtered into a taxonomy table like this:

| ID  | Phylum | Class  | Order | Family  | Genus | Species  | Similarity (%) |
| -----  | ----- | -----  | ----- | -----  | ----- | -----  | ----- |
| OTU_1  |  Chordata | Actinopteri  | Esociformes | Esocidae  | Esox | Esox lucius  | 100 |
| OTU_2  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  |  |  | 100 |
| OTU_3  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Rutilus | | 95 |
| OTU_4  |  Chordata | Actinopteri  | Cypriniformes | Leuciscidae  | Leuciscus | Leuciscus aspius  | 100 |


## Available databases for local BLAST

### Diat.barcode database
Available from here: https://www6.inrae.fr/carrtel-collection_eng/Barcoding-database/Database-download

Please download the latest .xlsx file!

### Midori2 database
Available from here: http://www.reference-midori.info/download.php#

Please download the latest .fasta file!

The PATH should be as follows: GenBank2xx/BLAST/longest/fasta/*.fasta.zip

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


<details><summary> The easiest way to construct a query is by using the the Genbank browser and copying the search details: </summary>

<img width="1175" alt="image" src="https://user-images.githubusercontent.com/48299746/173529849-ec54c414-dc98-4be9-9a1f-6ff5af922447.png">
  
</details>

The following command will download all 12S reference sequences for vertebrates:

`esearch -db nuccore -query '12S[All Fields] AND ("Vertebrata"[Organism] OR "Vertebrata"[Organism] OR Vertebrata[All Fields]) AND is_nuccore[filter]' | efetch -format gb > Desktop/vertebrate_sequences.gb`


### My database is missing!

Just let us know if there is need for further databases and we will try to add them.

## Citation

Please cite:
* [APSCALE](https://academic.oup.com/bioinformatics/article/38/20/4817/6677653)
* [VSEARCH](https://peerj.com/articles/2584/?report=reader)
* [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200)
* [LULU](https://www.nature.com/articles/s41467-017-01312-x)


