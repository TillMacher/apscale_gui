# APSCALE Graphical User Interface (GUI)
**Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data**

[![Downloads](https://static.pepy.tech/badge/apscale)](https://pepy.tech/project/apscale) - APSCALE  
[![Downloads](https://static.pepy.tech/badge/apscale-blast)](https://pepy.tech/project/apscale-blast) - APSCALE GUI

## Introduction
The **APSCALE Graphical User Interface (GUI)** is a tool for handling common tasks in DNA metabarcoding pipelines, such as:
- Paired-end merging
- Primer trimming
- Quality filtering
- Denoising
- Replicate merging
- Negative control handling
- Taxonomic assignment

APSCALE features an intuitive graphical interface and is configured via a single configuration file. It automatically utilizes the available resources on your machine, while allowing you to limit the usage if desired.

For more information on the pipeline running in the background, visit the [APSCALE GitHub repository](https://github.com/DominikBuchner/apscale).

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_1.png" alt="Startup window" width="500"/>

---

## Installation

APSCALE is compatible with Windows, Linux, and macOS. Below are two ways to install APSCALE GUI:

### 1. Via Conda (recommended):
The easiest way to install APSCALE GUI is by using our APSCALE4 conda environment. [Follow the installation guide here](https://github.com/TillMacher/apscale_installer).

### 2. Via pip:
To install via pip, run:

`pip install apscale_gui`

To update APSCALE GUI to the latest version, use:

`pip install --upgrade apscale_gui`

Note: APSCALE and APSCALE GUI requires python3.12 or higher!

---

## Tutorial

### Setup

1.	Create a new folder on your computer (e.g., on your desktop) and name it something like “APSCALE_projects.”

2.	Activate the metabarcoding conda environment by running:

`conda activate apscale4`

3.	To start the APSCALE GUI, type:

`apscale_gui`

4.	The streamlit graphical-user-interface will start.

5.	First, copy and paste the path to your "APSCALE_projects" folder. You can select to remember the current project folder.

6. *Important note: Use the "Refresh files and folders" button if directories or folders do not yet show up!*

7.	Next, press the button below to create the database and tagging scheme folders in your "APSCALE_projects".

### New projects

8.	To create a first project, simply type a name in the text field below and press the button.

9.	The folder structure of your newly create project looks as follows:

<pre>
/YOURPATH/APSCALE_projects/My_new_project/
├───01_raw_data
│   └───data
├───02_demultiplexing
│   └───data
├───03_PE_merging
│   └───data
├───04_primer_trimming
│   └───data
├───05_quality_filtering
│   └───data
├───06_dereplication
│   └───data
├───07_denoising
│   └───data
├───08_swarm_clustering
│   └───data
├───09_replicate_merging
│   └───data
├───10_nc_removal
│   └───data
├───11_read_table
│   └───data
├───12_analyze
│   └───data
Settings_My_new_project.xlsx
</pre>

9.	Now the project is set-up and the sequencing data can be imported:
* Non-demultiplexed raw data (.fastq.gz) can be copied to: 01_raw_data/data
* Demultiplexed raw data (.fastq.gz) can be copied to: 02_demultiplexing/data
* Ensure that the paired-end reads end with _R1.fastq.gz and _R2.fastq.gz.

## Demultiplexing

10. Demultiplexing can be performed using the [Demultiplexer2](https://github.com/DominikBuchner/demultiplexer2).

11. First, create a primer-set with the respective number of combinations required for the project.

12. Based on the primer-set, create a new tagging-scheme for the demultiplexing. Fill out the Excel sheet according to the tagging information of the sequenced library (the Excel sheet will open automatically).

13. When all information is provided, the demultiplexing can be initiated. After the demultiplexing is finished, the .fastq.gz files in "01_raw_data/data" will be moved to "01_raw_data/processed" and the demultiplexing option in APSCALE GUI will disappear.

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_2.png" alt="Main window" width="500"/>

## Raw Data Processing

14. Now, the raw data processing can be started. For more details on the available options, please refer to the [APSCALE manual](https://github.com/DominikBuchner/apscale/blob/main/manual/apscale_manual.pdf).

15. Most settings can be left at default. The primer sequences (forward and reverse) and the minimum and maximum lenght must be provided individually.

16. APSCALE can be executed in the following modes:
* Basic mode: Runs all modules except for "Replicate merging" and "Negative control removal".
* Complete mode: Runs all modules. Make sure that the the replicate delimiter and negative control prefix are adjusted.
* Individual commands: Each module can be run seperately.

17. The resulting read tables and .fasta files are written to the folder "11_read_table/data".

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_3.png" alt="Main window" width="500"/>

## Taxonomic assignment: BLASTn

18. The sequencing results can now be taxonomically assigned using [APSCALE-blast](https://github.com/TillMacher/apscale_blast).

19. APSCALE-blast is based on blast+ and performs the taxonomic assignments against [pre-compiled local databases](https://seafile.rlp.net/d/474b9682a5cb4193a6ad/). These can be indivudally downloaded or batch-downloaded and extracted using the "Download All Latest Databases" button. The APSCALE databases are stored in the folder "APSCALE_projects/APSCALE_databases".

20. To perform the taxonomic assignment, simply select the .fasta file and the database to begin the blast search.

21. For large datasets and large databases (e.g., COI), it is recommended to select the task "megablast", which decreases the runtimes immensely.

22. The blast results are written to the folder "11_read_table/data".

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_4.png" alt="Main window" width="500"/>

## Taxonomic assignment: BOLDigger3

23. Taxonomic assignment can additionally be performed against the [BOLDsystems database](https://boldsystems.org/), using [BOLDigger3](https://github.com/DominikBuchner/BOLDigger3).

24. To perform the taxonomic assignment, simply select the .fasta file and the database to begin the BOLDigger3 search.

25. The BOLDigger3 results are written to the folder "11_read_table/data".

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_5.png" alt="Main window" width="500"/>

### Citation

When using APSCALE, please cite the following:

* [APSCALE](https://academic.oup.com/bioinformatics/article/38/20/4817/6677653?login=false)
* [VSEARCH](https://peerj.com/articles/2584/)
* [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200)
* [BOLDigger](https://mbmg.pensoft.net/article/53535/)
* [The specific database for local BLASTn search](https://github.com/TillMacher/apscale_blast?tab=readme-ov-file#available-databases).
