# APSCALE Graphical User Interface (GUI)
**Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data**

[![Downloads](https://static.pepy.tech/badge/apscale)](https://pepy.tech/project/apscale) - APSCALE  
[![Downloads](https://static.pepy.tech/badge/apscale-blast)](https://pepy.tech/project/apscale-blast) - APSCALE GUI

## Introduction
The **APSCALE Graphical User Interface (GUI)** is a tool for handling common tasks in DNA metabarcoding pipelines, such as:
- Paired-end merging
- Primer trimming
- Quality filtering
- OTU clustering
- Denoising
- Taxonomic assignment

APSCALE features an intuitive graphical interface and is configured via a single configuration file. It automatically utilizes the available resources on your machine, while allowing you to limit the usage if desired.

For more information on the pipeline running in the background, visit the [APSCALE GitHub repository](https://github.com/DominikBuchner/apscale).

---

## Installation

APSCALE is compatible with Windows, Linux, and macOS. Below are two ways to install APSCALE GUI:

### 1. Via Conda (recommended):
The easiest way to install APSCALE GUI is by using our APSCALE conda environment. [Follow the installation guide here](https://github.com/TillMacher/apscale_installer).

### 2. Via pip:
To install via pip, run:

`pip install apscale_gui`

To update APSCALE GUI to the latest version, use:

`pip install --upgrade apscale_gui`

## Tutorial

Creating a New Project

1.	Create a new folder on your computer (e.g., on your desktop) and name it something like “APSCALE_projects.”
2.	Activate the metabarcoding conda environment by running:

`conda activate metabarcoding`

3.	To start the APSCALE GUI, type:

`apscale_gui`

4.	When prompted, select an output directory. It is recommended to create a new folder, such as “APSCALE_projects,” and select it using the “Browse” button.

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_1.png" alt="Startup window" width="500"/>

5.	After clicking “Continue,” the main menu will appear.

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_2.png" alt="Main window" width="500"/>

6.	If you are starting APSCALE GUI for the first time, click the “Create New Apscale Project” button to start a new project. You will be prompted to name your project, and a new project folder will be created with the following subdirectories:

<pre>
/YOUR_PROJECT_PATH/My_new_project/
├───1_raw_data
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
├───7_denoising
│   └───data
├───8_esv_table
Settings_My_new_project.xlsx
</pre>

6. Alternatively, you can load an existing project by clicking “Load Existing Apscale Project.”

## Raw Data Processing

### Input Data

APSCALE expects demultiplexed .fastq.gz files in the 2_demultiplexing/data folder (see project structure above).

Ensure that the paired-end reads end with _R1.fastq.gz and _R2.fastq.gz.

Optionally, you can place raw data in the 1_raw_data/data folder, but note that APSCALE does not handle demultiplexing directly. You can use tools like the [Demultiplexer tool](https://github.com/DominikBuchner/demultiplexer) for this step.

### Running APSCALE: All-in-One Analysis

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_3.png" alt="Startup window" width="500"/>

1. Adjust Settings: Before running APSCALE, make sure to adjust the necessary settings. You can do this:
* Directly in the GUI, or
* By manually editing the settings file.
2. While most settings can be left on default, you need to configure the following:
* Forward primer sequence (in 5’-3’ orientation)
* Reverse primer sequence (in 5’-3’ orientation)
* Length of the target fragment (after primer trimming)
9. Run the Pipeline
* You can choose to run the entire APSCALE pipeline or specific parts. To run the entire pipeline, select “Run Apscale” and click “Save & Submit.”

### Local BLASTn

The Local BLASTn module performs taxonomic assignments for your data. It uses the APSCALE-BLAST module.

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_4.png" alt="Startup window" width="500"/>

1. Download Database: [Download your database of choice from this server](https://seafile.rlp.net/d/754a25aa76e44b2381b6/).
2. Select Database: After unpacking the database into the “Apscale_databases” folder within your project directory, select the desired database for taxonomic assignment in the GUI.
3. Run the BLASTn search.
4. Save Results: The taxonomy table will automatically be saved to the 8_esv_table folder.

### BOLDigger2

The BOLDigger2 module is used for taxonomic assignments against the BOLD database, based on the [BOLDigger2 pipeline](https://github.com/DominikBuchner/BOLDigger2/tree/main).

<img src="https://github.com/TillMacher/apscale_gui/blob/master/_figures/Apscale_gui_5.png" alt="Startup window" width="500"/>

1. Enter Credentials: You will need to provide your BOLDsystems username and password to allow BOLDigger2 to perform the taxonomic assignment.
2. Save Results: The resulting taxonomy table will be written into the 8_esv_table folder.

### Citation

When using APSCALE, please cite the following:

* [APSCALE](https://academic.oup.com/bioinformatics/article/38/20/4817/6677653?login=false)
* [VSEARCH](https://peerj.com/articles/2584/)
* [Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200)
* [BOLDigger](https://mbmg.pensoft.net/article/53535/)
* [The specific database for local BLASTn search](https://github.com/TillMacher/apscale_blast?tab=readme-ov-file#available-databases).



