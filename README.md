# APSCALE graphical user interface
Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

## Introduction
The Apscale Graphical User Interface is a metabarcoding pipeline that handles the most common tasks in metabarcoding
pipelines like paired-end merging, primer trimming, quality filtering, otu clustering and
denoising. It uses a Graphical interface and is configured via a single configuration file.
It automatically uses the available ressources on the machine it runs on while still providing the option
to use less if desired.

For more information on the pipeline running in the background visit [apscale](https://github.com/DominikBuchner/apscale)

## Installation

Apscale can be installed on all common operating systems (Windows, Linux, MacOS).
Apscale requires Python 3.7 or higher and can be easily installed via pip in any command line:

`pip install apscale_gui`

To update apscale_gui run:

`pip install --upgrade apscale_gui`

### Further dependencies - vsearch

Apscale calls vsearch for multiple modules. It should be installed and be in PATH to be executed
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

Apscale also calls cutadapt with some modules. Cutadapt should be downloaded and installed
automatically with the Apscale installation. To check this, type:

`cutadapt --version`

and it should return the version number, for example:

`3.5`

## Tutorial
### Creating a new project
* Create a new folder (e.g. on your desktop) and name it e.g. apscale_project
* Now run the APSCALE GUI with: 
`python -m apscale_gui`
* You will be asked to select an output directory
* Select the folder you just created
* Now create a project using the GUI by typing the name of a new project
* Already existing project folders can be loaded from here
* A blank project folder will be created
<img src="https://github.com/TillMacher/apscale_gui/blob/master/_data/apscale_start.png" width="50%" height="50%">

### The interace
* When loading a project you will be greeted by the APSCALE home window
* From here a multiude of DNA metabarcoding related tools can be started
* Check out the Tutorual for a more detailed overview of the respective tools (COMING SOON) 
<img src="https://github.com/TillMacher/apscale_gui/blob/master/_data/apscale_gui.png" width="50%" height="50%">

### Running apscale: All-in-One Analysis
* The APSCALE pipeline can easily be started via the All-in-One window
* Simply select the steps to run and lean back!
<img src="https://github.com/TillMacher/apscale_gui/blob/master/_data/apscale_run.png" width="50%" height="50%">

## Output

Apscale will output an OTU table and an ESV table, as well as two .fasta files, which can be used for taxnomic assignment. For example, for COI sequences,
BOLDigger (https://github.com/DominikBuchner/BOLDigger) can be used directly with the output of Apscale to assign taxomoy to the OTUs / ESVs using the Barcode of Life Data system (BOLD) database.
