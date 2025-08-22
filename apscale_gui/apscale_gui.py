import os
import shutil
import time
import streamlit as st
from streamlit_file_browser import st_file_browser
import importlib
from update_checker import update_check
from pathlib import Path
import glob
import pandas as pd
import subprocess
import multiprocessing
import platform
from ete3 import NCBITaxa
from distutils.util import strtobool
import webbrowser
import importlib.util
import glob, datetime, os
import pandas as pd
import numpy as np
from playwright.sync_api import sync_playwright
import zipfile
import importlib.metadata
import subprocess

# demultiplexer2
from demultiplexer2.create_primerset import create_primerset
from demultiplexer2 import find_file_pairs
from demultiplexer2.create_tagging_scheme import main as create_tagging_scheme
from demultiplexer2.demultiplexing import main as demultiplexing

# Apscale
import apscale.a_create_project as a_create_project
import apscale.b_pe_merging as b_pe_merging
import apscale.c_primer_trimming as c_primer_trimming
import apscale.d_quality_filtering as d_quality_filtering
import apscale.e_dereplication as e_dereplication
import apscale.f_denoising as f_denoising
import apscale.g_swarm_clustering as g_swarm_clustering
import apscale.h_replicate_merging as h_replicate_merging
import apscale.i_nc_removal as i_nc_removal
import apscale.j_generate_read_table as j_generate_read_table

# Apscale blast
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
from apscale_blast.__main__ import organism_filter as organism_filter

# BOLDigger3
from boldigger3 import id_engine
from boldigger3 import metadata_download
from boldigger3 import add_metadata
from boldigger3 import select_top_hit

# help texts
n_cores_help = """
    All settings required by Apscale are configured within the project’s settings file. Each processing step has its
    own sheet in the document. The first tab you’ll see when opening the settings file is "0_general_settings". In this
    tab, you can define how many CPU cores Apscale should use for processing. By default, it uses the total number of
    available cores minus two.
  """

compression_level_help = """
    You can also set the gzip compression level here. Apscale compresses all output files to conserve disk space. The
    default compression level is 6, which is suitable for most use cases. If you need to save more space, you can
    increase it to 9—this will produce smaller files but may slow down processing.
  """

b_pe_merging_help = """ 
The first step performed by Apscale is merging paired-end reads using vsearch. The default settings are fairly relaxed
to merge the largest possible portion of reads, as quality filtering is handled in later steps.
                    """

c_primer_trimming_help = """ 
The next step performed by Apscale is primer trimming, which removes the primers used for target amplification since
they do not contain biologically relevant information.
                    """

d_quality_filtering_help = """ 
After primer trimming, Apscale performs quality filtering. This step filters out reads with an expected error higher
than the threshold defined in the settings, as well as sequences whose lengths fall outside the specified target range.
Typically, we use a tolerance of ±10 bases around the target length to allow for some biological variation while
removing artifacts such as primer dimers.
                    """

e_dereplication_help = """ 
Before running the modular workflow, the reads from all samples must be dereplicated. This step does not alter the data
itself but optimizes how it is stored. The output is a FASTA file containing all unique sequences with size annotations (e.g., >seq_1;size=100).
                    """

f_denoising_help = """ 
The denoising module performs sequence denoising using vsearch, processing each sample file individually. Pooling is
intentionally avoided to ensure that the resulting sequences remain independent of the overall dataset size. This
design choice guarantees that previously processed data remains unaffected when new samples are added and the project
is reanalyzed. During denoising, Apscale automatically assigns unique identifiers to each sequence using the SHA3-256
hashing algorithm.
                    """

f_denoising_threshold_help = """ 
Several threshold types are available to control which reads are considered for denoising.  By default, an absolute
threshold is applied (minsize = 4), meaning that only sequences with an abundance of four or more are retained —
effectively removing a substantial amount of low-abundance noise.  Alternatively, a relative threshold can be used to
retain only those sequences that represent a defined percentage of the sample’s total read count (e.g., 0.01%). Since
both absolute and relative thresholds are inherently arbitrary, we introduced a third option in version 4.0.0: power
law–based filtering. Read abundance distributions typically follow a power law, where a few sequences are highly
abundant (true biological signals) and many are rare (a mixture of real low-abundance taxa, sequencing noise, and
PCR artifacts). This filtering method fits a power law model to each sample’s read distribution and sets the threshold
at the point where the observed distribution deviates from the expected power law curve. The underlying assumption is
that this inflection marks a shift in the signal-to-noise ratio, with noise becoming dominant. This approach results
in a data-driven, rather than arbitrary, threshold for denoising.
                    """

g_swarm_clustering_help = """ 
Alternatively, files can be clustered individually using the Swarm algorithm with d=1 and the fastidious option enabled
by default (see https://github.com/torognes/swarm for details). We consider Swarm clustering as an alternative to denoising; when retaining
only the chosen centroid sequences, the results are generally quite similar. Additionally, the output from the
denoising module can be further clustered with Swarm if desired.
                    """

j_generate_read_table_help = """ 
This module generates the read table and performs threshold-based sequence grouping, similar to classical OTU
clustering. Apscale always outputs both sequences (ESVs) and sequence groups (OTUs). The read table is saved in
Parquet format and, if the dataset contains fewer than 1,000,000 distinct sequences, also in Excel format. Additionally,
this module creates a “read data store,” a DuckDB (https://duckdb.org/) database that contains comprehensive
information about sequences, groups, samples, and read counts. The read data store efficiently handles even very
large datasets—potentially billions of sequences—at high speed, without requiring the entire dataset to be loaded into
memory. This makes it especially useful for scaling up analyses.
                    """

def check_dependencies(tools=["cutadapt", "vsearch", "swarm", "blastn"]):
    missing = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        missing_tools = ', '.join(missing)
        if len(missing) == 1:
            st.error(f'WARNING: The following tool is missing: **{missing_tools}**')
        else:
            st.error(f'WARNING: The following tools are missing: **{missing_tools}**')
        st.warning('⚠️ Please install all required tools, either manually or using the "apscale_installer".')

def check_package_update(package_name):
    # Get the currently installed version
    try:
        installed_version = importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        print(f"{package_name} is not installed.")
        return

    # Check for updates
    res = update_check(package_name, installed_version)
    return res

def get_package_versions():
    packages = ["apscale", "apscale_gui", "apscale_blast", "boldigger3", "cutadapt"]

    for pkg in packages:
        try:
            version = importlib.metadata.version(pkg)
            st.write(f"{pkg}: {version}")
        except importlib.metadata.PackageNotFoundError:
            st.write(f"{pkg}: not installed")

    for tool in ['vsearch', 'swarm']:
        try:
            result = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True)
            version = result.stderr.strip().split('\n')[0]
            st.write(f"{tool}: {version}")
        except:
            st.write(f"{tool}: not installed")

    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
        version = result.stdout.strip().split('\n')[0]
        st.write(version)
    except:
        st.write(f"blastn: not installed")

def open_folder(folder_path):
    # Get the current operating system
    current_os = platform.system()

    # Open the folder based on the OS
    try:
        if current_os == "Windows":
            subprocess.Popen(f'explorer "{folder_path}"')
        elif current_os == "Darwin":  # macOS
            subprocess.Popen(['open', folder_path])
        else:  # Linux
            subprocess.Popen(['xdg-open', folder_path])
    except Exception as e:
        print(f"Failed to open folder: {e}")

def open_file(path: Path):
    path = str(path)
    if platform.system() == "Darwin":
        subprocess.run(["open", path])
    elif platform.system() == "Windows":
        os.startfile(path)
    else:
        subprocess.run(["xdg-open", path])

def collect_primerset_information(primerset_path):
    try:
        general_information = pd.read_excel(
            primerset_path, sheet_name="general_information"
        )
    except (FileNotFoundError, ValueError):
        st.error("❌ Please select a valid primerset file.")
        return None, None, None

    forward_primer = general_information["Forward primer (5' - 3')"].replace(np.nan, "").values[0]
    reverse_primer = general_information["Reverse primer (5' - 3')"].replace(np.nan, "").values[0]

    if not forward_primer or not reverse_primer:
        st.error("❌ The primerset is missing primer information.")
        return None, None, None

    forward_tags = pd.read_excel(
        primerset_path, sheet_name="forward_tags", index_col=0
    ).rename(columns={"name": "name_forward_tag", "sequence": "sequence_forward_tag"})
    reverse_tags = pd.read_excel(
        primerset_path, sheet_name="reverse_tags", index_col=0
    ).rename(columns={"name": "name_reverse_tag", "sequence": "sequence_reverse_tag"})

    tag_information = pd.concat([forward_tags, reverse_tags], axis=1)
    return forward_primer, reverse_primer, tag_information

def create_tagging_scheme(tagging_scheme_name, data_dir, primerset_path, path_to_projects):
    # find fastq.gz file pairs
    all_files = sorted(glob.glob(str(Path(data_dir).joinpath("*.fastq.gz"))))
    file_pairs, singles = find_file_pairs.find_file_pairs(all_files)

    if singles:
        st.warning("⚠️ Some files could not be paired:")
        for file in singles:
            st.text(file)
        return None

    if not file_pairs:
        st.error("❌ No file pairs were found in the selected folder.")
        return None

    # read primers + tags
    forward_primer, reverse_primer, tag_information = collect_primerset_information(primerset_path)
    if forward_primer is None:
        return None

    st.subheader("Primer Information")
    st.write(f"**Forward primer:** {forward_primer}")
    st.write(f"**Reverse primer:** {reverse_primer}")

    # let user select combinations
    st.subheader("Select Primer Tag Combinations")
    st.dataframe(tag_information)

    # Create all possible fwd-rev index combinations
    possible_combinations = [
        f"{i}- {j}" for i in tag_information.index for j in tag_information.index
    ]

    selected = st.multiselect(
        "Choose primer tag combinations (format: fwdIndex-revIndex):",
        options=possible_combinations,
    )

    if not selected:
        st.info("ℹ️ Please select at least one combination to continue.")
        return None

    # convert selection into tuples
    combinations_to_use = []
    for comb in selected:
        fwd_idx, rev_idx = comb.split("-")
        fwd_idx, rev_idx = int(fwd_idx.strip()), int(rev_idx.strip())
        fwd_tag_name = tag_information.loc[fwd_idx, "name_forward_tag"]
        rev_tag_name = tag_information.loc[rev_idx, "name_reverse_tag"]
        combinations_to_use.append((fwd_tag_name, rev_tag_name))

    # build tagging scheme dataframe
    file_pair_names = [(pair[0].name, pair[1].name) for pair in file_pairs]
    file_data = [(file_pair + file_pair_name) for file_pair, file_pair_name in zip(file_pairs, file_pair_names)]

    tagging_scheme = pd.DataFrame(
        file_data,
        columns=["forward file path", "reverse file path", "forward file name", "reverse file name"],
    )

    primer_columns = ["{}-{}".format(tag_fwd, tag_rev) for tag_fwd, tag_rev in combinations_to_use]
    primer_columns = pd.DataFrame([], columns=primer_columns)
    tagging_scheme = pd.concat([tagging_scheme, primer_columns], axis=1)

    # save
    savepath = Path(path_to_projects / 'APSCALE_tagging_schemes' / f"{tagging_scheme_name}_tagging_scheme.xlsx")

    if st.button('Save tagging scheme'):
        tagging_scheme.to_excel(savepath, index=False)
        st.success(f"✅ Tagging scheme saved at {savepath}")
        open_file(savepath)
        return savepath

def move_raw_files(project_folder):
    # move raw files from "data" to "processed" folder
    # this indicates that samples were processed
    raw_data_folder_data = project_folder / '01_raw_data' / 'data'
    raw_data_folder_processed = project_folder / '01_raw_data' / 'processed'
    os.makedirs(raw_data_folder_processed, exist_ok=True)

    # Iterate over FASTQ files
    for file in raw_data_folder_data.glob("*.fastq.gz"):
        new_file = raw_data_folder_processed / file.name
        shutil.move(str(file), str(new_file))
        print(f'Moved {file.name}.')

def read_settings_file(settings_xlsx, settings_dfs):
    for sheet_name, df in settings_dfs.items():
        for col in df.columns:
            if col not in st.session_state:  # only initialize if missing
                if len(df[col]) == 1 and pd.notna(df[col].iloc[0]):
                    val = df[col].iloc[0]

                    # Convert boolean-like strings to actual bool
                    if isinstance(val, str) and val.lower() in [True, False]:
                        val = val.lower() == 'true'
                    else:
                        # Everything else -> string
                        val = str(val)

                    st.session_state[col] = val

                else:
                    # Multiple values: replace NaN with '' and convert all to strings
                    st.session_state[col] = df[col].fillna('').astype(str).tolist()

def update_settings_file(settings_xlsx, settings_dfs):
    # Create a dict to hold updated DataFrames
    updated_dfs = {}

    for sheet_name, df in settings_dfs.items():
        # Make a copy
        sub_df = df.copy()

        # Update each column if it exists in session_state
        for col in sub_df.columns:
            if col in st.session_state:
                if st.session_state['run_apscale_mode'] == 'Run apscale (basic mode)' and col in ['perform nc removal', 'perform replicate merging']:
                    sub_df[col] = False
                else:
                    sub_df[col] = st.session_state[col]

        # Store updated DataFrame
        updated_dfs[sheet_name] = sub_df

    # Write all sheets back to Excel
    with pd.ExcelWriter(settings_xlsx, engine='openpyxl') as writer:
        for sheet_name, df in updated_dfs.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    st.success(f"Settings updated and saved to {settings_xlsx}")

def run_apscale(task, project_folder):

    st.info('Starting apscale analysis! Please refer to the terminal for live outputs!')
    print('')

    if task == "Run apscale (basic mode)":
        b_pe_merging.main(project_folder)
        c_primer_trimming.main(project_folder)
        d_quality_filtering.main(project_folder)
        e_dereplication.main(project_folder)
        f_denoising.main(project_folder)
        g_swarm_clustering.main(project_folder)
        j_generate_read_table.main(project_folder)
    if task == "Run apscale (complete mode)":
        b_pe_merging.main(project_folder)
        c_primer_trimming.main(project_folder)
        d_quality_filtering.main(project_folder)
        e_dereplication.main(project_folder)
        f_denoising.main(project_folder)
        g_swarm_clustering.main(project_folder)
        h_replicate_merging.main(project_folder)
        i_nc_removal.main(project_folder)
        j_generate_read_table.main(project_folder)
    if task == "PE-merging":
        b_pe_merging.main(project_folder)
    if task == "Primer-trimming":
        c_primer_trimming.main(project_folder)
    if task == "Quality-filtering":
        d_quality_filtering.main(project_folder)
    if task == "Dereplication":
        e_dereplication.main(project_folder)
    if task == "Denoising":
        f_denoising.main(project_folder)
    if task == "SWARM clustering":
        g_swarm_clustering.main(project_folder)
    if task == "Replicate merging":
        h_replicate_merging.main(project_folder)
    if task == "NC removal":
        i_nc_removal.main(project_folder)
    if task == "Generate read table":
        j_generate_read_table.main(project_folder)

    st.success('Finished apscale analysis!')
    print('')

def run_apscale_blast(project_folder, available_fasta_files, available_databases):

    st.info('Starting blastn! Please refer to the terminal for live outputs!')
    print('')

    if st.session_state['database'] != 'remote':
        database = str(available_databases[st.session_state['database']])
    else:
        database = 'remote'
    query_fasta = str(available_fasta_files[st.session_state['query_fasta']])

    output_folder = project_folder / '11_read_table' / 'data' / f'blastn_{Path(query_fasta).name.replace(".fasta", "")}'
    os.makedirs(output_folder, exist_ok=True)

    organism_mask = []
    if database == 'remote':
        for organism in st.session_state['organism_filter'].replace(' ', '').split(','):
            organism_mask.append(organism_filter(organism))

    if ' ' in str(database):
        st.error(
            "Invalid database path: it contains spaces.\n"
            "blastn cannot handle paths with whitespace.\n"
            "Please choose a path without spaces and try again.")
        return

    a_blastn(
        'blastn',
        query_fasta,
        database,
        output_folder,
        int(st.session_state['n_cores']),
        st.session_state['task'],
        int(st.session_state['subset_size']),
        int(st.session_state['max_target_seqs']),
        st.session_state['masking'],
        st.session_state['disable_headless'],  # headless
        organism_mask,  # organism mask
        st.session_state['include_uncultured']  # include uncultured
    )

    categories = ['t_species', 't_genus', 't_family', 't_order', 't_class']
    thresholds = ','.join([str(st.session_state[i]) for i in categories])
    b_filter(output_folder, database, thresholds, str(st.session_state['n_cores']))

    st.success('Finished blastn!')
    print('')

def run_update_taxids():
    st.info("Updating NCBI taxonomy database...")
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    st.success("Taxonomy database updated successfully.")
    if Path('./taxdump.tar.gz').exists():
        os.remove('./taxdump.tar.gz')
        st.info('Removed taxdmup.tar.gz')

def unzip_and_remove(zip_path, extract_to=None):
    zip_path = Path(zip_path)
    if extract_to is None:
        extract_to = zip_path.with_suffix("")  # folder named after zip
    extract_to = Path(extract_to)
    extract_to.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)

    # remove the zip after extraction
    zip_path.unlink()

    print(f"Extracted to {extract_to} and removed {zip_path.name}")
    return extract_to

def flatten_zip_files(output_dir):
    output_dir = Path(output_dir)
    # Find all .zip files inside subfolders
    zip_files = list(output_dir.rglob("*.zip"))

    for zip_file in zip_files:
        target = output_dir / zip_file.name
        # Move to output_dir, overwrite if exists
        shutil.move(str(zip_file), str(target))
        print(f"Moved {zip_file} -> {target}")

    # Optional: clean up empty subfolders
    for folder in sorted(output_dir.rglob("*"), reverse=True):
        if folder.is_dir() and not any(folder.iterdir()):
            folder.rmdir()
            print(f"Removed empty folder {folder}")

def download_seafile_zip(public_url, output_dir="downloads"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    st.info("Downloading all latest database ZIP files. This may take a few minutes…")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(public_url)

        # Trigger the download inside the expect_download block
        with page.expect_download() as download_info:
            page.click("button:has-text('ZIP')")
        download = download_info.value

        # Save the downloaded file
        download_zip = output_dir / download.suggested_filename
        download.save_as(str(download_zip))
        st.success(f"Saved {download.suggested_filename}")
        browser.close()

    st.info('Extracting databases.')
    if download_zip.exists():
        unzip_and_remove(download_zip)
        flatten_zip_files(output_dir)

        all_zip_files = list(output_dir.glob("*.zip"))
        for zip_file in all_zip_files:
            unzip_and_remove(zip_file)

    st.success(f"Finished the extraction of {len(list(output_dir.glob('*')))} databases!")

def run_boldigger3(available_fasta_files, bold_modes, bold_databases):

    st.info('Starting BOLDigger3! Please refer to the terminal for live outputs!')
    print('')

    # only use the threshold provided by the user replace the rest with defaults
    default_thresholds = [97, 95, 90, 85, 75]
    categories = ['bold_species', 'bold_genus', 'bold_family', 'bold_order', 'bold_class']
    user_thresholds = [int(st.session_state[i]) for i in categories]
    thresholds = []
    for i in range(5):
        try:
            thresholds.append(user_thresholds[i])
        except (IndexError, TypeError):
            thresholds.append(default_thresholds[i])

    # add an virtual treshold of 50 to the thresholds list, so hits with only phylum information can be handled
    thresholds.append(50)

    # Collect input values
    query_fasta = str(available_fasta_files[st.session_state['bold_query_fasta']])
    mode = bold_modes[st.session_state['bold_mode']]
    db = bold_databases[st.session_state['bold_database']]

    # download the current metadata from BOLD
    metadata_download.main()

    # run the id engine
    id_engine.main(
        query_fasta,
        db,
        mode,
    )

    # add additional data via the metadata
    add_metadata.main(query_fasta)

    # select the top hit
    select_top_hit.main(query_fasta, thresholds)

    st.success('Finished BOLDigger3!')
    print('')

def main():

    st.set_page_config(layout="wide")

    check_dependencies()

    # Sidebar inputs & outputs
    with st.sidebar:
        st.subheader("APSCALE projects")

        # read user_data.txt
        script_path = Path(__file__).resolve()
        user_data_txt = script_path.parent / '_user_data' / 'user_data.txt'

        default_value = ""
        if user_data_txt.exists():
            with open(user_data_txt, 'r', encoding='utf-8') as f:
                default_value = f.read().strip()  # strip removes newlines/spaces

        # Get user input
        path_to_projects = Path(st.text_input('Enter Path to APSCALE Projects', value=default_value))
        database_folder = path_to_projects / 'APSCALE_databases'
        tagging_scheme_folder = path_to_projects / 'APSCALE_tagging_schemes'

        if path_to_projects == Path('.'):
            st.write('Please select your APSCALE projects folder.')
        else:
            try:
                # Collect only folders that contain "_apscale" in their name
                if path_to_projects.exists() and path_to_projects.is_dir():
                    apscale_folders = [
                        f for f in path_to_projects.iterdir()
                        if f.is_dir() and "_apscale" in f.name
                    ]

                    if st.button(label='Remember project folder', key='remember_project_folder',use_container_width=True):
                        script_path = Path(__file__).resolve()
                        user_data_txt = script_path.parent / '_user_data' / 'user_data.txt'
                        user_data_txt.parent.mkdir(parents=True, exist_ok=True)
                        with open(user_data_txt, 'w', encoding='utf-8') as f:
                            f.write(str(path_to_projects))
                        st.success(f"Saved project folder!")

                    if apscale_folders:
                        project_folder = st.selectbox(
                            "Select an APSCALE project folder:",
                            sorted(apscale_folders, key=lambda p: p.name.lower()),  # sort by name only
                            format_func=lambda p: p.name
                        )
                        project_name = project_folder.name.replace('_apscale', '')
                        if st.button('Open project folder', use_container_width=True):
                            open_folder(project_folder)

                    else:
                        st.warning("No APSCALE project folders found here.")

                else:
                    st.error("The given path does not exist or is not a directory.")

                if not database_folder.exists():
                    st.error("Please first initialise the database and tagging scheme folder.")
                    if st.button('Create database and tagging scheme folders', use_container_width=True):
                        os.makedirs(database_folder, exist_ok=True)
                        os.makedirs(tagging_scheme_folder, exist_ok=True)
                        st.rerun()

                st.markdown("---")
                st.subheader('Create new project')
                st.text_input(label='Enter name of new project', key='new_project_name')
                new_project_path = path_to_projects / Path(st.session_state['new_project_name'])
                if new_project_path != path_to_projects:
                    if st.button('Create new project', use_container_width=True):
                        a_create_project.create_project(new_project_path)
                        st.success('Create new project folder!')

                st.markdown("---")
                st.subheader('APSCALE databases')
                n_databases = len(glob.glob(str(database_folder / '*')))
                st.info(f'{n_databases} databases available')

                if st.button('Open Database Hub', use_container_width=True):
                    webbrowser.open('https://seafile.rlp.net/d/474b9682a5cb4193a6ad/')
                if st.button('Download All Latest Databases', use_container_width=True):
                    public_url = 'https://seafile.rlp.net/d/474b9682a5cb4193a6ad/?p=%2FLatest&mode=list'
                    download_seafile_zip(public_url, database_folder)
                if st.button('Open APSCALE database folder', use_container_width=True):
                    open_folder(database_folder)

            except Exception as e:
                st.error(f"An error occurred: {e}")

        st.markdown("---")
        st.subheader("Refresh")

        if st.button("🔄 Refresh files and folders", use_container_width=True):
            st.success("Files and folders refreshed.")


    ####################################################################################################################
    if not path_to_projects.exists() and not path_to_projects.is_dir():
        st.info('Please select a project to continue!')
    elif 'project_folder' not in locals():
        st.info('Project folder is not selected yet!')
    else:
        ################################################################################################################
        # Raw data check
        st.header("APSCALE Raw Data Processing")

        raw_data = glob.glob(str(project_folder / '01_raw_data' / 'data' / '*.fastq.gz'))
        n_raw_data = len(raw_data)
        if n_raw_data == 0:
            st.info('No raw data files to demultiplex!')
        else:
            #############################################################################################################
            # Demultiplexer2
            st.subheader("Demultiplexer2")
            st.info(f"Found {n_raw_data} files for demultiplexing!")

            # Locate Demultiplexer2 package and primersets
            spec = importlib.util.find_spec("demultiplexer2")
            demultiplexer2_path = Path(spec.origin).parent
            primersets_path = demultiplexer2_path / "data" / "primersets"
            primersets = [f for f in primersets_path.iterdir() if f.is_file() and "_primerset" in f.name]

            with st.expander("📑 Manage primer sets", expanded=True):
                # Select existing primer set
                active_primerset = st.selectbox(
                    "Select from existing primer sets:",
                    sorted(primersets, key=lambda p: p.name.lower()),
                    format_func=lambda p: p.name,
                    key="active_primerset"
                )
                if st.button('Open primer set folder'):
                    open_folder(primersets_path)

                # Input fields to create new primer set
                st.markdown("**➕ Create new primer set**")
                col1, col2 = st.columns(2)
                with col1:
                    st.text_input("Primer set name", key="primerset_name")
                with col2:
                    st.text_input("Number of primers", key="n_primers")

                # Action buttons
                col3, col4 = st.columns(2)
                with col3:
                    if st.button("Create new primer set", use_container_width=True):
                        create_primerset(
                            st.session_state['primerset_name'],
                            int(st.session_state['n_primers'])
                        )
                        st.success(f"New primer set '{st.session_state['primerset_name']}' created")
                with col4:
                    if st.button("Open selected primer set", use_container_width=True):
                        open_file(st.session_state["active_primerset"])

            if active_primerset:
                with st.expander("🏷️ Manage tagging schemes", expanded=True):
                    # Select existing tagging scheme
                    tagging_schemes = [
                        f for f in tagging_scheme_folder.iterdir()
                        if f.is_file() and "_tagging_scheme" in f.name and not f.name.startswith(".")
                    ]
                    active_tagging_scheme = st.selectbox(
                        "Select from existing tagging schemes:",
                        sorted(tagging_schemes, key=lambda p: p.name.lower()),
                        format_func=lambda p: p.name,
                        key="active_tagging_scheme"
                    )

                    if st.button("Open selected tagging scheme", use_container_width=True):
                        open_file(st.session_state["active_tagging_scheme"])
                    create_new_tagging_scheme = st.toggle('Create new tagging scheme')
                    if create_new_tagging_scheme == True:
                        st.markdown("**➕ Create new tagging scheme**")
                        st.text_input(label='Enter name of new tagging scheme name:', key='tagging_scheme_name')
                        if st.session_state['tagging_scheme_name'] != '':
                            create_tagging_scheme(
                                st.session_state['tagging_scheme_name'],
                                project_folder / '01_raw_data' / 'data',
                                st.session_state['active_primerset'],
                                path_to_projects
                            )

            if active_primerset and active_tagging_scheme:
                with st.expander("🔀 Demultiplexing", expanded=True):
                    st.success('All required demultiplexing information is provided!')
                    st.write(f"Primer set: **{st.session_state['active_primerset'].name}**")
                    st.write(f"Tagging scheme: **{st.session_state['active_tagging_scheme'].name}**")
                    output_dir = project_folder / '02_demultiplexing' / 'data'
                    st.write(f"Output directory: **{output_dir}**")

                    move_files = st.checkbox("Move raw files to 'processed' after demultiplexing?", value=True)
                    if st.button('Begin demultiplexing', use_container_width=True):
                        st.info('Starting demultiplexing! Please refer to the terminal for live outputs!')
                        print('')
                        demultiplexing(st.session_state['active_primerset'], st.session_state['active_tagging_scheme'],output_dir)
                        if move_files:
                            move_raw_files(project_folder)
                            st.success('Finished demultiplexing and moved raw files!')
                        else:
                            st.success('Finished demultiplexing (raw files left in place).')
                        print('')
                        st.rerun()

        ################################################################################################################
        demultiplexing_data = glob.glob(str(project_folder / '02_demultiplexing' / 'data' / '*.fastq.gz'))
        n_demultiplexing_data = len(demultiplexing_data)
        if n_demultiplexing_data == 0:
            st.markdown("---")
            st.error(f'No input data files were found! Please copy your .fastq.gz files to "{project_folder}/02_demultiplexing/data"')

        else:
            st.info(f'Found {n_demultiplexing_data} files to proceed raw data analysis!')

            ################################################################################################################
            # Input
            settings_xlsx = project_folder / f'Settings_{project_name}.xlsx'

            if not settings_xlsx.exists():
                st.warning('⚠️ Could not find the settings file!')
            else:

                # Read Settings file
                settings_dfs = pd.read_excel(settings_xlsx, sheet_name=None)
                # Replace NaN with empty string in all sheets
                for sheet_name, df in settings_dfs.items():
                    settings_dfs[sheet_name] = df.fillna('')
                read_settings_file(settings_xlsx, settings_dfs)

                ############################################################################################################
                st.subheader('APSCALE4')
                # General settings
                with st.expander("⚙️ General settings", expanded=False):
                    col1, col2 = st.columns(2)
                with col1:
                    st.text_input(label='cores to use', key='cores to use', help=n_cores_help)
                with col2:
                    st.text_input(label='compression level', key='compression level', help=compression_level_help)
                # PE Merging
                with st.expander("🔗 PE Merging", expanded=False):
                    col1, col2, col3 = st.columns(3)
                with col1:
                    st.text_input(label='Maximum difference (%)', key='maxdiffpct', help=b_pe_merging_help)
                with col2:
                    st.text_input(label='Maximum differences', key='maxdiffs')
                with col3:
                    st.text_input(label='Minimum overlap length', key='minovlen')
                # Primer trimming
                with st.expander("🧬 Primer trimming", expanded=True):
                    col1, col2, col3 = st.columns(3)
                with col1:
                    st.text_input(label="P5 Primer (5' - 3')", key="P5 Primer (5' - 3')", help=c_primer_trimming_help)
                with col2:
                    st.text_input(label="P7 Primer (5' - 3')", key="P7 Primer (5' - 3')")
                with col3:
                    st.selectbox(label='Anchoring', options=[True, False], key='anchoring')
                # Quality filtering
                with st.expander("📉 Quality filtering", expanded=True):
                    col1, col2, col3 = st.columns(3)
                with col1:
                    st.text_input(label="Maximum Expected Error", key='maxEE', help=d_quality_filtering_help)
                with col2:
                    st.text_input(label="Minimum length", key='min length')
                with col3:
                    st.text_input(label="Maximum length", key='max length')
                # Dereplication
                with st.expander("🔁 Dereplication", expanded=False):
                    st.text_input(label='Minimum sequence abundance', key='minimum sequence abundance', help=e_dereplication_help)
                # Denoising
                with st.expander("✨ Denoising", expanded=False):
                    col1, col2, col3 = st.columns(3)
                with col1:
                    st.selectbox(label='Perform denoising', options=[True, False], key='perform denoising', help=f_denoising_help)
                with col2:
                    st.text_input(label='alpha', key='alpha')
                with col3:
                    st.selectbox(label='Threshold type', options=['absolute', 'relative', 'power law'], key='threshold type', help=f_denoising_threshold_help)
                    if st.session_state['threshold type'] == 'power law':
                        st.info('Power-law provides a data-driven threshold.')
                    else:
                        st.text_input(label='Size threshold [absolute nr / %]', key='size threshold [absolute nr / %]')
                # SWARM Clustering
                with st.expander("🐝 SWARM Clustering", expanded=False):
                    st.selectbox(label='Perform SWARM clustering', options=[True, False], key='perform swarm clustering', help=g_swarm_clustering_help)
                # Replicate Merging
                with st.expander("🧪 Replicate Merging", expanded=False):
                    col1, col2, col3 = st.columns(3)
                with col1:
                    st.selectbox(label='Perform replicate merging', options=[True, False], key='perform replicate merging')
                with col2:
                    st.text_input(label='Replicate delimiter', key='replicate delimiter')
                with col3:
                    st.text_input(label='Minimum replicate presence', key='minimum replicate presence')
                # NC removal
                with st.expander("🚫 NC removal", expanded=False):
                    col1, col2 = st.columns(2)
                with col1:
                    st.selectbox(label='Negative control removal', options=[True, False], key='perform nc removal')
                with col2:
                    st.text_input(label='Negative control prefix', key='negative control prefix')
                # Generate Read Table
                with st.expander("📊 Generate Read Table", expanded=False):
                    col1, col2 = st.columns(2)
                with col1:
                    st.selectbox(label='Generate read table', options=[True, False], key='generate read table', help=j_generate_read_table_help)
                with col2:
                    st.text_input(label='Sequence group threshold', key='sequence group threshold')

                ############################################################################################################
                st.subheader('Run apscale')
                options = ['Run apscale (basic mode)', 'Run apscale (complete mode)', 'PE-merging', 'Primer-trimming', 'Quality-filtering', 'Dereplication', 'Denoising', 'SWARM clustering', 'Replicate merging', 'NC removal', 'Generate read table']
                st.selectbox(label='Select module to run', options=options, index=0, key='run_apscale_mode')

                if st.session_state["P5 Primer (5' - 3')"] == '' or st.session_state["P7 Primer (5' - 3')"] == '' or st.session_state['min length'] == '' or st.session_state['max length'] =='':
                    st.error('Please fill out all required fields!')
                else:
                    if st.session_state['run_apscale_mode'] == 'Run apscale (basic mode)':
                        st.info('The "Basic mode" mode skips "Replicate merging" and "NC removal".')
                    if st.session_state['run_apscale_mode'] == 'Run apscale (complete mode)':
                        st.info('The "Complete mode" runs all modules (except specifically disabled above).')
                    if st.button('Start raw data analysis'):
                        update_settings_file(settings_xlsx, settings_dfs)
                        run_apscale(st.session_state['run_apscale_mode'], project_folder)


                ############################################################################################################
                st.markdown("---")
                st.header("APSCALE BLAST")
                st.subheader('Settings')
                # General settings
                with st.expander("⚙️ General settings", expanded=False):
                    col1, col2 = st.columns(2)
                with col1:
                    st.text_input(label='n_cores', key='n_cores', value=multiprocessing.cpu_count()-2)
                    st.text_input(label='subset_size', key='subset_size', value=100)
                with col2:
                    st.selectbox(label='Task', key='task', options=['blastn', 'megablast', 'dc-megablast'])
                    st.text_input(label='max_target_seqs', key='max_target_seqs', value=20)
                # Database & query
                available_databases = {Path(i).name:Path(i) for i in glob.glob(str(path_to_projects / 'APSCALE_databases' / '*'))}
                available_fasta_files = {Path(i).name:Path(i) for i in glob.glob(str(project_folder / '11_read_table' / 'data' / '*.fasta'))}
                with st.expander("🗄️ Database & Query", expanded=False):
                    st.selectbox(label='Database', key='database', options=list(available_databases.keys()) + ['remote'])
                    st.selectbox(label='Query FASTA', key='query_fasta', options=list(available_fasta_files.keys()))
                # Thresholds
                with st.expander("⚖️ Thresholds", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.text_input(label='Species (%)', key='t_species', value=97)
                        st.text_input(label='Family (%)', key='t_family', value=90)
                        st.text_input(label='Class (%)', key='t_class', value=85)
                    with col2:
                        st.text_input(label='Genus (%)', key='t_genus', value=95)
                        st.text_input(label='Order (%)', key='t_order', value=87)
                        st.selectbox(label='Masking', key='masking', options=[True, False], index=1)
                # Remote blastn
                with st.expander("🌐 Remote blastn settings", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.selectbox(label='Include "uncultured"', options=[False, True], key='include_uncultured')
                        st.selectbox(label='Disable headless mode', options=[True, False], key='disable_headless')
                    with col2:
                        st.text_input(label='Organism filter (Comma separated taxa)', key='organism_filter', value='Eukaryota')
                        if st.button('Update taxid database'):
                            run_update_taxids()

                st.subheader('Run apscale-blast')
                if st.button(f'Start taxonomic assignment'):
                    if st.session_state['task'] != 'megablast' and st.session_state['database'] == 'remote':
                        st.error('Please select task "megablast" to perform the remote blast.')
                    else:
                        run_apscale_blast(project_folder, available_fasta_files, available_databases)

                ############################################################################################################
                st.markdown("---")
                st.header("BOLDigger3")
                st.subheader('Settings')
                # General settings
                with st.expander("🗄️ Database & Operating Mode", expanded=False):
                    bold_databases = {
                                    'ANIMAL LIBRARY (PUBLIC)': 1,
                                    'ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)': 2,
                                    'ANIMAL LIBRARY (PUBLIC+PRIVATE)': 3,
                                    'VALIDATED CANADIAN ARTHROPOD LIBRARY': 4,
                                    'PLANT LIBRARY (PUBLIC)': 5,
                                    'FUNGI LIBRARY (PUBLIC)': 6,
                                    'ANIMAL SECONDARY MARKERS (PUBLIC)': 7,
                                    'VALIDATED ANIMAL RED LIST LIBRARY': 8
                                }
                    st.selectbox(label='BOLD database', key='bold_database', options=list(bold_databases.keys()))
                    bold_modes = {
                                    'Rapid Species Search': 1,
                                    'Genus and Species Search': 2,
                                    'Exhaustive Search': 3
                                }
                    st.selectbox(label='Operating mode', key='bold_mode', options=list(bold_modes.keys()))
                    st.selectbox(label='Query FASTA', key='bold_query_fasta', options=list(available_fasta_files.keys()))
                # Thresholds
                with st.expander("⚖️ Thresholds", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.text_input(label='Species (%)', key='bold_species', value=97)
                        st.text_input(label='Family (%)', key='bold_family', value=90)
                        st.text_input(label='Class (%)', key='bold_class', value=75)
                    with col2:
                        st.text_input(label='Genus (%)', key='bold_genus', value=95)
                        st.text_input(label='Order (%)', key='bold_order', value=85)

                st.subheader('Run BOLDigger3')
                if st.button(f'Start taxonomic assignment against BOLD'):
                    run_boldigger3(available_fasta_files, bold_modes, bold_databases)

                ############################################################################################################
                st.markdown("---")
                st.header("📚 Links and Tutorials")
                # GitHub Projects
                with st.expander("🔗 GitHub Repositories", expanded=False):
                    st.markdown("""
                    - [APSCALE](https://github.com/DominikBuchner/apscale)  
                    - [APSCALE-blast](https://github.com/TillMacher/apscale_blast)  
                    - [APSCALE-GUI](https://github.com/TillMacher/apscale_gui)  
                    - [TaxonTableTools2](https://github.com/TillMacher/TaxonTableTools2)  
                    - [SWARM](https://github.com/torognes/swarm)  
                    - [Demultiplexer2](https://github.com/DominikBuchner/demultiplexer2)  
                    """)
                # Tutorials
                with st.expander("🎥 Video Tutorials", expanded=False):
                    st.markdown("""
                    - [Raw data processing tutorial](https://www.youtube.com/watch?v=SV7EJ1w-0u4&t=27s)  
                    - [Installation tutorial](https://www.youtube.com/watch?v=SV7EJ1w-0u4&t=27s)  
                    """)
                # Manual
                with st.expander("📖 Documentation", expanded=False):
                    st.markdown("""
                    - [APSCALE Manual (PDF)](https://github.com/DominikBuchner/apscale/blob/main/manual/apscale_manual.pdf)  
                    """)
                # Citations
                with st.expander("📑 Citations", expanded=False):
                    st.markdown("""
                    
                    - Buchner, D., Macher, T.-H., & Leese, F. (2022). APSCALE: Advanced pipeline for simple yet comprehensive analyses of DNA metabarcoding data. Bioinformatics, btac588. https://doi.org/10.1093/bioinformatics/btac588
                    
                    - Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. Journal, 17(1), Article 1.
                    
                    - Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.

                    - Mahé, F., Czech, L., Stamatakis, A., Quince, C., de Vargas, C., Dunthorn, M., & Rognes, T. (2021). Swarm v3: Towards tera-scale amplicon clustering. Bioinformatics, 38(1), 267–269. https://doi.org/10.1093/bioinformatics/btab493

                    """)
                with st.expander("📦 Package versions", expanded=False):
                    get_package_versions()

if __name__ == "__main__":
    main()


