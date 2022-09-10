import glob, sys, os, pkgutil, ast, pkg_resources, datetime, psutil, subprocess
import PySimpleGUI as sg
from pathlib import Path
import pandas as pd
import webbrowser
from apscale.b_pe_merging import main as pe_merging_main
from apscale.c_primer_trimming import main as primer_trimming
from apscale.d_quality_filtering import main as quality_filtering
from apscale.e_dereplication_pooling import main as dereplication_pooling
from apscale.f_otu_clustering import main as otu_clustering
from apscale.g_denoising import main as g_denoising
from apscale.h_lulu_filtering import main as h_lulu_filtering
from apscale_gui.clean_up_data import clean_up_data
from apscale_gui.fastq_eestats import main as fastq_eestats
from apscale.a_create_project import create_project
from apscale_gui.blast_utilities import subset_fasta
from apscale_gui.blast_utilities import blast_xml_to_taxonomy
from apscale_gui.blast_utilities import create_database_diat_barcode
from apscale_gui.blast_utilities import create_database_midori2
from apscale_gui.blast_utilities import blastn
from apscale_gui.blast_utilities import filter_blast_results
from apscale_gui.blast_utilities import filter_blastn_results_NCBI
from apscale_gui.blast_utilities import create_database_NCBI
from apscale_gui.summary_stats import main as summary_stats
from apscale_gui.settings_file import load_settings
from apscale_gui.settings_file import apply_settings
from apscale_gui.settings_file import settings_integrity
from apscale_gui.rename_files import create_rename_sheet
from apscale_gui.rename_files import rename_files
from apscale_gui.export_results import export_results
from lastversion import lastversion

##########################################################################################################################
# update version here (will be displayed on the main layout)
# Support for: u = ubuntu, w = windows, m = macintosh
apscale_version = '1.1.6'

## check for the latest version of TTT
try:
    latest_version = lastversion.has_update(repo="https://pypi.org/project/apscale-gui/", at='pip', current_version=apscale_version)
    if latest_version:
        sg.Popup('Newer apscale GUI version is available: ' + str(latest_version) + "\nPlease update apscale GUI via pip!", title="Update available!")
except Exception:
    pass

##########################################################################################################################
# general functions

# slice function for lists to split up lists
def slices(list, slice):
    for i in range(0, len(list), slice):
        yield list[i : i + slice]

def open_file(file):
    if sys.platform == 'win32':
        os.startfile(file)
    else:
        opener = 'open' if sys.platform == 'darwin' else 'xdg-open'
        subprocess.call([opener, file])

def check_installations():

    out = subprocess.Popen('vsearch --version', shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess_return = out.stdout.read().decode("utf-8")
    if subprocess_return != '':
        # print(subprocess_return)
        pass
    else:
        text = 'Could not find VSEARCH!\nPlease make sure you have VSEARCH installed in PATH!'
        sg.PopupError(text, title='Warning')


    out = subprocess.Popen('blastn -version', shell=True, stdout=subprocess.PIPE)
    subprocess_return = out.stdout.read().decode("utf-8")
    if subprocess_return != '':
        # print(subprocess_return)
        pass
    else:
        text = 'Could not find blastn!\nPlease make sure you have blastn installed in PATH!'
        print(text)
        sg.PopupError(text, title='Warning')

    out = subprocess.Popen('cutadapt --version', shell=True, stdout=subprocess.PIPE)
    subprocess_return = out.stdout.read().decode("utf-8")
    if subprocess_return != '':
        # print('Cutadapt version {}'.format(subprocess_return))
        pass
    else:
        text = 'Could not find cutadapt!\nPlease make sure you have cutadapt installed in PATH!'
        print(text)
        sg.PopupError(text, title='Warning')

##########################################################################################################################
#########################################################################################################################
##########################################################################################################################


def main():

    sg.ChangeLookAndFeel('Reddit')

    ##########################################################################################################################
    # check installations

    check_installations()

    ##########################################################################################################################
    # start Popup window

    ##########################################################################################################################
    # start Popup window
    # assign picture paths

    ## define paths to directories
    primersets_path = Path(pkg_resources.resource_filename(__name__, 'user_data/primersets/'))
    tagging_templates_path = Path(pkg_resources.resource_filename(__name__, 'user_data/tagging_templates/'))

    ## define paths to files
    user_data_txt = Path(pkg_resources.resource_filename(__name__, '_user_data/user_data.txt'))
    primer_sheet_xlsx = Path(pkg_resources.resource_filename(__name__, 'user_data/primer_trimming/primer_sheet.xlsx/'))
    ESV_reference_fasta_path = Path(pkg_resources.resource_filename(__name__, 'user_data/ESV_references/'))

    ## load figures
    crash_png = Path(pkg_resources.resource_filename(__name__, '/_source/crash.png'))
    github_png = Path(pkg_resources.resource_filename(__name__, '/_source/github.png'))
    twitter_png = Path(pkg_resources.resource_filename(__name__, '/_source/twitter.png'))
    quality_control_png = Path(pkg_resources.resource_filename(__name__, '/_source/quality_control.png'))
    sample_renaming_png = Path(pkg_resources.resource_filename(__name__, '/_source/sample_renaming.png'))
    demultiplexing_png = Path(pkg_resources.resource_filename(__name__, '/_source/demultiplexing.png'))
    all_in_one_analysis_png = Path(pkg_resources.resource_filename(__name__, '/_source/all_in_one_analysis.png'))
    local_blast_png = Path(pkg_resources.resource_filename(__name__, '/_source/local_blast.png'))
    ncbi_blast_png = Path(pkg_resources.resource_filename(__name__, '/_source/ncbi_blast.png'))
    boldigger_png = Path(pkg_resources.resource_filename(__name__, '/_source/boldigger.png'))
    summary_statistics_png = Path(pkg_resources.resource_filename(__name__, '/_source/summary_statistics.png'))
    export_results_png = Path(pkg_resources.resource_filename(__name__, '/_source/export_results.png'))

    ## open user_data_txt to save the standard output put
    f = open(user_data_txt)
    projects_main_path = f.read()

    # fresh start: there is an empty user_data file
    # ask for user Input
    # stay open until a path was defined
    # then write it the user_data file to reload
    while projects_main_path == '':
        projects_main_path = sg.PopupGetFolder('Enter path to ouput directory:', title='Output directory')
        if projects_main_path == None:
            sys.exit()
        f = open(user_data_txt, 'w')
        f.write(projects_main_path)
        f.close()

    # create display text
    current_path = 'Current path: ' + str(projects_main_path)
    # load all available projects
    projects = glob.glob(str(projects_main_path) + '/*_apscale')
    projects_list = []

    for project in projects:
        projects_list.append(Path(project).stem)

    projects_radio = list(slices([sg.Radio(name, 'projects', default=True) for name in sorted(projects_list)], 3))

    start_window_layout = [
                [sg.Text('',size=(1,1))],
    			[sg.Text('Output directory', size=(50,1), font=('Arial', 11, 'bold'))],
                [sg.Text(current_path)],
                [sg.Text('Define new output directory:')],
                [sg.Input(key = 'new_projects_main_path', size=(40,1)), sg.FolderBrowse(), sg.Button('Refresh')],
                [sg.Text('',size=(1,1))],
    			[sg.Text('Project management', size=(50,1), font=('Arial', 11, 'bold'))],
                [sg.Text('Create new project folder:')],
                [sg.Input('', key='new_project_folder', size=(40,1)), sg.Button('Create new')],
                [sg.Text('Load existing project folder:')],
                [sg.Frame(layout = projects_radio, title = '')],
                [sg.Button('Load')],
                [sg.Text('',size=(1,1))],
                [sg.Button('Exit', button_color=('black', 'red'))],
                ]

    start_window = sg.Window('Projects', start_window_layout)
    event, values = start_window.read()

    while True:

        new_projects_main_path = values['new_projects_main_path']

        if event == 'Create new':
            if values['new_project_folder'] != '':
                ## define output path and project folder
                project_folder = projects_main_path + '/' + values['new_project_folder'].replace(' ', '_')
                ## create a new project folder
                create_project(project_folder)
                project_folder = project_folder + '_apscale'
                break
            else:
                project_folder = 'Default_project'
                ## create a new project folder
                create_project(project_folder)
                project_folder = project_folder + '_apscale'
                break

        if event == 'Load':
            project_folder = ''
            for key, value in values.items():
                if value == True:
                    project_folder = projects_main_path + '/' + sorted(projects_list)[key]
            if project_folder == '':
                project_folder = 'Default_project'
            break

        if event == 'Refresh':
            if new_projects_main_path == None:
                break
            f = open(user_data_txt, 'w')
            f.write(new_projects_main_path)
            f.close()
            sg.Popup('Please reload APSCALE to apply changes', title='Refresh output directory')
            sys.exit()

        if event == 'Exit':
            sys.exit()

    start_window.close()

    ##########################################################################################################################
    # check output directory
    try:
        if not os.path.exists(Path(str(projects_main_path))):
            os.mkdir(Path(str(projects_main_path) + '/Projects'))
    except Exception:
        sg.PopupError('The output directory does not exist anymore! Please refresh the output folder.')
        sys.exit()

    ##########################################################################################################################
    ## load standard output path
    path_to_outdirs = Path(project_folder + '/')
    settings_file = Path(path_to_outdirs).joinpath('Settings.xlsx')

    ##########################################################################################################################
    # create gui related folders
    try:
        os.mkdir(Path(path_to_outdirs).joinpath('0_statistics'))
    except FileExistsError:
        pass

    try:
        os.mkdir(Path(path_to_outdirs).joinpath('10_local_BLAST'))
    except FileExistsError:
        pass

    try:
        os.mkdir(Path(path_to_outdirs).joinpath('11_NCBI_BLAST'))
    except FileExistsError:
        pass

    try:
        os.mkdir(Path(path_to_outdirs).joinpath('11_NCBI_BLAST', 'XML2_files'))
    except FileExistsError:
        pass

    databases_path = Path(projects_main_path).joinpath('APSCALE_databases')
    try:
        os.mkdir(databases_path)
    except FileExistsError:
        pass

    ##########################################################################################################################
    ##########################################################################################################################
    ##########################################################################################################################

    layout_tools = [
                        # row 1
                        [sg.Button(key='open_quality_control', button_color=('white', 'white'), image_filename=quality_control_png),
                        sg.Button(key='open_sample_renaming', button_color=('white', 'white'), image_filename=sample_renaming_png),
                        sg.Button(key='open_demultiplexing', button_color=('white', 'white'), image_filename=demultiplexing_png)],
                        # row 2
                        [sg.Button(key='open_run_analyses', button_color=('white', 'white'), image_filename=all_in_one_analysis_png),
                        sg.Button(key='open_local_blast', button_color=('white', 'white'), image_filename=local_blast_png),
                        sg.Button(key='open_ncbi_blast', button_color=('white', 'white'), image_filename=ncbi_blast_png)],
                        # row 3
                        [sg.Button(key='open_boldigger', button_color=('white', 'white'), image_filename=boldigger_png),
                        sg.Button(key='open_analysis_statistics', button_color=('white', 'white'), image_filename=summary_statistics_png),
                        sg.Button(key='open_compress_results', button_color=('white', 'white'), image_filename=export_results_png)],
                        ]

    ## define variables for the main window
    welcome_text = datetime.datetime.now().strftime('%H:%M:%S') + ': Welcome to the APSCALE!\n'

    layout = [  [sg.Text('APSCALE', font=('Arial', 13, 'bold')),
                sg.Text('Advanced Pipeline for Simple yet Comprehensive Analyses of DNA metabarcoding data')],
                [sg.Text('Current project:', font=('Arial', 11)), sg.Text(project_folder, font=('Arial', 11))],
    			[sg.Text('',size=(1,1))],
    			[sg.TabGroup([[sg.Tab('Metabarcoding processing tools', layout_tools)]])], ## size = (w, h)
    			[sg.Text('',size=(2,2))],
    			[sg.Exit(button_color=('black', 'red')),
                sg.Text('', size=(25,1)),
                sg.Image(str(github_png)),
                sg.Button('apscale', button_color=('black', 'white'), key='open_apscale_github'),
                sg.Button('GUI', button_color=('black', 'white'), key='open_apscale_gui_github'),
                sg.Text('', size=(2,1)),
                sg.Image(str(twitter_png)),
                sg.Button('TM', button_color=('black', 'white'), key='open_TM_twitter'),
                sg.Button('DB', button_color=('black', 'white'), key='open_DB_twitter'),
                sg.Text('', size=(4,1)),
                sg.Text('Version: ' + str(apscale_version), font=('Arial', 8))]]

    # Create the Window
    MP_window = sg.Window('APSCALE', layout)
    win2_active=False

    ##########################################################################################################################

    while True:
        try:
            event, values = MP_window.Read()

            if event is None or event == 'Exit':
                break

            if event == 'open_demultiplexing':
                try:
                    import demultiplexer.__main__ as demultiplexer
                    MP_window.hide()
                    demultiplexer.main()
                except Exception:
                    sg.PopupError('You have to install the demultipexer first!')

            if event == 'open_boldigger':
                try:
                    import boldigger.__main__ as boldigger
                    MP_window.hide()
                    boldigger.main()
                except Exception:
                    sg.PopupError('You have to install BOLDigger first!')

            if event == 'open_compress_results':
                MP_window.hide()
                answer = sg.PopupOKCancel('Export and compress results?', title='Log file')
                if answer == 'OK':
                    export_results(path_to_outdirs)
                MP_window.UnHide()

            if event == 'open_run_analyses':

                MP_window.hide()

                settings_dict = load_settings(settings_file)

                modify_settings_sheet = [[
                    sg.Text('Load settings file and modify as required.'),
                    sg.Button('Open settings file', key='open_settings_file'),
                    sg.Button('Apply new settings', key='apply_new_settings', button_color=('black', 'lightgreen')),
                    ]]

                general_settings_layout = [[
                    sg.Text('Write to Excel:'), sg.Combo(['True', 'False'], key='settings_clustering_to_excel', default_value=settings_dict['to excel']),
                    sg.Text('Write to Parquet:'), sg.Combo(['True'], default_value='True'),
                    ]]

                pe_merging_layout = [[
                    sg.Text('maxdiffpct:'), sg.Spin([i for i in range(1,101)], initial_value=settings_dict['maxdiffpct'], key='settings_maxdiffpct', size=(6,1)),
                    sg.Text('maxdiffs:'), sg.Spin([i for i in range(1,200)], initial_value=settings_dict['maxdiffs'], key='settings_maxdiffs', size=(6,1)),
                    sg.Text('minovlen:'), sg.Spin([i for i in range(1,26)], initial_value=settings_dict['minovlen'], key='settings_minovlen', size=(6,1))
                    ]]

                primer_trimming_layout = [[
                    sg.Text("P5 Primer (5' - 3'):"), sg.Input(settings_dict["P5 Primer (5' - 3')"], key='settings_p5_primer', size=(20,1)),
                    sg.Text("P7 Primer (5' - 3'):"), sg.Input(settings_dict["P7 Primer (5' - 3')"], key='settings_p7_primer', size=(20,1)),
                    sg.Text('Anchoring:'), sg.Combo(['True', 'False'], key='settings_anchoring', default_value=settings_dict['anchoring'])
                    ]]

                read_filter_layout = [[
                    sg.Text('maxEE:'), sg.Input(settings_dict['maxEE'], key='settings_maxEE', size=(4,1)),
                    sg.Text('min length:'), sg.Input(settings_dict['min length'], key='settings_min_length', size=(6,1)),
                    sg.Text('max length:'), sg.Input(settings_dict['max length'], key='settings_max_length', size=(6,1)),
                    ]]

                pre_processing_layout = [[
                    sg.Text('min size to pool:'), sg.Spin([i for i in range(2,6)], initial_value=settings_dict['min size to pool'], key='settings_derep_minsize', size=(6,1)),
                    ]]

                clustering_layout = [[
                    sg.Text('pct id:'),
                    sg.Spin([i for i in range(80,101)], initial_value=settings_dict['pct id'], key='settings_pct_id', size=(6,1)),
                    ]]

                denoising_layout = [[
                    sg.Text('alpha:'), sg.Spin([i for i in range(1,11)], initial_value=settings_dict['alpha'], key='settings_alpha', size=(6,1)),
                    sg.Text('min size:'), sg.Spin([i for i in range(1,11)], initial_value=settings_dict['minsize'], key='settings_min_size', size=(6,1)),
                    ]]

                lulu_layout = [[
                    sg.Text('min similarity:'), sg.Spin([i for i in range(1,101)], initial_value=settings_dict['minimum similarity'], key='settings_min_sim', size=(6,1)),
                    sg.Text('min rel. co-occurrence:'), sg.Spin([i for i in range(1,101)], initial_value=settings_dict['minimum relative cooccurence'], key='settings_min_rel_coo', size=(6,1)),
                    sg.Text('min ratio:'), sg.Input(settings_dict['minimum ratio'], key='settings_min_ratio', size=(6,1)),
                    ]]

                clean_up_layout = [[
                    sg.Text('Remove all intermediate data (2.-4.) to safe storage space. OTUs and ESVs can still be generated.'),
                ]]

                layout_run_analyses = [
                					[sg.Text('All-in-One analysis', size=(50,1), font=('Arial', 12, 'bold'))],
                					[sg.Text('_'*115)],

                                    [sg.Text('', size=(4,1)), sg.Text('  Modify settings:', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=modify_settings_sheet, title='')],
                					[sg.Text('')],

                                    [sg.Text('', size=(4,1)), sg.Text('1. General settings:', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=general_settings_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_pe_merging'), sg.Text('2. Paired-end merging', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=pe_merging_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_primer_trimming'), sg.Text('3. Primer trimming', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=primer_trimming_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_quality_filtering'), sg.Text('4. Quality filtering', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=read_filter_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_dereplication_pooling'), sg.Text('5. Dereplication & pooling', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=pre_processing_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_otu_clustering'), sg.Text('6.1 OTU clustering', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=clustering_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_denoising'), sg.Text('6.2 Denoising', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=denoising_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=True, key='cb_lulu'), sg.Text('6.3 LULU filtering', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=lulu_layout, title='')],
                					[sg.Text('')],

                                    [sg.CB('', default=False, key='cb_clean_up'), sg.Text('Data clean-up', size=(22,1), font=('Arial', 11, 'bold')), sg.Frame(layout=clean_up_layout, title='')],
                					[sg.Text('')],

                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Run analysis', size=(10,2)), sg.CB('Minimize APSCALE', default=True, key='minimize_window')],
                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Exit', button_color=('black', 'red'))]
                					]

                # create the demultiplexing window
                tools_window = sg.Window('All-in-One analysis', layout_run_analyses, keep_on_top=False)

                while (True):
                    ######################################
                    event, values2 = tools_window.Read()

                    ######################################
                    if event in ('Exit', None):
                        break

                    if event == 'Run analysis':
                        ## check if settings from GUI and settings file match
                        answer = sg.PopupOKCancel('Please remember to apply your settings prior to running APSCALE!', title='Warning')

                        ## only proceed if settings have been applied from gui to sheet
                        if answer == 'OK':
                            ######################################
                            if values2['minimize_window'] == True:
                                tools_window.hide()
                                tools_window.refresh()

                            ## check the settings file for integrity
                            test = settings_integrity(settings_file)

                            if test != False and test != 'Cancel':
                                ## continue with pipeline only if settings file is legit
                                print('')

                                ############################################################################
                                if values2['cb_pe_merging'] == True:
                                    print('Starting paired-end merging')
                                    pe_merging_main(str(path_to_outdirs))
                                    print('')

                                if values2['cb_primer_trimming'] == True:
                                    print('Starting primer trimming')
                                    primer_trimming(str(path_to_outdirs))
                                    print('')

                                if values2['cb_quality_filtering'] == True:
                                    print('Starting quality filtering')
                                    quality_filtering(str(path_to_outdirs))
                                    print('')

                                if values2['cb_dereplication_pooling'] == True:
                                    print('Starting dereplication and pooling')
                                    dereplication_pooling(str(path_to_outdirs))
                                    print('')

                                if values2['cb_otu_clustering'] == True:
                                    print('Starting OTU clustering')
                                    otu_clustering(str(path_to_outdirs))
                                    print('')

                                if values2['cb_denoising'] == True:
                                    print('Starting denoising')
                                    g_denoising(str(path_to_outdirs))
                                    print('')

                                if values2['cb_lulu'] == True:
                                    print('Starting LULU filtering')
                                    h_lulu_filtering(str(path_to_outdirs))
                                    print('')

                                if values2['cb_clean_up'] == True:
                                    answer = sg.PopupOKCancel('Warning: You selected the data clean-up:\nThis will erase excess data permanently. Raw reads and pooled reads will be kept!', title='Warning')
                                    if answer == 'OK':
                                        print('Cleaning up the working directory.')
                                        clean_up_data(path_to_outdirs)
                                        print('')

                            ############################################################################

                            ######################################
                            ## finish the command chain
                            finished = sg.PopupYesNo('Jobs finished!\n\nContinue analyses?', title='')
                            if finished == 'No':
                                break
                                tools_window.Close()
                            elif values2['minimize_window'] == True:
                                tools_window.UnHide()

                            elif values2['minimize_window'] == True:
                                tools_window.UnHide()

                    if event == 'apply_new_settings':
                        settings = {i:str(values2[i]) for i in list(values2.keys()) if 'settings_' in str(i)}
                        answer = sg.PopupOKCancel('Warning: This will overwrite the settings file!', title='Warning')
                        if answer == 'OK':
                            apply_settings(settings_file, settings)
                            sg.Popup('New settings have been applied.', title='')

                    if event == 'open_settings_file':
                        open_file(settings_file)

                    if event == 'modify_primer_sheet':
                        open_file(primer_sheet_xlsx)

                        ##########################################################################################

                tools_window.Close()

            if event == 'open_ncbi_blast':

                MP_window.hide()

                layout_ncbi_blast = [
                					[sg.Text('NCBI BLAST', size=(50,1), font=('Arial', 12, 'bold'))],
                					[sg.Text('_'*115)],
                                    [sg.Frame(layout=[
                                    [sg.Text('1. BLAST your OTU or ESV fasta file against the NCBI database:'), sg.Button('NCBI blast website', key='open_ncbi_website')],
                                    [sg.Text('2. When results are available, click on \'Download all\' and select the \'Single-file XML2\' format.')],
                                    [sg.Text('3. Store the xml2 file in the \'9_NCBI_BLAST/xml2_files\' folder.')]
                                    ], title='BLAST search')],
                                    [sg.Frame(layout=[
                                    [sg.Text('If the fasta files contains too many entries and the NCBI BLAST keeps failing, create subsets here:')],
                                    [sg.Text('Sequence batch size:'), sg.Input('150', size=(5,1), key='batch_size'), sg.Button('Split fasta file', key='subset_fasta')],
                                    [sg.Text('BLAST all fasta files seperately and download the xml2 files as described above.')],
                                    ], title='Fasta subset')],
                                    [sg.Frame(layout=[
                                    [sg.Text('1. Select the original fasta file.')],
                                    [sg.Text('Fasta file:', size=(20,1)), sg.Input('', size=(30,1), key='fasta_file'), sg.FileBrowse('Browse', initial_folder = path_to_outdirs)],
                                    [sg.Text('2. Select the read table.')],
                                    [sg.Text('Read table:', size=(20,1)), sg.Input('', size=(30,1), key='read_table'), sg.FileBrowse('Browse', initial_folder = path_to_outdirs)],
                                    [sg.Text('3. Select all xml files downloaded from NCBI.')],
                                    [sg.Text('BLAST xml file(s):', size=(20,1)), sg.Input('', size=(30,1), key='xml_files'), sg.FilesBrowse('Browse', initial_folder = path_to_outdirs)],
                                    [sg.Text('4. Select a limit of how many hits should be considered.')],
                                    [sg.Text('BLAST hit limit:'), sg.Input('10', size=(5,1), key='limit')],
                                    ], title='Input files')],
                					[sg.Text('')],
                                    [sg.Button('Create taxonomy table', size=(20,2), key=('fetch_taxonomy')), sg.CB('Minimize APSCALE', default=True, key='minimize_window')],
                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Exit', button_color=('black', 'red'))]
                					]

                # create the demultiplexing window
                blast_window = sg.Window('NCBI Blast', layout_ncbi_blast, keep_on_top=False)
                while (True):
                    ######################################
                    event, values2 = blast_window.Read()
                    ######################################

                    fasta_file = values2['fasta_file']
                    read_table = values2['read_table']
                    limit = values2['limit']
                    batch_size = values2['batch_size']
                    xml_files = values2['xml_files'].split(';')

                    if event in ('Exit', None):
                        break

                    if event == 'subset_fasta':
                        if fasta_file == '':
                            sg.PopupError('Please provide a fasta file!', title='Error')
                        else:
                            print('')
                            subset_fasta(fasta_file, batch_size)
                            print('')

                    if event == 'fetch_taxonomy':
                        if fasta_file == '':
                            sg.PopupError('Please provide a fasta file!', title='Error')
                        elif read_table == '':
                            sg.PopupError('Please provide a read table!', title='Error')
                        else:
                            print('')
                            blast_xml_to_taxonomy(fasta_file, xml_files, read_table, limit)
                            print('')

                    if event == 'open_ncbi_website':
                        webbrowser.open('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome')

                blast_window.Close()

            if event == 'open_sample_renaming':

                MP_window.hide()

                layout_sample_renaming = [
                					[sg.Text('Sample renaming', size=(50,1), font=('Arial', 12, 'bold'))],
                					[sg.Text('_'*115)],
                                    [sg.Text('Select folder with files to rename:'), sg.Combo(['1_raw data', '2_demultiplexing'], default_value='1_raw data', key='rename_folder')],
                                    [sg.Text('Paired-end reads:'), sg.Combo(['True', 'False'], default_value='True', key='paired_end_reads')],
                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Create rename sheet', key='run_create_rename_sheet'), sg.Button('Open rename sheet  ', key='run_open_rename_sheet')],
                                    [sg.Button('Rename samples     ', key='run_rename_samples')],
                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Exit', button_color=('black', 'red'))]
                					]

                # create the demultiplexing window
                rename_files_window = sg.Window('Rename files', layout_sample_renaming, keep_on_top=False)
                while (True):
                    ######################################
                    event, values2 = rename_files_window.Read()
                    ######################################

                    if event in ('Exit', None):
                        break

                    if event == 'run_create_rename_sheet':
                        create_rename_sheet(values2['paired_end_reads'], values2['rename_folder'], path_to_outdirs)
                        answer = sg.PopupOKCancel('Finished.\nOpen rename sheet?', title='')
                        if answer == 'OK':
                            open_file(Path(str(path_to_outdirs) + '/rename_sheet.xlsx'))

                    if event == 'run_open_rename_sheet':
                        open_file(Path(str(path_to_outdirs) + '/rename_sheet.xlsx'))

                    if event == 'run_rename_samples':
                        answer = sg.PopupOKCancel('Warning this will overwrite ALL files from the rename sheet!\nThere is no going back!\nContinue?', title='Warning')
                        if answer == 'OK':
                            rename_files(values2['paired_end_reads'], path_to_outdirs)
                            sg.Popup('All files were renamed accordingly.')

                rename_files_window.Close()

            if event == 'open_quality_control':

                MP_window.hide()

                quality_control_layout = [
                					[sg.Text('Collect summary statistics for this project?')],
                                    [sg.Text('Folder:'), sg.DropDown(['1_raw data','2_demultiplexing', '3_PE_merging', '4_primer_trimming'], default_value='2_demultiplexing', key='raw_reads_analysis_folder')],
                                    [sg.Button('OK', key='run_analyse_raw_reads'), sg.Button('Cancel')]
                					]

                # create the demultiplexing window
                quality_control_window = sg.Window('Quality control', quality_control_layout, keep_on_top=False)
                while (True):
                    ######################################
                    event, values2 = quality_control_window.Read()
                    ######################################

                    if event in ('Cancel', None):
                        break

                    if event == 'run_analyse_raw_reads':
                        ## hide window
                        quality_control_window.hide()
                        quality_control_window.refresh()

                        print('Starting raw read analysis.')
                        fastq_eestats(values2['raw_reads_analysis_folder'], path_to_outdirs)
                        print('')

                        ## finish script
                        sg.Popup('Finished quality control\nAll plots are found in the 0_statistics folder.')

                        ## unhide window
                        quality_control_window.UnHide()

                quality_control_window.Close()

            if event == 'open_local_blast':
                MP_window.hide()

                reference_databases = ['Diat.barcode', 'Custom NCBI', 'Midori2']

                layout_local_blast = [
                					[sg.Text('Local BLAST', size=(50,1), font=('Arial', 12, 'bold'))],
                					[sg.Text('_'*90)],

                                    [sg.Frame(layout=[
                                    [sg.Text('1. Select your fasta file:', size=(25,1)), sg.Input('', size=(10,1)), sg.FileBrowse(key='local_blast_fasta_file', initial_folder = path_to_outdirs)],
                                    [sg.Text('2. Select your read table:', size=(25,1)), sg.Input('', size=(10,1)), sg.FileBrowse(key='local_blast_read_table', initial_folder = path_to_outdirs)],
                                    [sg.Text('', size=(90,1))],
                                    ], title='Load required files')],

                                    [sg.Frame(layout=[
                                    [sg.Text('1. Select a database format:', size=(25,1)), sg.Combo(reference_databases, default_value='Diat.barcode', key='database_selection'),
                                    sg.Text('Get more information about the database:'), sg.Button(' ? ', key='database_info')],
                                    [sg.Text('2. Select database source file:', size=(25,1)), sg.Input('', size=(10,1)), sg.FileBrowse(key='database_source_file', initial_folder = path_to_outdirs)],
                                    [sg.Text('3. Build a new database:', size=(25,1)), sg.Button('Go', key='database_builder')],
                                    [sg.Text('', size=(90,1))],
                                    ], title='Build a new database')],

                                    [sg.Frame(layout=[
                                    [sg.Text('1. Select a database', size=(25,1)), sg.Input('', size=(10,1)), sg.FolderBrowse(key='blast_database', initial_folder = databases_path)],
                                    [sg.Text('2. Run BLASTn', size=(25,1)), sg.Button('Go!', key='run_blast'), sg.Combo(['Highly similar sequences (megablast)', 'More dissimilar sequences (discontiguous megablast)', 'Somewhat similar sequences (blastn)'], default_value='Somewhat similar sequences (blastn)', key='blast_task')],
                                    [sg.Text('')],
                                    [sg.Text('3. Select BLAST results (.csv)', size=(25,1)), sg.Input('', size=(10,1)), sg.FileBrowse(key='blast_csv', initial_folder = Path(path_to_outdirs).joinpath('10_local_BLAST'))],
                                    [sg.Text('4. Filter BLAST results', size=(25,1)), sg.Button('Go!', key='run_blast_filter'), sg.Text('Select taxonomy source:'), sg.Combo(reference_databases, default_value=reference_databases[0], key='database_taxonomy_source')],
                                    [sg.Text('Similarity threshold (%):'),
                                    sg.Text('Species >='), sg.Spin([i for i in range(1,101)], initial_value=98, key='threshold_species'),
                                    sg.Text('Genus >='), sg.Spin([i for i in range(1,101)], initial_value=95, key='threshold_genus'),
                                    sg.Text('Family >='), sg.Spin([i for i in range(1,101)], initial_value=90, key='threshold_family'),
                                    sg.Text('Order >='), sg.Spin([i for i in range(1,101)], initial_value=85, key='threshold_order')],
                                    [sg.Text('E-value threshold:'), sg.Text('E-value >='), sg.Input('0.0005', size=(10,1), key='threshold_e_value')],
                                    [sg.Text('', size=(90,1))],
                                    ], title='Assign taxonomy using BLASTn')],

                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Exit', button_color=('black', 'red'))]
                					]

                # create the demultiplexing window
                blast_window = sg.Window('Local Blast', layout_local_blast, keep_on_top=False)
                while (True):
                    ######################################
                    event, values2 = blast_window.Read()
                    ######################################

                    if event in ('Exit', None):
                        blast_window.close()
                        break

                    if event == 'database_info':
                        webbrowser.open('https://github.com/TillMacher/apscale_gui#available-databases-for-local-blast')

                    if event == 'database_builder':

                        blast_window.Hide()

                        if values2['database_source_file'] == '':
                            sg.Popup('Please provide the database source file!', title='Error')

                        elif values2['database_selection'] == 'Diat.barcode':
                            create_database_diat_barcode(values2['database_source_file'], project_folder, databases_path)

                        elif values2['database_selection'] == 'Custom NCBI':
                            create_database_NCBI(values2['database_source_file'], project_folder, databases_path)

                        elif values2['database_selection'] == 'Midori2':
                            create_database_midori2(values2['database_source_file'], project_folder, databases_path)

                        blast_window.UnHide()

                    if event == 'run_blast':
                        if values2['local_blast_fasta_file'] == '' or values2['local_blast_read_table'] == '' or values2['blast_database'] == '':
                            sg.Popup('Please provide all files!', title='Error')
                        else:
                            blast_window.Hide()
                            sg.Popup('Starting BLAST!\nThis may take a while!', title='Finished')
                            blastn(values2['local_blast_fasta_file'], values2['blast_database'], project_folder, 4, values2['blast_task'])
                            sg.Popup('Finished BLAST!', title='Finished')
                            blast_window.UnHide()

                    if event == 'run_blast_filter':
                        if values2['blast_csv'] == '' or values2['local_blast_read_table'] == '':
                            sg.Popup('Please provide all files!', title='Error')
                        else:
                            blast_window.Hide()

                            ## ask if the correct database was selected.
                            answer = sg.PopupOKCancel('You selected \'{}\' as taxonomy source. Is that correct?'.format(values2['database_taxonomy_source']))
                            if answer == 'OK':
                                sg.Popup('Starting to filter BLAST results!\nThis may take a while!', title='Finished')

                                ## filter blast results
                                filter_tresholds = [int(values2['threshold_species']), int(values2['threshold_genus']), int(values2['threshold_family']), int(values2['threshold_order'])]
                                filter_blast_results(values2['blast_csv'], values2['local_blast_read_table'], project_folder, values2['blast_database'], filter_tresholds, values2['database_taxonomy_source'], float(values2['threshold_e_value']))

                                sg.Popup('Finished filtering BLAST results!', title='Finished')

                            blast_window.UnHide()

                blast_window.close()
                MP_window.UnHide()

            if event == 'open_analysis_statistics':
                MP_window.hide()
                answer = sg.PopupOKCancel('Collect summary statistics for this project?', title='Summary statistics')
                if answer == 'OK':
                    summary_stats(path_to_outdirs)
                    sg.Popup('Finished building summary statistics plots.\nAll plots are found in the 0_statistics folder.')
                    print('')
                    MP_window.UnHide()
                MP_window.UnHide()

            if event == 'open_TM_twitter':
                webbrowser.open('https://twitter.com/TillMacher')

            if event == 'open_DB_twitter':
                webbrowser.open('https://twitter.com/buchner_dominik')

            if event == 'open_apscale_gui_github':
                webbrowser.open('https://github.com/TillMacher/apscale_gui')

            if event == 'open_apscale_github':
                webbrowser.open('https://github.com/DominikBuchner/apscale')

        ###########################################################
            try:
                MP_window.UnHide()
            except Exception:
                pass
        ###########################################################

        # define exceptions
        # if there are unexpected errors print a message and continue the script!
        except Exception:
            ## unhide the main window if neccessary
            try:
                MP_window.UnHide()
            except Exception:
                pass
            ## close the secondary window if neccessary
            try:
                tools_window.Close()
            except Exception:
                pass

            ## create layout for the error message
            layout = [
                        [sg.Image(str(crash_png)), sg.Text(" You've been awarded with the gold medal in program crashing!")],
                        [sg.Text('', size=(1,1))],
                        [sg.Text('Unexpected error: ' + str(sys.exc_info()[0]))],
                        [sg.Text('', size=(1,1))],
                        [sg.Text('An unexpected error occured:')],
                        [sg.Text('Please refer to the manual.')],
                        [sg.Text('', size=(1,1))],
                        [sg.Button('Return'), sg.Button('Exit', button_color=('black', 'red'))]
                     ]

            window_error = sg.Window('Error', layout, keep_on_top=True)

            while (True):

                event, values = window_error.Read()

                ## return -> ignore error and go back to the main window
                if event in ('Return', None):
                    break
                ## exit -> raise error and go close the program
                if event in ('Exit', None):
                    raise

            window_error.Close()

        ###############################################################################################

    MP_window.Close()

## run only if called as toplevel script
if __name__ == '__main__':
    main()
