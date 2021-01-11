import glob, sys, os, pkgutil, ast, multiprocessing, pkg_resources, datetime, psutil, subprocess
import PySimpleGUI as sg
import matplotlib.pyplot as plt
from pathlib import Path
import webbrowser

##########################################################################################################################
# update version here (will be displayed on the main layout)
# Support for: u = ubuntu, w = windows, m = macintosh
meta_tools_version = "Version 0.3.beta"

##########################################################################################################################
# general functions

# slice function for lists to split up lists
def slices(list, slice):
    for i in range(0, len(list), slice):
        yield list[i : i + slice]

##########################################################################################################################
#########################################################################################################################
##########################################################################################################################

def main():

    sg.ChangeLookAndFeel('Reddit')

    ##########################################################################################################################
    # start Popup window

    ##########################################################################################################################
    # start Popup window
    # assign picture paths

    ## define paths to directories
    primersets_path = Path(pkg_resources.resource_filename(__name__, 'user_data/primersets/'))
    tagging_templates_path = Path(pkg_resources.resource_filename(__name__, 'user_data/tagging_templates/'))

    ## define paths to files
    user_data_txt = Path(pkg_resources.resource_filename(__name__, 'user_data/user_data.txt'))
    primer_sheet_xlsx = Path(pkg_resources.resource_filename(__name__, 'user_data/primer_trimming/primer_sheet.xlsx/'))

    ## load figures
    crash_png = Path(pkg_resources.resource_filename(__name__, '/_source/crash.png'))
    github_png = Path(pkg_resources.resource_filename(__name__, "/_source/github.png"))
    twitter_png = Path(pkg_resources.resource_filename(__name__, "/_source/twitter.png"))
    demultiplexing_png = Path(pkg_resources.resource_filename(__name__, '/_source/1_demultiplexing.png'))
    pe_merging_png = Path(pkg_resources.resource_filename(__name__, '/_source/2_pe_merging.png'))
    primer_trimming_png = Path(pkg_resources.resource_filename(__name__, '/_source/3_primer_trimming.png'))
    length_filtering_png = Path(pkg_resources.resource_filename(__name__, '/_source/4_length_filtering.png'))
    quality_filtering_png = Path(pkg_resources.resource_filename(__name__, '/_source/5_quality_filtering.png'))
    clustering_png = Path(pkg_resources.resource_filename(__name__, '/_source/6_clustering.png'))
    postprocessing_png = Path(pkg_resources.resource_filename(__name__, '/_source/7_postprocessing.png'))
    blast_png = Path(pkg_resources.resource_filename(__name__, '/_source/8_blast.png'))
    all_in_one_analysis_png = Path(pkg_resources.resource_filename(__name__, '/_source/9_all_in_one_analysis.png'))

    ## open user_data_txt to save the standard output put
    f = open(user_data_txt)
    projects_main_path = f.read()

    # fresh start: there is an empty user_data file
    # ask for user Input
    # stay open until a path was defined
    # then write it the user_data file to reload
    while projects_main_path == "":
        projects_main_path = sg.PopupGetFolder("Enter path to ouput directory:", title="Output directory")
        if projects_main_path == None:
            sys.exit()
        f = open(user_data_txt, "w")
        f.write(projects_main_path)
        f.close()

    # create display text
    current_path = "Current path: " + str(projects_main_path)
    # load all available projects
    projects = glob.glob(str(projects_main_path) + '/Projects/*')
    projects_list = []

    for project in projects:
        projects_list.append(Path(project).stem)

    projects_radio = list(slices([sg.Radio(name, "projects", default=True) for name in sorted(projects_list)], 3))

    start_window_layout = [
                [sg.Text('',size=(1,1))],
    			[sg.Text('Output directory', size=(50,1), font=('Arial', 11, "bold"))],
                [sg.Text(current_path)],
                [sg.Text('Define new output directory:')],
                [sg.Input(key = "new_projects_main_path", size=(40,1)), sg.FolderBrowse(), sg.Button("Refresh")],
                [sg.Text('',size=(1,1))],
    			[sg.Text('Project management', size=(50,1), font=('Arial', 11, "bold"))],
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

        new_projects_main_path = values["new_projects_main_path"]

        if event == 'Create new':
            if values["new_project_folder"] != '':
                project_folder = values["new_project_folder"].replace(" ", "_")
                break
            else:
                project_folder = "Default_project"
                break

        if event == 'Load':
            project_folder = ''
            for key, value in values.items():
                if value == True:
                    project_folder = sorted(projects_list)[key]
            if project_folder == '':
                project_folder = "Default_project"
            break

        if event == 'Refresh':
            if new_projects_main_path == None:
                break
            f = open(user_data_txt, "w")
            f.write(new_projects_main_path)
            f.close()
            sg.Popup("Please reload TaxonTableTools to apply changes", title="Refresh output directory")
            sys.exit()

        if event == 'Exit':
            sys.exit()

    start_window.close()

    ##########################################################################################################################
    # check folders and create new folders if neccessary
    try:
        if not os.path.exists(Path(str(projects_main_path) + "/Projects")):
            os.mkdir(Path(str(projects_main_path) + "/Projects"))
    except:
        sg.PopupError("The output directory does not exist anymore! Please refresh the output folder.")
        sys.exit()

    path_to_outdirs = Path(str(projects_main_path) + "/Projects/" + project_folder)
    if not os.path.exists(path_to_outdirs):
        os.mkdir(path_to_outdirs)

    directories_to_create = ["0_raw_data", "1_Quality_check", "2_Demultiplexing", "3_Paired-end_merging", "4_Primer_trimming", "5_length_filtering", "6_Quality_filtering", "7_Clustering", "8_Read_tables"]

    for directory in directories_to_create:
        dirName = Path(str(path_to_outdirs) + "/" + directory + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

        for sub_directory in ["_data", "_data_discarded" ,"_stats", "_log"]:
            dirName = Path(str(path_to_outdirs) + "/" + directory + "/" + sub_directory)
            if not os.path.exists(dirName):
                os.mkdir(Path(dirName))

    ##########################################################################################################################
    ##########################################################################################################################
    ##########################################################################################################################

    layout_getting_started = [
				        [sg.Text('',size=(1,1))],
                        [sg.Text("Before starting the analysis, some adjustments are required.", font=('Arial', 11, "bold"))],
    					[sg.Text('',size=(1,1))],
                        [sg.Text("1. Raw data format", font=('Arial', 10, "bold"))],
                        [sg.Text("Raw data files are required to be present as .fastq or .fastq.gz files.\nUncompressed files will be compressed.")],
    					[sg.Text('',size=(1,1))],
                        [sg.Text("2. Raw data suffix", font=('Arial', 10, "bold"))],
                        [sg.Text("Raw data files must be present as paired-end files with specific suffixes (e.g. \'_r1\' and \'_r2\').\nThis can any suffix as long as the files only differ in one character.")],
    					[sg.Text('',size=(1,1))],
                        [sg.Text("3. Negative controls.", font=('Arial', 10, "bold"))],
                        [sg.Text("The number of reads in the negative controls can be tracked.\nThe following ID will used to identify the negative controls."), sg.Input("nc_", size=(10,1), key="nc_ID")],
    					[sg.Text('',size=(1,1))],
                        [sg.Text("4. Validate your input before starting the analysis.", font=('Arial', 10, "bold"))],
                        [sg.Text("Verify files:", size=(20,1)), sg.Button("Check")],
    					]

    layout_quality_check = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Quality check', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('_'*100)],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('Enter path to raw data folder')],
                        [sg.Radio('Custom path', 'raw_data_folder'), sg.Input(), sg.FolderBrowse(key="raw_data_folder")],
    					[sg.Radio('0_raw_data', 'raw_data_folder', key="raw_data_folder_default", default=True)],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('FastQC', size=(10,1)), sg.Button("Run", key="run_fastqc"), sg.Text('',size=(2,1)),
                        sg.Frame(layout=[[sg.Text("Check the quality of the raw reads.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('MD5sum', size=(10,1)), sg.Button("Run", key="run_md5sum"), sg.Text('',size=(2,1)),
                        sg.Frame(layout=[[sg.Text("Check the md5 sum of the raw reads.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('Read plot', size=(10,1)), sg.Button("Run", key="run_raw_read_plot"), sg.Text('',size=(2,1)),
                        sg.Frame(layout=[[sg.Text("Plot the number of reads per sample.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					]

    layout_sample_renaming = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Sample renaming', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('_'*100)],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('Enter path to raw data folder')],
                        [sg.Radio('Custom path', 'sample_renaming_folder'), sg.Input(), sg.FolderBrowse(key="sample_renaming_folder")],
    					[sg.Radio('0_raw_data', 'sample_renaming_folder', key="sample_renaming_folder_default", default=True)],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('Create rename sheet', size=(30,1)), sg.Button("Run", key="run_create_rename_sheet"), sg.Text('',size=(2,1)),
                        sg.Frame(layout=[[sg.Text("Check the quality of the raw reads.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('Load rename sheet', size=(30,1)), sg.Input(key="sample_renaming_sheet"), sg.FileBrowse()],
    					[sg.Text('',size=(1,1))],
                        [sg.Text('Rename samples', size=(30,1)), sg.Button("Run", key="run_rename_samples"), sg.Text('',size=(2,1)),
                        sg.Frame(layout=[[sg.Text("Check the quality of the raw reads.")]], title="Information")],
    					]

    layout_tools = [
    					[sg.Text('',size=(1,1))],
                        [sg.Text("MetaProcessor tools", font=('Arial', 11, "bold"))],
    					[sg.Text('',size=(1,1))],
                        # row 1
                        [sg.Button(key="open_demultiplexing", button_color=('white', 'white'), image_filename=demultiplexing_png),
                        sg.Button(key="open_pe_merging", button_color=('white', 'white'), image_filename=pe_merging_png),
                        sg.Button(key="open_primer_trimming", button_color=('white', 'white'), image_filename=primer_trimming_png)],
                        # row 2
                        [sg.Button(key="open_length_filtering", button_color=('white', 'white'), image_filename=length_filtering_png),
                        sg.Button(key="open_quality_filtering", button_color=('white', 'white'), image_filename=quality_filtering_png),
                        sg.Button(key="open_clustering", button_color=('white', 'white'), image_filename=clustering_png)],
                        # row 3
                        [sg.Button(key="open_postprocessing", button_color=('white', 'white'), image_filename=postprocessing_png),
                        sg.Button(key="open_blast", button_color=('white', 'white'), image_filename=blast_png),
                        sg.Button(key="open_quick_analysis", button_color=('white', 'white'), image_filename=all_in_one_analysis_png)],
                        ]

    ## define variables for the main window
    num_cores = psutil.cpu_count(logical = False)
    welcome_text = datetime.datetime.now().strftime("%H:%M:%S") + ": Welcome to the MetaProcessor!\n"

    layout = [  [sg.Text('MetaProcessor'), sg.Text('Project:', font=('Arial', 12, "bold")), sg.Text(project_folder, font=('Arial', 12, "bold"))],
    			[sg.Text('',size=(1,1))],
    			[sg.TabGroup([[
                sg.Tab("Getting started", layout_getting_started),
                sg.Tab('Quality check', layout_quality_check),
                sg.Tab('Sample renaming', layout_sample_renaming),
                sg.Tab("Metabarcoding processing tools", layout_tools)
                ]]),
                sg.Multiline(welcome_text, size = (40, 32), key = '_OUTSTREAM_', autoscroll = True)], ## size = (w, h)
    			[sg.Text('',size=(1,1))],
    			[sg.Text('',size=(1,1))],
    			[sg.Exit(button_color=('black', 'red')), sg.Text("", size=(6,1)), sg.Button("Check installations", key="check_installations"), sg.Text("", size=(6,1)),
                sg.Text('Number of cores:'), sg.Slider(range=(1, num_cores), default_value=num_cores - 1, size=(15, 10), orientation='horizontal', key="num_cores_to_use"), sg.Text("", size=(6,1)), sg.Button(image_filename=github_png, image_size=(26,26), key='open_github', button_color=('black', 'white')),
                sg.Button(key='open_twitter', button_color=('white', 'white'), image_filename=twitter_png, image_size=(26,26)), sg.Text('', size=(1,1)), sg.Text(meta_tools_version, font=('Arial', 8))]]

    # Create the Window
    window = sg.Window('MetaProcessor', layout)
    win2_active=False

    ##########################################################################################################################

    while True:
        try:
            event, values = window.Read()

            if event is None or event == 'Exit':
                break

            i = 0
            # define variables
            ## TODO: DELETE AND DIRECTLY READ FROM VALUES
            raw_data_folder = values['raw_data_folder']
            raw_data_folder_default = values['raw_data_folder_default']
            sample_renaming_folder = values['sample_renaming_folder']
            sample_renaming_folder_default = values['sample_renaming_folder_default']
            sample_renaming_sheet = values['sample_renaming_sheet']
            num_cores_to_use = int(values["num_cores_to_use"])
            nc_ID = values["nc_ID"]

            # call events

            if event == 'run_fastqc':
                if raw_data_folder_default == True:
                    raw_data_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")

                if (raw_data_folder == ''):
                    sg.PopupError("Please provide a folder", keep_on_top=True)
                    print("Error: Please provide a folder")

                else:
                    from metaprocessor.fastqc_quality_check import fastqc_quality_check
                    fastqc_quality_check(raw_data_folder, path_to_outdirs, num_cores_to_use)

            if event == 'run_md5sum':
                if raw_data_folder_default == True:
                    raw_data_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")

                if (raw_data_folder == ''):
                    sg.PopupError("Please provide a folder", keep_on_top=True)
                    print("Error: Please provide a folder")

                else:
                    from metaprocessor.md5sum import md5sum
                    md5sum(raw_data_folder)

            if event == 'run_raw_read_plot':
                if raw_data_folder_default == True:
                    raw_data_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")

                if (raw_data_folder == ''):
                    sg.PopupError("Please provide a folder", keep_on_top=True)
                    print("Error: Please provide a folder")

                else:
                    from metaprocessor.raw_read_plot import raw_read_plot
                    raw_read_plot(raw_data_folder, path_to_outdirs, num_cores_to_use)

            if event == 'run_create_rename_sheet':
                if sample_renaming_folder_default == True:
                    sample_renaming_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")

                if (sample_renaming_folder == ''):
                    sg.PopupError("Please provide a folder", keep_on_top=True)
                    print("Error: Please provide a folder")

                else:
                    from metaprocessor.create_rename_sheet import create_rename_sheet
                    create_rename_sheet(sample_renaming_folder, path_to_outdirs)

            if event == 'run_rename_samples':

                if sample_renaming_folder_default == True:
                    sample_renaming_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")

                if (sample_renaming_folder == ''):
                    sg.PopupError("Please provide a folder", keep_on_top=True)
                    print("Error: Please provide a folder")

                security = sg.PopupOKCancel("Warning: You are about to rename all files in the directory." + "\n" + "Always create a backup of your files!", title="Warning")

                if security == 'Cancel':
                    break

                else:
                    from metaprocessor.rename_samples import rename_samples
                    rename_samples(sample_renaming_folder, sample_renaming_sheet, path_to_outdirs)

            if event == "run_create_demultiplexing_sheet":
                print("Create new demultiplexing sheet")
                from metaprocessor.create_demultiplexing_sheet import create_demultiplexing_sheet
                create_demultiplexing_sheet(path_to_outdirs)

            if event == 'run_demultiplexing':
                print("Demultiplexing of raw reads")

                if (r1_reads == ''):
                    sg.PopupError("Please provide the r1 reads file", keep_on_top=True)
                    print("Error: Please provide the r1 reads file")
                elif(r2_reads == ''):
                    sg.PopupError("Please provide the r2 reads file", keep_on_top=True)
                    print("Error: Please provide the r2 reads file")
                elif(demultiplexing_sheet == ''):
                    sg.PopupError("Please provide the demultiplexing sheet file", keep_on_top=True)
                    print("Error: Please provide the demultiplexing sheet file")
                else:
                    from metaprocessor.demultiplexing import demultiplexing
                    demultiplexing(r1_reads, r2_reads, demultiplexing_sheet, path_to_outdirs)

            if event == 'check_installations':
                print("Checking installations")
                from metaprocessor.check_installations import check_installations
                check_installations()

            if event == 'open_github':
                print("Open: https://github.com/thm93/TaxonTableTools")
                webbrowser.open('https://github.com/thm93/TaxonTableTools')

            if event == 'open_twitter':
                print("Open: https://twitter.com/THM_93")
                webbrowser.open('https://twitter.com/THM_93')

            if event == 'Check':
                from metaprocessor.check_raw_files import check_raw_files
                check_raw_files(nc_ID, path_to_outdirs, window["_OUTSTREAM_"], window)

###########################################################
###########################################################
            if event == 'open_demultiplexing':
                ## defines the layout of demultiplexer
                ## define the raw data folder path
                raw_data_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")
                layout_demultiplexing = [
                                [sg.Text('',size=(1,1))],
                                [sg.Text('Paired-end merging', size=(50,2), font=('Arial', 11, "bold"))],
                                [sg.Text('_'*115)],
                                [sg.Text('',size=(1,1))],

                                [sg.Text("1. Choose the path to your raw data.", font=('Arial', 10, "bold"))],
                                [sg.Text("Select all files to demultiplex:", size=(30, 1)), sg.InputText(size = (20, 1), key = '_FILES_'), sg.FilesBrowse(initial_folder = raw_data_folder), sg.Button('Load files')],
                                [sg.Text('',size=(1,1))],

                                [sg.Text("2. Create or load primer set:", font=('Arial', 10, "bold"))],
                                [sg.Text('a) Number of primers per direction:', size=(30, 1)), sg.Spin([i for i in range(1, 21)], size = (9, 1), key = '_NUM_PRIMERS'), sg.Button('Create new primer set')],
                                [sg.Text('b) Select existing primerset:', size=(30, 1)), sg.InputText(size=(10, 1), key='_PRIMERSET_'), sg.FileBrowse(initial_folder = primersets_path, file_types = (("Text Files", "*.txt"),))],
                                [sg.Text('',size=(1,1))],

                                [sg.Text("3. Create or load tagging scheme:", font=('Arial', 10, "bold"))],
                                [sg.Text('a) Number of primer combinations used:', size=(30, 1)), sg.Spin([i for i in range(1, 21)], size = (9, 1), key = '_USED_COMBS_'), sg.Button('Create new tagging scheme')],
                                [sg.Text('b) Select existing tagging scheme:', size=(30, 1)), sg.InputText(size=(10,1), key = '_TAGGING_SCHEME_'), sg.FileBrowse(initial_folder = tagging_templates_path, file_types = (("Worksheets", "*.xlsx"),)), sg.Button('Modify selected scheme')],
                                [sg.Text('',size=(1,1))],

                                [sg.Text("4. Raw data demultiplexing.", font=('Arial', 10, "bold"))],
                                [sg.Text("Demultipex selected files:", size=(30, 1)), sg.Button('Run', key='Start demultiplexing'), sg.Text(""), sg.Checkbox('Remove tags from sequence', key = '_TAG_REMOVAL_')],
                                [sg.Text('',size=(1,1))],

                                [sg.Button('Return', button_color=('black', 'red'))]
                                ]

                ## create the demultiplexing window
                tools_window = sg.Window('Demultiplexing', layout_demultiplexing, keep_on_top=False)
                primer_window_active = False
                scheme_gen_active = False

                ## import modules
                from metaprocessor import file_pairs
                from metaprocessor import save_primerset
                from metaprocessor.scheme_generator import scheme_layout
                from metaprocessor.scheme_generator import save_scheme
                from metaprocessor import demultiplexing

                ## define path to output files
                _OUTPUT_PATH_ = Path(str(path_to_outdirs) + "/2_Demultiplexing/_data/")

                ### main loop
                while True:
                    event, values = tools_window.read(timeout = 100)

                    ## code to close the program
                    if event == None or event == 'Return':
                        break

                    ## code to check the selected files
                    if event == 'Load files':
                        if not values['_FILES_']:
                            window['_OUTSTREAM_'].print('{}: Please select files first.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                        else:
                            files = values['_FILES_'].split(';')
                            pairs = [i for i in file_pairs.main(files) if len(i) == 2]
                            window['_OUTSTREAM_'].print('{}: {} files were loaded.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(files)))
                            window['_OUTSTREAM_'].print('{}: {} valid file pairs were found.'.format(datetime.datetime.now().strftime("%H:%M:%S"), len(pairs)))

                    ## opens the primerset generator, second part will prevent to windows open at the same time
                    if event == 'Create primer set' and not primer_window_active:
                        ## activate the primer generator window
                        primer_window_active = True

                        ## define the primer generator layout
                        fwd_col = [[sg.Input(size = (15, 1)), sg.Input(size = (10, 1))] for i in range(int(values['_NUM_PRIMERS']))]
                        rev_col = [[sg.Input(size = (15, 1)), sg.Input(size = (10, 1))] for i in range(int(values['_NUM_PRIMERS']))]

                        ## no better way to adjust the padding around the headers
                        ## layout definition for the primerset generator
                        primerset_layout = [[sg.Text('Name of primerset:'), sg.Input(size = (20, 1), key = '_PRIMERSETNAME_')],
                                           [sg.Frame(layout = [
                                           [sg.Text(' Name'), sg.Text('                 Sequence')],
                                           [sg.Column(fwd_col)]],
                                           title = 'Forward primers'),
                                           sg.Frame(layout = [
                                           [sg.Text(' Name'), sg.Text('                 Sequence')],
                                           [sg.Column(rev_col)]],
                                           title = 'Reverse primers')],
                                           [sg.Button('Close'), sg.Button('Save')]
                                           ]

                        primer_window = sg.Window('Primer generator', primerset_layout)

                    ## read the window only if active, handle all input to primer generator window
                    if primer_window_active:
                        primer_ev, primer_vals = primer_window.read(timeout = 100)

                        if primer_ev == None or primer_ev == 'Close':
                            primer_window_active = False
                            primer_window.close()

                        ## save the primerset in the modules data path
                        if primer_ev == 'Save':
                            ## code to save the primerset, only works if a name is in the primersets name Input Box
                            if primer_vals['_PRIMERSETNAME_'] != '':
                                save_primerset.save_primerset(primer_vals, primersets_path, window['_OUTSTREAM_'])
                            else:
                                window['_OUTSTREAM_'].print('{}: Please choose a valid name for your primerset.'.format(datetime.datetime.now().strftime("%H:%M:%S")))


                    ## opens the tagging scheme generator
                    ## all tagging scheme generator code can be found here
                    if event == 'Create tagging scheme' and values['_PRIMERSET_'] != '' and 'pairs' in locals():
                        ## activate the scheme generator window
                        scheme_gen_active = True

                        ## generate the layout at runtime
                        scheme_gen_layout = scheme_layout(values['_USED_COMBS_'], values['_PRIMERSET_'])

                        scheme_gen_window = sg.Window('Tagging scheme generator', scheme_gen_layout)

                    ## read the window
                    if scheme_gen_active:
                        scheme_ev, scheme_vals = scheme_gen_window.read(timeout = 100)

                        ## close the scheme gen window if X or Close is hit
                        if scheme_ev == None or scheme_ev == 'Close':
                            scheme_gen_active = False
                            scheme_gen_window.close()

                        if scheme_ev == 'Save':
                            if scheme_vals['_TAGGING_SCHEME_NAME_'] == '':
                                window['_OUTSTREAM_'].print('{}: Please choose a valid name for your tagging scheme.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                            else:
                                save_scheme(pairs, scheme_vals, tagging_templates_path, window['_OUTSTREAM_'])

                    ## check if files are selected and a primerset is selected
                    if event == 'Create tagging scheme':
                        if not 'pairs' in locals():
                            window['_OUTSTREAM_'].print('{}: Please load files to demultiplex.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                        if values['_PRIMERSET_'] == '':
                            window['_OUTSTREAM_'].print('{}: Please select a primerset.'.format(datetime.datetime.now().strftime("%H:%M:%S")))

                    ## code to modify the scheme
                    if event == 'Modify selected scheme':
                        if values['_TAGGING_SCHEME_'] == '':
                            window['_OUTSTREAM_'].print('{}: Please select a tagging scheme.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                        else:
                            if sys.platform == "win32":
                                os.startfile(values['_TAGGING_SCHEME_'])
                            else:
                                opener = "open" if sys.platform == 'darwin' else 'xdg-open'
                                subprocess.call([opener, values['_TAGGING_SCHEME_']])

                    if event == 'Start demultiplexing':
                        if not 'pairs' in locals():
                            window['_OUTSTREAM_'].print('{}: Please load files to demultiplex.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                        if values['_PRIMERSET_'] == '':
                            window['_OUTSTREAM_'].print('{}: Please select a primerset.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                        if values['_TAGGING_SCHEME_'] == '':
                            window['_OUTSTREAM_'].print('{}: Please select a tagging scheme.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
                        else:
                            from demultiplexer import demultiplexing
                            demultiplexing.main(values['_PRIMERSET_'], values['_TAGGING_SCHEME_'], _OUTPUT_PATH_, values['_TAG_REMOVAL_'], window['_OUTSTREAM_'], tools_window)

                tools_window.close()

###########################################################
            if event == 'open_pe_merging':

                layout_paired_end_merging = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('Paired-end merging', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*115)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("1. Choose the path to your raw data.", font=('Arial', 10, "bold"))],
                                    [sg.Radio('0_raw_data', 'paired_end_merging_folder', key="paired_end_merging_folder_0", default=True),
                                    sg.Radio('2_Demultiplexing', 'paired_end_merging_folder', key="paired_end_merging_folder_2"),
                                    sg.Radio('Custom path', 'paired_end_merging_folder'), sg.Input(size=(10,1)), sg.FolderBrowse(key="paired_end_merging_folder")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("2. Adjust merging option to your need.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Maxdiffpct:',size=(13,1)), sg.Input('25', size=(5,1), key="maxdiffpct"), sg.Text('Maxdiffs:',size=(8,1)), sg.Input('199', size=(5,1), key="maxdiffs"), sg.Text('Minovlen:',size=(8,1)), sg.Input('5', size=(5,1), key="minovlen"),],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("3. Start the paired-end merging.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Paired-end merging', size=(30,1)), sg.Button("Run", key="run_paired_end_merging"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("Merge paired-end reads using Vsearch")]], title="Information")],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('Paired-end merging', layout_paired_end_merging, keep_on_top=False)
                while (True):
                    ######################################
                    event, values = tools_window.Read()
                    paired_end_merging_folder = values['paired_end_merging_folder']
                    paired_end_merging_folder_0 = values['paired_end_merging_folder_0']
                    paired_end_merging_folder_2 = values['paired_end_merging_folder_2']
                    ######################################
                    if event == 'run_paired_end_merging':
                        if paired_end_merging_folder_0 == True:
                            paired_end_merging_folder = Path(str(path_to_outdirs) + "/0_raw_data/_data/")
                        elif paired_end_merging_folder_2 == True:
                            paired_end_merging_folder = Path(str(path_to_outdirs) + "/2_Demultiplexing/_data/")
                        if (paired_end_merging_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")
                        else:
                            from metaprocessor.v_merge import v_merge
                            v_merge(paired_end_merging_folder, path_to_outdirs, num_cores_to_use, values["minovlen"], values["maxdiffs"], values["maxdiffpct"], window["_OUTSTREAM_"], window)
                    ######################################
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_primer_trimming':

                layout_primer_trimming = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('Primer trimming', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*115)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("1. Choose the path to your data.", font=('Arial', 10, "bold"))],
                                    [sg.Radio('3_Paired-end_merging', 'cutadapt_folder',key="cutadapt_folder_default", default=True),
                                    sg.Radio('Custom path', 'cutadapt_folder'), sg.Input(size=(10,1)), sg.FolderBrowse(key="cutadapt_folder")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("2. Add or remove primers in your primer library.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Primer sheet', size=(30,1)), sg.Button('Modify primer sheet', key="modify_primer_sheet")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("3. Remove the primer sequences from your reads.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Trim primers', size=(30,1)), sg.Button("Run", key="run_cutadapt_primer_trimming"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("Trim primers using Cutadapt.")]], title="Information"), sg.Text('', size=(2,1)),
                                    sg.Text('Anchoring:'), sg.Radio('True', 'anchoring', key='anchored_true'), sg.Radio('False', 'anchoring', key='anchored_false', default=True)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('Primer trimming', layout_primer_trimming, keep_on_top=False)
                while (True):
                    ######################################
                    event, values = tools_window.Read()
                    cutadapt_folder = values['cutadapt_folder']
                    cutadapt_folder_default = values["cutadapt_folder_default"]
                    anchored_true = values["anchored_true"]
                    anchored_false = values["anchored_false"]
                    ######################################
                    if event == 'modify_primer_sheet':
                        if sys.platform == "win32":
                            os.startfile(primer_sheet_xlsx)
                        else:
                            opener = "open" if sys.platform == 'darwin' else 'xdg-open'
                            subprocess.call([opener, primer_sheet_xlsx])

                    if event == 'run_cutadapt_primer_trimming':
                        if cutadapt_folder_default == True:
                            cutadapt_folder = Path(str(path_to_outdirs) + "/3_Paired-end_merging/_data/")

                        if (cutadapt_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")

                        else:
                            if anchored_true == True:
                                anchoring = True
                            else:
                                anchoring = False

                            from metaprocessor.load_primers import load_primers
                            primers = load_primers(primer_sheet_xlsx)

                            if primers != "Back":
                                p5_forward_primer = primers[0]
                                p7_reverse_primer = primers[1]

                                from metaprocessor.cutadapt_trimming import cutadapt_trimming
                                cutadapt_trimming(cutadapt_folder, path_to_outdirs, p5_forward_primer, p7_reverse_primer, anchoring, num_cores_to_use, window["_OUTSTREAM_"], window)
                    ######################################
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_length_filtering':

                layout_length_filtering = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('Length filtering', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*115)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("1. Choose the path to your data.", font=('Arial', 10, "bold"))],
                                    [sg.Radio('4_Primer_trimming', 'length_filtering_folder', key="length_filtering_folder_default", default=True),
                                    sg.Radio('Custom path', 'length_filtering_folder'), sg.Input(size=(10,1)), sg.FolderBrowse(key="length_filtering_folder")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("2. Calculate the read lenghts of your data.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Get length distribution', size=(30,1)), sg.Button("Run", key="run_get_length_distribution"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("Get length distribution plot")]], title="Information")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("3. Enter the length Thresholds to accordingly filter your data.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Minimum read length:', size=(30,1)), sg.Input(size=(10,1), key="min_length"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("The lower threshold to cut.")]], title="Information")],
                                    [sg.Text('Maximum read length:', size=(30,1)), sg.Input(size=(10,1), key="max_length"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("The upper threshold to cut.")]], title="Information")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("4. Remove reads that do not meet the length requirements.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Trim lengths', size=(30,1)), sg.Button("Run", key="run_length_filtering"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("Filter reads using Vsearch.")]], title="Information")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('Lenght filtering', layout_length_filtering, keep_on_top=False)
                while (True):
                    ######################################
                    event, values = tools_window.Read()
                    length_filtering_folder = values['length_filtering_folder']
                    length_filtering_folder_default = values["length_filtering_folder_default"]
                    min_length = values['min_length']
                    max_length = values['max_length']
                    ######################################
                    if event == 'run_length_filtering':
                        if length_filtering_folder_default == True:
                            length_filtering_folder = Path(str(path_to_outdirs) + "/4_Primer_trimming/_data/")
                        if (length_filtering_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")
                        else:
                            from metaprocessor.v_length_filtering import v_length_filtering
                            v_length_filtering(length_filtering_folder, path_to_outdirs, min_length, max_length, num_cores_to_use, window["_OUTSTREAM_"], window)
                    ######################################
                    if event == 'run_get_length_distribution':
                        if length_filtering_folder_default == True:
                            length_filtering_folder = Path(str(path_to_outdirs) + "/4_Primer_trimming/_data/")
                        if (length_filtering_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")
                        else:
                            from metaprocessor.get_length_distribution import get_length_distribution
                            get_length_distribution(length_filtering_folder, path_to_outdirs, num_cores_to_use)
                    ######################################
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_quality_filtering':

                layout_quality_filtering = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('Quality filtering', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*115)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("1. Choose the path to your data.", font=('Arial', 10, "bold"))],
                                    [sg.Radio('5_length_filtering', 'maxee_filtering_folder', key="maxee_filtering_folder_default", default=True),
                                    sg.Radio('Custom path', 'maxee_filtering_folder'), sg.Input(size=(10,1)), sg.FolderBrowse(key="maxee_filtering_folder")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("2. Enter the maximum expected error threshold to accordingly filter your data.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Maxee Threshold:', size=(30,1)), sg.Input("1", size=(10,1), key="maxee_value"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("The maximum expected error threshold to cut.")]], title="Information")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("3. Remove reads that do not meet the quality expectations.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Quality fitering', size=(30,1)), sg.Button("Run", key="run_maxee_filtering"), sg.Text('',size=(2,1)),
                                    sg.Frame(layout=[[sg.Text("Filter reads using Vsearch.")]], title="Information")],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('Quality filtering', layout_quality_filtering, keep_on_top=False)
                while (True):
                    ######################################
                    event, values = tools_window.Read()
                    maxee_value = values['maxee_value']
                    maxee_filtering_folder = values['maxee_filtering_folder']
                    maxee_filtering_folder_default = values["maxee_filtering_folder_default"]
                    ######################################
                    if event == 'run_maxee_filtering':
                        if maxee_filtering_folder_default == True:
                            maxee_filtering_folder = Path(str(path_to_outdirs) + "/5_length_filtering/_data/")
                        if (maxee_filtering_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")
                        else:
                            from metaprocessor.v_maxee_filtering import v_maxee_filtering
                            v_maxee_filtering(maxee_filtering_folder, path_to_outdirs, maxee_value, num_cores_to_use, window["_OUTSTREAM_"], window)
                    ######################################
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_clustering':

                layout_clustering = [
        					[sg.Text('',size=(1,1))],
        					[sg.Text('Clustering', size=(50,2), font=('Arial', 11, "bold"))],
        					[sg.Text('_'*115)],
        					[sg.Text('',size=(1,1))],
                            [sg.Text("1. Choose the path to your data.", font=('Arial', 10, "bold"))],
                            [sg.Radio('6_Quality_filtering', 'clustering_folder', key="clustering_folder_default", default=True),
                            sg.Radio('Custom path', 'clustering_folder'), sg.Input(size=(10,1)), sg.FolderBrowse(key="clustering_folder")],
        					[sg.Text('',size=(1,1))],
                            [sg.Text("2. Before the clustering, the data needs to be dereplicated and both singletons and chimeras removed.", font=('Arial', 10, "bold"))],
                            [sg.Text('Process files:', size=(30,1)), sg.Button("Run", key="run_v_derep_singeltons_uchime"), sg.Text('',size=(2,1)),
                            sg.Frame(layout=[[sg.Text("Dereplication, singleton removal, chimera removal.")]], title="Information"),
                            sg.Text("Large files: "), sg.Radio('True', 'large_file_option', key='large_file_option_true'), sg.Radio('False', 'large_file_option', key='large_file_option_false', default=True)],
        					[sg.Text('',size=(1,1))],
        					[sg.Text('_'*115)],
        					[sg.Text('',size=(1,1))],
                            [sg.Text("3. Choose your clustering algorithm.", font=('Arial', 10, "bold"))],
                            [sg.Text("3.1 SWARM clustering is based on the small local linking threshold (d) clustering.", font=('Arial', 10, "bold"))],
                            [sg.Text('Swarm clustering', size=(30,1)), sg.Button("Run", key="run_swarm_clustering"), sg.Text("", size=(2,1)),
                            sg.Text("d min:"), sg.Input("1", size=(3,1), key="d_min"), sg.Text("d max:"), sg.Input("5", size=(3,1), key="d_max")],
        					[sg.Text('',size=(1,1))],
                            [sg.Text("3.2 Vsearch clustering is based on relative distances clustering.", font=('Arial', 10, "bold"))],
                            [sg.Text('Vsearch clustering', size=(30,1)), sg.Button("Run", key="run_v_clustering"), sg.Text("", size=(2,1)), sg.Text('Representative:'), sg.Combo(['Centroid', 'Consensus'], default_value='Centroid', key='clustering_representative'), sg.Text("", size=(2,1)),
                            sg.Text("Threshold:"), sg.Input("0.97", size=(5,1), key="cluster_threshold")],
        					[sg.Text('',size=(1,1))],
                            [sg.Text("3.3. Vsearch unoise will denoise your data into exact sequence variants (ESVs).", font=('Arial', 10, "bold"))],
                            [sg.Text('Vsearch unoise', size=(30,1)), sg.Button("Run", key="run_v_unoise"), sg.Text("", size=(2,1)), sg.Text('Representative:'), sg.Combo(['Centroid', 'Consensus'], default_value='Centroid', key='denoise_representative'),
                            sg.Text("Minimum size:"), sg.Input("8", size=(5,1), key="unoise_minsize")],
        					[sg.Text('',size=(1,1))],
                            [sg.Button('Return', button_color=('black', 'red'))]
        					]
                # create the demultiplexing window
                tools_window = sg.Window('Clustering', layout_clustering, keep_on_top=False)
                while (True):
                    ######################################
                    event, values = tools_window.Read()
                    clustering_folder = values['clustering_folder']
                    clustering_folder_default = values["clustering_folder_default"]
                    d_min = values['d_min']
                    d_max = values['d_max']
                    large_file_option_true = values["large_file_option_true"]
                    large_file_option_false = values["large_file_option_false"]
                    ######################################
                    if event == 'run_v_derep_singeltons_uchime':
                        if clustering_folder_default == True:
                            clustering_folder = Path(str(path_to_outdirs) + "/6_Quality_filtering/_data/")
                        if (clustering_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")
                        else:
                            if large_file_option_true == True:
                                large_file_option = True
                                answer = sg.PopupOKCancel("Large file option activated. This is not recommended unless the processing takes too long!", title="Warning")
                                if answer != "Cancel":
                                    from metaprocessor.v_derep_singletons_uchime import v_derep_singletons_uchime
                                    v_derep_singletons_uchime(clustering_folder, path_to_outdirs, num_cores_to_use, large_file_option, window["_OUTSTREAM_"], window)
                            else:
                                large_file_option = False
                                from metaprocessor.v_derep_singletons_uchime import v_derep_singletons_uchime
                                v_derep_singletons_uchime(clustering_folder, path_to_outdirs, num_cores_to_use, large_file_option, window["_OUTSTREAM_"], window)
                    ######################################
                    if event == 'run_v_clustering':
                        if clustering_folder_default == True:
                            clustering_folder = Path(str(path_to_outdirs) + "/4_Primer_trimming/_data/")

                        if (clustering_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")

                        else:
                            from metaprocessor.v_clustering import v_clustering
                            v_clustering(clustering_folder, values['cluster_threshold'], path_to_outdirs, values["clustering_representative"], window["_OUTSTREAM_"], window)
                    ######################################
                    if event == 'run_v_unoise':
                        if clustering_folder_default == True:
                            clustering_folder = Path(str(path_to_outdirs) + "/4_Primer_trimming/_data/")

                        if (clustering_folder == ''):
                            sg.PopupError("Please provide a folder", keep_on_top=True)
                            print("Error: Please provide a folder")

                        else:
                            from metaprocessor.v_unoise import v_unoise
                            v_unoise(clustering_folder, path_to_outdirs, values["denoise_representative"], values["unoise_minsize"], window["_OUTSTREAM_"], window)
                    ######################################
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_postprocessing':

                layout_postprocessing = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('Postprocessing', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*115)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("1. Filter your read table to exclude low-abundant reads.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Filter OTU tables', size=(30,1)), sg.Button("Run", key="run_otu_table_filtering"), sg.Text("", size=(2,1)), sg.Text("Threshold:"), sg.Input("0.01", size=(5,1), key="filter_threshold"), sg.Text("%")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("Calculate a heatmap to display the read distribution across your data.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Create heatmap', size=(30,1)), sg.Button("Run", key="run_otu_table_heatmap")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("Compare the number of OTUs for your different SWARM d values.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Compare SWARM clusters', size=(30,1)), sg.Button("Run", key="run_compare_swarm_clusters")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('Postprocessing', layout_postprocessing, keep_on_top=False)
                while (True):
                    ######################################
                    event, values = tools_window.Read()
                    filter_threshold = values['filter_threshold']
                    ######################################
                    if event == 'run_otu_table_filtering':
                        from metaprocessor.read_table_filtering import read_table_filtering
                        read_table_filtering(filter_threshold, path_to_outdirs, window["_OUTSTREAM_"], window)
                    if event == 'run_otu_table_heatmap':
                        from metaprocessor.read_table_heatmap import read_table_heatmap
                        read_table_heatmap(path_to_outdirs, window["_OUTSTREAM_"], window)
                    ######################################
                    if event == 'run_compare_swarm_clusters':
                        from metaprocessor.compare_swarm_clusters import compare_swarm_clusters
                        compare_swarm_clusters(path_to_outdirs, window["_OUTSTREAM_"], window)
                    ######################################
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_blast':

                layout_blast = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('BLAST', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*115)],
                                    [sg.Text("1. Enter the path to your BLAST results file.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Blast results:', size=(30,1)), sg.Input(), sg.FileBrowse()],
                                    [sg.Text("2. Enter the path to your clustering fasta file.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Fasta file:', size=(30,1)), sg.Input(), sg.FileBrowse()],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text("3. Process and filter your BLAST results to create a taxonomy table.", font=('Arial', 10, "bold"))],
                                    [sg.Text('Process BLAST hits', size=(30,1)), sg.Button("Run")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('BLAST', layout_blast, keep_on_top=True)
                while (True):
                    event, values = tools_window.Read()
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################
            if event == 'open_quick_analysis':

                layout_quick_analysis = [
                					[sg.Text('',size=(1,1))],
                					[sg.Text('Quick analysis', size=(50,2), font=('Arial', 11, "bold"))],
                					[sg.Text('_'*100)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Raw data folder:'), sg.Radio('0_raw_data', 'quick_analysis_folder', key="quick_analysis_folder_2", default=True), sg.Radio('2_Demultiplexing', 'quick_analysis_folder', key="quick_analysis_folder_3"), sg.Radio('Custom path', 'quick_analysis_folder'), sg.Input(), sg.FolderBrowse(key="quick_analysis_folder_1"),],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Paired-end merging', size=(20,1), font=('Arial', 10, "bold"))],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Primer trimming', size=(20,1), font=('Arial', 10, "bold")), sg.Button('Choose primers')],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Length filtering', size=(20,1), font=('Arial', 10, "bold")), sg.Text('Min. read length:', size=(15,1)), sg.Input(size=(10,1), key="qa_min_length"), sg.Text('Max. read length:', size=(15,1)), sg.Input(size=(10,1), key="qa_max_length"),],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Quality filtering', size=(20,1), font=('Arial', 10, "bold")), sg.Text('Maxee threshold:', size=(15,1)), sg.Input(size=(10,1), key="qa_maxee")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Clustering', size=(20,1), font=('Arial', 10, "bold")),
                                    sg.Radio('Vsearch', 'qa_clustering', key="qa_vsearch", default=True), sg.Input("0.97", size=(5,1), key="qa_cluster_threshold"), sg.Text('%'),
                                    sg.Text('',size=(1,1)),
                                    sg.Radio('Swarm', 'qa_clustering', key="qa_swarm"), sg.Input("1", size=(5,1), key="qa_swarm_dmin"), sg.Input("4", size=(5,1), key="qa_swarm_dmax"), sg.Text('d'),
                                    sg.Text('',size=(1,1)),
                                    sg.Radio('Unoise', 'qa_clustering', key="qa_unoise")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Read table filtering', size=(20,1), font=('Arial', 10, "bold")), sg.Text("Threshold:"), sg.Input("0.01", size=(5,1), key="qa_filter_threshold"), sg.Text("%")],
                					[sg.Text('',size=(1,1))],
                					[sg.Text('_'*115)],
                					[sg.Text('',size=(1,1))],
                                    [sg.Text('Start quick workflow', size=(20,1), font=('Arial', 10, "bold")), sg.Button("Run")],
                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Return', button_color=('black', 'red'))]
                					]
                # create the demultiplexing window
                tools_window = sg.Window('Quick analysis', layout_quick_analysis, keep_on_top=False)
                while (True):
                    event, values = tools_window.Read()
                    if event in ('Return', None):
                        break
                tools_window.Close()

###########################################################

        # define exceptions
        # if there are unexpected errors print a message and continue the script!
        except:
            ## unhide the main window if neccessary
            try:
                window.UnHide()
            except:
                pass
            ## close the secondary window if neccessary
            try:
                tools_window.Close()
            except:
                pass

            ## create layout for the error message
            layout = [
                        [sg.Image(crash_png), sg.Text(" You've been awarded with the gold medal in program crashing!")],
                        [sg.Text("", size=(1,1))],
                        [sg.Text("Unexpected error: " + str(sys.exc_info()[0]))],
                        [sg.Text("", size=(1,1))],
                        [sg.Text("An unexpected error occured:")],
                        [sg.Text("Please refer to the manual.")],
                        [sg.Text("", size=(1,1))],
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

    window.Close()

## run only if called as toplevel script
if __name__ == "__main__":
    main()
