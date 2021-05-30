import glob, sys, os, pkgutil, ast, multiprocessing, pkg_resources, datetime, psutil, subprocess
import PySimpleGUI as sg
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import webbrowser

##########################################################################################################################
# update version here (will be displayed on the main layout)
# Support for: u = ubuntu, w = windows, m = macintosh
meta_tools_version = "Version 0.4"

##########################################################################################################################
# general functions

# slice function for lists to split up lists
def slices(list, slice):
    for i in range(0, len(list), slice):
        yield list[i : i + slice]

def open_file(file):
    if sys.platform == "win32":
        os.startfile(file)
    else:
        opener = "open" if sys.platform == 'darwin' else 'xdg-open'
        subprocess.call([opener, file])

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
    quality_control_png = Path(pkg_resources.resource_filename(__name__, '/_source/quality_control.png'))
    sample_renaming_png = Path(pkg_resources.resource_filename(__name__, '/_source/sample_renaming.png'))
    demultiplexing_png = Path(pkg_resources.resource_filename(__name__, '/_source/demultiplexing.png'))
    all_in_one_analysis_png = Path(pkg_resources.resource_filename(__name__, '/_source/all_in_one_analysis.png'))
    postprocessing_png = Path(pkg_resources.resource_filename(__name__, '/_source/postprocessing.png'))
    ncbi_blast_png = Path(pkg_resources.resource_filename(__name__, '/_source/ncbi_blast.png'))
    boldigger_png = Path(pkg_resources.resource_filename(__name__, '/_source/boldigger.png'))
    analysis_statistics_png = Path(pkg_resources.resource_filename(__name__, '/_source/analysis_statistics.png'))
    log_file_png = Path(pkg_resources.resource_filename(__name__, '/_source/log_file.png'))

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

    directories_to_create = ["1_raw_data"]

    for directory in directories_to_create:
        dirName = Path(str(path_to_outdirs) + "/" + directory + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

    ## load the log file
    main_log_file = Path(str(path_to_outdirs) + "/log.xlsx")

    ## load the available primers
    df = pd.read_excel(Path(primer_sheet_xlsx))
    available_fw_primers_dict = {}
    available_rv_primers_dict = {}
    for primer in df.values.tolist():
        direction = primer[0]
        ID = primer[1]
        sequence = primer[2]

        if direction == 'forward':
            available_fw_primers_dict[ID] = sequence
        elif direction == 'reverse':
            available_rv_primers_dict[ID] = sequence

    available_fw_primers_list = list(available_fw_primers_dict.keys())
    available_rv_primers_list = list(available_rv_primers_dict.keys())

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

    layout_tools = [
    					[sg.Text('',size=(1,1))],
                        [sg.Text("MetaProcessor tools", font=('Arial', 11, "bold"))],
    					[sg.Text('',size=(1,1))],
                        # row 1
                        [sg.Button(key="open_quality_control", button_color=('white', 'white'), image_filename=quality_control_png),
                        sg.Button(key="open_sample_renaming", button_color=('white', 'white'), image_filename=sample_renaming_png),
                        sg.Button(key="open_demultiplexing", button_color=('white', 'white'), image_filename=demultiplexing_png)],
                        # row 2
                        [sg.Button(key="open_run_analyses", button_color=('white', 'white'), image_filename=all_in_one_analysis_png),
                        sg.Button(key="open_postprocessing", button_color=('white', 'white'), image_filename=postprocessing_png),
                        sg.Button(key="open_ncbi_blast", button_color=('white', 'white'), image_filename=ncbi_blast_png)],
                        # row 3
                        [sg.Button(key="open_boldigger", button_color=('white', 'white'), image_filename=boldigger_png),
                        sg.Button(key="open_analysis_statistics", button_color=('white', 'white'), image_filename=analysis_statistics_png),
                        sg.Button(key="open_log_file", button_color=('white', 'white'), image_filename=log_file_png)],
                        ]

    ## define variables for the main window
    num_cores = psutil.cpu_count(logical = False)
    welcome_text = datetime.datetime.now().strftime("%H:%M:%S") + ": Welcome to the MetaProcessor!\n"

    layout = [  [sg.Text('MetaProcessor'), sg.Text('Project:', font=('Arial', 12, "bold")), sg.Text(project_folder, font=('Arial', 12, "bold"))],
    			[sg.Text('',size=(1,1))],
    			[sg.TabGroup([[
                sg.Tab("Getting started", layout_getting_started),
                sg.Tab("Metabarcoding processing tools", layout_tools)
                ]]),
                sg.Multiline(welcome_text, size = (40, 32), key = '_OUTSTREAM_', autoscroll = True)], ## size = (w, h)
    			[sg.Text('',size=(1,1))],
    			[sg.Text('',size=(1,1))],
    			[sg.Exit(button_color=('black', 'red')), sg.Text("", size=(6,1)), sg.Button("Check installations", key="check_installations"), sg.Text("", size=(6,1)),
                sg.Text('Number of cores:'), sg.Slider(range=(1, num_cores), default_value=num_cores, size=(15, 10), orientation='horizontal', key="num_cores_to_use"), sg.Text("", size=(6,1)), sg.Button(image_filename=github_png, image_size=(26,26), key='open_github', button_color=('black', 'white')),
                sg.Button(key='open_twitter', button_color=('white', 'white'), image_filename=twitter_png, image_size=(26,26)), sg.Text('', size=(1,1)), sg.Text(meta_tools_version, font=('Arial', 8))]]

    # Create the Window
    MP_window = sg.Window('MetaProcessor', layout)
    win2_active=False

    ##########################################################################################################################

    default_options_dict = {
                        "Maxdiffpct":25,
                        "Maxdiffs":199,
                        "Minovlen":5,
                        "f_primers":["fwhF2", "teleo2_f", "BF2", "BF1"],
                        "r_primers":["fwhR2n", "teleo2_r", "BR2", "BR1"],
                        "min_len":130,
                        "max_len":210,
                        "maxee":1,
                        "unoise_minsize":8,
                        "unoise_alpha":2,
                        "clustering_threshold":0.97
                            }

    ##########################################################################################################################

    while True:
        try:
            event, values = MP_window.Read()

            if event is None or event == 'Exit':
                break

            if event == "open_demultiplexing":
                try:
                    import demultiplexer.__main__ as demultiplexer
                    MP_window.hide()
                    demultiplexer.main()
                except:
                    sg.PopupError("You have to install the demultipexer first!")

            if event == "open_boldigger":
                try:
                    import boldigger.__main__ as boldigger
                    MP_window.hide()
                    boldigger.main()
                except:
                    sg.PopupError("You have to install BOLDigger first!")

            if event == "open_log_file":
                try:
                    open_file(main_log_file)
                except:
                    sg.PopupError("Log file not found!")

            if event == "open_run_analyses":

                MP_window.hide()

                raw_reads_layout = [[ sg.Text("Input folder:"), sg.Input("", size=(30,1), key="raw_data_folder"), sg.FolderBrowse("Browse", initial_folder = path_to_outdirs) ]]

                pe_merging_layout = [[ sg.Text('Maxdiffpct:'), sg.Input(default_options_dict["Maxdiffpct"], size=(5,1), key="maxdiffpct"),
                                    sg.Text('Maxdiffs:'), sg.Input(default_options_dict["Maxdiffs"], size=(5,1), key="maxdiffs"),
                                    sg.Text('Minovlen:'), sg.Input(default_options_dict["Minovlen"], size=(5,1), key="minovlen") ]]

                primer_trimming_layout = [[ sg.Text("Forward primer:"), sg.Combo(available_fw_primers_list, size=(10,1), default_value=available_fw_primers_list[0], key="fw_primer"),
                                            sg.Text("Reverse primer:"), sg.Combo(available_rv_primers_list, size=(10,1), default_value=available_rv_primers_list[0], key="rv_primer"),
                                            sg.CB("Anchoring", default=False, key="anchoring"),
                                            sg.Button("Modify primers", key="modify_primer_sheet") ]]

                read_filter_layout = [[     sg.Text('Min. length:'), sg.Input(default_options_dict["min_len"], size=(5,1), key="max_length"),
                                            sg.Text('Max. length:'), sg.Input(default_options_dict["max_len"], size=(5,1), key="min_length"),
                                            sg.Text('Maxee threshold:'), sg.Input(default_options_dict["maxee"], size=(5,1), key="maxee_value"),
                                        ]]

                pre_processing_layout = [[ sg.Text("Dereplication, singleton removal, and chimera removal.") ]]

                denoising_layout = [[ sg.Text("Min. size:"), sg.Input(default_options_dict["unoise_minsize"], key="unoise_minsize", size=(5,1)),
                                      sg.Text("Alpha:"), sg.Input(default_options_dict["unoise_alpha"], key="unoise_alpha", size=(5,1)) ]]

                clustering_layout = [[  sg.Text('Cluster threshold:'), sg.Input(default_options_dict["clustering_threshold"], key="clustering_threshold", size=(5,1)),
                                        sg.Text('Input data:'), sg.Combo(['Denoised reads', 'Processed reads'], default_value='Denoised reads', key="clustering_input_data") ]]

                layout_run_analyses = [
                					[sg.Text('All-in-One analysis', size=(50,1), font=('Arial', 12, "bold"))],
                					[sg.Text('_'*115)],

                                    [sg.CB("", default=True, key="cb_raw_reads"), sg.Text("1. Raw reads:", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=raw_reads_layout, title="")],
                					[sg.Text('')],

                                    [sg.CB("", default=True, key="cb_pe_merging"), sg.Text("2. Paired-end merging", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=pe_merging_layout, title="")],
                					[sg.Text('')],

                                    [sg.CB("", default=True, key="cb_primer_trimming"), sg.Text("3. Primer trimming", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=primer_trimming_layout, title="")],
                					[sg.Text('')],

                                    [sg.CB("", default=True, key="cb_read_filtering"), sg.Text("4. Read filtering", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=read_filter_layout, title="")],
                					[sg.Text('')],

                                    [sg.CB("", default=True, key="cb_preprocessing"), sg.Text("5. Pre-processing", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=pre_processing_layout, title="")],
                					[sg.Text('')],

                                    [sg.CB("", default=True, key="cb_denoising"), sg.Text("6.1 Denoising", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=denoising_layout, title="")],
                					[sg.Text('')],

                                    [sg.CB("", default=True, key="cb_clustering"), sg.Text("6.2 OTU clustering", size=(18,1), font=('Arial', 11, "bold")), sg.Frame(layout=clustering_layout, title="")],

                					[sg.Text('',size=(1,1))],
                                    [sg.Button('Run analysis', size=(10,2)), sg.CB("Minimize MetaProcessor", default=True, key="minimize_window")],
                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Exit', button_color=('black', 'red'))]
                					]

                # create the demultiplexing window
                tools_window = sg.Window('All-in-One analysis', layout_run_analyses, keep_on_top=False)
                while (True):
                    ######################################
                    event, values2 = tools_window.Read()
                    ######################################
                    ## load command chain to execute
                    command_chain = [str(k) for k,v in values2.items() if "cb_" in str(k) and v == True]
                    ## load log file to check for possible overwriting
                    if os.path.isfile(main_log_file):
                        main_log_file_df = pd.read_excel(main_log_file)
                        performed_steps = main_log_file_df.columns.tolist()[1:]
                    else:
                        main_log_file_df = pd.DataFrame()
                        performed_steps = []

                    ######################################
                    if event in ('Exit', None):
                        break

                    if event == "Run analysis":

                        if values2["minimize_window"] == True:
                            tools_window.hide()

                        ## check if steps will be overwritten
                        performed_steps_dict = {"cb_raw_reads":"Raw reads", "cb_pe_merging":"Merged reads", "cb_primer_trimming":"Primer trimming", "cb_read_filtering":"Filtered reads", "cb_preprocessing":"Chimera removal"}
                        ## add special cases to the dict for denoising and clustering
                        performed_steps_dict["cb_denoising"] = "ESVs (id=" + str(values2["unoise_alpha"]) + ")"
                        if values2["clustering_input_data"] == "Processed reads":
                            performed_steps_dict["cb_clustering"] = "OTUs (id=" + str(values2["clustering_threshold"]) + ")"
                        else:
                            performed_steps_dict["cb_clustering"] = "OTUs (id=" + str(values2["clustering_threshold"]) + ";a=" + str(values2["unoise_alpha"]) + ")"

                        ## check if the command chain was already performed
                        ask = False
                        answer = "OK"
                        for command in command_chain:
                            if performed_steps_dict[command] in performed_steps:
                                ask = True
                        ## otherwise ask to overwrite
                        if ask == True:
                            answer = sg.PopupOKCancel("One or more tools were already run. Re-running the analysis will overwrite the results.\n\nContinue?")

                        ##########################################################################################
                        ## RUN ANALYSES ##

                        if answer == "OK":

                            break_test = ""

                            ## use the standard raw data folder if empty
                            if values2["raw_data_folder"] == "":
                                values2["raw_data_folder"] = Path(str(path_to_outdirs) + "/1_raw_data/_data/")

                            ## check if the samples are present
                            ## otherwise run the verify raw data script
                            if ("Sample ID" not in main_log_file_df.columns.tolist() or not os.path.isfile(main_log_file)):
                                if "cb_raw_reads" not in command_chain:
                                    print("")
                                    sg.Popup("Warning: You need to count the raw reads first!")
                                    from metaprocessor.verify_raw_data import verify_raw_data
                                    verify_raw_data(values2["raw_data_folder"], path_to_outdirs, main_log_file)
                                    print("")

                            ## 1 ##
                            ## verify raw reads
                            if ("cb_raw_reads" in command_chain and break_test != "OK"):
                                print("")
                                from metaprocessor.verify_raw_data import verify_raw_data
                                verify_raw_data(values2["raw_data_folder"], path_to_outdirs, main_log_file)
                                print("")
                                ## ask to continue
                                break_test = sg.popup_auto_close('Automatically continuing analysis...\n\nPress OK to abort!', auto_close_duration=4)

                            ## 2 ##
                            if ("cb_pe_merging" in command_chain and break_test != "OK"):
                                print("")
                                from metaprocessor.v_merge import v_merge
                                v_merge(values2["raw_data_folder"], path_to_outdirs, main_log_file, values["num_cores_to_use"], str(values2["minovlen"]), str(values2["maxdiffs"]), str(values2["maxdiffpct"]))
                                print("")
                                ## ask to continue
                                break_test = sg.popup_auto_close('Automatically continuing analysis...\n\nPress OK to abort!', auto_close_duration=4)

                            ## 3 ##
                            if ("cb_primer_trimming" in command_chain and break_test != "OK"):
                                print("")
                                p5_forward_primer = available_fw_primers_dict[values2["fw_primer"]]
                                p7_reverse_primer = available_rv_primers_dict[values2["rv_primer"]]
                                anchoring = values2["anchoring"]
                                from metaprocessor.cutadapt_trimming import cutadapt_trimming
                                cutadapt_trimming(path_to_outdirs, p5_forward_primer, p7_reverse_primer, anchoring, values["num_cores_to_use"], main_log_file)
                                print("")
                                ## ask to continue
                                break_test = sg.popup_auto_close('Automatically continuing analysis...\n\nPress OK to abort!', auto_close_duration=4)

                            ## 4 ##
                            if ("cb_read_filtering" in command_chain and break_test != "OK"):
                                print("")
                                from metaprocessor.v_read_filtering import v_read_filtering
                                v_read_filtering(path_to_outdirs, values2["min_length"], values2["max_length"], values2["maxee_value"], values["num_cores_to_use"], main_log_file)
                                print("")
                                ## ask to continue
                                break_test = sg.popup_auto_close('Automatically continuing analysis...\n\nPress OK to abort!', auto_close_duration=4)

                            ## 5 ##
                            if ("cb_preprocessing" in command_chain and break_test != "OK"):
                                print("")
                                from metaprocessor.v_derep_singletons_uchime import v_derep_singletons_uchime
                                v_derep_singletons_uchime(path_to_outdirs, values["num_cores_to_use"], main_log_file)
                                print("")
                                ## ask to continue
                                break_test = sg.popup_auto_close('Automatically continuing analysis...\n\nPress OK to abort!', auto_close_duration=4)

                            ## 6 ##
                            if ("cb_denoising" in command_chain and break_test != "OK"):
                                print("")
                                from metaprocessor.v_unoise import v_unoise
                                v_unoise(path_to_outdirs, values2["unoise_minsize"], values2["unoise_alpha"], main_log_file)
                                print("")
                                ## ask to continue
                                break_test = sg.popup_auto_close('Automatically continuing analysis...\n\nPress OK to abort!', auto_close_duration=4)

                            ## 7 ##
                            if ("cb_clustering" in command_chain and break_test != "OK"):
                                if values2["clustering_input_data"] == "Processed reads":
                                    print("")
                                    from metaprocessor.v_clustering import v_clustering_raw
                                    v_clustering_raw(path_to_outdirs, values2["clustering_threshold"],main_log_file)
                                    print("")
                                else:
                                    print("")
                                    from metaprocessor.v_clustering import v_clustering_denoised
                                    v_clustering_denoised(path_to_outdirs, values2["clustering_threshold"], values2["unoise_alpha"], main_log_file)
                                    print("")

                            ## finish the command chain
                            finished = sg.PopupYesNo("Jobs finished!\n\nContinue analyses?", title="")
                            if finished == "No":
                                break
                                tools_window.Close()
                            elif values2["minimize_window"] == True:
                                tools_window.UnHide()

                        elif values2["minimize_window"] == True:
                            tools_window.UnHide()

                    if event == "modify_primer_sheet":
                        open_file(primer_sheet_xlsx)

                        ##########################################################################################

                tools_window.Close()

            if event == "open_ncbi_blast":

                MP_window.hide()

                layout_ncbi_blast = [
                					[sg.Text('NCBI BLAST', size=(50,1), font=('Arial', 12, "bold"))],
                					[sg.Text('_'*115)],
                                    [sg.Frame(layout=[
                                    [sg.Text("Fasta file:", size=(20,1)), sg.Input("", size=(30,1), key="fasta_file"), sg.FileBrowse("Browse", initial_folder = path_to_outdirs)],
                                    [sg.Text("Read table:", size=(20,1)), sg.Input("", size=(30,1), key="read_table"), sg.FileBrowse("Browse", initial_folder = path_to_outdirs)],
                                    [sg.Text("BLAST xml file(s):", size=(20,1)), sg.Input("", size=(30,1), key="xml_files"), sg.FilesBrowse("Browse", initial_folder = path_to_outdirs)],
                                    [sg.Text("BLAST hit limit:"), sg.Input("10", size=(5,1), key="limit")],
                                    ], title="Input files")],
                                    [sg.Frame(layout=[
                                    [sg.Text('If the fasta files contains to many entries for an NCBI BLAST, create subsets here:')],
                                    [sg.Text("Batch size:"), sg.Input("150", size=(5,1), key="batch_size"), sg.Button("Run", key="subset_fasta")]
                                    ], title="Fasta subset")],
                					[sg.Text('')],
                                    [sg.Button('Create taxonomy table', size=(20,2), key=("fetch_taxonomy")), sg.CB("Minimize MetaProcessor", default=True, key="minimize_window")],
                                    [sg.Text('',size=(1,1))],
                                    [sg.Button('Exit', button_color=('black', 'red'))]
                					]

                # create the demultiplexing window
                blast_window = sg.Window('NCBI Blast', layout_ncbi_blast, keep_on_top=False)
                while (True):
                    ######################################
                    event, values2 = blast_window.Read()
                    ######################################

                    fasta_file = values2["fasta_file"]
                    read_table = values2["read_table"]
                    limit = values2["limit"]
                    batch_size = values2["batch_size"]
                    xml_files = values2["xml_files"].split(";")

                    if event in ('Exit', None):
                        break

                    if event == "subset_fasta":
                        if fasta_file == "":
                            sg.PopupError("Please provide a fasta file!", title="Error")
                        else:
                            print("")
                            from metaprocessor.ncbi_blast import subset_fasta
                            subset_fasta(fasta_file, batch_size)
                            print("")

                    if event == "fetch_taxonomy":
                        if fasta_file == "":
                            sg.PopupError("Please provide a fasta file!", title="Error")
                        elif read_table == "":
                            sg.PopupError("Please provide a read table!", title="Error")
                        else:
                            print("")
                            from metaprocessor.ncbi_blast import blast_xml_to_taxonomy
                            blast_xml_to_taxonomy(fasta_file, xml_files, read_table, limit)
                            print("")

                blast_window.Close()

        ###########################################################
            try:
                MP_window.UnHide()
            except:
                pass
        ###########################################################

        # define exceptions
        # if there are unexpected errors print a message and continue the script!
        except:
            ## unhide the main window if neccessary
            try:
                MP_window.UnHide()
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

    MP_window.Close()

## run only if called as toplevel script
if __name__ == "__main__":
    main()
