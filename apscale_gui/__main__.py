import os.path
import tkinter as tk
from tkinter import messagebox, filedialog, simpledialog
from pathlib import Path
import subprocess, glob, sys
from tkinter import font
import pandas as pd
import apscale.b_pe_merging as b_pe_merging
import apscale.c_primer_trimming as c_primer_trimming
import apscale.d_quality_filtering as d_quality_filtering
import apscale.e_dereplication as e_dereplication
import apscale.f_denoising as f_denoising
import apscale.g_generate_esv_table as g_generate_esv_table
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
import multiprocessing
import boldigger2.id_engine_coi as id_engine_coi
import importlib.metadata
from update_checker import update_check
import importlib.metadata
import platform

## check for updates
def check_package_update(package_name):
    # Get the currently installed version
    try:
        installed_version = importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        print(f"{package_name} is not installed.")
        return

    # Check for updates
    res = update_check(package_name, installed_version)

check_package_update('apscale_gui')
check_package_update('apscale')
check_package_update('apscale_blast')
check_package_update('boldigger2')

## initialize a pressed button
def button_pressed(button_name, root):
    global result
    result = button_name
    root.destroy()

## Return to main menu
def cancel(project_folder, current_project, root):
    root.destroy()
    apscale_gui(project_folder, current_project)

## restart script function
def return_to_start(project_folder, current_project):
    apscale_gui(project_folder, current_project)

## provide a window to select a folder from
def get_working_folder():
    selected_folder = None

    # Create the main window
    root = tk.Tk()
    root.title("Open Folder")

    # Create a button to open the folder dialog
    open_button = tk.Button(root, text="Open Folder", command=open_folder)
    open_button.pack(pady=20)

    # Run the application
    root.mainloop()

    # Return the selected folder path
    return selected_folder

## import the user default project location
def import_user_data():
    ## open user_data_txt to save the standard output put
    try:
        user_data_txt = Path(__file__).resolve().parent.joinpath('_user_data', 'user_data.txt')
        f = open(user_data_txt)
        projects_main_path = f.read()
        return Path(projects_main_path)
    except:
        folder = Path(__file__).resolve().parent.joinpath('_user_data')
        file = Path(__file__).resolve().parent.joinpath('_user_data', 'user_data.txt')
        os.mkdir(folder)
        f = open(file, 'w')
        f.close()
        return Path('./')

## if apscale gui is started for the first time, initialize a window that asks for the apscale projects folder
def create_initialize_window():
    # Set a default project folder
    default_project_folder = import_user_data()
    project_folder = default_project_folder

    def select_project_folder():
        nonlocal project_folder
        selected_folder = filedialog.askdirectory()
        if selected_folder:  # Check if a folder was selected
            project_folder = Path(selected_folder)
            folder_label.config(text=f"Selected Folder: {project_folder}", fg="green")
        else:
            folder_label.config(text=f"Default Folder: {default_project_folder} (Not Selected)", fg="gray")

    def submit():
        if project_folder.exists():  # Ensure a valid folder has been selected
            user_data_txt = Path(__file__).resolve().parent.joinpath('_user_data', 'user_data.txt')
            f = open(user_data_txt, 'w')
            f.write(str(project_folder))
            f.close()
            root.destroy()
        else:
            messagebox.showwarning("Warning", "Please select a valid folder before proceeding.")

    # Create the main window
    root = tk.Tk()
    root.title("Initialize Apscale-GUI")
    root.geometry("500x300")  # Adjust window size for better UI

    # Add welcome text
    label = tk.Label(root, text="Welcome to Apscale-GUI", font=font.Font(weight="bold", size=14))
    label.pack(pady=10)

    # Horizontal separator
    spacer1 = tk.Frame(root, height=2, bd=1, relief=tk.SUNKEN)
    spacer1.pack(fill=tk.X, padx=5, pady=5)

    # Welcome message with more detail
    welcome_label = tk.Label(root, text="Please choose a folder where your Apscale projects will be stored.\nIf you don't select one, the default folder will be used.")
    welcome_label.pack(pady=10)

    # Button to open the folder dialog
    folder_button = tk.Button(root, text="Browse", command=select_project_folder, width=15)
    folder_button.pack(pady=10)

    # Label to display the selected folder, with clear text for default and selected states
    folder_label = tk.Label(root, text=f"Default Folder: {default_project_folder}", fg="gray")
    folder_label.pack(pady=10)

    # Submit button with a wider area for better clickability
    submit_button = tk.Button(root, text="Continue", command=submit, width=15)
    submit_button.pack(pady=20)

    # Add padding for better visual spacing
    for widget in [label, welcome_label, folder_button, folder_label, submit_button]:
        widget.pack_configure(padx=20)

    # Create apscale_database folder if not existing, but delay creation until submit is clicked
    def create_project_folders():
        apscale_databases = project_folder.joinpath('apscale_databases')
        os.makedirs(apscale_databases, exist_ok=True)

    # Ensure project folders are created only after the user confirms the selection
    submit_button.config(command=lambda: [create_project_folders(), submit()])

    # Run the application
    root.mainloop()

    return Path(project_folder)

## collect the existing projects
def get_existing_projects(project_folder):
    # Collect all projects
    projects = sorted(glob.glob(str(project_folder) + "/*_apscale"))
    if projects == []:
        projects = ['Not selected']
    return projects

## open the project folder in a file browser
def open_project_folder(folder_path):
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

## create the start window
def create_start_window(project_folder, current_project):
    # Create the main window
    root = tk.Tk()
    root.title("Apscale-GUI")

    # Set up a grid layout
    root.columnconfigure(0, weight=1)

    # Welcome Section
    welcome_frame = tk.Frame(root)
    welcome_frame.pack(pady=10)

    welcome_label = tk.Label(welcome_frame, text="Welcome to Apscale-GUI", font=font.Font(weight="bold", size=14))
    welcome_label.pack()

    # Current Project Section
    project_frame = tk.Frame(root)
    project_frame.pack(pady=10)

    project_info_label = tk.Label(project_frame, text=f"Current project folder: {project_folder}")
    project_info_label.pack(pady=5)

    current_project_name = Path(current_project).name
    active_project_label = tk.Label(project_frame, text=f"Active project: {current_project_name}")
    active_project_label.pack(pady=5)

    # Button to open the current project folder
    open_folder_button = tk.Button(project_frame, text="Open Current Project Folder",
                                    command=lambda: open_project_folder(project_folder))
    open_folder_button.pack(pady=5)

    # Project Actions Section
    action_frame = tk.Frame(root)
    action_frame.pack(pady=10)

    action_label = tk.Label(action_frame, text="Create or load a project", font=font.Font(weight="bold", size=12))
    action_label.pack(pady=5)

    button_a = tk.Button(action_frame, text="Create New Apscale Project",
                         command=lambda: button_pressed("Create new apscale project", root))
    button_a.pack(pady=5, padx=10)

    button_b = tk.Button(action_frame, text="Load Existing Apscale Project",
                         command=lambda: button_pressed("Load existing apscale project", root))
    button_b.pack(pady=5, padx=10)

    # Task Section
    task_frame = tk.Frame(root)
    task_frame.pack(pady=10)

    task_label = tk.Label(task_frame, text="Apscale Pipeline", font=font.Font(weight="bold", size=12))
    task_label.pack(pady=5)

    button_c = tk.Button(task_frame, text="Raw Data Processing",
                         command=lambda: button_pressed("Raw data processing", root))
    button_c.pack(pady=5, padx=10)

    button_d = tk.Button(task_frame, text="Local BLASTn", command=lambda: button_pressed("Local blastn", root))
    button_d.pack(pady=5, padx=10)

    button_e = tk.Button(task_frame, text="BOLDigger2", command=lambda: button_pressed("BOLDigger2", root))
    button_e.pack(pady=5, padx=10)

    # Exit Button Section
    exit_frame = tk.Frame(root)
    exit_frame.pack(pady=10)

    exit_button = tk.Button(exit_frame, text="Exit", command=sys.exit, bg='red', fg='black')
    exit_button.pack(pady=5)

    # Run the application
    root.mainloop()

## window to provide a project name
def get_project_name():
    # Prompt the user for the project name
    project_name = simpledialog.askstring("Input", "Enter the project name:")

    ## check if project name is empty
    if project_name == '':
        project_name = 'default'

    # Return the project name
    return project_name

## window to create a new project
def run_create_new_project(project_folder):
    ## tkinter input window for a project name
    try:
        project_name = get_project_name()
        current_project = project_folder.joinpath(project_name)
    except:
        current_project = project_folder.joinpath('default')

    # Example command
    command = ["apscale", "--create_project", current_project]
    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)
    # Print the output
    print(result.stdout)

    current_project = Path(str(current_project) + '_apscale')

    ## return to main
    return_to_start(project_folder, current_project)

## dropdown menu to select a project
def select_project(projects, project_folder):
    # Create a window with a dropdown menu
    root = tk.Tk()
    root.title("Select Project")

    # Add welcome text
    label = tk.Label(root, text="Select an existing project:", font=font.Font(weight="bold"))
    label.pack(pady=10)

    # Variable to store the selected project
    selected_project = tk.StringVar(root)
    selected_project.set(projects[0])  # Set default value

    # Create dropdown menu
    dropdown = tk.OptionMenu(root, selected_project, *projects)
    dropdown.pack(pady=20)

    # Button to confirm selection
    def on_confirm():
        root.destroy()

    confirm_button = tk.Button(root, text="Confirm", command=on_confirm)
    confirm_button.pack(pady=10)

    root.mainloop()
    return selected_project.get()

## load projects
def run_load_existing_project(project_folder):
    # Collect projects
    projects = get_existing_projects(project_folder)

    # Select project
    project_name = select_project(projects, project_folder)

    # Provide current project
    if project_name != 'Not selected':
        current_project = Path(project_folder).joinpath(project_name)
    else:
        current_project = 'Not selected'

    # Return to start
    return_to_start(project_folder, current_project)

## load settings file
def load_settings_file(project_folder, current_project):
    settings_xlsx = 'Settings_' + Path(current_project).name.replace('_apscale', '') + '.xlsx'
    settings_xlsx_path = Path(current_project).joinpath(settings_xlsx)

    data = {}
    for sheet in ['0_general_settings', '3_PE_merging', '4_primer_trimming', '5_quality_filtering', '6_denoising']:
        settings_df = pd.read_excel(settings_xlsx_path, sheet_name=sheet).fillna('')
        res = {i:str(settings_df[i].values.tolist()[0]) for i in settings_df.columns}
        data = data | res

    return data

## run apscale pipeline
def run_apscale(project_folder, current_project, updated_data):
    ## Create new settings file with updated values
    settings_xlsx = 'Settings_' + Path(current_project).name.replace('_apscale', '') + '.xlsx'
    settings_xlsx_path = Path(current_project).joinpath(settings_xlsx)

    # Create a Pandas Excel writer using XlsxWriter as the engine
    with pd.ExcelWriter(settings_xlsx_path, engine='xlsxwriter') as writer:
        # general settings
        data = [[int(updated_data['cores to use']), int(updated_data['compression level'])]]
        cols = ['cores to use', 'compression level']
        df = pd.DataFrame(data, columns=cols)
        df.to_excel(writer, sheet_name='0_general_settings', index=False)

        # PE merging
        data = [[int(updated_data['maxdiffpct']), int(updated_data['maxdiffs']), int(updated_data['minovlen'])]]
        cols = ['maxdiffpct', 'maxdiffs', 'minovlen']
        df = pd.DataFrame(data, columns=cols)
        df.to_excel(writer, sheet_name='3_PE_merging', index=False)

        # Primer trimming
        data = [[updated_data["P5 Primer (5' - 3')"], updated_data["P7 Primer (5' - 3')"], updated_data['anchoring']]]
        cols = ["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", 'anchoring']
        df = pd.DataFrame(data, columns=cols)
        df.to_excel(writer, sheet_name='4_primer_trimming', index=False)

        # Quality filtering
        data = [[float(updated_data["maxEE"]), int(updated_data["min length"]), int(updated_data['max length'])]]
        cols = ["maxEE", "min length", 'max length']
        df = pd.DataFrame(data, columns=cols)
        df.to_excel(writer, sheet_name='5_quality_filtering', index=False)

        # Denoising
        data = [[float(updated_data['alpha']), int(updated_data['minsize'])]]
        cols = ['alpha', 'minsize']
        df = pd.DataFrame(data, columns=cols)
        df.to_excel(writer, sheet_name='6_denoising', index=False)

    ## also collect the task to perform
    task = updated_data['dropdown']

    ## Run apscale
    if task == "Run apscale":
        b_pe_merging.main(current_project)
        c_primer_trimming.main(current_project)
        d_quality_filtering.main(current_project)
        e_dereplication.main(current_project)
        f_denoising.main(current_project)
        g_generate_esv_table.main(current_project)
    if task == "PE merging":
        b_pe_merging.main(current_project)
    if task == "primer trimming":
        c_primer_trimming.main(current_project)
    if task == "quality filtering":
        d_quality_filtering.main(current_project)
    if task == "dereplication":
        e_dereplication.main(current_project)
    if task == "denoising":
        f_denoising.main(current_project)
    if task == "generate ESV table":
        g_generate_esv_table.main(current_project)

## apscale main window
def create_apscale_main_window(project_folder, current_project):

    def save_and_submit(project_folder, current_project):
        updated_data = {
            "cores to use": entry1.get(),
            "compression level": entry2.get(),
            "maxdiffpct": entry3.get(),
            "maxdiffs": entry4.get(),
            "minovlen": entry5.get(),
            "P5 Primer (5' - 3')": entry6.get(),
            "P7 Primer (5' - 3')": entry7.get(),
            "anchoring": entry8.get(),
            "maxEE": entry9.get(),
            "min length": entry10.get(),
            "max length": entry11.get(),
            "alpha": entry12.get(),
            "minsize": entry13.get(),
            "dropdown": dropdown_var.get()
        }
        run_apscale(project_folder, current_project, updated_data)

    root = tk.Tk()
    root.title("Apscale")
    data = load_settings_file(project_folder, current_project)

    # Create labels and entry fields for each category
    labels = [
        "cores to use", "compression level",
        "maxdiffpct", "maxdiffs", "minovlen",
        "P5 Primer (5' - 3')", "P7 Primer (5' - 3')", "anchoring",
        "maxEE", "min length", "max length",
        "alpha", "minsize"
    ]
    default_values = list(data.values())
    entries = []
    entry_vars = []

    # Add welcome text
    label = tk.Label(root, text="Raw data processing", font=font.Font(weight="bold"))
    label.grid(row=0, column=0, columnspan=2, pady=10)

    # Headers and input fields with reduced spacing
    # General settings (category 1)
    tk.Label(root, text="General settings", font=font.Font(weight="bold")).grid(row=2, column=0, columnspan=2, pady=(10, 5))
    for i, (label, default_value) in enumerate(zip(labels[:2], default_values[:2]), start=3):
        tk.Label(root, text=label).grid(row=i, column=0, padx=10, pady=(2, 2), sticky="w")
        entry_var = tk.StringVar(value=default_value)
        entry = tk.Entry(root, textvariable=entry_var)
        entry.grid(row=i, column=1, padx=10, pady=(2, 2), sticky="w")
        entries.append(entry)

    # PE merging (category 2)
    tk.Label(root, text="PE merging", font=font.Font(weight="bold")).grid(row=5, column=0, columnspan=2, pady=(10, 5))
    for i, (label, default_value) in enumerate(zip(labels[2:5], default_values[2:5]), start=6):
        tk.Label(root, text=label).grid(row=i, column=0, padx=10, pady=(2, 2), sticky="w")
        entry_var = tk.StringVar(value=default_value)
        entry = tk.Entry(root, textvariable=entry_var)
        entry.grid(row=i, column=1, padx=10, pady=(2, 2), sticky="w")
        entries.append(entry)

    # Primer trimming (category 3)
    tk.Label(root, text="Primer trimming", font=font.Font(weight="bold")).grid(row=9, column=0, columnspan=2, pady=(10, 5))
    for i, (label, default_value) in enumerate(zip(labels[5:8], default_values[5:8]), start=10):
        tk.Label(root, text=label).grid(row=i, column=0, padx=10, pady=(2, 2), sticky="w")
        entry_var = tk.StringVar(value=default_value)
        entry = tk.Entry(root, textvariable=entry_var)
        entry.grid(row=i, column=1, padx=10, pady=(2, 2), sticky="w")
        entries.append(entry)

    # Quality filtering (category 4)
    tk.Label(root, text="Quality filtering", font=font.Font(weight="bold")).grid(row=13, column=0, columnspan=2, pady=(10, 5))
    for i, (label, default_value) in enumerate(zip(labels[8:11], default_values[8:11]), start=14):
        tk.Label(root, text=label).grid(row=i, column=0, padx=10, pady=(2, 2), sticky="w")
        entry_var = tk.StringVar(value=default_value)
        entry = tk.Entry(root, textvariable=entry_var)
        entry.grid(row=i, column=1, padx=10, pady=(2, 2), sticky="w")
        entries.append(entry)

    # Denoising (category 5)
    tk.Label(root, text="Denoising", font=font.Font(weight="bold")).grid(row=18, column=0, columnspan=2, pady=(10, 5))
    for i, (label, default_value) in enumerate(zip(labels[11:], default_values[11:]), start=19):
        tk.Label(root, text=label).grid(row=i, column=0, padx=10, pady=(2, 2), sticky="w")
        entry_var = tk.StringVar(value=default_value)
        entry = tk.Entry(root, textvariable=entry_var)
        entry.grid(row=i, column=1, padx=10, pady=(2, 2), sticky="w")
        entries.append(entry)

    entry1, entry2, entry3, entry4, entry5, entry6, entry7, entry8, entry9, entry10, entry11, entry12, entry13 = entries

    # Create a dropdown menu
    tk.Label(root, text="Select task:", font=font.Font(weight="bold")).grid(row=23, column=0, pady=10)
    dropdown_var = tk.StringVar(value="Run apscale")
    dropdown_options = ["Run apscale", "PE merging", "primer trimming", "quality filtering",
                        "dereplication", "denoising", "generate ESV table"]
    dropdown_menu = tk.OptionMenu(root, dropdown_var, *dropdown_options)
    dropdown_menu.grid(row=23, column=1, pady=10, sticky="w")

    # Submit and return buttons
    submit_button = tk.Button(root, text="Save & submit", command=lambda: save_and_submit(project_folder, current_project))
    submit_button.grid(row=24, column=1, columnspan=2, pady=10, sticky="w")
    return_button = tk.Button(root, text="Return", command=lambda: cancel(project_folder, current_project, root))
    return_button.grid(row=25, column=1, columnspan=2, pady=10, sticky="w")

    root.mainloop()

def create_local_blastn_window(project_folder, current_project):
    ####################################################################################################################
    ## run apscale pipeline
    def save_and_submit(project_folder, current_project):
        # Capture the updated values from input fields
        updated_data = {
            'cores': int(cores_entry.get()),
            'subset_size': int(subset_size_entry.get()),
            'max_target_seqs': int(max_target_seqs_entry.get()),
            'Species': int(species_entry.get()),
            'Genus': int(genus_entry.get()),
            'Family': int(family_entry.get()),
            'Order': int(order_entry.get()),
            'Class': int(class_entry.get()),
            'task': dropdown_var_task.get(),
            'query_fasta': str(Path(current_project).joinpath('8_esv_table', dropdown_var_fasta.get())),
            'db_folder': str(Path(project_folder).joinpath('apscale_databases', dropdown_var_dbs.get())),
            'out': str(Path(current_project).joinpath('8_esv_table', dropdown_var_fasta.get().replace('.fasta', '')))
        }

        a_blastn('blastn',
                 updated_data['query_fasta'],
                 updated_data['db_folder'],
                 updated_data['out'],
                 updated_data['cores'],
                 updated_data['task'],
                 updated_data['subset_size'],
                 updated_data['max_target_seqs'],
                 'nan')

        # create thresholds
        thresholds = [updated_data['Species'], updated_data['Genus'], updated_data['Family'], updated_data['Order'], updated_data['Class']]
        b_filter(updated_data['out'],
                 updated_data['db_folder'],
                 ','.join([str(i) for i in thresholds]),
                 updated_data['cores'],
                 )

    ####################################################################################################################

    # Initialize data dictionary
    data = {
        'cores': multiprocessing.cpu_count() - 2,
        'subset_size': 100,
        'max_target_seqs': 20,
        'Species': 97,
        'Genus': 95,
        'Family': 90,
        'Order': 87,
        'Class': 85
    }

    # Create the main window
    root = tk.Tk()
    root.title("Local blastn")

    # Add welcome text
    label = tk.Label(root, text="Blast identification engine", font=font.Font(weight="bold"))
    label.grid(row=0, column=0, columnspan=2, pady=10)

    # Create input fields for each data entry
    tk.Label(root, text="Cores:").grid(row=2, column=0, sticky="e")
    cores_entry = tk.Entry(root)
    cores_entry.insert(0, data['cores'])
    cores_entry.grid(row=2, column=1)

    tk.Label(root, text="Subset Size:").grid(row=3, column=0, sticky="e")
    subset_size_entry = tk.Entry(root)
    subset_size_entry.insert(0, data['subset_size'])
    subset_size_entry.grid(row=3, column=1)

    tk.Label(root, text="Max Target Seqs:").grid(row=4, column=0, sticky="e")
    max_target_seqs_entry = tk.Entry(root)
    max_target_seqs_entry.insert(0, data['max_target_seqs'])
    max_target_seqs_entry.grid(row=4, column=1)

    tk.Label(root, text="Species:").grid(row=5, column=0, sticky="e")
    species_entry = tk.Entry(root)
    species_entry.insert(0, data['Species'])
    species_entry.grid(row=5, column=1)

    tk.Label(root, text="Genus:").grid(row=6, column=0, sticky="e")
    genus_entry = tk.Entry(root)
    genus_entry.insert(0, data['Genus'])
    genus_entry.grid(row=6, column=1)

    tk.Label(root, text="Family:").grid(row=7, column=0, sticky="e")
    family_entry = tk.Entry(root)
    family_entry.insert(0, data['Family'])
    family_entry.grid(row=7, column=1)

    tk.Label(root, text="Order:").grid(row=8, column=0, sticky="e")
    order_entry = tk.Entry(root)
    order_entry.insert(0, data['Order'])
    order_entry.grid(row=8, column=1)

    tk.Label(root, text="Class:").grid(row=9, column=0, sticky="e")
    class_entry = tk.Entry(root)
    class_entry.insert(0, data['Class'])
    class_entry.grid(row=9, column=1)

    # Create a folder selection button
    label = tk.Label(root, text="Select task:", font=font.Font(weight="bold"))
    label.grid(row=16, column=0, pady=10)
    dropdown_var_task = tk.StringVar(value="blastn")
    dropdown_options = ["blastn", "megablast", "dc-megablast"]
    dropdown_menu = tk.OptionMenu(root, dropdown_var_task, *dropdown_options)
    dropdown_menu.grid(row=16, column=1, pady=10, sticky="w")

    # Select fasta file
    label = tk.Label(root, text="Select query fasta:", font=font.Font(weight="bold"))
    label.grid(row=17, column=0, pady=10)
    all_fasta_files = glob.glob(str(Path(current_project).joinpath('8_esv_table/*.fasta')))
    if all_fasta_files == []:
        fasta_names = ['No files available.']
    else:
        fasta_names = [Path(i).name for i in all_fasta_files]
    dropdown_var_fasta = tk.StringVar(value=fasta_names[0])
    dropdown_options = fasta_names
    dropdown_menu = tk.OptionMenu(root, dropdown_var_fasta, *dropdown_options)
    dropdown_menu.grid(row=17, column=1, pady=10, sticky="w")

    # Select fasta file
    label = tk.Label(root, text="Select database:", font=font.Font(weight="bold"))
    label.grid(row=18, column=0, pady=10)
    all_dbs = glob.glob(str(Path(project_folder).joinpath('apscale_databases/*')))
    if all_dbs == []:
        db_names = ['No databases available.']
    else:
        db_names = [Path(i).name for i in all_dbs]
    dropdown_var_dbs = tk.StringVar(value=db_names[0])
    dropdown_options = db_names
    dropdown_menu = tk.OptionMenu(root, dropdown_var_dbs, *dropdown_options)
    dropdown_menu.grid(row=18, column=1, pady=10, sticky="w")

    # Create a submit button
    submit_button = tk.Button(root, text="Save & Submit", command=lambda: save_and_submit(project_folder, current_project))
    submit_button.grid(row=19, column=1, columnspan=2, pady=10, sticky="w")

    # Create a return button
    return_button = tk.Button(root, text="Return", command=lambda: cancel(project_folder, current_project, root))
    return_button.grid(row=20, column=1, columnspan=2, pady=10, sticky="w")

    # Run the application
    root.mainloop()

def create_boldigger2_window(project_folder, current_project):
    ####################################################################################################################
    ## run apscale pipeline
    def save_and_submit(project_folder, current_project):
        # Capture the updated values from input fields
        updated_data = {
            'username': username_entry.get(),
            'password': pw_entry.get(),
            'Species': int(species_entry.get()),
            'Genus': int(genus_entry.get()),
            'Family': int(family_entry.get()),
            'Order': int(order_entry.get()),
            'Class': int(class_entry.get()),
            'query_fasta': str(Path(current_project).joinpath('8_esv_table', dropdown_var_fasta.get())),
        }

        ## collect thresholds
        thresholds = [updated_data['Species'], updated_data['Genus'], updated_data['Family'], updated_data['Order'], updated_data['Class']]

        ## start boldigger2
        id_engine_coi.main(
            updated_data['query_fasta'],
            updated_data['username'],
            updated_data['password'],
            thresholds,
        )

    ####################################################################################################################

    # Initialize data dictionary
    data = {
        'username': '',
        'password': '',
        'Species': 97,
        'Genus': 95,
        'Family': 90,
        'Order': 85,
        'Class': 50
    }

    # Create the main window
    root = tk.Tk()
    root.title("Boldigger2")

    # Add welcome text
    label = tk.Label(root, text="Boldigger2 identification engine", font=font.Font(weight="bold"))
    label.grid(row=0, column=0, columnspan=2, pady=10)

    # Create input fields for each data entry
    tk.Label(root, text="Username:").grid(row=2, column=0, sticky="e")
    username_entry = tk.Entry(root)
    username_entry.insert(0, data['username'])
    username_entry.grid(row=2, column=1)

    tk.Label(root, text="Password:").grid(row=3, column=0, sticky="e")
    pw_entry = tk.Entry(root, show="*")  # Hide the password by showing '*' for each character
    pw_entry.insert(0, data['password'])
    pw_entry.grid(row=3, column=1)

    tk.Label(root, text="Species:").grid(row=5, column=0, sticky="e")
    species_entry = tk.Entry(root)
    species_entry.insert(0, data['Species'])
    species_entry.grid(row=5, column=1)

    tk.Label(root, text="Genus:").grid(row=6, column=0, sticky="e")
    genus_entry = tk.Entry(root)
    genus_entry.insert(0, data['Genus'])
    genus_entry.grid(row=6, column=1)

    tk.Label(root, text="Family:").grid(row=7, column=0, sticky="e")
    family_entry = tk.Entry(root)
    family_entry.insert(0, data['Family'])
    family_entry.grid(row=7, column=1)

    tk.Label(root, text="Order:").grid(row=8, column=0, sticky="e")
    order_entry = tk.Entry(root)
    order_entry.insert(0, data['Order'])
    order_entry.grid(row=8, column=1)

    tk.Label(root, text="Class:").grid(row=9, column=0, sticky="e")
    class_entry = tk.Entry(root)
    class_entry.insert(0, data['Class'])
    class_entry.grid(row=9, column=1)

    # Select fasta file
    label = tk.Label(root, text="Select query fasta:", font=font.Font(weight="bold"))
    label.grid(row=17, column=0, pady=10)
    all_fasta_files = glob.glob(str(Path(current_project).joinpath('8_esv_table/*.fasta')))
    if all_fasta_files == []:
        fasta_names = ['No files available.']
    else:
        fasta_names = [Path(i).name for i in all_fasta_files]
    dropdown_var_fasta = tk.StringVar(value=fasta_names[0])
    dropdown_var_fasta = tk.StringVar(value=fasta_names[0])
    dropdown_options = fasta_names
    dropdown_menu = tk.OptionMenu(root, dropdown_var_fasta, *dropdown_options)
    dropdown_menu.grid(row=17, column=1, pady=10, sticky="w")

    # Create a submit button
    submit_button = tk.Button(root, text="Save & Submit", command=lambda: save_and_submit(project_folder, current_project))
    submit_button.grid(row=19, column=1, columnspan=2, pady=10, sticky="w")

    # Create a return button
    return_button = tk.Button(root, text="Return", command=lambda: cancel(project_folder, current_project, root))
    return_button.grid(row=20, column=1, columnspan=2, pady=10, sticky="w")

    # Run the application
    root.mainloop()

## error message
def show_error_window(message):
    error_root = tk.Tk()
    error_root.title("Error")
    tk.Label(error_root, text=message, padx=20, pady=20).pack()
    tk.Button(error_root, text="OK", command=error_root.destroy).pack(pady=10)
    error_root.mainloop()

## main function (needs to be seperate to be called again when exiting certain windows
def apscale_gui(project_folder, current_project):
    global result
    result = None

    ## create current_project dummy for display
    if current_project is None:
        current_project = 'Not selected'

    ## create the apscale main window
    create_start_window(project_folder, current_project)

    if result == "Create new apscale project":
        current_project = run_create_new_project(project_folder)

    if result == "Load existing apscale project":
        current_project = run_load_existing_project(project_folder)

    if result == "Raw data processing":
        if current_project != 'Not selected':
            create_apscale_main_window(project_folder, current_project)
        else:
            show_error_window("No project folder selected. Please create or load a project first.")
            return_to_start(project_folder, current_project)

    if result == "Local blastn":
        if current_project != 'Not selected':
            create_local_blastn_window(project_folder, current_project)
        else:
            show_error_window("No project folder selected. Please create or load a project first.")
            return_to_start(project_folder, current_project)

    if result == "BOLDigger2":
        if current_project != 'Not selected':
            create_boldigger2_window(project_folder, current_project)
        else:
            show_error_window("No project folder selected. Please create or load a project first.")
            return_to_start(project_folder, current_project)

def main():
    ## create the initialize window
    project_folder = create_initialize_window()
    project_folder = Path(project_folder)
    current_project = None

    apscale_gui(project_folder, current_project)

## run only if called as toplevel script
if __name__ == "__main__":
    main()
