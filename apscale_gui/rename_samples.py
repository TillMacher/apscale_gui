def create_rename_sheet(data_folder, path_to_outdirs):

    import os, glob
    import pandas as pd
    from pathlib import Path

    # collect the files from the inout folder
    data_folder = Path(data_folder)
    input_files = sorted(glob.glob(str(data_folder) + "/*.fastq.gz"))

    ## create a rename sheet
    rename_sheet_list = []
    for i, file in enumerate(input_files):
        if (i % 2) == 0:
            replicate = "_r1.fastq.gz"
        else:
            replicate = "_r2.fastq.gz"
        file_name = str(Path(file).stem).replace(".fastq", "")
        file_path = str(Path(file))
        new_name = ""
        new_name_function = "=C" + str(i+2) + "&D" + str(i+2)
        rename_sheet_list = rename_sheet_list + [[file_path, file_name, new_name, replicate, new_name_function]]

    rename_sheet = str(path_to_outdirs) + "/rename_sheet.xlsx"
    df = pd.DataFrame(rename_sheet_list, columns=["File path", "File name", "ENTER NEW NAME", "Replicate suffix", "New name"])
    df.to_excel(rename_sheet, index=False)

def rename_samples(rename_sheet):
    import os
    import pandas as pd
    from pathlib import Path
    import PySimpleGUI as sg

    # rename_sheet = "/Volumes/Rex/rename_sheet_spring_fisch.xlsx"

    df = pd.read_excel(rename_sheet)
    rename_list = df.values.tolist()

    ## check table
    n_files = len(df["New name"].values.tolist())
    n_pairs = n_files / 2
    ## i) for duplicates
    if len(set(df["New name"].values.tolist())) != n_pairs*2:
        sg.PopupError("Error: Duplicates found!\nPlease check your rename sheet!")
    else:
        answer = sg.PopupYesNo("Warning: Files will be overwritten!\nMake sure to create a backup first!\n\nContinue?")
        if answer == "Yes":
            for file in rename_list:
                old_file = file[0]
                new_file = str(Path(old_file).parent) + "/" + file[-1]
                os.rename(old_file, new_file)








#
