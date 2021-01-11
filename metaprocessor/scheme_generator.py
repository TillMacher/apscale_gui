import datetime, openpyxl
import PySimpleGUI as sg
from pathlib import Path

## funtion to generate the layout for the scheme generator at runtime
## combinations is an integer, primerset a path to the currently selected primerset
def scheme_layout(combinations, primerset):

    ## read data from the primerset
    data = [line.strip().split(',') for line in open(primerset, 'r')]
    data = {primer[0]: primer[2] for primer in data}

    ## extract primer names per direction
    fwd_names = [name for name in data.keys() if data[name] == 'fwd']
    rev_names = [name for name in data.keys() if data[name] == 'rev']

    ## define the layout for the combo elements to select the primers
    layout = [[sg.Text('Combination {}'.format(comb + 1)), sg.Combo(fwd_names, size = (15, 1)), sg.Combo(rev_names, size = (15, 1))] for comb in range(combinations)]

    ## insert the missing elements e.g. Input for name, cancal and save button
    layout.insert(0, [sg.Text('Name of tagging scheme:'), sg.Input(size = (26, 1), key = '_TAGGING_SCHEME_NAME_')])
    layout.append([sg.Button('Close'), sg.Button('Save')])

    return layout

## function to save the tagging scheme
## gets the paired file pairs for the loaded files to demultiplex
## the selected primerpairs, still have to be parsed from the scheme dict
## and the print handle to stream the output to the main window
def save_scheme(file_pairs, combinations, savepath, print_handle):

    # save the savename for later, remove after so keys are less cluttered
    savename =  combinations['_TAGGING_SCHEME_NAME_']
    del combinations['_TAGGING_SCHEME_NAME_']

    ## check if filename already exist before saving
    if Path(savepath).joinpath('{}.xlsx'.format(savename)).is_file():
        print_handle.print('{}: This tagging scheme already exists. Please choose another name.'.format(datetime.datetime.now().strftime("%H:%M:%S")))

    ## check if the combination input is valid e.g. all combinations were set
    elif not all([bool(combinations[key]) for key in combinations.keys()]):
        print_handle.print('{}: There are missing values in the combinations. Please fill all combination fields'.format(datetime.datetime.now().strftime("%H:%M:%S")))

    ## if all criteria are met save the tagging template
    else:
        ## put all primer combinations in a list as pairs, also generate the header of the tagging template
        combolist = [' - '.join((combinations[i * 2], combinations[i * 2 + 1])) for i in range(int(len(combinations) / 2))]
        header = ['File path forward', 'File path reverse', 'File name forward', 'File name reverse']

        ## excel save code
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = savename

        ## add the header, freeze column and header, so you only scroll thorugh the samples
        ws.append(header + combolist)
        ws.freeze_panes = 'E2'

        ## add the pathes to the files as well as the file names
        data = [[Path(path[0]).__str__(), Path(path[1]).__str__(), Path(path[0]).name, Path(path[1]).name] for path in file_pairs]

        ## write to the template
        for line in data:
            ws.append(line)

        ## asave the tagging scheme
        wb.save(Path(savepath).joinpath('{}.xlsx'.format(savename)))

        print_handle.print('{}: Tagging scheme was saved at: \n {}'.format(datetime.datetime.now().strftime("%H:%M:%S"), Path(savepath).joinpath('{}.xlsx'.format(savename))))
