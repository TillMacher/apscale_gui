import pandas as pd
import PySimpleGUI as sg

def settings_integrity(file):
    " Check the settings file for integrity and errors "

    ## open file
    xl = pd.ExcelFile(file)

    ## store errors
    lst = []

    ## loop through sheets
    for sheet in xl.sheet_names:
        ## check each sheet for errors
        if 'nan' in pd.read_excel(file, sheet_name=sheet).fillna('nan').values.tolist()[0]:
            lst.append(sheet)

    if lst != []:
        sg.PopupOK('Found missing data following sheet(s):', '; '.join(lst), 'Please add the missing information to continue!', title='Warning!')
        return False
