import pandas as pd
import PySimpleGUI as sg

def settings_integrity(file):
    " Check the settings file for integrity and errors "

    ## open file
    xl = pd.ExcelFile(file)

    ## store errors
    error_list = []

    ## store settings
    settings_list = ['Your settings:\n']

    ## loop through sheets
    for sheet in xl.sheet_names:
        ## check each sheet for errors
        df = pd.read_excel(file, sheet_name=sheet).fillna('nan')
        if 'nan' in df.values.tolist()[0]:
            error_list.append(sheet)

        settings_list = settings_list + [str(i) + ' = ' + str(df[i][0]) + '\n' for i in df.columns.values.tolist()]

    if error_list != []:
        sg.PopupOK('Found missing data following sheet(s):', '; '.join(lst), 'Please add the missing information to continue!', title='Warning!')
        return False

    return sg.PopupOKCancel(' '.join(settings_list), title='Settings')
