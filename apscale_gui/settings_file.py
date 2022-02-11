import pandas as pd
import PySimpleGUI as sg

def load_settings(file):
    " Load the settings from the project's settings file to add them to the GUI "

    ## open file
    xl = pd.ExcelFile(file)

    ## store settings
    settings = []

    ## loop through sheets
    for sheet in xl.sheet_names:
        ## check each sheet for errors
        df = pd.read_excel(file, sheet_name=sheet).fillna('')
        settings = settings + [str(df[i][0]) for i in df.columns.values.tolist()]

    return settings

def apply_settings(file, settings):
    " Write settings from GUI to settings file "

    ## convert to in whenver possible
    settings_converted = []
    for i in settings:
        try:
            settings_converted.append(int(i))
        except:
            settings_converted.append(i)

    ## collect settings from the GUI
    df0 = pd.read_excel(file, sheet_name='0_general_settings')
    df1 = pd.DataFrame([[settings_converted[0], settings_converted[1], settings_converted[2]]], columns=['maxdiffpct', 'maxdiffs', 'minovlen'])
    df2 = pd.DataFrame([[settings_converted[3], settings_converted[4], settings_converted[5]]], columns=["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", 'anchoring'])
    df3 = pd.DataFrame([[settings_converted[6], settings_converted[7], settings_converted[8]]], columns=['maxEE', 'min length', 'max length'])
    df4 = pd.DataFrame([[settings_converted[9]]], columns=['min size to pool'])
    df5 = pd.DataFrame([[settings_converted[10], settings_converted[11], settings_converted[12]]], columns=['pct id', 'to excel', 'to parquet'])
    df6 = pd.DataFrame([[settings_converted[13], settings_converted[14], settings_converted[15], settings_converted[16]]], columns=['alpha', 'minsize', 'to excel', 'to parquet'])


    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(file, engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    df0.to_excel(writer, sheet_name='0_general_settings', index=False)
    df1.to_excel(writer, sheet_name='3_PE_merging', index=False)
    df2.to_excel(writer, sheet_name='4_primer_trimming', index=False)
    df3.to_excel(writer, sheet_name='5_quality_filtering', index=False)
    df4.to_excel(writer, sheet_name='6_dereplication_pooling', index=False)
    df5.to_excel(writer, sheet_name='7_otu_clustering', index=False)
    df6.to_excel(writer, sheet_name='8_denoising', index=False)


    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

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
        sg.PopupOK('Found missing data following sheet(s):', '; '.join(error_list), 'Please add the missing information to continue!', title='Warning!')
        return False

    return sg.PopupOKCancel(' '.join(settings_list), title='Settings')
