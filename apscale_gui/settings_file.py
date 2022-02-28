import pandas as pd
import PySimpleGUI as sg

def load_settings(file):
    " Load the settings from the project's settings file to add them to the GUI "

    ## open file
    xl = pd.ExcelFile(file)

    ## store settings
    settings_dict = {}

    ## loop through sheets
    for sheet in xl.sheet_names:
        ## check each sheet for errors
        df = pd.read_excel(file, sheet_name=sheet).fillna('')
        for value in df.columns.tolist():
            settings_dict[value] = str(df[value].tolist()[0])

    return settings_dict

def apply_settings(file, settings):
    " Write settings from GUI to settings file "

    settings_keys = ['settings_clustering_to_excel',
                     'settings_maxdiffpct',
                     'settings_maxdiffs',
                     'settings_minovlen',
                     'settings_p5_primer',
                     'settings_p7_primer',
                     'settings_anchoring',
                     'settings_maxEE',
                     'settings_min_length',
                     'settings_max_length',
                     'settings_derep_minsize',
                     'settings_pct_id',
                     'settings_alpha',
                     'settings_min_size',
                     'settings_min_sim',
                     'settings_min_rel_coo',
                     'settings_min_ratio']

    ## convert to in whenver possible
    settings_converted = {}
    for key,value in settings.items():
        try:
            settings_converted[key] = int(value)
        except Exception:
            settings_converted[key] = value

    ## collect settings from the GUI
    df0 = pd.read_excel(file, sheet_name='0_general_settings')
    df1 = pd.DataFrame([[settings_converted['settings_maxdiffpct'], settings_converted['settings_maxdiffs'], settings_converted['settings_minovlen']]], columns=['maxdiffpct', 'maxdiffs', 'minovlen'])
    df2 = pd.DataFrame([[settings_converted['settings_p5_primer'], settings_converted['settings_p7_primer'], settings_converted['settings_anchoring']]], columns=["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", 'anchoring'])
    df3 = pd.DataFrame([[settings_converted['settings_maxEE'], settings_converted['settings_min_length'], settings_converted['settings_max_length']]], columns=['maxEE', 'min length', 'max length'])
    df4 = pd.DataFrame([[settings_converted['settings_derep_minsize']]], columns=['min size to pool'])
    df5 = pd.DataFrame([[settings_converted['settings_pct_id'], settings_converted['settings_clustering_to_excel']]], columns=['pct id', 'to excel'])
    df6 = pd.DataFrame([[settings_converted['settings_alpha'], settings_converted['settings_min_size'], settings_converted['settings_clustering_to_excel']]], columns=['alpha', 'minsize', 'to excel'])
    df7 = pd.DataFrame([[settings_converted['settings_min_sim'], settings_converted['settings_min_rel_coo'], settings_converted['settings_min_ratio'], settings_converted['settings_clustering_to_excel']]], columns=['minimum similarity', 'minimum relative cooccurence', 'minimum ratio', 'to excel'])


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
    df7.to_excel(writer, sheet_name='9_lulu_filtering', index=False)


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
