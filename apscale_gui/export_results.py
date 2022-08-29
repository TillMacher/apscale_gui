import os, glob, datetime
from pathlib import Path
import shutil
import PySimpleGUI as sg

def export_results(path_to_outdirs):
    " Export all relevant files to a gzip archive "

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': Compressing all result files and folders.')

    ## collect all files and folders to export
    xlsx_files = glob.glob('{}/*.xlsx'.format(str(path_to_outdirs)))
    folders = [Path('{}/{}'.format(path_to_outdirs,i)) for i in ['7_otu_clustering', '8_denoising', '9_lulu_filtering', '10_local_BLAST', '11_NCBI_BLAST']]

    ## copy them to a new folder
    project_name = path_to_outdirs.stem
    directory = Path('{}/{}_export'.format(path_to_outdirs, project_name))
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

    for f in xlsx_files:
        try:
            f_new = Path('{}/{}'.format(directory, Path(f).name))
            shutil.copyfile(f, f_new)
        except FileExistsError:
            pass

    for f in folders:
        try:
            f_new = Path('{}/{}'.format(directory, Path(f).name))
            shutil.copytree(f, f_new)
        except FileExistsError:
            pass

    ## compress the folder
    shutil.make_archive(directory, 'zip', directory)

    ## remove export folder
    shutil.rmtree(directory)

    print(datetime.datetime.now().strftime('%H:%M:%S') + ': All results were exported to a new .zip archive!.')
    sg.Popup('All results were exported to a new .zip archive!', title='Finished')