import glob, os
import pandas as pd
import PySimpleGUI as sg
from demultiplexer import file_pairs
from pathlib import Path

def create_rename_sheet(pe_reads, folder, project):
    ## collect all files
    input_folder = Path(project).joinpath(folder, 'data', '*.fastq.gz')
    input = glob.glob(str(input_folder))

    if pe_reads == 'True':
        ## option 1: paired-end reads
        pairs = file_pairs.main(input)
        lines = []
        for pair in pairs:
            name = '_'.join(Path(pair[0]).name.split('_')[:-1])
            line = pair + [name]
            lines.append(line)

        df = pd.DataFrame(lines, columns = ['R1 file', 'R2 file', 'New name'])
        df.to_excel(Path(project).joinpath('rename_sheet.xlsx'), index=False)
    else:
        ## option 2: single files
        lines = []
        for file in input:
            name = '_'.join(Path(file).name.split('_')[:-1])
            line = [file] + [name]
            lines.append(line)

        df = pd.DataFrame(lines, columns = ['File', 'New name'])
        df.to_excel(Path(project).joinpath('rename_sheet.xlsx'), index=False)


def rename_files(pe_reads, project):
    " Rename paired-end files "

    rename_sheet = Path(project).joinpath('rename_sheet.xlsx')
    df = pd.read_excel(rename_sheet)

    if pe_reads == 'True':
        ## option 1: paired-end reads
        for line in df.values.tolist():
            ## collect values
            r1_file = line[0]
            r2_file = line[1]
            new_r1_file = line[2] + '_R1.fastq.gz'
            new_r2_file = line[2] + '_R2.fastq.gz'
            rename_r1 = Path(r1_file).parents[0].joinpath(new_r1_file)
            rename_r2 = Path(r1_file).parents[0].joinpath(new_r2_file)

            ## rename files
            os.rename(r1_file, rename_r1)
            os.rename(r2_file, rename_r2)
    else:
        ## option 2: single files
        for line in df.values.tolist():
            ## collect values
            file = line[0]
            new_file = line[1] + '.fastq.gz'
            rename = Path(file).parents[0].joinpath(new_file)
            ## rename files
            os.rename(file, rename)





















#
