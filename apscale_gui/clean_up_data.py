import os, glob
from pathlib import Path
from datetime import datetime

def clean_up_data(path_to_outdirs):
    " Clean-up the directory to safe storage space "

    ## calcualte the storage before the clean-up
    size_before = sum(f.stat().st_size for f in path_to_outdirs.glob('**/*') if f.is_file())
    ## only delete data in excess folders
    folders = ['3_PE_merging', '4_primer_trimming', '5_quality_filtering']

    ## remove data
    for folder in folders:
        Path
        [os.remove(i) for i in glob.glob(str(path_to_outdirs) + '/' + folder + '/data/*.gz')]
        [os.remove(i) for i in glob.glob(str(path_to_outdirs) + '/' + folder + '/data/*.fasta')]
        [os.remove(i) for i in glob.glob(str(path_to_outdirs) + '/' + folder + '/data/*.fastq')]

    ## calculate safed storage space
    size_after = sum(f.stat().st_size for f in path_to_outdirs.glob('**/*') if f.is_file())
    size_safed = round((size_before - size_after) / 1000000000, 2)

    ## collect time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")

    ## print results
    if size_safed < 1:
        print(current_time + ': Safed ' + str(size_safed*1000) + ' MB.')
    else:
        print(current_time + ': Safed ' + str(size_safed) + ' GB.')
