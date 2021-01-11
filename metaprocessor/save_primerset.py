import datetime
from pathlib import Path

## function to save a primerset for the demultiplexer module
def save_primerset(generator_dict, savepath, print_handle):

    ## save the savename for later, remove from dict to loop through the integer keys
    savename = generator_dict['_PRIMERSETNAME_']
    del generator_dict['_PRIMERSETNAME_']

    ## check if a primerset with the chosen name already exists
    ## if so let the user chose another name
    if Path(savepath.joinpath('{}.txt'.format(savename))).is_file():
        print_handle.print('{}: This primerset already exists. Please choose another name.'.format(datetime.datetime.now().strftime("%H:%M:%S")))
    else:
        with open(savepath.joinpath('{}.txt'.format(savename)), 'w+') as output:

            ## count the valid primers for output
            valid_primers = 0
            ## since the data is stored in pairs only len(dict) / 2 iterations are needed
            for i in sorted(range(int(len(generator_dict) / 2))):
                ## the first half are always the forward primers
                ## the second half are then the reverse primers
                ## direction needed for later
                if i < len(generator_dict) / 4:
                    dir = 'fwd'
                else:
                    dir = 'rev'

                ## only write the primer files if valid input in both input fields
                if generator_dict[i * 2] != '' and generator_dict[i * 2 + 1] != '':
                    valid_primers += 1
                    output.write('{},{},{}\n'.format(generator_dict[i * 2], generator_dict[i * 2 + 1].upper(), dir))

        ## user output
        print_handle.print('{}: {} valid primers were collected from the input.'.format(datetime.datetime.now().strftime("%H:%M:%S"), valid_primers))
        print_handle.print('{}: Primerset was saved at:\n {}.'.format(datetime.datetime.now().strftime("%H:%M:%S"), savepath.joinpath('{}.txt'.format(savename))))
