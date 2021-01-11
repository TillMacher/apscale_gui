from glob import glob
from difflib import SequenceMatcher
from pathlib import Path

## function that accepts a list of path strings and finds matching file names
## e.g. ('1_r1.fastq.gz', 1_r2.fastq.gz)
## the list with both files / paths removed will be returned to improve runtime
## recursion is run through the main function until all pairs are found
def path_matcher(path_list):
    ## number of elements to cut off from the back of the file name
    cut_off = 1

    ## collect SequenceMatcher ratios here. If there are only two items left
    ## in this list two items have the same (maximum) ratio --> match
    max_ratios = []

    ## list to return at the end of the function. e.g. full list with the searched
    ## pathes removed
    return_list = path_list.copy()

    ## generate Path object from the path strings and remove the file extension
    ## e.g. '1_r1.fastq.gz' --> '1_r1'
    path_list = [Path(path).with_suffix('').with_suffix('').stem for path in path_list]

    ## compute the max cut_off that will stop the loop
    ## if all strings are empty they all will have the same ratio --> infinite loop
    max_cut_off = max([len(filename) for filename in path_list])

    ## calculate the matching ratio of all file names against each other
    while len(max_ratios) != 2:
        ## cut of a letter from the file name each iteration until a matching pair is found
        ## if max cut off is reached there is no match for the first item
        if cut_off < max_cut_off:
            file_names = [(filename[:-cut_off], SequenceMatcher(None, path_list[0][:-cut_off], filename[:-cut_off]).ratio()) for filename in path_list]
            max_ratios = [element[0] for element in file_names if element[1] == max([element[1] for element in file_names])]

            ## find the indices of the max ratio items
            indices = [i for i, x in enumerate(file_names) if x[0] == list(set(max_ratios))[0]]

            ## return the original pathes of the matching files
            ## sort to have the forward reads first in the list
            matched = sorted([return_list[i] for i in indices])
            remainder = [i for i in return_list if i not in matched]

            cut_off += 1
        else:
            ## if no match can be generated return the search item without match
            ## also return the remainder of the patch list
            matched = path_list[0]
            remainder = [i for i in return_list if i != matched]
            cut_off = 1
        break

    ## return match and the remainder of the path list
    return matched, remainder

## main function to run the algorithm
def main(path_list):

    ## collect result here
    result = []

    ## run until all matches are found
    while path_list:
        match, path_list = path_matcher(path_list)
        result.append(match)

    ## return a list of all matching file pairs. Non matching pairs are just
    ## appended to the list and have to be checked later
    return result
