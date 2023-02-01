# File Name: main.py
# Author: Hannah Nguyen
# Date: 9/22/22
# Description: Runs the main smd pipeline

# importing classes into main
from class_input import theInput
from class_orthologs import theOrtholog
import sys

if __name__ == "__main__":

    # gets json file name
    if len(sys.argv) != 2:
        raise Exception("The file name was not added as an argument. Please try again.")
    json_name = sys.argv[1]

    # inputs information from json file
    the_input = theInput(json_name)

    # loads json file
    the_input.load_json()

    # do blast search (either regular, clustered, or hierarchical)
    hits = []
    if the_input.blast_param["hierarchical_taxon_level"] != "None":
        hits = the_input.hierarchical_BLAST_search()
    else:
        hits = the_input.reg_BLAST_search()

    # treats each hit from the blast search as a potential ortholog and adds to list of orthologs
    orthologs = []
    for i in range(len(hits)):
        for j in range(len(hits[i])):
            orthologs.append(theOrtholog(the_input, hits[i][j]))






