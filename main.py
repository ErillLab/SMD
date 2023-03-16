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

    '''
    # do blast search (either regular, clustered, or hierarchical)
    hits = []
    if the_input.blast_param["hierarchical_taxon_level"] != "None":
        hits = the_input.hierarchical_BLAST_search()
    else:
        hits = the_input.reg_BLAST_search()

    '''
    #hits = [({'prot_id': 'WP_011079176'}, 'WP_011079176.1')]
    hits = [({'prot_id': 'WP_011079176'}, 'WP_272923659.1')]
    #hits = [({'prot_id': 'WP_011079176'}, 'WP_039469417.1')]
    #hits = [({'prot_id': 'WP_011079176'}, 'WP_053542922.1')]

    # treats each hit from the blast search as a potential ortholog and adds to list of orthologs
    orthologs = []
    for hit in hits:
        orthologs.append(theOrtholog(the_input, hit))

    # prepares an ortholog dictionary for later
    orthologs_dict = {}
    #goes through orthologs list to add to ortho dictionary and get records/sequences
    for potential_olog in orthologs:
        orthologs_dict[potential_olog.hit] = ""
        potential_olog.get_nuc_rec() #gets nucleotide record for each ortholog
        potential_olog.get_nuc_seq() #get nucleotide sequence from the record for each ortholog

    # goes through list of orthologs and does popping algorithm that marks sequences that are too similar
    for first_ortholog in range(len(orthologs)):
        compared_ortholog = first_ortholog+1 # compares first ortholog with ortholog after it
        while compared_ortholog < len(orthologs):
            # gets percent similarity and marks in dictionary if above the threshold
            per_similarity = orthologs[first_ortholog].percent_similarity(orthologs[compared_ortholog])
            if per_similarity < 8: #checks if ortholog is less
                orthologs_dict[compared_ortholog.hit] = "DELETE"
            compared_ortholog += 1

    # goes through ortholog dictionary to finalize which orthologs can stay
    for ortholog_key in orthologs_dict:
        if orthologs_dict[ortholog_key] == "DELETE":
            for ortholog in orthologs:
                if (ortholog.hit == ortholog_key):
                    orthologs.remove(ortholog)
