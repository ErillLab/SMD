# File Name: class_input.py
# Author: Hannah Nguyen
# Date: 9/22/22
# Description: Input class takes in information from json file

from Bio import Seq, Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from ete3 import NCBITaxa
from os.path import exists
import json
import time

class theInput:
    def __init__(self, json_file):
        # content from json file input
        self.json_file = json_file
        self.ref_proteins = []
        self.target_clade = 0
        self.blast_param = {}
        self.entrez_param = {}


    def load_json(self):
        print("Loading json file... ")
        # checks if file actually exists
        file_exists = exists(self.json_file)
        if file_exists == False:
            raise Exception("The file you inputted does not exist. Please try again.")
        
        # opens and loads json file as class variables
        open_file = open(self.json_file)
        data = json.load(open_file)
        self.ref_proteins = data["reference_proteins"]
        self.target_clade = data["target_clade"]
        self.blast_param = data["BLAST_parameters"]
        self.entrez_param = data["entrez_parameters"]
        Entrez.email = self.entrez_param["email"]
        Entrez.apikey = self.entrez_param["NCBI_API_key"]


    def reg_BLAST_search(self):
        # check the length of the input_records.
        if len(self.ref_proteins) < 1:
            raise Exception("Need at least one protein record to conduct BLAST search.")

        # The final list which we will add the BLAST hits to
        hits = []

        # Gets the accession numbers for all the hits in the BLAST search and adds to hits[] with every unique record.
        for ref_protein in self.ref_proteins:
            # Fetches the protein record based off the accession number
            print("Getting protein record...")

            print('Hello!')
            for i in range(self.entrez_param["retry_number"]):

                try:
                    handle = Entrez.efetch("protein", id=ref_protein["prot_id"], rettype="fasta",
                                           retmode="text")
                    time.sleep(self.entrez_param["sleep_time"])
                    break

                except:
                    print("\tNCBI exception raised on attempt " + str(i) + "\n\treattempting now...")
                    if i == (self.entrez_param["retry_number"] - 1):
                        print("\tCould not download record after " + str(self.entrez_param["retry_number"]) + " attempts")

            # Fetches the protein sequence (from record) to be used in the BLAST search
            print("Getting protein sequence...\n")
            input_seq = (SeqIO.read(handle, "fasta")).seq

            # Performs the BLAST using parameters from json file
            print("Performing BLAST search: " + str(ref_protein["prot_id"]) + ", this may take a while...")
            if self.target_clade is not None:
                # Holds the parameter for the target clade
                taxon = "txid" + str(self.target_clade) + "[orgn]"

                for i in range(self.entrez_param["retry_number"]):

                    try:
                        result_handle = NCBIWWW.qblast("blastp", "nr", input_seq,
                                                       entrez_query=taxon, expect=self.blast_param["e_value"],
                                                       hitlist_size=self.blast_param["max_results"])
                        # Parses the resulting hits as a list
                        print("\tGetting records")
                        blast_records = list(NCBIXML.parse(result_handle))
                        time.sleep(self.entrez_param["sleep_time"])
                        break

                    except:
                        print("\tNCBI exception raised on attempt " + str(i) + "\n\treattempting now...")
                        if i == (self.entrez_param["retry_number"] - 1):
                            print("\tCould not download record after " + str(self.entrez_param["retry_number"]) + " attempts")

            else:

                for i in range(self.entrez_param["retry_number"]):

                    try:
                        result_handle = NCBIWWW.qblast(program="blastp", database=self.blast_param["database"],
                                                       sequence=input_seq, expect=self.blast_param["e_value"],
                                                       hitlist_size=self.blast_param["max_results"])
                        time.sleep(self.entrez_param["sleep_time"])
                        # Parses the resulting hits as a list
                        print("\tGetting records")
                        blast_records = list(NCBIXML.parse(result_handle))
                        break

                    except:
                        print("\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                        if i == (self.entrez_param["retry_number"] - 1):
                            print("\tCould not download record after " + str(self.entrez_param["retry_number"]) + " attempts")

            # Adds each unique accession number to hits[]
            for record in blast_records[0].alignments:
                curr_hit_rec = record.hit_id.split('|')[-2]
                print("\t\tAnalyzing hit " + str(curr_hit_rec))

                for hit in record.hsps:

                    # Checks if hit meets the minimum coverage if provided
                    if self.blast_param["query_coverage"]:
                        cov = (hit.query_end - hit.query_start + 1) / (len(input_seq))

                        if (cov >= self.blast_param["query_coverage"]):

                            # Check if the hit is already in the return list
                            if len(hits) == 0:
                                print("\t\tAdding first hit (Coverage = " + str(cov) + "): " + str(curr_hit_rec))
                                hits.append((ref_protein, curr_hit_rec))

                            elif (not (curr_hit_rec in list(zip(*hits))[1])):
                                print("\t\tAdding hit (Coverage = " + str(cov) + "): " + str(curr_hit_rec))
                                hits.append((ref_protein, curr_hit_rec))

                        # Prints error if the minimum coverage is not met
                        else:
                            print("\t\tHit did not meet coverage (Coverage = " + str(cov) + ") requirement: " + str(
                                curr_hit_rec))
                    else:
                        # Check if the hit is already in the return list
                        if len(hits) == 0:
                            print("\t\tAdding first hit: " + str(curr_hit_rec))
                            hits.append((ref_protein, curr_hit_rec))

                        elif (not (curr_hit_rec in list(zip(*hits))[1])):
                            print("\t\tAdding hit: " + str(curr_hit_rec))
                            hits.append((ref_protein, curr_hit_rec))

        print("\tReturning " + str(len(hits)) + " unique hits\n")
        return hits


    def clustered_BLAST_search(self):
        save_database = self.blast_param["database"]
        self.blast_param["database"] = "clustered_nr"

        hits = self.reg_BLAST_search()
        self.blast_param["database"] = save_database

        return hits


    def hierarchical_BLAST_search(self):
        # convert target clade to name
        clade_dict = NCBITaxa().get_taxid_translator([self.target_clade])
        clade_name = clade_dict[self.target_clade]

        # get all descendants of the reference protein
        descendants = NCBITaxa().get_descendant_taxa(clade_name, intermediate_nodes=True)
        desired_descendants = []

        # getting all descendants under taxa and putting them into list
        for i in range(len(descendants)):
            rank = ""
            rank_dict = NCBITaxa().get_rank([descendants[i]])

            for key in rank_dict:
                rank = rank_dict[key]

            if str(rank) == (self.blast_param["hierarchical_taxon_level"]):
                desired_descendants.append(descendants[i])

        # saving original parameters of the search
        target_clade = self.target_clade

        # BLAST search for each group in the desired descendants
        original_hits = []
        for i in range(len(desired_descendants)):
            self.target_clade = desired_descendants[i]
            original_hits.append(self.reg_BLAST_search())

        # make the target clade back to what it was originally (before search)
        self.target_clade = target_clade

        # deleting duplicates from list
        hits = []
        for i in range(len(original_hits)):
            for j in range(len(original_hits[i])):
                if original_hits[i][j] not in hits:
                    hits.append(original_hits[i][j])

        return hits