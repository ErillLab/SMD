# File Name: class_input.py
# Author: Hannah Nguyen
# Date: 9/22/22
# Description: orthologs class

from Bio import Seq, Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import time

class theOrtholog:
    def __init__(self, the_input, hit):
        # content from the hit
        self.hit = hit
        # content from the input
        self.json_file = the_input.json_file
        self.ref_proteins = the_input.ref_proteins
        self.target_clade = the_input.target_clade
        self.blast_param = the_input.blast_param
        self.entrez_param = the_input.entrez_param
        # result from get_nuc_rec
        self.nuc_rec = []
        # result from get_nec_seq
        self.nuc_seq = []


    def get_nuc_rec(self):
        print("Getting IPG nucelotide records for " + str(self.hit[1]) + "...")
        print("\tFetching the protein records...")
        # Fetches protein records for each hit
        for i in range(self.entrez_param["retry_number"]):
            try:
                handle = Entrez.efetch(db="protein", id=self.hit[1], rettype="ipg",
                                                    retmode="xml")
                records = Entrez.read(handle)
                time.sleep(self.entrez_param["sleep_time"])
                break
            except:
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                time.sleep(self.entrez_param["sleep_time"])
                if i == (self.entrez_param["retry_number"] - 1):
                    print("\t\tCould not download record after " + str(self.entrez_param["retry_number"]) + " attempts")
                    print("\tNo gene records found for " + str(self.hit[1]))
                    return None

        # The priority scores for the types of gene records available
        # NC and ac = ref seq sequences (prefer those bc complete genomes)
        # AE and CP = complete not as good
        # everything else not as good
        p_scores = {"NC_": 7, "AC_": 7,
                    "AE": 6, "CP": 6, "CY": 6,
                    "NZ_": 5, "NT_": 5, "NW_": 5,
                    "AAAA-AZZZ": 4,
                    "U": 3, "AF": 3, "AY": 3, "DQ": 3}

        # Gene records are appended to this list
        genome_list = []
        print("\tRecording the gene records with priorities...")
        if 'ProteinList' in records['IPGReport'].keys():
            for idprotein in records['IPGReport']['ProteinList']:
                # focuses on coding regions in the record
                if 'CDSList' in idprotein.keys():
                    for cds in idprotein['CDSList']:
                        cds_acc = cds.attributes['accver']
                        cds_start = cds.attributes['start']
                        cds_stop = cds.attributes['stop']
                        cds_strand = cds.attributes['strand']
                        cds_scr = 0
                        # assign priority
                        for key in p_scores:
                            if cds_acc.startswith(key):
                                cds_scr = p_scores[key]
                        # create and append record
                        cds_rec = {'acc': cds_acc, 'start': cds_start,
                                   'stop': cds_stop, 'strand': cds_strand,
                                   'p_score': cds_scr}
                        genome_list.append(cds_rec)

        else:
            print("\tNo gene records found for " + self.hit[1])
            return None

        # Finds the genome with the max p-score
        if len(genome_list) > 0:
            max_p_record = genome_list[0]
            for genome in genome_list:
                if genome["p_score"] > max_p_record["p_score"]:
                    max_p_record = genome

            print("\tReturning gene record for " + self.hit[1] + ". p-score: " +
                  str(max_p_record["p_score"]))
            self.nuc_rec = max_p_record
            return max_p_record

        else:
            print("\tNo gene records found for " + self.hit[1])
            return None


    def get_nuc_seq(self):
        if self.nuc_rec is None:
            print(self.nuc_rec)
            print("\tRecord is not valid")
            return None

        print("Get nucleotide sequence for " + str(self.nuc_rec['acc']) + "...")
        # grabs upstream region some parameter base pairs away from start to check for operons
        if self.nuc_rec['strand'] == '+':
            s_start = int(self.nuc_rec['start']) - self.blast_param["base_pairs_to_grab"]
            if s_start <= 0:
                s_start = 1
            s_stop = int(self.nuc_rec['stop'])
            s_strand = 1
        else:
            s_stop = int(self.nuc_rec['stop'])
            s_start = int(self.nuc_rec['start']) + self.blast_param["base_pairs_to_grab"]
            s_strand = 2

        print("\t\tGetting genbank record...")

        # Fetch and read the annotated GenBank record
        for i in range(self.entrez_param["retry_number"]):
            try:
                handle = Entrez.efetch(db="nuccore", id=self.nuc_rec['acc'], strand=s_strand,
                                       seq_start=s_start, seq_stop=s_stop,
                                       rettype='gbwithparts', retmode="text")
                time.sleep(self.entrez_param["sleep_time"])
                genome_record = SeqIO.read(handle, "genbank")
                break
            except:
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                time.sleep(self.entrez_param["sleep_time"])
                if i == (self.entrez_param["retry_number"] - 1):
                    print("\t\tCould not download record after " + str(self.entrez_param["retry_number"]) + " attempts")

        sequence = genome_record.seq

        print(genome_record.features)

'''
        print("\t\tParsing intervals for coding regions...")
        # Find all coding regions in the returned GenBank sequence.
        coding_intervals = []

        sequence = genome_record[0]['GBSeq_sequence']

        for feature in genome_record[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'gene':
                if feature['GBFeature_location'].startswith("complement"):
                    strand = 2
                else:
                    strand = 1
                if "GBInterval_from" in feature['GBFeature_intervals'][0]:
                    coding_start = feature['GBFeature_intervals'][0]['GBInterval_from']
                    coding_end = feature['GBFeature_intervals'][0]['GBInterval_to']
                    coding_intervals.append((int(coding_start), int(coding_end), strand))

        print(coding_intervals)

        # setting up return id
        return_id = "p_" + str(self.nuc_rec['acc'])

        #Setting up the description for the FASTA record
        return_description = "Original_Query_Protein " + str(self.hit[0]) + " BLAST_Hit_Accession " + str(self.hit[1])

'''
        # If there is only one coding region in the selected sequence, then
        # the sequence is returned unmodified.
        if len(coding_intervals) == 1:
            # Appends information to record description
            print("\t\tOnly one coding regions found. Returning full sequence...")
            return SeqRecord(Seq.Seq(sequence), id=return_id, description=return_description)

        # If no coding intervals are indentified, None is returned.
        elif len(coding_intervals) == 0:
            print("\t\t|~> No coding intervals found for record: " + str(self.nuc_rec['acc']) + ".")
            return None

        the_end = 0 #for now
        upstream_regions = []
        upstream_sequences = []

        # LATER make sure to implement a way to add 10,000 more when you dont reach the end
        for i in range(len(coding_intervals)):
            if i > 0:
                print("poop ye " + str(coding_intervals[i][0] - coding_intervals[i-1][1]))
                if coding_intervals[i][0] - coding_intervals[i-1][1] <= self.blast_param["min_intergenic_distance"]:
                    if coding_intervals[i][0] - coding_intervals[i-1][1] > 0:
                        upstream_regions.append((coding_intervals[i-1][1], coding_intervals[i][0]))
                elif coding_intervals[i][2] == coding_intervals[i-1][2]:
                    if coding_intervals[i][1] - coding_intervals[i-1][0] > 0:
                        upstream_regions.append((coding_intervals[i-1][1], coding_intervals[i][0]))
                else:
                    the_end = coding_intervals[i][0]
                    break

        print(upstream_regions)
        # get sequences of each upstream region
        for i in range(len(upstream_regions)):
            print("poop")
            upstream_sequence = sequence[upstream_regions[i][0]:upstream_regions[i][1]]
            upstream_sequences.append(SeqRecord(Seq.Seq(upstream_sequence), id=return_id, description=return_description))

        print(upstream_sequences)
        return upstream_sequences