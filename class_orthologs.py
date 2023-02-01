# File Name: class_input.py
# Author: Hannah Nguyen
# Date: 9/22/22
# Description: Input class  , in information from json file

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
        self.nuc_rec
        # result from get_nec_seq
        self.nuc_seq


    def get_nuc_rec(self):
        print("Getting IPG nucelotide records for " + str(self.hit[0]['prot_id']) + "...")
        # Fetch the protein record
        print("Fetching the protein records...")

        for i in range(self.entrez_param["retry_number"]):

            try:
                records = Entrez.read(Entrez.efetch(db="protein", id=self.hit[0]['prot_id'], rettype="ipg",
                                                    retmode="xml"))
                time.sleep(self.entrez_param["sleep_time"])
                break

            except:
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                if i == (self.entrez_param["retry_number"] - 1):
                    print("\t\tCould not download record after " + str(self.entrez_param["retry_number"]) + " attempts")
                    print("\tNo gene records found for " + str(self.hit[0]['prot_id']))
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

        print("\t|~> Recording the gene records with priorities")
        if 'ProteinList' in records['IPGReport'].keys():
            for idprotein in records['IPGReport']['ProteinList']:
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
                    continue
        else:
            print("\tNo gene records found for " + self.hit[0]['prot_id'])
            return None

        # Finds the genome with the max p-score
        if len(genome_list) > 0:
            max_p_record = genome_list[0]
            for genome in genome_list:
                if genome["p_score"] > max_p_record["p_score"]:
                    max_p_record = genome

            print("\tReturning gene record for " + self.hit[0]['prot_id'] + ". p-score: " +
                  str(max_p_record["p_score"]))
            self.nuc_rec = max_p_record
            return max_p_record
        else:
            print("\tNo gene records found for " + self.hit[0]['prot_id'])
            return None


        def get_nuc_seq(self):
            if self.nuc_rec is None:
                print("\tRecord is not valid...")
                return None

            print("Get nucleotide sequence for " + str(nuc_rec['acc']) + "...")

            print("\tAdjusting start and stop positions...")

            if nuc_rec['strand'] == '+':
                s_start = int(nuc_rec['start']) - start_adj
                s_stop = int(nuc_rec['start']) + stop_adj
                s_strand = 1
            else:
                s_stop = int(nuc_rec['stop']) + start_adj
                s_start = int(nuc_rec['stop']) - stop_adj
                s_strand = 2

            if isolate_promoters:

                print("\t|~> Getting genbank record")

                # Fetch and read the annotated GenBank record

                for i in range(REQUEST_LIMIT):

                    try:

                        handle = Entrez.efetch(db="nuccore", id=nuc_rec['acc'], strand=s_strand,
                                               seq_start=s_start, seq_stop=s_stop,
                                               rettype='gbwithparts', retmode="XML")

                        genome_record = Entrez.read(handle, "xml")

                        time.sleep(SLEEP_TIME)
                        break

                    except:

                        print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                        if i == (REQUEST_LIMIT - 1):
                            print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

                print("\t|~> Parsing intervals for coding regions")
                # Find all coding regions in the returned GenBank sequence.
                coding_intervals = []

                sequence = genome_record[0]['GBSeq_sequence']

                for feature in genome_record[0]['GBSeq_feature-table']:
                    if feature['GBFeature_key'] == 'gene':
                        if "GBInterval_from" in feature['GBFeature_intervals'][0]:
                            coding_start = feature['GBFeature_intervals'][0]['GBInterval_from']
                            coding_end = feature['GBFeature_intervals'][0]['GBInterval_to']
                            coding_intervals.append((coding_start, coding_end))

                # The FASTA ID for the promoter sequence is in the following format:
                # p_NucleotideRecord
                print("\t|~> Returning promoter sequence")
                return_id = "p_" + str(nuc_rec['acc'])

                # Setting up the description for the FASTA record
                return_description = "Original_Query_Protein " + str(nuc_rec_in[0]) + " BLAST_Hit_Accession " + str(
                    nuc_rec_in[1])

                # If there is only one coding region in the selected sequence, then
                # the sequence is returned unmodified.
                if len(coding_intervals) == 1:
                    # Appends information to record description
                    print("\t\t|~>No-additional-coding-regions-found-Returning-full-sequence")
                    return SeqRecord(Seq.Seq(sequence), id=return_id, description=return_description)

                # If no coding intervals are indentified, None is returned.
                elif len(coding_intervals) == 0:
                    print("\t\t|~> No coding intervals found for record: " + str(nuc_rec['acc']) + ".")
                    return None

                # The start of the promoter is set to the start/end of the upstream gene
                # based on the directionality. ( --> --> or <-- -->)
                # If there was no downstream adjustment, then the last record in the list is upstream from the feature of interest.
                # If there was a downstream adjustment, then the second to last record in the list is upstream from the feature of interest.

                if stop_adj > 0:
                    promoter_start = max(int(coding_intervals[-2][0]), int(coding_intervals[-2][1]))
                else:
                    promoter_start = max(int(coding_intervals[-1][0]), int(coding_intervals[-1][1]))

                # Everything upstream of the promoter start is clipped off the
                # sequence and the substring is returned.
                return_seq = str(sequence[promoter_start:])
                # Appends information to record description
                print("\t\t|~>Successfully-clipped-off-upstream-coding-regions")
                return SeqRecord(Seq.Seq(return_seq), id=return_id, description=return_description)

            # If promoters aren't being isolated
            else:

                print("\t|~> Getting FASTA record")
                # Fetch the requested nucleotide sequence and return without any
                # modification.

                for i in range(REQUEST_LIMIT):

                    try:

                        handle = Entrez.efetch(db="nuccore", id=nuc_rec['acc'], strand=s_strand,
                                               seq_start=s_start, seq_stop=s_stop,
                                               rettype='fasta', retmode="txt")

                        time.sleep(SLEEP_TIME)
                        break

                    except:

                        print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                        if i == (REQUEST_LIMIT - 1):
                            print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

                print("\t|~> Returnig sequence")
                return SeqIO.read(handle, "fasta")
