from Bio import Align

def align_sequence():
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
    alignments = aligner.align("GTACACACTC", "ATACAATATAT")
    alignment = alignments[0]
    matches = alignment.counts().identities
    mismatches = alignment.counts().mismatches
    gaps = alignment.counts().gaps
    print(alignment)
    print(matches)
    print(mismatches)
    print(gaps)

    '''
    aligner2 = Align.PairwiseAligner()
    aligner2.match_score = 0
    aligner2.mismatch_score = 1
    aligner2.gap_score = 0
    alignments = aligner2.align("GTACACACTC", "ATACAATATAT")
    alignment = alignments[0]
    print(alignment)
    print(alignment.score)
    '''


if __name__ == "__main__":
    align_sequence()