#!/usr/bin/env python3

from collections import namedtuple, defaultdict

BASIC_NT = ['A','C','G','T']
UNIQUE_NT = ['M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N']
ALPHABET = {
    "A" : ['A'],
    "C" : ['C'],
    "T" : ['T'],
    "G" : ['G'],
    "M" : ['A','C'],
    "R" : ['A','G'],
    "W" : ['A','T'],
    "S" : ['C','G'],
    "Y" : ['C','T'],
    "K" : ['G','T'],
    "V" : ['A','C','G'],
    "H" : ['A','C','T'],
    "D" : ['A','G','T'],
    "B" : ['C','G','T'],
    "N" : ['A','C','G','T']
}

ALPHABETS = list(ALPHABET.keys())
FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s' 

## Default Settings
MULTIPLIER = 2
KMER_LENGTH = 3 * MULTIPLIER    
MISMATCH_PER_KMER = 1 * MULTIPLIER   
SEED_LENGTH = 12 
MAX_MM = 8

def create3MerList():
    import itertools
    perm3mer = [''.join(p) for p in itertools.product(BASIC_NT, repeat=KMER_LENGTH)]
    perm3mer.sort()
    return perm3mer

KMER_LIST = create3MerList()

def myers(P,T,k):
    B = defaultdict(int)
    i = 1
    for c in P:
        B[c] = B[c] | i 
        i = i << 1
    
    m = len(P)
    VP = (1 << m) - 1 #set m bits
    VN = 0
    score = m

    for i, c in enumerate(T):
        Eq = B[c]
        Xv = Eq | VN 
        Xh = (((Eq & VP) + VP) ^ VP) | Eq
        HP = VN | ~(Xh | VP)
        HN = VP & Xh
        if HP & (1 <<(m-1)):
            score += 1
        elif HN & (1 <<(m-1)):
            score -= 1
        if score < k:
            return 0

        HP = HP << 1 
        VP = (HN << 1) | ~(Xv|HP)
        VN = HP & Xv
    
    return 1
    #return 0 if score < k else 1


##-------------------------------------------------------------------------------------------------
## Read FASTQ | FASTA
def readfq(fp, gzipped = False): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                l = l.decode() if gzipped else l
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            l = l.decode() if gzipped else l
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                l = l.decode() if gzipped else l
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

##-------------------------------------------------------------------------------------------------
## Multiprocessing Utils

def partition(lst, n):
    division = len(lst) / n
    return [lst[round(division * i):round(division * (i + 1))] for i in range(n)]

def ranges(N, nb):
    if isinstance(N, list):
        tmp = len(N)
    else:
        N
    temp_list = [(r.start, r.stop) for r in partition(range(N), nb)]
    for i in temp_list: 
        yield i[0],i[1]

##-------------------------------------------------------------------------------------------------
## General Utils

def reverseC(seq):
    RT = {'A':'T','C':'G','G':'C','T':'A'}
    reverseComplement = ''
    for i in seq:
        nt = RT.get(i)
        reverseComplement += nt
    return reverseComplement[::-1]