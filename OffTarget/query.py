#!/usr/bin/env python3

import os
import argparse
import math
import numpy as np
import logging
import logging.handlers
import ctypes
import sys
import time
import multiprocessing as mp
from bitarray.util import ba2int

from OffTarget.wavelet import *
from OffTarget.automata import *

AUTOMATA = Matcher(KMER_LIST)
PYTHON_FILE = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)
alive = mp.Value(ctypes.c_bool, lock=False)
alive.value = True
MAX_LENGTH = 50

def openReference(referenceList, fasta):   
    with open(fasta, 'r') as fin:
        count1 = 0
        count2 = 0
        for name, seq, _ in readfq(fin):
            length = len(seq)
            if length >= (1<<16):
                count1 += 1
                fragment_c = (length // (1 << 16)) + 1
                logger.debug(f"Fragmented: {count1}, {name}, {length}, {fragment_c}")
                for c, (start, stop) in enumerate(ranges(length, fragment_c)):
                    if start != 0: start -= MAX_LENGTH
                    referenceList.append(Reference(name, seq[start:stop], start, c))
                    Reference.totalCount += 1
            else:
                count2 += 1
                referenceList.append(Reference(name, seq))
                Reference.totalCount += 1

class Query:
    totalCount = 0

    def __init__(self, ID, seq):
        self.id = ID
        self.length = len(seq)
        self.forward = QuerySubclass(seq)
        self.reverse = QuerySubclass(reverseC(seq), reverse = True)
        self.unique_NT = False #self._unique_NT(seq.upper())
        self.processed = False

    def preprocess(self):
        if not self.processed:
            self.forward.prepareQuery2()
            self.reverse.prepareQuery2()
            self.processed = True

    def updateQueryID(self):
        for f in self.forward.queryFoundList:
            f.query_id = self.id

        for r in self.reverse.queryFoundList:
            r.query_id = self.id

    def preprocessDestructor(self):
        self.forward.destructor2()
        self.reverse.destructor2()
        self.processed = False

    def write_empty(self, file):
        file.write(f'{self.id}\t4\t*\t0\t0\t*\t*\t0\t0\t{self.forward.seq}\t{"".join(["I" for i in range(self.forward.length)])}\tXM:i:0\n')

seed = namedtuple('seed', 'pos bitarray adjustedKmerIndexes')
class QuerySubclass:

    def __init__(self, seq, reverse = False):
        self.seq = seq
        self.length = len(seq)
        self.rev = reverse
        self.kmerListLen = -1
        self.queryFoundList = []
        self.seeds = [] 
        self.tmp = {}

    def destructor2(self):
        del self.seeds

    def prepareQuery2(self, kmerLen = KMER_LENGTH):
        '''
        Every Seed contain: bitarray, adjustedKmerIndexes
        '''
        self.kmerListLen = self.length - kmerLen + 1
        seedsCount = math.ceil(self.length/SEED_LENGTH)
        for i in range(seedsCount):
            i *= SEED_LENGTH
            end = i+SEED_LENGTH-1 if i+SEED_LENGTH-1 < self.length else self.length
            calculatedRange = range(i, end)
            bitArray = bitarray(self.length - kmerLen + 1)
            bitArray.setall(0)
            adjustedKmerList = []

            if end == self.length:
                lastMer = self.seq[-kmerLen:]
                bitArray[-1] = 1
                adjustedKmerList.append(lastMer)

            for j in range(self.length - kmerLen):
                if j not in calculatedRange: continue
                if j % kmerLen == 0:
                    seq = self.seq[j:j+kmerLen]
                    bitArray[j] = 1
                    adjustedKmerList.append(seq)

            adjustedKmerIndexes = self._levenshteinAutomata2(list(set(adjustedKmerList)))
            tempSeed = seed(i, ba2int(bitArray), adjustedKmerIndexes)
            self.seeds.append(tempSeed)

    def _levenshteinAutomata2(self, adjustedKmerList):
        tempList = []
        adjustedKmerIndexes = []
        for kmer in adjustedKmerList:
            a = list(find_all_matches(kmer, MISMATCH_PER_KMER, AUTOMATA))
            tempList += a

        tempList = list(set(tempList))
        tempList.sort()

        for i in tempList:
            adjustedKmerIndexes.append(KMER_LIST.index(i))
        return adjustedKmerIndexes

    def _referenceBitArray3(self, r, adjustedKmerIndexes, kmerLen = KMER_LENGTH):
        diff = kmerLen - 1
        finalBitArray = 1 << (r.wavelet.length-diff)
        for index in adjustedKmerIndexes:
            finalBitArray |= r.BFbitArray[index]
        return finalBitArray

    def BloomFilterSearch4(self, reference: Reference, mm: int = MAX_MM, maxLen = 1 << 16, processNum = 100):
        
        def allOnes(n):
            return ((n+1) & n == 0) and (n!=0)

        def set_bit(value, bit_index):
            return value | (1 << bit_index)

        def get_normalized_bit(value, bit_index):
            return (value >> bit_index) & 1

        referenceEnd = reference.wavelet.length
        cache = 1 << referenceEnd ## Store in Reverse lol
        potentialList = None
        reference_fragment = 1 if reference.adjusted else 0
        for s in self.seeds:
            bitArray = self._referenceBitArray3(reference, s.adjustedKmerIndexes)
            count = 0
            interestBA = s.bitarray
            while True:
                if count >= referenceEnd: break
                count += 1
                bitArray >>= 1
                
                if interestBA != interestBA & bitArray: 
                    continue

                #Adjust Seed Length
                stop = referenceEnd - count
                start = stop - SEED_LENGTH  - s.pos + 1
                stop = start + self.length
                #start = stop - self.length
                #if reference_fragment: 
                #    if start < MAX_LENGTH:
                #        continue

                if start < 0: continue
                if stop > referenceEnd: continue
                if stop - start < self.length: continue

                cacheValue = 1 << self.length
                for i in range(start, stop):
                    if get_normalized_bit(cache, i):
                        cacheValue = set_bit(cacheValue, i-start)

                if allOnes(cacheValue): continue
                for i in range(start, stop): cache = set_bit(cache, i)
                
                sequence = reference.wavelet.ReconstructSequence(start, stop)
                if myers(self.seq, sequence, mm): 
                    if myers(sequence, self.seq, mm):
                        continue

                if potentialList is None: 
                    if maxLen < (1<<16):
                        potentialList = np.array([start], dtype = np.uint16)
                    else:
                        potentialList = np.array([start], dtype = np.uint32)
                else:
                    potentialList = np.append(potentialList, start)

        
        if potentialList is None: 
            return

        self.queryFoundList += [QueryFound(reference.id, s+reference.adjusted, reference.wavelet.ReconstructSequence(s, s+self.length), self.seq, reverse = self.rev) for s in potentialList]

    @staticmethod
    def commonPrefixLen(str1, str2):
        count = 0
        n1 = len(str1)
        n2 = len(str2)
     
        # Compare str1 and str2
        i = 0
        j = 0
        while i <= n1 - 1 and j <= n2 - 1:
            if (str1[i] != str2[j]):
                break
                 
            count += 1
            i += 1
            j += 1
     
        return count

class QueryFound():

    def __init__(self, reference_ID, reference_pos, seqReference, query_seq, reverse = False, get_backtrace = True):
        self.query_id = None 
        self.reference_id = reference_ID
        
        self.reference_pos = reference_pos - 1
        self.reference = seqReference 
        self.query = query_seq 
        self.length = len(query_seq)
        self.reverse = reverse 
        self.mapping_quality = 255

        self.cigar = None
        self.editDistance = self.length     
        self.alignment = []
        self.ignore = False
        self.nm = float('inf')

    def write_result1(self, file, offset = 0):
        ## Write to SAM output 
        reverse_score = 16 if self.reverse else 0
        reference_pos = self.reference_pos + offset
        file.write(f"{self.query_id}\t{reverse_score}\t{self.reference_id}\t{self.reference_pos + 1}\t")
        file.write(f"{self.mapping_quality}\t{self.length}M\t*\t0\t0\t")
        #file.write(f"{decode_DNA(ba2int(self.reference),self.length)}\t{decode_DNA(ba2int(self.query), self.length)}\t{''.join(['I' for i in range(self.length)])}\tXA:i:{self.editDistance}\tMD:Z:{self.cigar}\tNM:i:{self.editDistance}\n")
        #file.write(f"{self.reference}\t{self.query}\t{''.join(['I' for i in range(self.length)])}\tXA:i:{self.editDistance}\tMD:Z:{self.cigar}\tNM:i:{self.editDistance}\n")
        file.write(f"{self.query}\t{''.join(['I' for i in range(self.length)])}\tXA:i:{self.editDistance}\tMD:Z:{self.cigar}\tNM:i:{self.nm}\n")


    def _edit_distance(self, matrix_a, matrix_b, skip = 0, get_backtrace=True, mm=MAX_MM):
        
        def __get_min(arr, validate = True):
            minimum = np.min(arr)
            minIndex = np.min(np.where(arr == minimum))
            if validate:
                assert(arr[minIndex] == arr.min())
            return minIndex

        string_x = self.reference #Y-axis
        string_y = self.query     #X-axis (Always the same)

        for i in range(1, len(string_x)+1):
            if i<skip: continue
            matrix_a[i, 0] = i
            cache = 0
            c = 0
            for j in range(1, len(string_y)+1):
                if abs(i-j) > mm: 
                    matrix_a[i, j] = mm+1
                    continue
                
                c += 1
                deletion = matrix_a[i - 1, j] + 1
                insertion = matrix_a[i, j - 1] + 1
                if (string_x[i - 1] == string_y[j - 1]):
                    substitution = matrix_a[i - 1, j - 1]
                else:
                    substitution = matrix_a[i - 1, j - 1] + 1
                
                score = min(([deletion, insertion, substitution]))  
                matrix_a[i, j] = score
                if get_backtrace:
                    temp = 0
                    if score == deletion:
                        temp |= 4
                    if score == substitution:
                        temp |= 2
                    if score == insertion:
                        temp |= 1    

                    matrix_b[i, j] = temp

                if score > mm:
                    cache += 1

            if cache == c:
                return 0

        i, j = len(string_x) + 1, len(string_y) + 1
        tmp = matrix_a[-1, -1]

        '''
        adjusted = 0
        skip = 0
        if tmp > mm:
            last_row = matrix_a[-1, :]
            last_row_index = __get_min(last_row)
            if last_row[last_row_index] < tmp:
                tmp = last_row[last_row_index]
                j = last_row_index #+ 1
                skip = True
                adjusted = True
            
            if not skip:
                last_col = matrix_a[:, -1]
                last_col_index = __get_min(last_col)
                if last_col[last_col_index] < tmp:
                    tmp = last_col[last_col_index]
                    i = last_col_index #+ 1
                    adjusted = True
        
        '''
        self.editDistance = tmp
        if tmp > mm:
            return 0

        if get_backtrace: 
            self._align_gRNA(self._naive_backtrace(matrix_b, i-1, j-1, matrix_a, skip))

        return 1

    def _naive_backtrace(self, B_matrix, i, j, A_matrix, skip):
        og = (i,j)
        backtrace_idxs = [(i, j)]

        while (i, j) != (0, 0):
            if B_matrix[i,j] & 2:
                i, j = i-1, j-1
            elif B_matrix[i,j] & 1:
                i, j = i, j-1
            elif B_matrix[i,j] & 4:
                i, j = i-1, j
            
            if len(backtrace_idxs) > 50: 
                print (i, j, B_matrix[i,j], len(backtrace_idxs), skip, og)
                print (B_matrix)
                print (A_matrix)
                raise ("Error")
 
            backtrace_idxs.append((i,j))

        return backtrace_idxs

    def _align_gRNA(self, bt):
        '''
        string_x - reference 
        string_y - query
        return cigar/ get position of Indels
        '''
        string_x = self.reference
        string_y = self.query
        op = []
        mutated_pos = [] #relative to reference
        backtrace = bt[::-1]
        
        mismatch_count = 0
        for k in range(len(backtrace)-1): 
            x0, y0 = backtrace[k]
            x1, y1 = backtrace[k+1]

            if x1 > x0 and y1 > y0:
                if string_x[x0] == string_y[y0]:
                    op.append('M')
                else:
                    mutated_pos.append(x0)
                    op.append(string_x[x0])
                    mismatch_count += 1
            elif x0 == x1:
                mutated_pos.append(y0)
                op.append('I')
            elif y0 == y1:
                mutated_pos.append(x0)
                op.append('D')

        operations = ['M', 'I', 'D']
        current_op = None
        count = 0 
        cigar = ''
        for char in op:
            if char in operations:
                if char != current_op:
                    if count >= 1: 
                        if current_op == 'M':
                            cigar += f'{count}{current_op}'
                        else:
                            cigar += (current_op * count)
                    current_op = char
                    count = 1
                else:
                    count += 1
            else:
                if current_op:
                    if count >= 1: 
                        if current_op == 'M':
                            cigar += f'{count}{current_op}'
                        else:
                            cigar += (current_op * count)
                    count = 0 
                    current_op = None
                cigar += f'{char}'
        else:
            if count >= 1 and current_op: 
                if current_op == 'M':
                    cigar += f'{count}{current_op}'
                else:
                    cigar += (current_op * count)

        if cigar[0] in BASIC_NT or cigar[0] in ['I', 'D']:
            cigar = '0M' + cigar

        self.nm = mismatch_count
        self.cigar = cigar.replace('M', '')
        if 'D' * KMER_LENGTH in self.cigar:
            self.ignore = True

def queueSearch4(queue, notsurewhy):
    while True:
        item = queue.get()
        if item == 'Done':
            return

def reset_matrix(matrix_a, matrix_0, length):
    matrix_a = np.zeros([length+1, length+1], dtype = np.uint8)
    matrix_a[0, :] = np.arange(0, length+1, 1)

    matrix_b = [[0 if x > 0 else 4 for x in range(length + 1)] for y in range(length + 1)]
    matrix_b[0] = [0 if x == 0 else 1 for x in range(length + 1)]
    matrix_b = np.array(matrix_b, dtype = np.uint8)

def refine_search1(query, mm = MAX_MM):
    length = query.length
    query.queryFoundList.sort(key=lambda x: x.reference)

    prev = None
    success = 1

    matrix_a = np.zeros([length+1, length+1], dtype = np.uint8)
    matrix_a[0, :] = np.arange(0, length+1, 1)
    matrix_a[:, 0] = np.arange(0, length+1, 1)

    matrix_b = [[0 if x > 0 else 4 for x in range(length + 1)] for y in range(length + 1)]
    matrix_b[0] = [0 if x == 0 else 1 for x in range(length + 1)]
    matrix_b = np.array(matrix_b, dtype = np.uint8)

    for q in query.queryFoundList:
        if not success: 
            reset_matrix(matrix_a, matrix_b, length)
            skip = 0
        else:
            skip = 0 if not prev else QuerySubclass.commonPrefixLen(prev, q.reference)

        success = q._edit_distance(matrix_a, matrix_b, skip = skip, mm=mm)
        prev = q.reference
            

def processSearch6(args, queryList, Reference, queue, processNum, lock, mm = MAX_MM):

    def check_alive(queue):
        queued = False
        while not queued:
            try:
                queued = True
            except:
                pass
            par_alive = mp.parent_process().is_alive()
            if not (par_alive and alive.value):
                alive.value = False
                queue.close()
                sys.exit(1)

    def get_maxLen(referenceList):
        maxLen = 0
        for r in referenceList:
            if r.wavelet.length > maxLen:
                maxLen = r.wavelet.length
        return maxLen

    maxLen = get_maxLen(Reference)

    #Preprocess (Lev Automata for the ends and multiple of 3s)
    for q in queryList:
        q.preprocess()

    for r in Reference:
        if r.wavelet.length > 1<<16: continue

        r.createBloomFilters()
        for q in queryList:
            if q.forward.length > r.wavelet.length: break
            q.forward.BloomFilterSearch4(r, mm, maxLen, processNum)
            check_alive(queue)
            q.reverse.BloomFilterSearch4(r, mm, maxLen, processNum)
            check_alive(queue)
        r.BFdestructor()

    logger.debug(f"{processNum} Completed Query Search!")
    del Reference
    for c, q in enumerate(queryList):
        length = q.length
        q.preprocessDestructor()

        refine_search1(q.forward, mm = mm)
        check_alive(queue)

        refine_search1(q.reverse, mm = mm)
        check_alive(queue)

        with lock:
            #print (processNum, c, 1, f"Length of Query List {len(q.forward.queryFoundList) + len(q.reverse.queryFoundList)}")
            with open(f'{args.prefix}.sam', 'a') as f:
                for qf in q.forward.queryFoundList:
                    if qf.ignore: continue
                    if qf.editDistance > mm: continue
                    qf.query_id = q.id
                    qf.write_result1(f)
                for qf in q.reverse.queryFoundList:
                    if qf.ignore: continue
                    if qf.editDistance > mm: continue
                    qf.query_id = q.id
                    qf.write_result1(f)
            #print (processNum, c, "Done")

        del q

    logger.debug(f"{processNum} Completed Refined Search!")
    queue.put("Complete")

def multiSearch(args, query, referenceList):

    def sort_reference(newLst, referenceList, count):
        # Sort by length
        c = 0
        direction = 1
        for r in referenceList:
            if r is None: continue
            p = c if direction else count - c - 1
            newLst[p].append(r)
            c += 1
            if c > count - 1:
                c = 0
                direction = not direction 
        return

    mm = args.mismatch  
    queue = mp.Queue()
    manager = mp.Manager()
    lock = mp.Lock()
    current_processes = []

    referenceList.sort(key = lambda x:x.wavelet.length)
    adjustedReferenceList = [[] for a in range(args.CPU)]
    sort_reference(adjustedReferenceList, referenceList, args.CPU)
    del referenceList

    for i in range(args.CPU):
        length_list = [i.wavelet.length for i in adjustedReferenceList[i]]
        logger.debug(f"Process {i}: Total   :{len(adjustedReferenceList[i])} References")
        logger.debug(f"Process {i}: Max Len :{max(length_list)} References")
        logger.debug(f"Process {i}: Med Len :{np.median(length_list)} References")
        logger.debug(f"Process {i}: Sum Len :{sum(length_list)} References")
        pp = mp.Process(target = processSearch6, args = (args, query, adjustedReferenceList[i], queue, i, lock, mm))
        current_processes.append(pp)

    for pp in current_processes:
        pp.start()

    sink = mp.Process(target=queueSearch4, args=(queue, True))
    sink.start()

    for pp in current_processes:
        pp.join()

    queue.put('Done')
    sink.join()

    return None

def write_output1(args, reference, minLen):
    prefix = args.prefix
    with open(f'{prefix}.sam', 'w') as f:
        f.write('@HD\tVN:1.0\tSO:unsorted\n')
        for r in reference:
            if r.wavelet.length < minLen: continue
            f.write(f'@SQ\tSN:{r.id}\tLN:{r.wavelet.length}\n')

        command_line = f"{PYTHON_FILE}/query.py --reference {','.join(args.reference)} --query {args.query} -p {args.prefix} -m {args.mismatch}"
        f.write(f'@PG\tID:Off-Target\tVN:0.0.1\tCL:"{command_line}"\n')

def open_query(args, checkUnique = True):
    file = args.query
    queryDict = {}
    adjusted = False
    minLen = float('inf')
    with open(file, 'r') as f: 
        seq = None
        for line in f:
            if line.startswith('>'): 
                if seq:
                    Query.totalCount += 1
                    if curr_id not in list(queryDict.keys()):
                        queryDict[curr_id] = Query(curr_id, seq)
                    else:
                        combi = ''.join(choice(ALPHABET_LIST) for i in range(4))
                        logger.info(f"{curr_id} exist already, will append {combi} to make unique!")
                        curr_id += combi
                        queryDict[curr_id] = Query(curr_id, seq)
                        adjusted = True
                curr_id = line.strip('\n').replace('>','')
                curr_id = curr_id.split(' ')[0]
                seq = ''
                continue
            line = line.strip('\n')
            seq += line.upper()
            if len(seq) > MAX_LENGTH: 
                logger.info(f"{seq} length > 50, Will Ignore!")
                continue
            if len(seq) < minLen:
                minLen = len(seq)
        
        if seq:
            Query.totalCount += 1
            if curr_id not in list(queryDict.keys()):
                queryDict[curr_id] = Query(curr_id, seq)
            else:
                combi = ''.join(choice(ALPHABET_LIST) for i in range(4))
                logger.info(f"{curr_id} exist already, will append {combi} to make unique!")
                curr_id += combi
                queryDict[curr_id] = Query(curr_id, seq)
                adjusted = True

    if adjusted:    
        prefix = file.replace('fasta', '') if 'fasta' in file else file.replace('fa','')
        with open(f'{prefix}.adjusted.fasta', 'w') as f:
            for k, v in queryDict.items():
                f.write(f'>{k}\n{v}\n')
        args.query = f'{prefix}.adjusted.fasta'
    
    return [v for k,v in queryDict.items()], minLen, args

__version__ = 3.0
PROG = f"Off-Target Pythonize Version: {__version__}"
DESCRIPTION = '''
Off-Target maps queries to hit with high mismatches (default: 1mm per 3bp) 
''' 

def get_args():
    parser = argparse.ArgumentParser(
        prog = PROG,
        description = DESCRIPTION
    )

    parser.add_argument('-r', '--reference', dest = 'reference', metavar="fasta", nargs = '+', required=True,
                        help="Reference Sequences")
    parser.add_argument('-q','--query', dest = "query", metavar="fasta", required=True,
                        help="Query Sequences")
    parser.add_argument('-m', '--mismatch', default=5, type=int,  choices=range(0,10), dest ='mismatch',
                        help='Number of mismatches allowed for Consensus Search (Default: %(default)s)')
    parser.add_argument('-p', '--prefix', default = 'Offtarget-Test', 
                        help='file name prefix to your file names (Default: %(default)s)')  
    parser.add_argument('--threads', default=2, type=int, dest = 'CPU', 
                        help="Number of Core to use (Default: %(default)s)")
    parser.add_argument('--mode', default='ALL', choices=['ALL', 'BEST'],
                        help="BEST searchers for the least mismatches (slighlty faster); ALL searches for the all mismatches")
    parser.add_argument('-d', '--debug',
                        help='Print lots of debugging statements',
                        action="store_const",dest="loglevel",const=logging.DEBUG,
                        default=logging.INFO)

    args = parser.parse_args()
    return args


def main():
    args = get_args()
    logging.basicConfig(level=args.loglevel, format=FORMAT)

    logger.info("*"*25)
    logger.info(f"{PROG} starting!")
    
    logger.info(f"Reference File     : {','.join(args.reference)}")
    logger.info(f"Query File         : {args.query}")
    logger.info(f"Writing to Output  : {args.prefix}.sam")
    logger.info(f"Number of threads  : {args.CPU}")
    logger.info(f"Number of mismatch : {args.mismatch}")
    logger.info("*"*25)

    ### START
    logger.info(f"Processing Reference File!")
    referenceList = []
    for f in args.reference:
        openReference(referenceList, f)

    logger.info(f'Number of Reference Processed : {Reference.totalCount}')
    logger.info(f'Median Length of Reference    : {np.median([i.wavelet.length for i in referenceList])}')
    
    logger.info("*"*25)
    logger.info(f"Processing Query File!")
    queryList, minLen, args = open_query(args)
    logger.info(f'Number of Query Processed: {Query.totalCount}')

    write_output1(args, referenceList, minLen)
    
    start_time = time.time()
    logger.info("*"*25)
    logger.info("Starting Off-Target Search!")
    queryList = multiSearch(args, queryList, referenceList)
    
    end_time = time.time()
    total_time = round(end_time - start_time, 2)
    logger.info(f"{PROG} has Ended in {total_time} seconds!")
    logger.info("*"*25)
    ### END

if __name__ == '__main__':
    main()