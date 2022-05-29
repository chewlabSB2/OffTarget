# Off-Target

Off-Target is a short-string (<50bp) search query tool which maps queries to hits with high mismatches or Insertion Deletions (InDels) with a default setting of 1 mismatch per 3 nucleotides. 

Whilst the algorithm is constantly being updated and a C version is under production, a brief explanation behind Off-Target algorithm will be updated soon! (Stay Tuned!)

Installation 
------------

Installation:
```bash 
git clone https://github.com/chewlabSB2/OffTarget.git
cd OffTarget 
python setup.py install
```

Simple Usage
------------

```bash
OffTarget -r References.fasta -q Query.fasta -p Test -m 8 --threads 2 -d
```

Output
------

Off-Target output is according to SAM file conventions. A sample of the results can be found [here](https://github.com/chewlabSB2/OffTarget/blob/main/sample/Test.sam).