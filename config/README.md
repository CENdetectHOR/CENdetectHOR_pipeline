# config options for CENdetectHOR pipeline

All the values used here are the default ones.

Path of the folder where the pipeline will be deployed. Here it is expected the `fasta` folder with the input files. 
```
workdir: test/
```

File with the consensus sequence of the centromeric satellite for the analyzed species. Leave blank if executing without it. 
```
consFile: human_cons_alphasat.txt
```

Size of the windows for the periodicity step.
```
windowsize: 1000
```

Size of the kmer for the periodicity step.
```
kmerSize: 8
```

Maximum allowed gap (bp) between monomers in HORs.
```
maxHORgap: 10
```

Maximum HOR length in terms of monomer families.
```
maxHORlength: 50
```

Minumum number of contiguous repetitions to detect a HOR
```
minHORrep: 3
```

Method to construct the HOR tree, if discrete levels are defined by height from leaves. 
```
discreteTree: True
```

Number of threads to use for building the distance matrix, this is correlated with the size of the chunks.
```
HORdet_threads: 16
```
