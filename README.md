# CENdetectHOR pipeline

This Snakemake pipeline has been created to deeply analyze centromeric region, by taking as input single-chromosome (or centromere, if you already know its location) fasta sequences and optionally the consensus sequence of the centromeric satellite for the analysed organism. 
Without needin any other information, CENdetectHOR is able to:
- precisely locate the higly repetitive region of the input sequence;
- infer the size of the repetitive unit;
- extract all the monomeric sequences of the higly repetitive region;
- cluster monomers to phylogenetically identify Higher-Order Repeats (HORs) and monomer faimies.

![Workflow](/img/Workflow.pdf)

CENdetectHOR has been implemented with [cen_detect_hor library](https://github.com/CENdetectHOR/CENdetectHOR_lib) for the HOR identification andd with the [PhyloTreeGUI](https://github.com/CENdetectHOR/PhyloTreeGUI) to visualize the detected HOR, select the most informative ones, and export output files 

## Running CENdetectHOR pipeline

### Setting up the environment

To run CENdetectHOR, we raccomend creating a new `conda env` with a current version of `snakemake` to run this workflow. To get `snakemake` please referr the [install instructions on their website](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

Alternatively, you can create a new env with:
```
conda create -n cendetecthor
```

Once created, you can activate the `conda env` env and download the repository. 
```
conda activate cendetecthor
git clone https://CENdetectHOR/CENdetectHOR_pipeline.git
conda install --file requirements.txt
```

### Running

To run CENdetectHOR, place all the fasta files to be analyzed in a directory and specify the path in the [config.yaml](config/config.yaml) file accordingly.

Once this is done and you have activated your `conda env` you can run the pipeline like so:
```
snakemake --cores 24
```
Or do a dry run of the pipeline:
```
snakemake --cores 24 -n
```
All parameters are described above and you can modify any of them by directly editing the [config.yaml](config/config.yaml) file or via the command line. For example, to change the sample identifier and fasta options do:
```
snakemake --cores 24 --config fasta=/some/fasta/path/
```

## Running by to submitting jobs to a HTCondor cluster. 

The pipeline has been developed and tested by submitting jobs to a HTCondor cluster, by using [HTCondor Snakemake profile](https://github.com/Snakemake-Profiles/htcondor). 
The required libraries are already present in the `requirements.txt`, but you will need to deploy the profile as specified in the [repo](https://github.com/Snakemake-Profiles/htcondor). 
Then, you can run the pipleine by using:
```
snakemake --profile htcondor
```


