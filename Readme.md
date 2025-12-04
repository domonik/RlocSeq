# Minimal documentation

This repository exists for reproducibility of the corresponding publication.

*Spatial organization of a cyanobacterial transcriptome *



## System requirements

### Hardware requirements

This pipeline requires a standard computer with at ~ 16 GB of RAM.

### Software requirements

#### OS
This pipeline requires a Linux system since some dependencies are not available for Windows or macOS it was tested using the following setup:

- Linux: Ubuntu 22.04 and Ubuntu 20.04

####  Dependencies

Dependencies are stored in the [environment.yml](./environment.yml) and in the yaml files in the [envs](./envs) directory. 
The specific versions used for the publication are stored in the `*.pin.txt` files


## How to run

To run this pipeline you must install the required packages from the environment.yml file using conda:

```shell
conda env create -f environment.yml
conda activate RlocSeq
```

afterwards you can run the pipeline with the following command. Please adjust the cores parameter according to your system.


```shell
snakemake -s pipeline.smk --cores 4 --use-conda --configfile configFile.yml --rerun-incomplete
```

This will run the whole pipeline and place generated files in a Pipeline directory. 
