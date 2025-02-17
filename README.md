# The Precise Basecalling of Short-Read Nanopore Sequencing

## Introduction

Nucleotide modifications deviate nanopore sequencing readouts, therefore generating artifacts during the basecalling of sequence backbones. Here, we present an iterative approach to polish modification-disturbed basecalling results. We show such an approach is able to promote the basecalling accuracy of both artificially-synthesized and real-world molecules. With demonstrated efficacy and reliability, we exploit the approach to precisely basecall therapeutic RNAs consisting of artificial or natural modifications, as the basis for quantifying the purity and integrity of vaccine mRNAs which are transcribed in vitro, and for determining modification hotspots of novel therapeutic RNA interference (RNAi) molecules which are bioengineered (BioRNA) in vivo.

## Major Contribution: 3-step sampling

Our study shows that compromised basecalling can be improved through an iterative workflow. To enhance polishing at the 3’ and 5’ ends, which is crucial for short reads, we developed a 3-step sampling strategy. Reads are sampled from the 5’ end, the full molecule, and the 3’ end, ensuring even coverage and better basecalling at both termini.

# Installation

## Pre-request

+ [singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
+ [git](https://git-scm.com/) 

## Step 1: git clone

```bash
git clone https://github.com/wangziyuan66/iterative-labeling-toolkit-bonito
```

## Step 2: build sif file

```bash
singularity build bonito.sif bonito.recipe
```

Create a standalone envrionment to run the iterative-labeling-bonito.

# Usage

```bash
# raw=$1

# reference=$2

# bonito=$3

# basecall=$4

bash scripts/sr_iterative_labelling.sh raw reference ./ ./scripts/sr_basecall.py
```

+ **raw** : The path to folder containing raw pod5 files.

+ **reference** : Reference genome path.

+ **bonito** : The path to the bonito singularity image.

+ **basecall** : The path to the "basecall.py" file.

# Miscellaneous

If the sequencing kit is RNA002, we recommend you to use [iterative-labeling-toolkit-taiyaki(https://github.com/wangziyuan66/iterative-labeling-toolkit-taiyaki). Currently, for the first round of basecalling we are using RNA004 hac 5.0.0 model.

# Data availability

Sample **raw pod5 files** are provided for **BioRNA-Leu** **BioRNA-Ser** **ChemoRNA-Leu** and **ChemoRAN-Ser** which you can downloaded in [here](example/bioRNA). In sra, only bam can be uploaded. If some need more rawdata, contact us. 

# Contact

Ziyuan Wang princezwang@arizona.edu
