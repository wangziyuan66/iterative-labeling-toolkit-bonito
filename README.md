# The Precise Basecalling of Short-Read Nanopore Sequencing

## Introduction

Nucleotide modifications deviate nanopore sequencing readouts, therefore generating artifacts during the basecalling of sequence backbones. Here, we present an iterative approach to polish modification-disturbed basecalling results. We show such an approach is able to promote the basecalling accuracy of both artificially-synthesized and real-world molecules. With demonstrated efficacy and reliability, we exploit the approach to precisely basecall therapeutic RNAs consisting of artificial or natural modifications, as the basis for quantifying the purity and integrity of vaccine mRNAs which are transcribed in vitro, and for determining modification hotspots of novel therapeutic RNA interference (RNAi) molecules which are bioengineered (BioRNA) in vivo.

## Major Contribution: 3-step sampling

Our recent study demonstrated that compromised basecalling can be polished via an iterative workflow. Our workflow builds upon the community consensus that Bonito can achieve acceptable accuracy for generic basecalling tasks. We then polish the yielded sketch sequences by aligning them to the ground-truth reference. We next take polished sequences, and their nanopore sequencing readouts, as training data to update Bonito. To make our workflow more effective in polishing 3’ and 5’-ends, we further develop a 3-step sampling strategy for preparing training data. Specifically, we sample reads that are mapped to the 5’-end, the entire molecule and the 3’-end subsequently. Our strategy guarantees the candidate sequence to be evenly covered by training reads, and therefore promises to improve the basecalling of 3’ and 5’-ends

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

# Data availability

Sample raw pod5 files are provided for **BioRNA-Leu** **BioRNA-Ser** **ChemoRNA-Leu** and **ChemoRAN-Ser** which you can downloaded in ![here](example/bioRNA). In sra, only bam can be uploaded. If some need more rawdata, contact us. 

# Contact

Ziyuan Wang princezwang@arizona.edu
