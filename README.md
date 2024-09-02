# Introduction

Nucleotide modifications deviate nanopore sequencing readouts, therefore generating artifacts during the basecalling of sequence backbones. Here, we present an iterative approach to polish modification-disturbed basecalling results. We show such an approach is able to promote the basecalling accuracy of both artificially-synthesized and real-world molecules. With demonstrated efficacy and reliability, we exploit the approach to precisely basecall therapeutic RNAs consisting of artificial or natural modifications, as the basis for quantifying the purity and integrity of vaccine mRNAs which are transcribed in vitro, and for determining modification hotspots of novel therapeutic RNA interference (RNAi) molecules which are bioengineered (BioRNA) in vivo.

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

# Contact

Ziyuan Wang princezwang@arizona.edu
