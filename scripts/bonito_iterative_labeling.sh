#basecall=/athena/chenlab/scratch/ziw4007/bioRNA/analysis/scripts/sr_basecall.py
#bonito=/athena/chenlab/scratch/ziw4007/RNA004/software/bonito

raw=$1

reference=$2

bonito=$3

basecall=$4

chunks=5000
num_reads=20000
lr=5e-4
iters=9

echo "raw data : " $1
echo "Reference Genome : " $2
echo "bonito Software : " $2
echo "sr_basecall.py : " $2
source ~/.bashrc 
conda activate bonito_env
ml samtools


mkdir -p round0/test

singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller rna004_130bps_hac@v5.0.0 $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > round0/basecalls.bam
singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./round0/ ./round0/fine-tuned-model --pretrain rna004_130bps_hac@v5.0.0 --chunks $chunks

singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./round0/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > round0/test/basecalls.bam
samtools view -b round0/test/basecalls.bam > round0/test/test.bam
samtools sort round0/test/test.bam -o round0/test/basecalls.sorted.bam
samtools index round0/test/basecalls.sorted.bam
# ############ Round 0 of Iterative Labeling ##############

round=round1
model=round0
mkdir -p $round/test

singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > $round/basecalls.bam
singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/fine-tuned-model --pretrain $model/fine-tuned-model  --chunks $chunks

singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$round/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > $round/test/basecalls.bam
samtools view -b $round/test/basecalls.bam > $round/test/test.bam
samtools sort $round/test/test.bam -o $round/test/basecalls.sorted.bam
samtools index $round/test/basecalls.sorted.bam
# ############ Round1 #############

round=round2
model=round1
mkdir -p $round/test

singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > $round/basecalls.bam
singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/fine-tuned-model --pretrain $model/fine-tuned-model  --chunks $chunks


singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$round/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > $round/test/basecalls.bam
samtools view -b $round/test/basecalls.bam > $round/test/test.bam
samtools sort $round/test/test.bam -o $round/test/basecalls.sorted.bam
samtools index $round/test/basecalls.sorted.bam

# ############ Round2 #############

round=round3
model=round2
mkdir -p $round/test

singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > $round/basecalls.bam
singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/fine-tuned-model --pretrain $model/fine-tuned-model  --chunks $chunks


singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$round/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > $round/test/basecalls.bam
samtools view -b $round/test/basecalls.bam > $round/test/test.bam
samtools sort $round/test/test.bam -o $round/test/basecalls.sorted.bam
samtools index $round/test/basecalls.sorted.bam
############ Round3 #############

model=round3
mkdir -p all
singularity exec --bind /athena:/athena --nv $bonito/bonito.sif bonito basecaller ./$model/fine-tuned-model $validation --reference $reference --rna --recursive  --max-reads 200000 --chunksize 8000 > all/basecalls.bam
samtools view all/basecalls.bam -b > all/test.bam; samtools sort all/test.bam -o all/basecalls.sorted.bam; samtools index all/basecalls.sorted.bam
############ 3 round of Iterative Labeling ##############