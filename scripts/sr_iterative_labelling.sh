raw=$1

reference=$2

basecall=$3

bonito=$4

chunks=5000
num_reads=20000
lr=5e-4
iters=3

echo "raw data : " $1
echo "Reference Genome : " $2

source ~/.bashrc 
conda activate bonito_env
ml samtools

mkdir -p default/
singularity exec  --nv $bonito/bonito.sif python3 $basecall rna004_130bps_hac@v5.0.0 $raw --reference $reference --rna --recursive  --max-reads 200000 --chunksize 8000 > default/basecalls.bam
samtools view -b default/basecalls.bam > default/test.bam
samtools sort default/test.bam -o default/basecalls.sorted.bam
samtools index default/basecalls.sorted.bam

mkdir -p round0/test

singularity exec  --nv $bonito/bonito.sif python3 $basecall rna004_130bps_hac@v5.0.0 $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --chemo --min-accuracy-save-ctc 0.8 > round0/front.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./round0/ ./round0/front --pretrain rna004_130bps_hac@v5.0.0 --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall rna004_130bps_hac@v5.0.0 $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > round0/general.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./round0/ ./round0/general --pretrain ./round0/front --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall rna004_130bps_hac@v5.0.0 $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --sr --min-accuracy-save-ctc 0.8 > round0/tail.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./round0/ ./round0/fine-tuned-model --pretrain ./round0/general --chunks $chunks

singularity exec  --nv $bonito/bonito.sif python3 $basecall ./round0/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > round0/test/basecalls.bam
samtools view -b round0/test/basecalls.bam > round0/test/test.bam
samtools sort round0/test/test.bam -o round0/test/basecalls.sorted.bam
samtools index round0/test/basecalls.sorted.bam
# ############ Round 0 of Iterative Labeling ##############

round=round1
model=round0
mkdir -p $round/test

singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --chemo --min-accuracy-save-ctc 0.8 > $round/front.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/front --pretrain $model/fine-tuned-model --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > $round/general.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/general --pretrain ./$round/front --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --sr --min-accuracy-save-ctc 0.8 > $round/tail.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/fine-tuned-model --pretrain ./$round/general --chunks $chunks

singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$round/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > $round/test/basecalls.bam
samtools view -b $round/test/basecalls.bam > $round/test/test.bam
samtools sort $round/test/test.bam -o $round/test/basecalls.sorted.bam
samtools index $round/test/basecalls.sorted.bam
# ############ Round1 #############

round=round2
model=round1
mkdir -p $round/test

singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --chemo --min-accuracy-save-ctc 0.8 > $round/front.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/front --pretrain $model/fine-tuned-model  --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > $round/general.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/general --pretrain ./$round/front  --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --sr --min-accuracy-save-ctc 0.8 > $round/tail.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/fine-tuned-model --pretrain ./$round/general  --chunks $chunks


singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$round/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > $round/test/basecalls.bam
samtools view -b $round/test/basecalls.bam > $round/test/test.bam
samtools sort $round/test/test.bam -o $round/test/basecalls.sorted.bam
samtools index $round/test/basecalls.sorted.bam

# ############ Round2 #############

round=round3
model=round2
mkdir -p $round/test

singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --chemo --min-accuracy-save-ctc 0.8 > $round/front.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/front --pretrain $model/fine-tuned-model --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --min-accuracy-save-ctc 0.8 > $round/general.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/general --pretrain ./$round/front  --chunks $chunks
singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$model/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads $num_reads --chunksize 5000 --save-ctc --sr --min-accuracy-save-ctc 0.8 > $round/tail.bam
singularity exec  --nv $bonito/bonito.sif bonito train --epochs $iters --lr $lr --directory ./$round/ ./$round/fine-tuned-model --pretrain ./$round/general --chunks $chunks


singularity exec  --nv $bonito/bonito.sif python3 $basecall ./$round/fine-tuned-model $raw --reference $reference --rna --recursive  --max-reads 50000 --chunksize 8000 > $round/test/basecalls.bam
samtools view -b $round/test/basecalls.bam > $round/test/test.bam
samtools sort $round/test/test.bam -o $round/test/basecalls.sorted.bam
samtools index $round/test/basecalls.sorted.bam
############ Round3 #############

