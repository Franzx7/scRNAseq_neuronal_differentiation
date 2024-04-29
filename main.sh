
nextflow run -resume ~/FAKE/SRCs/SCALPEL/scalpel.nf \
--sequencing dropseq \
--reads ~/FAKE/DATAS/NEURODIFF/FASTQs/ \
--samples ~/FAKE/PROJECTS/scRNAseq_neuronal_differentiation/datas/ \
--gtf ~/FAKE/DATAS/REFERENCES/gencode.v41/gencode.v41.annotation.gtf \
--transcriptome ~/FAKE/DATAS/REFERENCES/gencode.v41/gencode.v41.transcripts.fa \
--ipdb /home/fake/FAKE/DATAS/REFERENCES/gencode.v41/refdata-gex-GRCh38-2020-A.polyA.RESTRICTED_chr_corrected.track \
--cpus 40 \
-with-conda /home/fake/.conda/envs/scalpel_env
