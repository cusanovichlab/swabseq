# swabseq
Code for analyzing Swab-Seq data from the Cusanovich lab.

The basic steps of the data processing workflow are as follows:

1. Convert a 'plate map' .xls to a 'sample sheet' .csv using the `platemap2samp.py` script from https://github.com/octantbio/SwabSeq (we provide this script in our repo as well for convenience). A 'plate map' is a specially formatted Excel spread sheet that establishes the barcodes used and the experimental conditions present in a Swab-Seq experiment. Please see the Octant repo for further details. We provide an example here as well.
2. Convert the 'sample sheet' to a metadata file using the `parseSS.py` script available from the Pachter lab (https://caltech.box.com/shared/static/m6t4ok1bqwuhy3f6tut9sqantufvtiro.gz) as described in their notebook on processing Swab-Seq data (https://github.com/pachterlab/BLCSBGLKP_2020/blob/master/notebooks/swabseq.ipynb).
3. Convert the metadata to a 'whitelist' using the following command (from the Pachter lab workflow): `cat [metadata.txt] | awk '{print $1}' | tail -n +2 > [whitelist.txt]`.
4. Generate fastq files. In our case, we used the `fastq_barcode_correcter_reformatter_w_exclusion_list.py` script to demultiplex reads before further processing. An example of usage is: `python fastq_barcode_correcter_reformatter_w_exclusion_list.py [R1 fastq] [I1 fastq] [I2 fastq] [Sample sheet] [Minus list] [Out prefix]`. The `Sample sheet` is the same one that will be used to create the metadata file in subsequent steps. The `Minus list` is a tab-separated list of barcode combinations ('i7\ti5') that were included on the sequencing run, but should be excluded from demultiplexing (this was employed to guard against any erroneous assignments of barcodes due to tolerance of mismatches). The `Out prefix` should specify the directory and the name of the output files (suffixes are appended by the script).
5. Index the custom reference transcriptome available from the Pachter lab tarball downloaded in Step 2: `kallisto index -i colab/index.idx -k 11 colab/trunc_transcriptome_11.fa`.
6. Map reads to the Swab-Seq-specific reference. Again, following the Pachter lab workflow (https://github.com/pachterlab/BLCSBGLKP_2020/blob/master/notebooks/swabseq.ipynb), except that we used the SwabSeq10 configuration of kallisto (available on the `covid` branch, which we have forked for convenience - https://github.com/cusanovichlab/kallisto): `kallisto bus -x SwabSeq10 -o [Outdir] -t 2 -i colab/index.idx [I1 fastq] [I2 fastq] [R1 fastq]`.
7. Sort the kallisto output with bustools: `colab/bustools sort -o sort.bus output.bus`.
8. Correct barcodes with bustools: `colab/bustools correct -d dump.txt -w whitelist.txt -o sort.correct.bus sort.bus`.
9. Sort the barcode-corrected bus file: `colab/bustools sort -o sort.correct.sort.bus sort.correct.bus`.
10. Generate a text file of read counts mapping to target genes for each sample barcode: `colab/bustools text -p sort.correct.sort.bus > data.txt`.
11. Read data.txt file and metadata.txt files into R and generate appropriate plots with R scripts provided in https://github.com/cusanovichlab/swabseq/tree/main/analysis,
