# Cell Ranger

FASTQ files

- Read 1 - cell barcode & UMI (unique molecule identifier)
- Read 2 - cDNA

Reference transcriptome

# Workflow

1. Cell preparation

   - In droplets containing single cells, tag mRNA with 10x oligonucleotides
   - Reverse transcription of tagged transcripts yield cDNA
   - cDNA amplified with PCR

2. Library Prep

   - Sequencing adapters ligated to prepare for Illumina sequencing
   - **Adapters** allow for **paired-end sequencing**
   - R1 reads cell barcode + UMI
   - R2 reads cDNA

   FASTA structure:

   ```txt
   @<read_id>
   <sequence>
   +
   <quality_scores>
   ```

   R1:

   ```txt
   @A00471:69:H3Y23DRXX:1:1101:1043:1000 1:N:0:GGCTAC	# readID, should match R2 for same read
   AAACCTGAGACGTTGCCTATGAACAC	# cell barcode 1st 16bp, UMI next 10bp
   +
   FFFFFFFFFFFFFFFFFFFFFFFFFF	# quality score
   ```

   R2:

   ```txt
   @A00471:69:H3Y23DRXX:1:1101:1043:1000 2:N:0:GGCTAC	# readID, should match R1 for same read
   GAGTCCAGGAGGAGGCAGAGGACAGAGGAGGAGGAAGAGGAGGAAGGAGGAAG	# cDNA, 5' -> 3'
   +
   FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	# quality score
   ```

   

3. Extract cell barcodes & UMI from R1
   - Identify which cell & which RNA molecule a read came from

4. Align cDNA to reference transcriptome (**STAR** RNA-seq alignment (Spliced Transcript Alignment to Reference))
   - From alignment, maps cDNA to gene

5. Filter + UMI Counting

   - PCR duplicates cDNA with oligonucleotide (cell barcode + UMI)
   - Double counting PCR duplicates gives erroneous measure of gene expression
   - Only count unique UMI that are mapped to a gene after alignment
   - Result is a **UMI Count Matrix**

   | Gene  | Cell A | Cell B | Cell C |
   | ----- | ------ | ------ | ------ |
   | ACTB  | 12     | 3      | 0      |
   | GAPDH | 7      | 10     | 4      |

   

**Reverse Transcription**

- Reads mRNA 3' to 5', therefore synthesizes cDNA 5' to 3'

**Oligonucleotide** `5' cell barcode, UMI, & poly-dT tail 3'`.

- poly-dT tail hybridizes to 3' of mRNA poly-A tail.
- oligonucleotide acts as primer for RT
- oligonucleotide gets copied in cDNA



## Cell Preparation

10x uses **microfluid droplet system**. Each droplet has:

- single cell
- barcoded gel bead with millions of oligonucleotides; **same cell barcode**, **unique UMI**
  - there are multiple mRNA transcripts, so UMI helps unique tag them
  - reads with identical cell barcode + UMI implies PCR duplicates

**Steps:**

1. Cell lysed inside droplet; mRNA is released
2. oligonucleotide hybridizes to mRNA via 3' poly-A tail & 3' poly-dT
3. Cell lysed & PCR



# Code

Cell Ranger v8.0 requires `--create-bam` parameter for `cellranger count`. This replaces `--no-bam` option.

`cellranger count` aligns sequencing reads in FASTQ to reference transcriptome.

```bash
#!/bin/bash
#$ -cwd
#$ -o /u/scratch/c/cbui/cellranger_tutorial/job-logs/joblog.$JOB_ID
#$ -j y
### CHANGE RESOURCES AS NEEDED:
#$ -l h_rt=3:00:00,h_data=8G,exclusive
#$ -pe shared 4
### CHANGE NAME OF JOB AS NEEDED:
## $ -N NAMEOFJOB
### EMAIL ADDRESS TO NOTIFY:
#$ -M $USER@mail
### NOTIFY WHEN
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load environment
. /u/local/Modules/default/init/modules.sh

module load cellranger/6.1.1
echo " "

################################################################

ID=/u/scratch/c/cbui/cellranger_tutorial/run_count_1kpbmc
TRANSCRIPTOME=/u/scratch/c/cbui/cellranger_tutorial/refdata-gex-GRCh38-2020-A
FASTQS=/u/scratch/c/cbui/cellranger_tutorial/pbmc_1k_v3_fastqs
SAMPLE=pbmc_1k_v3
JOBMODE=sge

# cellranger count
cellranger count \
    --disable-ui \
    --id=$ID \
    --transcriptome=$TRANSCRIPTOME \
    --fastqs=$FASTQS \
    --sample=$SAMPLE \
    --jobmode=$JOBMODE \
    --jobinterval=10000

################################################################

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
```



# Misc.

**Library:** Collection of DNA / cDNA prepared (hybridized with oligonucleotides containing Adapter, Cell Barcode, UMI) for sequencing.

| Library Type          | Purpose                         | Input Material  |
| --------------------- | ------------------------------- | --------------- |
| Gene Expression (GEX) | Measure mRNA transcript levels  | Single-cell RNA |
| V(D)J                 | Study immune receptor sequences | BCR/TCR mRNA    |
| ATAC                  | Measure chromatin accessibility | Nuclei          |
| Multiome              | GEX + ATAC from same cell       | Nuclei          |

**Multiplex:** Putting libraries from multiple samples into the same flow cell for sequencing. Use cell barcode to demultiplex (group libraries into correct sample) for consistent analysis.



# Space Ranger

Similar to cell ranger, but uses spatial barcode to label which cell (& therefore position on image) that a specific column's gene expression is associated with.

| **Gene**   | AAAGATGGTCCGAAAG | AAAGATGGTGAGTGAC | AAAGCAATCACCTTAC |
| ---------- | ---------------- | ---------------- | ---------------- |
| **MALAT1** | 34               | 21               | 11               |
| **GAPDH**  | 102              | 87               | 96               |
| **ACTB**   | 56               | 49               | 60               |
| **COL1A1** | 8                | 0                | 2                |
| **FN1**    | 0                | 1                | 0                |
| **MYH11**  | 0                | 0                | 6                |



There is a separate `tissue_positions_list.csv` that provides details for each spatial barcode.

| **Barcode**      | **In_tissue** | **Array_row** | **Array_col** | **Pixel_row** | **Pixel_col** |
| ---------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| AAAGATGGTCCGAAAG | 1             | 35            | 43            | 3357          | 2894          |
| AAAGATGGTGAGTGAC | 1             | 36            | 43            | 3412          | 2890          |
| AAAGCAATCACCTTAC | 0             | 45            | 52            | 4022          | 3480          |

Details:

| Column            | Description                                         |
| ----------------- | --------------------------------------------------- |
| **Barcode**       | Spot barcode (matches columns in expression matrix) |
| **In_tissue**     | `1` = spot overlaps tissue; `0` = background        |
| **Array_row/col** | Grid layout of the spot on the slide                |
| **Pixel_row/col** | Spot's position on the full-resolution tissue image |