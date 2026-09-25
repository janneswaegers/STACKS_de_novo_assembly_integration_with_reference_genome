# Stacks: de novo assembly integrated with a reference genome

A step-by-step guide for running the Stacks pipeline de novo and then placing the de novo
loci on a reference genome. I wrote it for PhD students in our group who work with RAD-seq
data from species with a draft genome assembly, based on my work on *Ischnura* damselflies.

## Why combine de novo and reference?

Stacks can build RAD loci in two ways:

- **De novo** (`denovo_map.pl`) builds loci directly from the reads. It needs no genome
  and uses all reads, including those from regions that are missing or fragmented in a
  draft assembly.
- **Reference-based** (`ref_map.pl`) builds loci from reads aligned to a genome, which gives
  every locus a genomic position.

Integrating the two keeps the de novo loci but adds their position on the genome. Genomic
positions are needed for analyses along the genome, such as sliding-window statistics
(`populations -k`), splitting SNPs by chromosome (for example sex chromosomes versus
autosomes) or finding genes near outlier loci.

## Workflow

```
decloned reads
     │
     ├─ 1. Optimise M and n on a subset of samples      (denovo_map.pl, -r 0.8)
     │
     ├─ 2. De novo assembly of all samples              (denovo_map.pl)
     │
     ├─ 3. Align the de novo catalog to the genome      (bwa mem)
     │       └─ check alignment quality                 (samstat)
     │
     ├─ 4. Add genomic positions to the catalog         (stacks-integrate-alignments)
     │
     └─ 5. Filter SNPs and export                       (populations)
```

## Before you start

You need:

- **Cleaned, decloned reads.** Demultiplexed with `process_radtags` and PCR duplicates
  removed with `clone_filter`, named `<sample>.1.fq.gz` and `<sample>.2.fq.gz`.
- **A popmap.** A tab-separated file with one line per sample: `sample_name<TAB>population`.
- **A reference genome** in fasta format. With a fragmented draft assembly it can help to
  keep only the longer scaffolds (we used scaffolds > 50 kb).
- **Stacks 2**, BWA and SAMtools. `stacks-integrate-alignments` is part of Stacks 2.

The examples below use SLURM job scripts as on UPPMAX. Replace `<account>` with your own
project, and adjust module names and resources for your cluster.

To follow a running job:

```bash
scontrol show jobid -dd <job_id>
```

---

## 1. Optimise the de novo parameters M and n

The two main de novo parameters are:

- **M**: the number of mismatches allowed between stacks *within* one individual when
  merging them into a locus.
- **n**: the number of mismatches allowed between loci of *different* individuals when
  building the catalog.

Values that are too low split one real locus into several; values that are too high merge
different loci (paralogs) into one. Following Paris et al. (2017, *Methods in Ecology and
Evolution*), choose the values that give the highest number of **r80 loci**: loci present
in at least 80% of the samples of a population.

Use a subset of around 20 samples that covers the variation in your data set, for example
one individual per population (`popmap_1eachpop`). Then run M = 1 to 6, with n = M:

```bash
for M in 1 2 3 4 5 6; do
    denovo_map.pl -M $M -n $M -T 16 \
        -o ./opt/M${M} \
        --popmap ./popmap_1eachpop \
        --samples ../decloned_reads \
        --paired \
        -X "populations: -r 0.8"
done
```

Count the number of r80 loci for each run (the first column of `populations.hapstats.tsv`
is the locus ID):

```bash
for dir in ./opt/M*; do
    echo "$dir: $(awk '{print $1}' $dir/populations.hapstats.tsv | sort | uniq | wc -l) loci"
done
```

For the M values with the most loci, also try n = M − 1 and n = M + 1. In our data M = 2 and
n = 3 gave the best result.

## 2. De novo assembly of all samples

Run `denovo_map.pl` on all samples with the chosen parameters:

```bash
#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J stacks_denovo

module load bioinfo-tools Stacks

denovo_map.pl -M 2 -n 3 -T 16 \
    -o ./stacks_outputM2n3 \
    --popmap ./all_popmap \
    --samples ../decloned_reads \
    --paired
```

The catalog of all loci is written to `stacks_outputM2n3/catalog.fa.gz`. This is the file
we align to the genome in the next step.

## 3. Align the catalog to the reference genome

In our tests, BWA gave better results than Bowtie2 for aligning the catalog. Index the
genome once:

```bash
#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -J bwa_index

module load bioinfo-tools bwa

bwa index ./index/assembly_50k.fasta
```

Then align the catalog loci:

```bash
#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J bwa_catalog

module load bioinfo-tools bwa

bwa mem -t 16 ./index/assembly_50k.fasta ./stacks_outputM2n3/catalog.fa.gz \
    > ./stacks_outputM2n3/catalog.sam
```

### Check the alignment

SAMStat makes an HTML report with the proportion of aligned loci and their mapping quality:

```bash
module load bioinfo-tools samtools samstat
samstat ./stacks_outputM2n3/catalog.sam
```

To open the report, copy it to your own computer. Run this in a terminal on your computer,
not on the cluster:

```bash
scp <user>@rackham.uppmax.uu.se:/path/to/stacks_outputM2n3/catalog.sam.samstat.html .
```

A low alignment rate is not unusual with a fragmented draft genome; loci that do not align
simply get no genomic position.

## 4. Add genomic positions to the de novo catalog

Convert the alignment to BAM:

```bash
module load bioinfo-tools samtools

mkdir -p ./bwa
samtools view -b ./stacks_outputM2n3/catalog.sam > ./bwa/catalog.bam
```

Then integrate the alignments into the de novo catalog:

```bash
module load bioinfo-tools Stacks python3

stacks-integrate-alignments \
    -P ./stacks_outputM2n3 \
    -B ./bwa/catalog.bam \
    -O ./integrate
```

`stacks-integrate-alignments` writes new catalog files with genomic coordinates to
`./integrate`. Copy these into the de novo output folder (`stacks_outputM2n3`), replacing
the originals. Keep a backup of the original files first. The folder is then ready for
`populations`.

## 5. Filter SNPs with populations

```bash
#!/bin/bash -l
#SBATCH -A <account>
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J stacks_populations

module load bioinfo-tools Stacks

populations -P ./stacks_outputM2n3 \
    --popmap ./all_popmap \
    -O ./pops.maf05.r75.het70 \
    -r 0.75 \
    --min_maf 0.05 \
    --max_obs_het 0.7 \
    --fstats -f p_value -k \
    --vcf --plink --genepop --structure
```

What the options do:

| Option | Meaning |
|---|---|
| `-r 0.75` | Keep a locus only if it is present in at least 75% of the samples in a population |
| `--min_maf 0.05` | Remove SNPs with a minor allele frequency below 0.05 (likely sequencing errors) |
| `--max_obs_het 0.7` | Remove SNPs with observed heterozygosity above 0.7 (likely merged paralogs) |
| `--fstats` | Calculate F-statistics and divergence between populations |
| `-f p_value` | Apply a p-value correction to the FST values |
| `-k` | Kernel-smoothed statistics along the genome; only possible because loci now have genomic positions |
| `--vcf --plink --genepop --structure` | Output formats for downstream programs |

Other useful options are `-p` (minimum number of populations a locus must occur in),
`--write_random_snp` (one random SNP per locus, for analyses that assume unlinked markers)
and `-W` (a whitelist of loci to keep).

## References

- Catchen J, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013). Stacks: an analysis tool
  set for population genomics. *Molecular Ecology* 22: 3124–3140.
- Paris JR, Stevens JR, Catchen JM (2017). Lost in parameter space: a road map for Stacks.
  *Methods in Ecology and Evolution* 8: 1360–1373.
- Stacks manual: https://catchenlab.life.illinois.edu/stacks/manual/

## Author

Janne Swaegers
