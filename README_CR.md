# Unraveling theGenetics of Maternal Sex Determination in the Blowfly *Chrysomya rufifacies*

---

## Objectives;

1. to identify any **genetic differences between the two different types** of **females** and understand the genetic basis of maternal sex determination in this species,

2. to identify **differentially expressed genes** that might be involved in sex determination,

3. to determine if the putative sex determining region is found in an **area of reduced recombination** (as has been suggested for other monogenic taxa) 

4. to explore any **potential genetic differences between the sexes**

--- 

# <u>Chapter 0</u> Prepartion

---

## 0.1 Quality Control &Trimming

> The **main aim** of this step is to "produce" the best possible quality of reads from the raw reads, in order to use them for dowstream analysis

**0.1.1 Illumina Reads**

In this step we used;

- `fastqc` *v.0.12.1* ([Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), to assess the quality of the raw illumina data before (pre-trimming) and after the treaming (post-trimming). 

- `fastp` ([GitHub - OpenGene/fastp: An ultra-fast all-in-one FASTQ preprocessor (QC/adapters/trimming/filtering/splitting/merging...)](https://github.com/OpenGene/fastp)) *v. 0.22.0* to trimm, quality control and read merge the illumina raw reads. 

##### Directories

| Type of Directory | Directory                                                                             | Description              |
| ----------------- | ------------------------------------------------------------------------------------- | ------------------------ |
| Input             | `/data/ross/sequencing/raw/blowflies/picard_lab_illumina`                             | raw reads                |
|                   | `/data/ross/flies/raw/Chrysomya_rufifacies/illumina/`                                 | trimmed reads            |
| Scripts           | `/data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/scripts/`       |                          |
| Output            | `/data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/outputs/01_QC/` | QC (.html files)         |
|                   | `/data/ross/flies/raw/Chrysomya_rufifacies/illumina/`                                 | trimmed reads (fastq.gz) |

##### Analysis

Both `fastp `and `fastqc `versions used in this study were downloaded by using *anaconda* (https://anaconda.org/). 

**Pre-trim QC**

Perform a QC analysis with `fastqc`BEFORE trimming takes place, in order in the future (by performing a QC AFTER trimming) to meassure the effect of trimming in quality of the reads.

```bash
fastqc -t 4 *.fastq.gz 
```

    - * is refering to all the fastq.qz files present in the directory

**Trimming etc**. 

Illumina sequencing reads often contain adapters, with low-quality bases at the 3' end of the read, which can lead to reduced mapping efficiency, assembly errors, and decreased accuracy in downstream analysis. To avoid that we **trimmed** the Illumina reads by using `fastp`. 

*Furthermore*,  `fastp` program performs quality control and read merging. **Quality control** includes filtering out reads with low-quality scores. **Read merging** involves merging the paired-end reads into longer sequences if they overlap, which can increase the accuracy and length of the final sequences. In our cases we do not expect any read merging, as the DNA fragments that were sequenced were >300bp.

```bash
fastp -i <forward reads>.fastq.gz -I <reverse reads>.fastq.gz -o <forward reads>.trimmed.fastq.gz -O <reverse reads>.trimmed.fastq.gz
```

     -i and - I specify input files (forward and reverse)

    -o and -O specifies the output files (forward and reverse)

**Post-trimm QC** 

Use `fastqc` again on the trimmed Illumina reads. (see pre trim QC for why)

```bash
fastqc -t 4 *.trimmed.fastq.gz
```

    -  * is refering to all the .**trimmed**.fastq.qz files present in the directory

<u>Output;</u> 

- The plots belows are representative of all the available, so for saving space, we only depict the plots from the forward reads of one individual (TF19 R1)

<img title="" src="file:///C:/Users/Elpida/AppData/Roaming/marktext/images/2023-04-21-16-31-50-image.png" alt="" width="525">

Per sequence GC content plot | We notice the pressence of three peaks insdead of one. That suggests the presence of **condamination**, so we decided to move forward examining this scenarion by using `bloptools`.

<img title="" src="file:///C:/Users/Elpida/AppData/Roaming/marktext/images/2023-04-21-16-50-54-image.png" alt="" width="371">

Sequence Duplication Levels | We notice that there is a 10% of the total sequence with >10k duplicate sequences. That could suggest the presence of a PCR step before the illumna sequence step. The presence of a PCR step could be visible in the kmer plots, as a small hill (due to the 10% part) at the far end of the x axis (due to the >10k part). 

<img title="" src="file:///C:/Users/Elpida/AppData/Roaming/marktext/images/2023-04-21-16-57-04-image.png" alt="" width="380">

Adapter content | the adapter remove part was successful for the reads. As in all cases in the trimmed reads the ONLY visible lines are the ones that represent the polyA and polyG.

## 0.2 Genome Assembly

> The **aim** of this step was to create two assembled genomes, one for the AF individuals and one for the TF. For the TF individuals we used 2 illumina libraries to create the assemblies.

**0.2.1 Illumina Reads**

In this step we used; 

- `spades` *v. 3.15.5* ([GitHub - ablab/spades: SPAdes Genome Assembler](https://github.com/ablab/spades)), in order to assembly the Illumina trimmed reads.

**Directories**

| Type of Directory | Directory                                                                  | Discription                                          |
| ----------------- | -------------------------------------------------------------------------- | ---------------------------------------------------- |
| Input             | `/data/ross/flies/raw/Chrysomya_rufifacies/illumina/`                      | trimmed (.fastq.gz)                                  |
| Scripts           | `/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/scripts` | 2 scripts; 1 for TF (female prod) & 1 AF (male prod) |
| Output            | `/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs` | 2 (AF and TF) folders with spades outputs            |

##### Analysis

I did not perform any correction, that were suggested in SPAdes website because of memmory problems, plus illumina is high accurate allready.

Importantly, `spades` was difficult to run in the cluster, as it demanded a lot of memmory. I manged to accomodate its needs, thus make it run ONLY  when using this command to run the script in ths cluster. 

```bash
qsub -pe smp64 20 -N spades_AF -l mem_free=500g, h_vmem=500g spades_AF.sh
```

the special thing in this command are the `mem_free=500g, h_vmem=500g ` flags, which request 500 GB of free memory and 500 GB of virtual memory. Thus the job is not going to start running unless this amount of memmory is available and "booked" for the job. In that way spades, will not run out of memmory (a really frequent problem). 

For the memory problem to be fully resolved you need to specify how much memory should be used by `spades`. You do it by using the `-m` flag as done below.

```bash
spades.py --pe1-1 ./${TF19}R1_001.trimmed.fastq.gz --pe1-2 ./${TF19}R2_001.trimmed.fastq.gz --pe2-1 ./${TF11}R1_001.trimmed.fastq.gz --pe2-2 ./${TF11}R2_001.trimmed.fastq.gz -o . --threads 16 --only-assembler -m 400000
```

The above command was used to assembly the genome for the TF individuals, for which two libraries (TF19 and TF11) were used. A similar command was used for the AF assembly.

- `spades.py`: This is the command to run SPAdes.
- `--pe1-1 ./${TF19}R1_001.trimmed.fastq.gz`: This option specifies the path to the forward paired-end read file for the first sample (`${TF19}` is a variable that  contains the sample name).
- `--pe1-2 ./${TF19}R2_001.trimmed.fastq.gz`: This option specifies the path to the reversed paired-end read file for the first sample.
- `--pe2-1 ./${TF11}R1_001.trimmed.fastq.gz`: This option specifies the path to the forward paired-end read file for the second sample (`${TF11}` is a variable that contains the sample name).
- `--pe2-2 ./${TF11}R2_001.trimmed.fastq.gz`: This option specifies the path to the second paired-end read file for the second sample.
- `-o .`: This option specifies the output directory for the assembly results. In this case, the output will be saved in the current directory (`.`).
- `--threads 16`: This option specifies the number of CPU threads to use for the assembly process. In this case, 16 threads will be used.
- `--only-assembler`: This option specifies that only the assembly process should be run, without any additional steps such as error correction or scaffolding.
- `-m 400000`: This option specifies the amount of memory to be used by the assembly process. In this case, 400,000 MB of memory (i.e., 400 GB) will be allocated.

Check the quality of the sequences using `Assembly Stats` *v.1.0.1* ([GitHub - sanger-pathogens/assembly-stats: Get assembly statistics from FASTA and FASTQ files](https://github.com/sanger-pathogens/assembly-stats))

*AF assembly*

```bash
$assembly-stats contigs.fasta

stats for contigs.fasta
sum = 661265216, n = 2073607, ave = 318.90, largest = 785956
N50 = 780, n = 154289
N60 = 412, n = 273699
N70 = 245, n = 488677
N80 = 155, n = 829525
N90 = 108, n = 1297573
N100 = 78, n = 2073607
N_count = 0
Gaps = 0
```

There are 661265216 nucleotides in total **(3x10^8 more than the genome size in A. Andere et al. 2020 and 2x10^8 from the expected)**. Could be another clue for condamination 

In general really small contigs. Expected for the Illumina data.

*TF assembly*

```bash
$assembly-stats contigs.fasta
stats for contigs.fasta
sum = 638926939, n = 1987717, ave = 321.44, largest = 667008
N50 = 982, n = 139964
N60 = 556, n = 226790
N70 = 332, n = 378571
N80 = 155, n = 656371
N90 = 90, n = 1201263
N100 = 78, n = 1987717
N_count = 0
Gaps = 0
```

similar results with the one above 

## 0.3 Condamination detection

> From the QC chapter we got the suspision for the presence of DNA condamition (GC % plot). As a result we decided to go forward with examing this possibility by using `blobtools`.  `blobtools` (https://blobtools.readme.io/docs/my-first-blobplot) need two type of input data ( 1. coverage file(s) and 2. hit(s) file(s)). 

In this step we used; 

- `bowtie2` v. 2.2.5 ([GitHub - BenLangmead/bowtie2: A fast and sensitive gapped read aligner](https://github.com/BenLangmead/bowtie2)) an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
- `samtools` v. 1.6 ([GitHub - samtools/samtools: Tools (written in C using htslib) for manipulating next-generation sequencing data](https://github.com/samtools/samtools)) tools for handling SAM, BAM, CRAM
- `blastn` v. 2.12.0 ([Command Line BLAST &#8211; A Primer for Computational Biology](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/)) 
- **blobtools**
  - conda enviroment which was created by executing **Option A** as depicted in [GitHub - DRL/blobtools: Modular command-line solution for visualisation, quality control and taxonomic partitioning of genome datasets](https://github.com/DRL/blobtools) + conda download samstools (which have to be the same version as the one used to create the .bam files in the coverage files step)
    - some important packages used `pysam` v. 0.15.3 and `python` v. 3.7.16 (at least 3.7, so `pysam` works)
  - `/ceph/users/eskarlou/blobtools` downloaded by using 
    `git clone https://github.com/DRL/blobtools.git` on 25th April 2023

**Directories**

| Type of Directory | Directory                                                                                                                                                                                                                                                                                                                                                                                                                         | Description                                       |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------- |
| Input             | **for coverage files**<br/>`/data/ross/flies/raw/Chrysomya_rufifacies/illumina/`<br/>`/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/`<br/>**for hits files** <br/>`/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs/`<br/>**for blobplots**<br/> - coverage & hits files (see above)<br/> - NCBI taxonomy names + nodes files `/data/ross/flies/analyses/blowflies/05_NCBI_taxdump` |                                                   |
| Script            | `/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/scripts/AF`<br/> `/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/scripts/TF`                                                                                                                                                                                                                                                            | 3 scripts (2 for input and one for the blobtools) |
| Output            | **for coverage files** <br/>`/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out or TF_out/coverage_file` <br/> **for hits files** <br/>`/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out or TF_out/hits_file`<br/>**for blobplots**<br/> `/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out or TF_out/pre_partition/`          |                                                   |

Workflow A was used in this project ![4476](https://files.readme.io/21ef680-figure_1.png)

## Input files for blobplots

### 1st Coverage files

The file (e.g. `mapping_1.bam`) contains information regarding the base/read coverage of each sequence in an [assembly file](https://blobtools.readme.io/docs/assembly-file). Assuming an unbiased sequencing process, the base/read coverage reflects the **molarity**of the DNA molecule (represented by the sequence in the assembly) that went into the sequencing reaction.

*The example here is for AF but the same code was used for TF*

Firstly, index the asssembly genome (e.g. `contig.fasta`) by using bowtie2. 

```bash
bowtie2-build contigs.fasta AFcontigs_index
```

Then, align the reads used for the de-novo assembly (i.e. `contig.fasta`) to the assembly itself.  The output is the `.sam` file.

```bash
bowtie2 -x AFcontigs_index \
-1 AF7_Chrysomya-rufifacies_S2_R1_001.trimmed.fastq.gz \
-2 AF7_Chrysomya-rufifacies_S2_R2_001.trimmed.fastq.gz \
-S AFmappings.sam
```

convert the `.sam` file to a `.bam` one by using `samtools`

```bash
samtools view -b -o AFmappings.bam AFmappings.sam
```

The next two final commands take the `.bam` file and create a sorted and indexed BAM file in order to create the coverage file "mapping.sorted.bam".

```bash
samtools sort -o AFmapping.sorted.bam AFmappings.bam
samtools index AFmapping.sorted.bam
```

### 2nd Hit files

(https://blobtools.readme.io/docs/taxonomy-file)

For `blobtools` to assign taxonomies and thereafter "represent" the contamination, it uses **hit files** (`blast.out`) as an input.

To produce those hits files (in TSV format; more in https://blobtools.readme.io/docs/taxonomy-file) we performed a `blastn`. 

- Important coding tip; `blast` default version was 2.5.0 but it didn't work with the database available in `/ceph/software/databases/ncbi` (which we used for our `blastn`). To update `blastn` run `conda install -c bioconda blast==2.12.0`,  the important part is the `==` instead of `=`.

```bash
blastn -query contigs.fasta -db /ceph/software/databases/ncbi/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 4
```

- `blastn`: This is the name of the BLASTN program that is being run.

- `-query contigs.fasta`: This specifies the name of the query file, which contains the sequences to be searched against the database.

- `-db /ceph/software/databases/ncbi/nt`: This specifies the name and location of the database that will be searched, which is the NCBI nt nucleotide database available/downloaded in Ashworth cluster.

- `-outfmt '6 qseqid staxids bitscore std'`: This specifies the format of the output file. In this case, the output will be in tab-delimited format, and will include the query sequence ID (qseqid), the taxonomic IDs of the hit sequences (staxids), the bitscore, and the standard deviation of the bitscore (std).

- `-max_target_seqs 1`: This specifies the maximum number of hits to return for each query sequence. In this case, only the top hit will be returned.

- `-max_hsps 1`: This specifies the maximum number of high-scoring segment pairs (HSPs) to return for each hit. In this case, only the top HSP will be returned.

- `-evalue 1e-25`: This specifies the e-value threshold for reporting hits. Hits with e-values greater than this threshold will not be reported.

- `-num_threads 4`: This specifies the number of threads or processors to use for the search. In this case, 4 threads will be used.

### 3rd NCBI taxonomy names + nodes files

- `--names`: This flag is used to specify the NCBI taxonomy names file (`names.dmp`) that contains the scientific names of the taxonomic groups.
- `--nodes`: This flag is used to specify the NCBI taxonomy nodes file (`nodes.dmp`) that provides the hierarchical structure of the taxonomic groups.

In order to get the files needed by flags `--names` and `--nodes` in `blobtools create`

we downloaded the `taxdump.tar.gz`

```bash
# paste the link on your internet browser
//ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# upload taxdump.tar.gz in the cluster 

# Extract files 
tar -xvzf taxdump.tar.gz
```

The `names.dmp` and `nodes.dmp` files will be located in the directory used 

## Blobtools; Pre-partition

Based on the **workflow A** we executed *step 1* and *2*. 

***Step 1***

```bash
./blobtools create -i assembly.fasta -b sorted.bam -t blast.out -o my_blobplot --names names.dmp --nodes nodes.dmp
```

- `-i assembly.fasta`: Use the `assembly.fasta` file as the input assembly file.
- `-b sorted.bam`: Use the `sorted.bam` file as the input sorted alignment file in BAM format. The .bam.bai file should be in the same directory, when the command is executed
- `-t blast.out`: Use the `blast.out` file as the input file containing blast results.
- `-o my_blobplot`: Save the output files to the directory named `my_blobplot`.
- `--names names.dmp`: Use the `names.dmp` file as the taxonomic names file.
- `--nodes nodes.dmp`: Use the `nodes.dmp` file as the taxonomic nodes file.

Output; 

- `myblobplot.blobDB.json`  is a JSON (JavaScript Object Notation) file that contains a database of information about the genomic data analyzed by `blobtools`. This database contains information on the sequences in the genome assembly, including their length, taxonomic classification, and read coverage.

***Step 2***

```bash
# create Blobplot and Cov plot
./blobtools plot -i my_blobplot.blobDB.json -o ./
# create table
./blobtools view -i my_blobplot.blobDB.json -o ./
```

Output (./blobtools plot);

- `.blobplot.bam0.png` blopplot

- `blobplot.read_cov.bam0.png` readcov plot

- `.cov` cov plot

- `.blobplot.stats.txt` summary info of the above

**AF Blobplot + Cov read**

![](C:\Users\Elpida\AppData\Roaming\marktext\images\2023-05-11-15-37-34-image.png)

![](C:\Users\Elpida\AppData\Roaming\marktext\images\2023-05-11-15-38-07-image.png)

**TF Blobplot + read plot**

![](C:\Users\Elpida\AppData\Roaming\marktext\images\2023-05-11-18-16-12-image.png)

![](C:\Users\Elpida\AppData\Roaming\marktext\images\2023-05-11-18-17-17-image.png)

In both of them, more than half of the assembly is a no-hit. But they have the same GC content as the arthropoda. So it should be that those contigs fly dna but the contigs are so small (N50 ~1000 nt) to be taxonomized.

Question rised; 

- the biased GC content find previously was biased on the left not the right like here...

- the contamination does not seem to be that much (3%)

- the problem is the no-hit which I guess has to do with the fact that we are talking about an illumina assembly...

-  should focus on the cleaning reads in order to have a descent kmer plot

Output (./blobtools view);

- `blobDB.table.txt`

# <u>Chapter 1</u>

# Differences between the two types of females

---

## 1.1 Kmer analysis

> the main aim of this part was to examine the heterozygosity levels of the females and check our hypothesis of the presence of a inversion in the female producing females

In this step we used

- `kmc`  *v. 3.2.1* ([GitHub - refresh-bio/KMC: Fast and frugal disk based k-mer counter](https://github.com/refresh-bio/KMC)). KMC is a disk-based program for counting k-mers from (possibly gzipped) FASTQ/FASTA files.

- `genomescope.R` ([GitHub - tbenavi1/genomescope2.0: Reference-free profiling of polyploid genomes](https://github.com/tbenavi1/genomescope2.0)). GenomeScope 2.0 uses the k-mer count distribution, e.g. from KMC or Jellyfish, and produces a report and several informative plots describing the genome properties.

#### 1.1.1 Use raw reads (contamination + PCR)

**Directories**

| Type of Directory | Directory                                                                | Description                    |
| ----------------- | ------------------------------------------------------------------------ | ------------------------------ |
| Input             | `/data/ross/sequencing/raw/blowflies/picard_lab_illumina`                | raw illumina reads             |
| Script            | `/data/ross/flies/analyses/blowflies/03_kmer/01_illumina_reads/scripts/` | script                         |
| Output            | `/data/ross/flies/analyses/blowflies/03_kmer/01_illumina_reads/outputs`  | output from raw illumina reads |

##### Analysis

With this analysis we wanted to explore the structure of the genome of the two types of females. Using k-mer distributiond plots (provided by `genomescope`) we can make inference regarding the aploidy and the heteroplasmy of the organism. Each peak represent a level of ploidy e.g. two peaks = diploid organism. In the case of diploid organisms, the frequency hight of each peak can give us information about the leve of heterogenity. So in our case, if there is an inversion that includes the sex determing area,  that inversion should be more diverge/heterozygous that the rest of the non-inversed genome. Thus, we expect that the first peak (heterozygosity peak) will be higher in individuals that posses the inversion (TF) than the ones that do not (AF). See figure below.

<img title="" src="file:///C:/Users/Elpida/AppData/Roaming/marktext/images/2023-04-27-14-35-25-image.png" alt="" width="547">

Firstly, we needed to compute the histogram of k-mer frequencies. For that we used `kmc` inside the below `for` loop.

- ```bash
  # A for loop, for all the 
  for file in $(ls *1_001.fastq.gz)
  do
  # In order to keep the part of the name of each file that diferentiates 
  # them from the rest. 
  base=$(basename $file "1_001.fastq.gz")
  cat ${base}1_001.fastq.gz ${base}2_001.fastq.gz > ${base}_files.fastq.gz
  # count  kmers
  kmc -k21 -t10 -m64 -ci1 -cs10000 -fq1 ${base}_files.fastq.gz ${base}_kmer_counts .  # use fq1 for fastq.gz input files
  kmc_tools transform ${base}_kmer_counts histogram ${base}_kmer_k21.histo -cx100000 
  done
  ```

- `for file in $(ls *1_001.fastq.gz)`: sets up a loop that iterates over all the files in the current directory that end with "1_001.fastq.gz". The `ls` command lists all the files that match the pattern, and the `$(...)` syntax captures the output of the command and uses it as the list for the loop.

- `base=$(basename $file "1_001.fastq.gz")`: extracts the base file name from the current file by removing the "1_001.fastq.gz" suffix. The `basename` command takes the file name as its argument and the suffix to be removed as the second argument.

- `cat ${base}1_001.fastq.gz ${base}2_001.fastq.gz > ${base}_files.fastq.gz`: concatenates the paired-end FASTQ files for the current sample into a single file with a name derived from the base file name and a suffix of "_files.fastq.gz". The output file contains all the reads from both paired-end files.

- `kmc -k21 -t10 -m64 -ci1 -cs10000 -fq1 ${base}_files.fastq.gz ${base}_kmer_counts`: runs KMC on the concatenated FASTQ file for the current sample to count 21-mers. The output file name is derived from the base file name and a suffix of "_kmer_counts".
  
  - "We recommend using a **k-mer length of 21** for most genomes, as this length is sufficiently long that most k-mers are not repetitive and is short enough that the analysis will be more robust to sequencing errors."

- `kmc_tools transform ${base}_kmer_counts histogram ${base}_kmer_k21.histo -cx100000`: processes the output of KMC using `kmc_tools` to generate a histogram of k-mer frequencies. The histogram file name is derived from the base file name and a suffix of "_kmer_k21.histo". The `-cx100000` option specifies the cutoff for the maximum frequency to be included in the histogram.

The we created a `GenomScope` plot for each library by executing the code below.

```bash
/ceph/users/eskarlou/miniconda3/envs/for_genomescope/bin/Rscript genomescope.R -i TF11_Chrysomya-rufifacies_S3_R_kmer_k21.histo -o ./TF11_plots -k 21 -l 34
```

- `/ceph/users/eskarlou/miniconda3/envs/for_genomescope/bin/Rscript`, in order to use Rscript. The path before before `/Rscript` is the directory of the enviroment that has the Rscript. 

- `genomescope.R`, the R script that is used to create all the plots and outputs of GenomeScope. In order to be used should be present in the directory were the script is executed.

- `-i TF11_Chrysomya-rufifacies_S3_R_kmer_k21.histo` input .histo file

- `-o ./TF11_plots` output directory

- `-k 21` k-mer lenght

- `-l 34`  initial guess for the average k-mer coverage of the sequencing

**TF individual**

<img src="file:///C:/Users/Elpida/AppData/Roaming/marktext/images/2023-05-11-11-21-08-image.png" title="" alt="" width="335">

**AF individuals**

<img title="" src="file:///C:/Users/Elpida/AppData/Roaming/marktext/images/2023-05-11-11-20-34-image.png" alt="" width="344">

- Regarding the two peaks right of the main one; (from [UCD Bioinformatics Core Workshop](https://ucdavis-bioinformatics-training.github.io/2020-Genome_Assembly_Workshop/kmers/kmers)) ->  *are the duplicated heterozygous regions and duplicated homozygous regions and correspond to two smaller peaks. The shape of these peaks are affected by the sequencing errors and sequencing duplicates.*
- The assumption of the presence of PCR does not allow us to make any inferences. As PCR can induce non-biological uneven coverage.
