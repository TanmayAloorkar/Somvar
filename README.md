# Somvar

A Nextflow pipeline for somatic variant analysis of matched tumor–normal whole-exome sequencing (WES) data, with a Python terminal progress display and automatic resume.

**Why the name?** The author comes from a Marathi family, where **सोमवार (Somvar)** means **Monday**, an auspicious day. The name also combines **Som**atic and **Var**iants.

## Workflow

```text
Raw FastQC → fastp → Clean FastQC → BWA-MEM → BAM processing
→ BQSR (optional) → Mutect2 → Contamination estimation
→ Variant filtering → PASS variants → VEP → Summaries → MultiQC → Final report
```

Stages run in order. Normal and tumor tasks can run together when resources permit; Mutect2 analyzes both samples in one task. Supply `--known_sites` for BQSR and both `--germline_resource` and `--pon` for somatic calling. Without the latter two, the workflow stops after BAM processing or BQSR.

## Tested dataset

Tested end to end on this matched human WES pair from NCBI SRA:

| Sample | Alias | Accession |
| --- | --- | --- |
| Normal | FTA05N | [SRR3546766](https://www.ncbi.nlm.nih.gov/sra/?term=srr3546766) |
| Tumor | FTA05T | [SRR3546767](https://www.ncbi.nlm.nih.gov/sra/?term=srr3546767) |

The completed run used Broad GRCh38, BQSR, Mutect2, VEP 116, MultiQC 1.35, and 10 threads. It produced **114 PASS records**, including **89 represented in VEP output**. All 114 remain in the variant tables. These counts describe this run; variant-calling accuracy was not benchmarked.

## System requirements

- **Platform:** 64-bit Linux; the conda installer targets Linux x86_64, including compatible WSL2 installations.
- **RAM:** at least **16 GB**; additional headroom is recommended for Mutect2.
- **Storage:** SSD, roughly **150–200 GB** for the tested example, plus VEP cache storage. Its archive is approximately 26 GB; allow space for extraction.

| Threads per task | CPUs for two concurrent tasks plus headroom | Suggested RAM |
| --- | --- | --- |
| 4 (default) | 9 | 16-20 GB |
| 8 | 17 | 32 GB |
| 12 | 25 | 48 or more GB |

These are planning estimates; 8- and 12-thread configurations have not been benchmarked. `--threads` applies per task. Most tasks request 4 GB; BWA, BAM processing, and BQSR request 10 GB, Mutect2 14 GB, and VEP 8 GB. VEP uses at most four forks; reporting tasks use one CPU.

Nextflow schedules tasks against detected host resources. On a 16 GB machine, the larger normal and tumor tasks run sequentially. Local memory requests are scheduling estimates, not hard limits; Mutect2's 14 GB Java heap leaves limited headroom. RAM is configured in `scripts/main.nf`, with no launcher RAM option.

## Software setup

### Conda setup

Install [Miniforge](https://github.com/conda-forge/miniforge) or Miniconda, then run from the project root:

```bash
bash setup_conda.sh
```

The opening prompt asks whether to proceed. Enter `y` to create `somvar`, or press Enter to cancel. After successful installation:

```bash
conda activate somvar
```

[setup_conda.sh](setup_conda.sh) installs and checks the tools listed in [environment.yml](environment.yml): Python, Nextflow, Java, FastQC, fastp, BWA, samtools, GATK, bcftools, HTSlib, VEP, MultiQC, and download/shell utilities. References, FASTQs, and the VEP cache are prepared separately.

For another environment name, use `bash setup_conda.sh --name somvar_custom`. Existing environments are protected. Successful setup saves exact package builds to `conda-<name>-linux-64.lock.txt`; recreate them with `conda create --name somvar_copy --file conda-somvar-linux-64.lock.txt`.

Package resolution was checked with a conda dry run. A complete pipeline run in the conda environment has not yet been tested; the dataset above used local tool installations.

### Manual installation

Alternatively, install the tools on `PATH`. The launcher/reporting helper need Python 3.9+. Use the [Nextflow installation guide](https://docs.seqera.io/nextflow/install), [GATK releases](https://github.com/broadinstitute/gatk/releases), [VEP repository](https://github.com/Ensembl/ensembl-vep) (branch `release/116`, [installer](https://github.com/Ensembl/ensembl-vep/blob/release/116/INSTALL.pl)), and [MultiQC installation guide](https://docs.seqera.io/multiqc/getting_started/installation). Match VEP's software release to its cache.

## Reference preparation and download sources

Use Broad's **hg38/v0** reference and resources with matching contig names and coordinates. Save each file in the destination below; download both the VCF and its index. Sources: [Broad reference bundle](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0), GATK's somatic-hg38 resource bucket, and Ensembl release 116.

| File | Download | Destination |
| --- | --- | --- |
| `Homo_sapiens_assembly38.fasta` | [FASTA](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) | `references/hg38_broad_v0/` |
| `Homo_sapiens_assembly38.fasta.fai` | [FAI](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai) | Same folder |
| `Homo_sapiens_assembly38.dict` | [Dictionary](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict) | Same folder |
| `Homo_sapiens_assembly38.dbsnp138.vcf.gz` | [VCF](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz) · [TBI](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi) | `references/hg38_broad_v0/known_sites/` |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | [VCF](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz) · [TBI](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi) | Same known-sites folder |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | [VCF](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz) · [TBI](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi) | Same known-sites folder |
| `af-only-gnomad.hg38.vcf.gz` | [VCF](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz) · [TBI](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi) | `references/hg38_broad_v0/mutect2/` |
| `1000g_pon.hg38.vcf.gz` | [VCF](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz) · [TBI](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi) | Same Mutect2 folder |
| `small_exac_common_3.hg38.vcf.gz` | [VCF](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz) · [TBI](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi) | Same Mutect2 folder |
| `homo_sapiens_vep_116_GRCh38.tar.gz` | [Cache archive](https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz) | `~/references/vep/` |

Build the five BWA indexes once, then extract the cache:

```bash
bwa index -a bwtsw references/hg38_broad_v0/Homo_sapiens_assembly38.fasta
tar -xzf "$HOME/references/vep/homo_sapiens_vep_116_GRCh38.tar.gz" \
    -C "$HOME/references/vep"
```

Keep `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa` beside the FASTA. The cache root must contain `homo_sapiens/116_GRCh38/info.txt`; VEP reuses the alignment FASTA. The common-SNP VCF is supplied as the pileup intervals and requires biallelic SNPs with `INFO/AF`.

<details>
<summary>Recorded download sizes, checksums, and provenance</summary>

| File | Size (bytes) | MD5 |
| --- | --- | --- |
| `Homo_sapiens_assembly38.dbsnp138.vcf.gz` | 1,560,889,937 | `110aaeb8130bf4edb544a72d2c7829f7` |
| `Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi` | 2,321,422 | `d77e60d4520f27493a1a2c13c88b7626` |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | 61,692,306 | `14cc588a271951ac1806f9be895fb51f` |
| `Homo_sapiens_assembly38.known_indels.vcf.gz.tbi` | 1,567,886 | `1a55fdfa6533ae5cbc70e8188e779229` |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | 20,685,880 | `2e02696032dcfe95ff0324f4a13508e3` |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi` | 1,500,013 | `4c807e2cbe0752c0c44ac82ff3b52025` |
| `small_exac_common_3.hg38.vcf.gz` | 1,297,183 | `4c75c1755a45c64e8af7784db7fde009` |
| `small_exac_common_3.hg38.vcf.gz.tbi` | 242,095 | `f650d1dda6bd68cba65d77f131147985` |

MD5 values were recorded from Google Cloud Storage metadata. Broad known sites were downloaded on **28 September 2026**; common-SNP URLs/checksums and VEP cache metadata were verified on **1 October 2026**. The extracted cache identifies `homo_sapiens`, GRCh38, and source assembly GRCh38.p14. On the test machine, the archive and extracted cache were under `/home/tanmay/references/vep/`.

Additional source records: [common-SNP VCF metadata](https://storage.googleapis.com/storage/v1/b/gatk-best-practices/o/somatic-hg38%2Fsmall_exac_common_3.hg38.vcf.gz), [index metadata](https://storage.googleapis.com/storage/v1/b/gatk-best-practices/o/somatic-hg38%2Fsmall_exac_common_3.hg38.vcf.gz.tbi), [Ensembl cache directory](https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/), and [Ensembl checksum listing](https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/CHECKSUMS).

</details>

## Run

From the project root, after preparing the resources:

```bash
python3 run_somvar.py \
    --normalR1 tests/FTA05N_R1.fastq.gz \
    --normalR2 tests/FTA05N_R2.fastq.gz \
    --tumorR1 tests/FTA05T_R1.fastq.gz \
    --tumorR2 tests/FTA05T_R2.fastq.gz \
    --threads 4 \
    --known_sites references/hg38_broad_v0/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.gz,references/hg38_broad_v0/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz,references/hg38_broad_v0/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --germline_resource references/hg38_broad_v0/mutect2/af-only-gnomad.hg38.vcf.gz \
    --pon references/hg38_broad_v0/mutect2/1000g_pon.hg38.vcf.gz
```

Replace the FASTQ paths for other pairs and use a separate `--outdir` for each pair. The example uses four threads; retain `--threads 10` when resuming the previously completed test. Keep the comma-separated known-sites value free of spaces.

| Option | Default / purpose |
| --- | --- |
| `--threads` | `4`, per analysis task |
| `--ref` | Project `references/hg38_broad_v0/Homo_sapiens_assembly38.fasta` |
| `--outdir` | Project `results/` |
| `--known_sites` | Optional indexed VCFs for BQSR |
| `--germline_resource`, `--pon` | Both required for somatic calling |
| `--intervals` | Optional compatible capture-target BED/interval list; recommended for WES |
| `--common_sites` | Project `references/hg38_broad_v0/mutect2/small_exac_common_3.hg38.vcf.gz` |
| `--vep_cache`, `--vep_cache_version` | `~/references/vep`, `116` |
| `--skip_vep` | Stop after PASS extraction |
| `--no-resume` | Disable automatic resume |
| `--multiqc_config` | Project `conf/multiqc_config.yaml` |

Explicit arguments override defaults. Run `python3 run_somvar.py --help` for launcher help. Additional Nextflow options are forwarded, including `-with-trace`; direct execution uses `nextflow run scripts/main.nf -resume` with the same input/resource arguments.

## Results and progress

Start with **`results/reports/index.html`**. Keep the entire `reports/` folder together when sharing it.

| Output folder | Contents |
| --- | --- |
| `fastqc_raw/`, `fastp/`, `fastqc_clean/` | Raw/clean QC and trimmed reads, separated into `normal/` and `tumor/` |
| `bwa/`, `bam/`, `bqsr/` | Alignments, duplicate-marked BAMs, and optional recalibrated BAMs |
| `mutect2/`, `contamination/`, `filtered/`, `pass/` | Calling, contamination tables, filtered and PASS VCFs |
| `vep/` | Annotated VCF, summary, warnings, and skipped variants |
| `summary/` | Variant TSV/CSV, transcript consequences, gene/consequence counts, and JSON summary |
| `multiqc/` | MultiQC HTML and data exports |
| `reports/` | Portable report bundle, downloadable VCFs/tables, and run provenance |

Counts refer to VCF records; multiallelic records count once. Every PASS record remains in the summary, including those omitted by VEP. Allele frequencies and depths are copied from the VCF. HTML search covers up to 1,000 displayed records; downloads contain the full table. Orientation-bias priors are not currently included in filtering.

The launcher timestamps stage changes and shows Mutect2's latest `ProgressMeter` line. Most bars advance when tasks finish, so long steps can appear paused. Driver logs and settings are saved under `logs/<run-id>/`; task logs remain under `work/`. Keep **`work/` and `.nextflow/`** for resume. Changing inputs, tool commands, resources, or threads can invalidate cached tasks.

For failures, check the latest `.nextflow.log` and the failed task's `.command.err`. Confirm tools are on `PATH`, VCFs have neighboring `.tbi`/`.idx` indexes, and the VEP cache matches its software release. Investigate memory usage for `Killed`/exit 137 errors.

**Project files:** [scripts/main.nf](scripts/main.nf), [run_somvar.py](run_somvar.py), [scripts/reporting.py](scripts/reporting.py), [setup_conda.sh](setup_conda.sh), [environment.yml](environment.yml), and [conf/multiqc_config.yaml](conf/multiqc_config.yaml). Nextflow closures use explicitly named parameters to avoid implicit-`it` deprecation warnings.
