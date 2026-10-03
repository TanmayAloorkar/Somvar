params.normalR1 = null
params.normalR2 = null
params.tumorR1 = null
params.tumorR2 = null
params.threads = 4
params.ref = "${projectDir.parent}/references/hg38_broad_v0/Homo_sapiens_assembly38.fasta"
params.outdir = "${projectDir.parent}/results"
params.known_sites = null
params.germline_resource = null
params.pon = null
params.intervals = null
params.common_sites = "${projectDir.parent}/references/hg38_broad_v0/mutect2/small_exac_common_3.hg38.vcf.gz"
params.vep_cache = "${System.getProperty('user.home')}/references/vep"
params.vep_cache_version = 116
params.skip_vep = false
params.multiqc_config = "${projectDir.parent}/conf/multiqc_config.yaml"

process FASTQC {
    cpus (params.threads as int)
    memory '4 GB'
    publishDir { "${params.outdir}/fastqc_raw/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path '*', emit: reports

    script:
    """
    fastqc --threads ${task.cpus} ${reads.join(' ')}
    """
}

process FASTP {
    cpus (params.threads as int)
    memory '4 GB'
    publishDir { "${params.outdir}/fastp/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val raw_qc_done

    output:
    tuple val(sample), path("${reads[0].simpleName}.trimmed.fastq.gz"), path("${reads[1].simpleName}.trimmed.fastq.gz"), emit: trimmed_reads
    path '*_fastp.html'
    path '*_fastp.json'

    script:
    """
    fastp --thread ${task.cpus} \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${reads[0].simpleName}.trimmed.fastq.gz \
        -O ${reads[1].simpleName}.trimmed.fastq.gz \
        --html ${sample}_fastp.html \
        --json ${sample}_fastp.json
    """
}

process FASTQC_CLEAN {
    cpus (params.threads as int)
    memory '4 GB'
    publishDir { "${params.outdir}/fastqc_clean/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)
    val fastp_done

    output:
    path '*_fastqc.html', emit: html
    path '*_fastqc.zip'

    script:
    """
    fastqc --threads ${task.cpus} ${read1} ${read2}
    """
}

process BWA_MEM {
    tag "${sample}"
    cpus (params.threads as int)
    memory '10 GB'
    publishDir { "${params.outdir}/bwa/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)
    path ref
    path index_files
    val raw_qc_done
    val clean_qc_done

    output:
    tuple val(sample), path("${sample}.sam"), emit: alignments

    script:
    """
    bwa mem -t ${task.cpus} -K 10000000 \
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA' \
        ${ref} ${read1} ${read2} > ${sample}.sam
    """
}

process BAM_PROCESSING {
    tag "${sample}"
    cpus (params.threads as int)
    memory '10 GB'
    publishDir { "${params.outdir}/bam/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.markdup.bam"), path("${sample}.markdup.bam.bai"), emit: indexed_bams
    tuple val(sample), path("${sample}.duplicate_metrics.txt"), emit: duplicate_metrics

    script:
    def extra_threads = Math.max(0, (task.cpus as int) - 1)
    """
    samtools sort \
        -@ ${extra_threads} \
        -m 512M \
        -T ${sample}.sort_tmp \
        -o ${sample}.sorted.bam \
        ${sam}

    mkdir -p markdup_tmp
    gatk --java-options "-Xmx3g -XX:ActiveProcessorCount=${task.cpus}" MarkDuplicates \
        -I ${sample}.sorted.bam \
        -O ${sample}.markdup.bam \
        -M ${sample}.duplicate_metrics.txt \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX false \
        --TMP_DIR markdup_tmp

    samtools index -@ ${extra_threads} ${sample}.markdup.bam
    samtools quickcheck -v ${sample}.markdup.bam
    """
}

process BQSR {
    tag "${sample}"
    cpus (params.threads as int)
    memory '10 GB'
    publishDir { "${params.outdir}/bqsr/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path ref_fai
    path ref_dict
    path known_vcfs
    path known_indexes

    output:
    tuple val(sample), path("${sample}.bqsr.bam"), path("${sample}.bqsr.bam.bai"), emit: recalibrated_bams
    path "${sample}.recal.table", emit: tables

    script:
    def extra_threads = Math.max(0, (task.cpus as int) - 1)
    def sites = (known_vcfs instanceof List ? known_vcfs : [known_vcfs])
        .collect { site_file -> "--known-sites ${site_file}" }.join(' ')
    """
    gatk --java-options "-Xmx6g -XX:ActiveProcessorCount=${task.cpus}" BaseRecalibrator \
        -R ${ref} \
        -I ${bam} \
        ${sites} \
        -O ${sample}.recal.table

    gatk --java-options "-Xmx6g -XX:ActiveProcessorCount=${task.cpus}" ApplyBQSR \
        -R ${ref} \
        -I ${bam} \
        --bqsr-recal-file ${sample}.recal.table \
        --create-output-bam-index false \
        -O ${sample}.bqsr.bam

    samtools index -@ ${extra_threads} ${sample}.bqsr.bam
    samtools quickcheck -v ${sample}.bqsr.bam
    """
}

process MUTECT2 {
    cpus (params.threads as int)
    memory '14 GB'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    tuple path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path ref
    path ref_fai
    path ref_dict
    path germline_vcf
    path germline_index
    path pon_vcf
    path pon_index
    path interval_files

    output:
    path 'somatic.unfiltered.vcf.gz', emit: variants
    path 'somatic.unfiltered.vcf.gz.stats', emit: stats
    path 'somatic.unfiltered.vcf.gz.tbi', optional: true, emit: index

    script:
    def target_intervals = interval_files ? "-L ${interval_files[0]}" : ''
    """
    gatk --java-options "-Xmx14g -XX:ActiveProcessorCount=${task.cpus}" Mutect2 \
        -R ${ref} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -normal normal \
        --germline-resource ${germline_vcf} \
        --panel-of-normals ${pon_vcf} \
        ${target_intervals} \
        --native-pair-hmm-threads ${task.cpus} \
        -O somatic.unfiltered.vcf.gz
    """
}

process GET_PILEUP_SUMMARIES {
    tag "${sample}"
    cpus (params.threads as int)
    memory '4 GB'
    publishDir "${params.outdir}/contamination", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref
    path ref_fai
    path ref_dict
    path common_vcf
    path common_index
    val mutect2_done

    output:
    tuple val(sample), path("${sample}.pileups.table"), emit: pileups

    script:
    """
    gatk --java-options "-Xmx2g -XX:ActiveProcessorCount=${task.cpus}" GetPileupSummaries \
        -R ${ref} \
        -I ${bam} \
        -V ${common_vcf} \
        -L ${common_vcf} \
        -O ${sample}.pileups.table
    """
}

process CALCULATE_CONTAMINATION {
    cpus (params.threads as int)
    memory '4 GB'
    publishDir "${params.outdir}/contamination", mode: 'copy'

    input:
    tuple path(normal_pileups), path(tumor_pileups)

    output:
    path 'tumor.contamination.table', emit: contamination
    path 'tumor.segments.table', emit: segments

    script:
    """
    gatk --java-options "-Xmx2g -XX:ActiveProcessorCount=${task.cpus}" CalculateContamination \
        -I ${tumor_pileups} \
        --matched-normal ${normal_pileups} \
        --tumor-segmentation tumor.segments.table \
        -O tumor.contamination.table
    """
}

process FILTER_MUTECT_CALLS {
    cpus (params.threads as int)
    memory '4 GB'
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path unfiltered_vcf
    path unfiltered_indexes
    path mutect_stats
    path contamination
    path segments
    path ref
    path ref_fai
    path ref_dict

    output:
    tuple path('somatic.filtered.vcf.gz'), path('somatic.filtered.vcf.gz.tbi'), emit: variants
    path 'somatic.filtering_stats.tsv', emit: stats

    script:
    """
    gatk --java-options "-Xmx2g -XX:ActiveProcessorCount=${task.cpus}" FilterMutectCalls \
        -R "${ref}" \
        -V "${unfiltered_vcf}" \
        --stats "${mutect_stats}" \
        --contamination-table "${contamination}" \
        --tumor-segmentation "${segments}" \
        --filtering-stats somatic.filtering_stats.tsv \
        --create-output-variant-index true \
        -O somatic.filtered.vcf.gz
    """
}

process PASS_VARIANTS {
    cpus (params.threads as int)
    memory '4 GB'
    publishDir "${params.outdir}/pass", mode: 'copy'

    input:
    tuple path(filtered_vcf), path(filtered_index)

    output:
    tuple path('somatic.pass.vcf.gz'), path('somatic.pass.vcf.gz.tbi'), emit: variants

    script:
    def extra_threads = Math.max(0, (task.cpus as int) - 1)
    """
    bcftools view --threads ${extra_threads} \
        -f PASS -Oz -o somatic.pass.vcf.gz "${filtered_vcf}"
    bcftools index --threads ${extra_threads} -t somatic.pass.vcf.gz
    """
}

process VEP {
    cpus (params.threads as int)
    memory '8 GB'
    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    tuple path(pass_vcf), path(pass_index)
    path ref
    path ref_fai
    path vep_cache
    val cache_version

    output:
    tuple path('somatic.pass.vep.vcf.gz'), path('somatic.pass.vep.vcf.gz.tbi'), emit: variants
    path 'somatic.pass.vep.summary.html', emit: summary
    path 'vep.warnings.txt', emit: warnings
    path 'vep.skipped_variants.txt', emit: skipped

    script:
    // Limit annotation forks to four for the fixed 8 GB memory request.
    def annotation_forks = Math.min(task.cpus as int, 4)
    """
    touch vep.warnings.txt vep.skipped_variants.txt
    if [ \$(bcftools index --nrecords "${pass_vcf}") -eq 0 ]; then
        bcftools view -Oz -o somatic.pass.vep.vcf.gz "${pass_vcf}"
        printf '%s\\n' '<!doctype html><html><head><title>Somvar VEP summary</title></head><body><h1>VEP annotation</h1><p>No PASS variants were available for annotation.</p></body></html>' > somatic.pass.vep.summary.html
    else
        vep --input_file "${pass_vcf}" \
            --output_file somatic.pass.vep.vcf.gz \
            --format vcf --vcf --compress_output bgzip \
            --cache --offline --species homo_sapiens --assembly GRCh38 \
            --dir_cache "${vep_cache}" --cache_version ${cache_version} \
            --fasta "${ref}" --fork ${annotation_forks} --buffer_size 500 \
            --symbol --canonical --mane --biotype --protein --hgvs \
            --sift b --polyphen b --af --af_gnomade --af_gnomadg --variant_class \
            --stats_file somatic.pass.vep.summary.html \
            --warning_file vep.warnings.txt --skipped_variants_file vep.skipped_variants.txt \
            --force_overwrite
    fi
    bcftools index -t somatic.pass.vep.vcf.gz
    """
}

process VARIANT_SUMMARY {
    cpus 1
    memory '4 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(pass_vcf), path(pass_index)
    tuple path(vep_vcf), path(vep_index)
    path vep_warnings
    path vep_skipped
    path contamination
    path reporting_script

    output:
    path 'summary', emit: summary

    script:
    """
    python3 "${reporting_script}" summary \
        --pass-vcf "${pass_vcf}" --vep-vcf "${vep_vcf}" \
        --warnings "${vep_warnings}" --skipped "${vep_skipped}" \
        --contamination "${contamination}" --outdir summary
    """
}

process MULTIQC {
    cpus 1
    memory '4 GB'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path raw_qc, stageAs: 'qc/raw/*'
    path clean_qc, stageAs: 'qc/clean/*'
    path fastp_json, stageAs: 'qc/fastp/*'
    path duplicate_metrics, stageAs: 'qc/duplicates/*'
    path vep_summary, stageAs: 'qc/vep/somatic_vep.html'
    path summary_dir
    path qc_config

    output:
    tuple path('multiqc_report.html'), path('multiqc_data'), path('multiqc.version.txt'), emit: qc

    script:
    """
    multiqc --version > multiqc.version.txt
    multiqc qc "${summary_dir}" --config "${qc_config}" \
        --module fastqc --module fastp --module picard --module vep --module custom_content \
        --force --filename multiqc_report.html --outdir .
    """
}

process FINAL_REPORT {
    cpus 1
    memory '4 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path summary_dir
    tuple path(multiqc_html), path(multiqc_data), path(multiqc_version)
    tuple path(pass_vcf), path(pass_index)
    tuple path(vep_vcf), path(vep_index)
    path vep_summary
    path vep_warnings
    path vep_skipped
    path fastp_json, stageAs: 'qc/fastp/*'
    path duplicate_metrics, stageAs: 'qc/duplicates/*'
    path reporting_script
    val manifest_base64

    output:
    path 'reports', emit: report

    script:
    """
    python3 -c 'import base64,sys; sys.stdout.write(base64.b64decode(sys.argv[1]).decode())' '${manifest_base64}' > run_manifest.json
    python3 "${reporting_script}" report \
        --summary-dir "${summary_dir}" \
        --multiqc-report "${multiqc_html}" --multiqc-data "${multiqc_data}" \
        --multiqc-version "${multiqc_version}" \
        --pass-vcf "${pass_vcf}" --pass-index "${pass_index}" \
        --vep-vcf "${vep_vcf}" --vep-index "${vep_index}" \
        --vep-summary "${vep_summary}" --warnings "${vep_warnings}" --skipped "${vep_skipped}" \
        --fastp-dir qc/fastp --duplicate-dir qc/duplicates \
        --manifest run_manifest.json --outdir reports
    """
}

workflow {
    // Nextflow CLI values may be strings, including "false".
    def skip_vep = params.skip_vep.toString().toBoolean()

    normal_reads_ch = channel.of(
        tuple('normal', [file(params.normalR1), file(params.normalR2)])
    )

    tumor_reads_ch = channel.of(
        tuple('tumor', [file(params.tumorR1), file(params.tumorR2)])
    )

    ref_ch = channel.value(file(params.ref, checkIfExists: true))
    bwa_index_ch = channel.value(
        ['.amb', '.ann', '.bwt', '.pac', '.sa'].collect { suffix ->
            file("${params.ref}${suffix}", checkIfExists: true)
        }
    )

    paired_reads_ch = normal_reads_ch.mix(tumor_reads_ch)

    // Normal and tumor may run together; each stage waits for both to finish.
    FASTQC(paired_reads_ch)
    raw_qc_done = FASTQC.out.reports.collect()

    FASTP(paired_reads_ch, raw_qc_done)
    fastp_done = FASTP.out.trimmed_reads.collect()

    FASTQC_CLEAN(FASTP.out.trimmed_reads, fastp_done)
    clean_qc_done = FASTQC_CLEAN.out.html.collect()

    BWA_MEM(FASTP.out.trimmed_reads, ref_ch, bwa_index_ch, raw_qc_done, clean_qc_done)

    alignments_for_bam = BWA_MEM.out.alignments
        .collect(flat: false)
        .flatMap { alignments -> alignments }

    BAM_PROCESSING(alignments_for_bam)

    // BQSR is optional and needs known variants matched to the reference.
    bams_for_calling = BAM_PROCESSING.out.indexed_bams
    if (params.known_sites) {
        known_site_names = (params.known_sites instanceof List
            ? params.known_sites : params.known_sites.toString().split(','))
            .collect { site -> site.toString().trim() }.findAll { site -> site }
        if (known_site_names.isEmpty()) {
            error 'Provide at least one VCF with --known_sites'
        }
        known_vcfs_ch = channel.value(known_site_names.collect { site -> file(site, checkIfExists: true) })
        known_indexes_ch = channel.value(known_site_names.collect { site ->
            def index_name = ["${site}.tbi", "${site}.idx"].find { index_path -> file(index_path).exists() }
            if (!index_name) {
                error "Missing .tbi or .idx index for known-sites VCF: ${site}"
            }
            file(index_name, checkIfExists: true)
        })
        ref_fai_ch = channel.value(file("${params.ref}.fai", checkIfExists: true))
        ref_dict_name = params.ref.toString().replaceFirst(/(?i)\.(fasta|fna|fa)$/, '') + '.dict'
        ref_dict_ch = channel.value(file(ref_dict_name, checkIfExists: true))

        indexed_bams_for_bqsr = BAM_PROCESSING.out.indexed_bams
            .collect(flat: false)
            .flatMap { bams -> bams }
        BQSR(indexed_bams_for_bqsr, ref_ch, ref_fai_ch, ref_dict_ch,
            known_vcfs_ch, known_indexes_ch)
        bams_for_calling = BQSR.out.recalibrated_bams
    }

    // Require reference-matched VCFs; capture intervals may be supplied when known.
    if (params.germline_resource || params.pon || params.intervals) {
        if (!(params.germline_resource && params.pon)) {
            error 'Mutect2 requires --germline_resource and --pon together'
        }
        ref_fai_ch = channel.value(file("${params.ref}.fai", checkIfExists: true))
        ref_dict_name = params.ref.toString().replaceFirst(/(?i)\.(fasta|fna|fa)$/, '') + '.dict'
        ref_dict_ch = channel.value(file(ref_dict_name, checkIfExists: true))
        germline_ch = channel.value(file(params.germline_resource, checkIfExists: true))
        pon_ch = channel.value(file(params.pon, checkIfExists: true))
        intervals_ch = channel.value(params.intervals
            ? [file(params.intervals, checkIfExists: true)] : [])
        germline_index_name = ["${params.germline_resource}.tbi", "${params.germline_resource}.idx"]
            .find { index_path -> file(index_path).exists() }
        pon_index_name = ["${params.pon}.tbi", "${params.pon}.idx"]
            .find { index_path -> file(index_path).exists() }
        if (!germline_index_name || !pon_index_name) {
            error 'Mutect2 germline resource and panel of normals each need a .tbi or .idx index'
        }
        germline_index_ch = channel.value(file(germline_index_name, checkIfExists: true))
        pon_index_ch = channel.value(file(pon_index_name, checkIfExists: true))
        common_sites_ch = channel.value(file(params.common_sites, checkIfExists: true))
        common_index_name = ["${params.common_sites}.tbi", "${params.common_sites}.idx"]
            .find { index_path -> file(index_path).exists() }
        if (!common_index_name) {
            error "Missing .tbi or .idx index for common-SNP VCF: ${params.common_sites}"
        }
        common_index_ch = channel.value(file(common_index_name, checkIfExists: true))
        if (!skip_vep) {
            cache_info = file("${params.vep_cache}/homo_sapiens/${params.vep_cache_version}_GRCh38/info.txt")
            if (!cache_info.exists()) {
                error "Missing GRCh38 VEP cache: ${cache_info}. See README.md (VEP cache) or use --skip_vep."
            }
        }

        paired_bams_ch = bams_for_calling.collect(flat: false).map { bams ->
            def normal = bams.find { bam -> bam[0] == 'normal' }
            def tumor = bams.find { bam -> bam[0] == 'tumor' }
            if (!normal || !tumor || bams.size() != 2) {
                error 'Mutect2 needs one normal BAM and one tumor BAM'
            }
            tuple(normal[1], normal[2], tumor[1], tumor[2])
        }
        MUTECT2(paired_bams_ch, ref_ch, ref_fai_ch, ref_dict_ch,
            germline_ch, germline_index_ch, pon_ch, pon_index_ch, intervals_ch)

        // Finish variant calling before estimating contamination in both samples.
        mutect2_done = MUTECT2.out.stats.collect()
        bams_for_pileups = bams_for_calling.collect(flat: false).flatMap { bams -> bams }
        GET_PILEUP_SUMMARIES(bams_for_pileups, ref_ch, ref_fai_ch, ref_dict_ch,
            common_sites_ch, common_index_ch, mutect2_done)

        paired_pileups_ch = GET_PILEUP_SUMMARIES.out.pileups.collect(flat: false).map { tables ->
            def normal = tables.find { pileup -> pileup[0] == 'normal' }
            def tumor = tables.find { pileup -> pileup[0] == 'tumor' }
            if (!normal || !tumor || tables.size() != 2) {
                error 'Contamination estimation needs one normal and one tumor pileup table'
            }
            tuple(normal[1], tumor[1])
        }
        CALCULATE_CONTAMINATION(paired_pileups_ch)
        FILTER_MUTECT_CALLS(MUTECT2.out.variants, MUTECT2.out.index.toList(), MUTECT2.out.stats,
            CALCULATE_CONTAMINATION.out.contamination, CALCULATE_CONTAMINATION.out.segments,
            ref_ch, ref_fai_ch, ref_dict_ch)
        PASS_VARIANTS(FILTER_MUTECT_CALLS.out.variants)

        if (!skip_vep) {
            vep_cache_ch = channel.value(file(params.vep_cache, checkIfExists: true))
            VEP(PASS_VARIANTS.out.variants, ref_ch, ref_fai_ch, vep_cache_ch,
                channel.value(params.vep_cache_version as int))

            reporting_script_ch = channel.value(file("${projectDir}/reporting.py", checkIfExists: true))
            qc_config_ch = channel.value(file(params.multiqc_config, checkIfExists: true))
            VARIANT_SUMMARY(PASS_VARIANTS.out.variants, VEP.out.variants,
                VEP.out.warnings, VEP.out.skipped, CALCULATE_CONTAMINATION.out.contamination,
                reporting_script_ch)

            // Use emitted files from this run rather than scanning published results.
            raw_qc_reports = FASTQC.out.reports.flatten().filter { report -> report.name.endsWith('.zip') }.toList()
            clean_qc_reports = FASTQC_CLEAN.out[1].flatten().toList()
            fastp_reports = FASTP.out[2].toList()
            duplicate_reports = BAM_PROCESSING.out.duplicate_metrics.map { sample, metrics -> metrics }.toList()
            MULTIQC(raw_qc_reports, clean_qc_reports, fastp_reports, duplicate_reports,
                VEP.out.summary, VARIANT_SUMMARY.out.summary, qc_config_ch)

            run_manifest = [
                run_name: workflow.runName,
                session_id: workflow.sessionId.toString(),
                nextflow_version: workflow.nextflow.version.toString(),
                command: workflow.commandLine,
                reference: params.ref.toString(),
                normal_reads: [params.normalR1, params.normalR2],
                tumor_reads: [params.tumorR1, params.tumorR2],
                threads_per_task: params.threads as int,
                known_sites: params.known_sites,
                germline_resource: params.germline_resource,
                panel_of_normals: params.pon,
                common_sites: params.common_sites,
                intervals: params.intervals,
                vep_cache: params.vep_cache.toString(),
                vep_cache_version: params.vep_cache_version as int,
                outdir: params.outdir.toString(),
                main_nf_sha256: java.security.MessageDigest.getInstance('SHA-256')
                    .digest(file("${projectDir}/main.nf").bytes).encodeHex().toString()
            ]
            manifest_base64 = groovy.json.JsonOutput.toJson(run_manifest).bytes.encodeBase64().toString()
            FINAL_REPORT(VARIANT_SUMMARY.out.summary, MULTIQC.out.qc,
                PASS_VARIANTS.out.variants, VEP.out.variants, VEP.out.summary,
                VEP.out.warnings, VEP.out.skipped, fastp_reports, duplicate_reports,
                reporting_script_ch, channel.value(manifest_base64))
        }
    }
}
