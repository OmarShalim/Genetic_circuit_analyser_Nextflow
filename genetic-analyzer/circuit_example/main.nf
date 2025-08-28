#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// -------------------- Params --------------------
params.settings      = "./data/settings.txt"
params.bin_path      = "../bin"
params.feature       = "gene"
params.attribute     = "Name"
params.strand_opt    = "reverse"
params.chrom_list    = "0x58v50"
params.only_samples  = null

// -------------------- Processes --------------------

process Setup {
    input:
    val sample

    output:
    path "logs"
    path "results/${sample}"
    path "tmp/${sample}"

    script:
    """
    mkdir -p logs results tmp results/${sample} tmp/${sample}
    """
}

process MapReads {
    input:
    tuple val(sample), path(settings_file), path(bin_dir)

    output:
    path "tmp/${sample}/*"

    publishDir "./", mode: 'copy'

    script:
    """
    python3 ${bin_dir}/map_reads.py \\
        -settings ${settings_file} \\
        -samples ${sample}
    """
}

process CountReads {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_file), path(count_script), path(settings_file)

    output:
    // Emit per-sample tuple containing all files the cohort step will need
    tuple val(sample_id),
          path("${sample_id}.counts.txt"),
          path("${sample_id}.mapped.reads.txt"),
          path("${sample_id}.gene.lengths.txt")

    publishDir "results/${sample_id}", mode: 'copy'

    script:
    """
    python3 ${count_script} \\
        -settings ${settings_file} \\
        -samples ${sample_id} \\
        -feature ${params.feature} \\
        -attribute ${params.attribute} \\
        -strand_opt ${params.strand_opt}
    """
}

process FragmentDist {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fragdist_script), path(settings_file)

    output:
    tuple val(sample_id), path("${sample_id}.fragment.distribution.txt")

    publishDir "results/${sample_id}", mode: 'copy'

    script:
    """
    python3 ${fragdist_script} \\
        -settings ${settings_file} \\
        -samples ${sample_id}
    """
}

/*
 * Run ONCE for the whole cohort after CountReads finishes.
 */
process ReadAnalysis {
    tag "AllSamples"

    input:
    path counts_files
    path mapped_files
    path lengths_files
    path settings_file
    path bin_path
    path norm_r_script

    output:
    path("counts.matrix.txt"),          emit: counts_matrix
    path("mapped.reads.matrix.txt"),    emit: mapped_matrix
    path("gene.lengths.matrix.txt"),    emit: lengths_matrix
    path("fpkm.normed.matrix.txt"),     emit: fpkm_normed
    path("norm.factors.matrix.txt"),    emit: norm_factors

    publishDir "results", mode: 'copy'

    script:
    """
    ls -l
    python3 ${bin_path}/read_analysis.py \\
        -settings ${settings_file} \\
        -bin_path ${bin_path}
    """
}

/*
 * Fan-out the single norm factors file so each sample gets a copy
 * under results/<sample>/norm.factors.matrix.txt.
 * NOTE: Nextflow stages the input with the same basename; avoid copying onto itself.
 */
process DistributeNormFactors {
    tag "$sample_id"

    input:
    val sample_id
    path norm_factors

    output:
    tuple val(sample_id), path("norm.factors.matrix.txt")

    publishDir "results/${sample_id}", mode: 'copy', overwrite: true

    script:
    """
    # If the staged file already exists with this name, do nothing; otherwise copy it.
    if [ ! -e norm.factors.matrix.txt ]; then
        cp "${norm_factors}" norm.factors.matrix.txt
    fi
    """
}

process TranscriptionProfile {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(profile_script), path(settings_file), val(chroms)
    tuple val(sample_id), path(norm_factors_file)
    tuple val(sample_id), path(frag_dist_file)

    output:
    path("${sample_id}.rev.profiles.txt")
    path("${sample_id}.fwd.profiles.txt")
    path("${sample_id}.rev.norm.profiles.txt")
    path("${sample_id}.fwd.norm.profiles.txt")

    publishDir "results/${sample_id}", mode: 'copy', overwrite: true

    script:
    """
    python3 ${profile_script} \\
        -settings ${settings_file} \\
        -samples ${sample_id} \\
        -chroms ${chroms}
    """
}

// -------------------- Workflow --------------------

workflow {
    def count_script     = file("${params.bin_path}/count_reads.py")
    def fragdist_script  = file("${params.bin_path}/fragment_distributions.py")
    def profile_script   = file("${params.bin_path}/transcription_profile.py")
    def norm_r_script    = file("${params.bin_path}/binnorm_fpkm.r")
    def settings_file    = file(params.settings)
    def bin_path         = file(params.bin_path)
    def chroms           = params.chrom_list

    // Parse sample names from settings; optionally restrict via --only_samples
    all_samples_ch = Channel.fromPath(params.settings)
        .map { f ->
            def lines = f.text.readLines()
            def names = lines.findAll { ln ->
                ln && !ln.startsWith("sample") && !ln.startsWith("None")
            }.collect { it.split('\\t')[0].trim() }
            if (params.only_samples) {
                def wanted = (params.only_samples as String).split(',')*.trim() as Set
                names = names.findAll { it in wanted }
            }
            names.unique()
        }
        .flatten()

    sample_info = all_samples_ch.map { s -> tuple(s, file(params.settings), file(params.bin_path)) }

    // Folders
    all_samples_ch | Setup

    // Map reads
    mapreads_ch = sample_info | MapReads

    // Count reads
    countreads_ch = mapreads_ch | flatMap { files ->
        files.findAll { it.name.endsWith('.bam') }.collect { bam ->
            def sample_id = bam.getBaseName().tokenize('.')[0]
            tuple(sample_id, bam, count_script, settings_file)
        }
    } | CountReads

    // Fragment distributions per sample
    fragdist_ch = countreads_ch \
        | map { tuple(it[0], fragdist_script, settings_file) } \
        | FragmentDist

    // ----- Cohort aggregation for ReadAnalysis -----
    counts_list_ch  = countreads_ch.map { it[1] }.collect()
    mapped_list_ch  = countreads_ch.map { it[2] }.collect()
    lengths_list_ch = countreads_ch.map { it[3] }.collect()

    read_analysis_out = ReadAnalysis(counts_list_ch, mapped_list_ch, lengths_list_ch,
                                     settings_file, bin_path, norm_r_script)

    // Fan-out norm factors to each sample
    sample_ids_ch = countreads_ch.map { it[0] }.unique()
    norm_factors_ch = DistributeNormFactors(sample_ids_ch, read_analysis_out.norm_factors)

    // Inputs for TranscriptionProfile
    profile_input_ch = countreads_ch.map {
        tuple(it[0], profile_script, settings_file, chroms)
    }

    // Wire up profile: (per sample) × (norm factors per sample) × (frag dist per sample)
    TranscriptionProfile(profile_input_ch, norm_factors_ch, fragdist_ch)
}
