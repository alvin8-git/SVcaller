// params.delly_inv_max_time -> whole seconds; 0 = no cap. Accepts a bare number or digit string
// (seconds, so 7200 == '7200s'), or a Nextflow duration ('2h', '90m', '7200s'). 0, '0', '', null and
// false mean no cap. Anything that parses to under 1 s (but not 0), or is negative or unparsable, is an error.
// Called from main.nf at workflow start (fail fast) and by sv_calling.nf, which passes it to the INV task only.
def delly_inv_cap_seconds(v) {
    if (v == null || v == false) return 0L
    def str = v.toString().trim()
    if (str == '' || str.toLowerCase() == 'false') return 0L
    long secs
    try {
        secs = (str ==~ /\d+/) ? str.toLong() : nextflow.util.Duration.of(str).toSeconds()
    } catch (Exception e) {
        error "ERROR: --delly_inv_max_time '${str}' is not a duration (use e.g. 2h, 90m, 7200s, a number of seconds, or 0 for no cap)"
    }
    if (secs == 0 && str ==~ /\d+/) return 0L
    if (secs < 1) error "ERROR: --delly_inv_max_time '${str}' is under 1 second; use 0 for no cap"
    return secs
}

// DELLY_CALL_SVTYPE runs one SV type per process invocation.
// sv_calling.nf fans out 5 types in parallel then collects into DELLY_MERGE.
process DELLY_CALL_SVTYPE {
    tag "${meta.id}:${svtype}"
    label 'process_medium'
    // INV only writes this record: header only = delly INV completed, one row = it hit
    // params.delly_inv_max_time and the sample continued WITHOUT delly inversion calls.
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.sv_not_assessed.tsv"

    input:
    tuple val(meta), path(bam), path(bai), val(svtype), val(cap_seconds)   // cap_seconds: 0 for every type but INV
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.delly.${svtype}.vcf"), emit: vcf
    tuple val(meta), path("${meta.id}.sv_not_assessed.tsv"), optional: true, emit: not_assessed
    path "versions.yml",                                       emit: versions

    script:
    // The cap arrives as an input (not read from params here) so only the INV task's hash depends on it.
    """
    delly_call_capped.sh ${meta.id} ${svtype} ${cap_seconds} ${fasta} ${fai} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | grep "DELLY" | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
