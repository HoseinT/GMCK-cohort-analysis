# Step 5: Annotate filtered VCFs with VEP v108 and CADD pre-scored files.

REF = config["reference"]["genome"]


rule vep_annotate:
    """Annotate with Ensembl VEP v108 (cache, GRCh37)."""
    input:
        vcf = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_filtered.vcf.gz",
    output:
        vcf = f"{RESULTS}/{{patient}}/annotated/{{tumor}}_vep.vcf.gz",
        tbi = f"{RESULTS}/{{patient}}/annotated/{{tumor}}_vep.vcf.gz.tbi",
        summary = f"{RESULTS}/{{patient}}/annotated/{{tumor}}_vep_summary.html",
    params:
        vep       = config["tools"]["vep"],
        cache_dir = config["vep"]["cache_dir"],
        assembly  = config["vep"]["assembly"],
        version   = config["vep"]["version"],
        extra     = config["vep"]["extra_flags"],
        ref       = REF,
    threads: config["threads"]["general"]
    log: f"{RESULTS}/{{patient}}/logs/vep_{{tumor}}.log"
    shell:
        """
        {params.vep} \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --stats_file {output.summary} \
            --vcf --compress_output bgzip \
            --fork {threads} \
            --cache --dir_cache {params.cache_dir} \
            --cache_version {params.version} \
            --assembly {params.assembly} \
            --fasta {params.ref} \
            --offline \
            {params.extra} \
            --force_overwrite \
        2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule vcf_to_table:
    """Convert VEP-annotated VCF to a flat TSV for downstream filtering."""
    input:
        vcf = f"{RESULTS}/{{patient}}/annotated/{{tumor}}_vep.vcf.gz",
    output:
        tsv = temp(f"{RESULTS}/{{patient}}/annotated/{{tumor}}_vep.tsv"),
    params:
        gatk   = config["tools"]["gatk"],
        ref    = REF,
    log: f"{RESULTS}/{{patient}}/logs/vcf_to_table_{{tumor}}.log"
    shell:
        """
        {params.gatk} VariantsToTable \
            -R {params.ref} \
            -V {input.vcf} \
            -F CHROM -F POS -F REF -F ALT -F FILTER \
            -F CSQ \
            -GF AD -GF AF -GF DP \
            --show-filtered \
            -O {output.tsv} \
        2> {log}
        """


rule add_cadd:
    """Look up pre-scored CADD scores (tabix) and append to TSV."""
    input:
        tsv = f"{RESULTS}/{{patient}}/annotated/{{tumor}}_vep.tsv",
    output:
        tsv = f"{RESULTS}/{{patient}}/annotated/{{tumor}}_cadd.tsv",
    params:
        script    = "../../scripts/add_cadd_scores.py",
        snv_file  = config["cadd"]["snv_file"],
        indel_file = config["cadd"]["indel_file"],
    log: f"{RESULTS}/{{patient}}/logs/cadd_{{tumor}}.log"
    shell:
        """
        python {params.script} \
            --input {input.tsv} \
            --snv-file {params.snv_file} \
            --indel-file {params.indel_file} \
            --output {output.tsv} \
        2> {log}
        """
