# Step 3: UMI-aware duplicate marking then base quality score recalibration (BQSR).

REF     = config["reference"]["genome"]
DBSNP   = config["reference"]["dbsnp"]
MILLS   = config["reference"]["mills"]
KG1     = config["reference"]["known_1000g"]
INTERVALS = config["reference"]["intervals"]


rule umi_mark_duplicates:
    """UMI-aware duplicate marking using GATK UmiAwareMarkDuplicatesWithMateCigar."""
    input:
        bam = f"{RESULTS}/{{patient}}/merged/{{stype}}_merged.bam",
    output:
        bam     = temp(f"{RESULTS}/{{patient}}/markdup/{{stype}}_markdup.bam"),
        metrics = f"{RESULTS}/{{patient}}/markdup/{{stype}}_markdup_metrics.txt",
    params:
        gatk = config["tools"]["gatk"],
    threads: config["threads"]["gatk"]
    log: f"{RESULTS}/{{patient}}/logs/markdup_{{stype}}.log"
    shell:
        """
        {params.gatk} UmiAwareMarkDuplicatesWithMateCigar \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --UMI_TAG_NAME RX \
            --DUPLEX_UMI false \
            --ASSUME_SORT_ORDER coordinate \
        2> {log}
        """


rule base_recalibrator:
    """Build BQSR recalibration table."""
    input:
        bam = f"{RESULTS}/{{patient}}/markdup/{{stype}}_markdup.bam",
    output:
        table = f"{RESULTS}/{{patient}}/bqsr/{{stype}}_recal.table",
    params:
        gatk      = config["tools"]["gatk"],
        ref       = REF,
        dbsnp     = DBSNP,
        mills     = MILLS,
        kg1       = KG1,
        intervals = INTERVALS,
    threads: config["threads"]["gatk"]
    log: f"{RESULTS}/{{patient}}/logs/bqsr_recal_{{stype}}.log"
    shell:
        """
        {params.gatk} BaseRecalibrator \
            -I {input.bam} \
            -R {params.ref} \
            --known-sites {params.dbsnp} \
            --known-sites {params.mills} \
            --known-sites {params.kg1} \
            -L {params.intervals} \
            -O {output.table} \
        2> {log}
        """


rule apply_bqsr:
    """Apply BQSR recalibration to produce the analysis-ready BAM."""
    input:
        bam   = f"{RESULTS}/{{patient}}/markdup/{{stype}}_markdup.bam",
        table = f"{RESULTS}/{{patient}}/bqsr/{{stype}}_recal.table",
    output:
        bam = f"{RESULTS}/{{patient}}/bqsr/{{stype}}_bqsr.bam",
        bai = f"{RESULTS}/{{patient}}/bqsr/{{stype}}_bqsr.bam.bai",
    params:
        gatk      = config["tools"]["gatk"],
        ref       = REF,
        intervals = INTERVALS,
    threads: config["threads"]["gatk"]
    log: f"{RESULTS}/{{patient}}/logs/apply_bqsr_{{stype}}.log"
    shell:
        """
        {params.gatk} ApplyBQSR \
            -I {input.bam} \
            -R {params.ref} \
            --bqsr-recal-file {input.table} \
            -L {params.intervals} \
            -O {output.bam} \
        2> {log}
        """
