# Step 4: Somatic variant calling with Mutect2 (tumor-vs-normal), contamination
# estimation, and hard filtering. Runs independently for primary and relapse.

REF       = config["reference"]["genome"]
GNOMAD    = config["reference"]["gnomad"]
PON       = config["reference"]["pon"]
INTERVALS = config["reference"]["intervals"]
DBSNP     = config["reference"]["dbsnp"]


rule mutect2:
    """Call somatic SNVs and indels (tumor vs germline) with Mutect2."""
    input:
        tumor_bam    = f"{RESULTS}/{{patient}}/bqsr/{{tumor}}_bqsr.bam",
        germline_bam = f"{RESULTS}/{{patient}}/bqsr/germline_bqsr.bam",
    output:
        vcf   = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_raw.vcf.gz",
        stats = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_raw.vcf.gz.stats",
        f1r2  = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_f1r2.tar.gz",
    params:
        gatk        = config["tools"]["gatk"],
        ref         = REF,
        gnomad      = GNOMAD,
        pon         = PON,
        intervals   = INTERVALS,
        tumor_name  = lambda wc: f"{wc.patient}_{wc.tumor}",
        normal_name = lambda wc: f"{wc.patient}_germline",
    threads: config["threads"]["gatk"]
    log: f"{RESULTS}/{{patient}}/logs/mutect2_{{tumor}}.log"
    shell:
        """
        {params.gatk} Mutect2 \
            -R {params.ref} \
            -I {input.tumor_bam} \
            -I {input.germline_bam} \
            -tumor {params.tumor_name} \
            -normal {params.normal_name} \
            --germline-resource {params.gnomad} \
            --panel-of-normals {params.pon} \
            -L {params.intervals} \
            --f1r2-tar-gz {output.f1r2} \
            -O {output.vcf} \
        2> {log}
        """


rule learn_read_orientation_model:
    """Estimate read orientation artifacts (FFPE, oxidative damage)."""
    input:
        f1r2 = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_f1r2.tar.gz",
    output:
        model = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_read_orientation.tar.gz",
    params:
        gatk = config["tools"]["gatk"],
    log: f"{RESULTS}/{{patient}}/logs/orientation_model_{{tumor}}.log"
    shell:
        """
        {params.gatk} LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.model} \
        2> {log}
        """


rule get_pileup_summaries:
    """Pileup at common germline sites for contamination estimation."""
    input:
        bam = f"{RESULTS}/{{patient}}/bqsr/{{tumor}}_bqsr.bam",
    output:
        pileup = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_pileup.table",
    params:
        gatk   = config["tools"]["gatk"],
        gnomad = GNOMAD,
        intervals = INTERVALS,
    threads: config["threads"]["gatk"]
    log: f"{RESULTS}/{{patient}}/logs/pileup_{{tumor}}.log"
    shell:
        """
        {params.gatk} GetPileupSummaries \
            -I {input.bam} \
            -V {params.gnomad} \
            -L {params.intervals} \
            -O {output.pileup} \
        2> {log}
        """


rule calculate_contamination:
    """Estimate cross-sample contamination fraction."""
    input:
        pileup = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_pileup.table",
    output:
        contamination = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_contamination.table",
        segments      = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_segments.table",
    params:
        gatk = config["tools"]["gatk"],
    log: f"{RESULTS}/{{patient}}/logs/contamination_{{tumor}}.log"
    shell:
        """
        {params.gatk} CalculateContamination \
            -I {input.pileup} \
            --tumor-segmentation {output.segments} \
            -O {output.contamination} \
        2> {log}
        """


rule filter_mutect_calls:
    """Apply Mutect2 artifact filters and contamination estimates."""
    input:
        vcf           = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_raw.vcf.gz",
        stats         = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_raw.vcf.gz.stats",
        orientation   = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_read_orientation.tar.gz",
        contamination = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_contamination.table",
        segments      = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_segments.table",
    output:
        vcf = f"{RESULTS}/{{patient}}/mutect2/{{tumor}}_filtered.vcf.gz",
    params:
        gatk = config["tools"]["gatk"],
        ref  = REF,
    log: f"{RESULTS}/{{patient}}/logs/filter_mutect_{{tumor}}.log"
    shell:
        """
        {params.gatk} FilterMutectCalls \
            -R {params.ref} \
            -V {input.vcf} \
            --tumor-segmentation {input.segments} \
            --contamination-table {input.contamination} \
            --ob-priors {input.orientation} \
            --min-allele-fraction 0.0 \
            -O {output.vcf} \
        2> {log}
        """
