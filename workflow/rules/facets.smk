# Step 7: Copy number analysis with FACETS v0.6.2.
# snp-pileup generates allele counts at common SNP positions;
# the R script runs FACETS and writes CNV segments + purity/ploidy.

SNP_VCF   = config["reference"]["snp_pileup_vcf"]
MAX_DEPTH = config["params"]["facets_max_depth"]


rule snp_pileup:
    """Generate snp-pileup counts (tumor + germline) for FACETS."""
    input:
        tumor_bam    = f"{RESULTS}/{{patient}}/bqsr/{{tumor}}_bqsr.bam",
        germline_bam = f"{RESULTS}/{{patient}}/bqsr/germline_bqsr.bam",
    output:
        pileup = f"{RESULTS}/{{patient}}/facets/{{tumor}}_pileup.csv.gz",
    params:
        snp_pileup = config["tools"]["snp_pileup"],
        snp_vcf    = SNP_VCF,
        max_depth  = MAX_DEPTH,
    log: f"{RESULTS}/{{patient}}/logs/snp_pileup_{{tumor}}.log"
    shell:
        """
        {params.snp_pileup} \
            --count-orphans \
            --max-depth {params.max_depth} \
            --min-map-quality 15 \
            --min-base-quality 20 \
            --gzip \
            {params.snp_vcf} \
            {output.pileup} \
            {input.germline_bam} \
            {input.tumor_bam} \
        2> {log}
        """


rule run_facets:
    """Run FACETS CNV analysis and emit segment TSV + purity/ploidy."""
    input:
        pileup = f"{RESULTS}/{{patient}}/facets/{{tumor}}_pileup.csv.gz",
    output:
        cnv_tsv = f"{RESULTS}/{{patient}}/facets/{{tumor}}_cnv.tsv",
        purity  = f"{RESULTS}/{{patient}}/facets/{{tumor}}_purity.txt",
        qc      = f"{RESULTS}/{{patient}}/facets/{{tumor}}_facets_qc.txt",
    params:
        script  = "../../scripts/run_facets.R",
        Rscript = config["tools"]["Rscript"],
        outdir  = lambda wc: f"{RESULTS}/{wc.patient}/facets",
        sample  = lambda wc: f"{wc.patient}_{wc.tumor}",
    log: f"{RESULTS}/{{patient}}/logs/facets_{{tumor}}.log"
    shell:
        """
        {params.Rscript} {params.script} \
            --pileup {input.pileup} \
            --outdir {params.outdir} \
            --sample {params.sample} \
        2> {log}
        """
