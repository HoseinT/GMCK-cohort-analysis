# Step 6: Apply post-annotation filters and merge primary + relapse variant tables.
# The Python script handles: PASS-only, depth >= 50, VAF > 5% (union),
# consequence exclusion, CADD >= 0.20, and VAF normalization by FACETS purity.


rule filter_and_merge_variants:
    """Merge primary/relapse CADD-annotated TSVs, apply all filters."""
    input:
        primary = f"{RESULTS}/{{patient}}/annotated/primary_cadd.tsv",
        relapse = f"{RESULTS}/{{patient}}/annotated/relapse_cadd.tsv",
        primary_purity = f"{RESULTS}/{{patient}}/facets/primary_purity.txt",
        relapse_purity = f"{RESULTS}/{{patient}}/facets/relapse_purity.txt",
    output:
        tsv = f"{RESULTS}/{{patient}}/final_variants/{{patient}}_final.tsv",
    params:
        script       = "../../scripts/filter_variants.py",
        min_depth    = config["params"]["min_depth"],
        min_vaf      = config["params"]["min_vaf"],
        cadd_cutoff  = config["params"]["cadd_cutoff"],
        excluded_csq = " ".join(config["excluded_consequences"]),
        patient      = lambda wc: wc.patient,
    log: f"{RESULTS}/{{patient}}/logs/filter_variants.log"
    shell:
        """
        python {params.script} \
            --primary {input.primary} \
            --relapse {input.relapse} \
            --primary-purity {input.primary_purity} \
            --relapse-purity {input.relapse_purity} \
            --patient-id {params.patient} \
            --min-depth {params.min_depth} \
            --min-vaf {params.min_vaf} \
            --cadd-cutoff {params.cadd_cutoff} \
            --exclude-consequences {params.excluded_csq} \
            --output {output.tsv} \
        2> {log}
        """
