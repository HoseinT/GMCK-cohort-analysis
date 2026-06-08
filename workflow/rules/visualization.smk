# Step 8: Generate Plotly interactive HTML reports for SNVs and CNVs.


rule plot_snv:
    """Lollipop plot comparing primary vs relapse VAF per patient."""
    input:
        tsv = f"{RESULTS}/{{patient}}/final_variants/{{patient}}_final.tsv",
    output:
        html = f"{RESULTS}/{{patient}}/plots/{{patient}}_snv.html",
    params:
        script  = "../../scripts/visualize.py",
        patient = lambda wc: wc.patient,
    log: f"{RESULTS}/{{patient}}/logs/plot_snv_{{patient}}.log"
    shell:
        """
        python {params.script} snv \
            --variants {input.tsv} \
            --patient-id {params.patient} \
            --output {output.html} \
        2> {log}
        """


rule plot_cnv:
    """Genome-wide CNV plot (primary + relapse) for one patient."""
    input:
        primary_cnv = f"{RESULTS}/{{patient}}/facets/primary_cnv.tsv",
        relapse_cnv = f"{RESULTS}/{{patient}}/facets/relapse_cnv.tsv",
    output:
        html = f"{RESULTS}/{{patient}}/plots/{{patient}}_cnv.html",
    params:
        script  = "../../scripts/visualize.py",
        patient = lambda wc: wc.patient,
    log: f"{RESULTS}/{{patient}}/logs/plot_cnv_{{patient}}.log"
    shell:
        """
        python {params.script} cnv \
            --primary-cnv {input.primary_cnv} \
            --relapse-cnv {input.relapse_cnv} \
            --patient-id {params.patient} \
            --output {output.html} \
        2> {log}
        """
