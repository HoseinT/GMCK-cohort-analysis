# Cohort-level analyses — run once all per-patient outputs are available.
# Rules depend on the full expand() of per-patient outputs so Snakemake
# schedules them only after every sample has been processed.

COHORT = f"{RESULTS}/cohort"
RSCRIPT = config["tools"]["Rscript"]
GENE_BED = config["cohort"]["gene_bed"]


rule cohort_tmb:
    """Calculate TMB (SNVs / panel Mb) for every primary and relapse sample."""
    input:
        variants = expand(
            f"{RESULTS}/{{patient}}/final_variants/{{patient}}_final.tsv",
            patient=PATIENTS,
        ),
    output:
        tsv = f"{COHORT}/tmb_summary.tsv",
    params:
        script       = "../../scripts/tmb.py",
        panel_size   = config["cohort"]["panel_size_mb"],
        results_dir  = RESULTS,
        patient_list = " ".join(PATIENTS),
    log: f"{COHORT}/logs/tmb.log"
    shell:
        """
        python {params.script} \
            --results-dir {params.results_dir} \
            --patients {params.patient_list} \
            --panel-size-mb {params.panel_size} \
            --output {output.tsv} \
        2> {log}
        """


rule discover_test:
    """Mutual exclusivity / co-occurrence test with the DISCOVER Python package."""
    input:
        variants = expand(
            f"{RESULTS}/{{patient}}/final_variants/{{patient}}_final.tsv",
            patient=PATIENTS,
        ),
    output:
        tsv = f"{COHORT}/mutual_exclusivity.tsv",
    params:
        script       = "../../scripts/mutual_exclusivity.py",
        results_dir  = RESULTS,
        patient_list = " ".join(PATIENTS),
        min_samples  = config["cohort"]["discover_min_samples"],
    log: f"{COHORT}/logs/discover.log"
    shell:
        """
        python {params.script} \
            --results-dir {params.results_dir} \
            --patients {params.patient_list} \
            --min-samples {params.min_samples} \
            --output {output.tsv} \
        2> {log}
        """


rule cohort_cnv:
    """Ploidy-adjusted CNV analysis: burden, recurrence, paired Wilcoxon+BH."""
    input:
        cnv = expand(
            f"{RESULTS}/{{patient}}/facets/{{tumor}}_cnv.tsv",
            patient=PATIENTS, tumor=TUMOR_TYPES,
        ),
        purity = expand(
            f"{RESULTS}/{{patient}}/facets/{{tumor}}_purity.txt",
            patient=PATIENTS, tumor=TUMOR_TYPES,
        ),
    output:
        burden    = f"{COHORT}/cnv_burden.tsv",
        recurrent = f"{COHORT}/cnv_recurrent.tsv",
        wilcoxon  = f"{COHORT}/cnv_wilcoxon.tsv",
        matrix    = f"{COHORT}/cnv_adjusted_matrix.tsv",
    params:
        script            = "../../scripts/cohort_cnv_analysis.py",
        results_dir       = RESULTS,
        patient_list      = " ".join(PATIENTS),
        gene_bed          = GENE_BED,
        amp_per_sample    = config["cohort"]["cnv_amp_per_sample"],
        del_per_sample    = config["cohort"]["cnv_del_per_sample"],
        avg_amp           = config["cohort"]["cnv_avg_amp"],
        avg_del           = config["cohort"]["cnv_avg_del"],
        burden_abs        = config["cohort"]["cnv_burden_abs"],
        wilcoxon_min_pat  = config["cohort"]["wilcoxon_min_patients"],
        outdir            = COHORT,
    log: f"{COHORT}/logs/cohort_cnv.log"
    shell:
        """
        python {params.script} \
            --results-dir {params.results_dir} \
            --patients {params.patient_list} \
            --gene-bed {params.gene_bed} \
            --amp-threshold {params.amp_per_sample} \
            --del-threshold {params.del_per_sample} \
            --avg-amp-threshold {params.avg_amp} \
            --avg-del-threshold {params.avg_del} \
            --burden-abs {params.burden_abs} \
            --wilcoxon-min-patients {params.wilcoxon_min_pat} \
            --outdir {params.outdir} \
        2> {log}
        """


rule cohort_plots:
    """Generate all cohort-level Plotly HTML reports."""
    input:
        tmb       = f"{COHORT}/tmb_summary.tsv",
        burden    = f"{COHORT}/cnv_burden.tsv",
        recurrent = f"{COHORT}/cnv_recurrent.tsv",
        wilcoxon  = f"{COHORT}/cnv_wilcoxon.tsv",
        matrix    = f"{COHORT}/cnv_adjusted_matrix.tsv",
        mutex     = f"{COHORT}/mutual_exclusivity.tsv",
    output:
        expand(
            f"{COHORT}/plots/{{plot}}.html",
            plot=["tmb", "cnv_burden", "cnv_heatmap",
                  "mutual_exclusivity", "cnv_wilcoxon_volcano"],
        ),
    params:
        script = "../../scripts/cohort_plots.py",
        outdir = f"{COHORT}/plots",
    log: f"{COHORT}/logs/cohort_plots.log"
    shell:
        """
        python {params.script} \
            --tmb        {input.tmb} \
            --burden     {input.burden} \
            --recurrent  {input.recurrent} \
            --wilcoxon   {input.wilcoxon} \
            --matrix     {input.matrix} \
            --mutex      {input.mutex} \
            --outdir     {params.outdir} \
        2> {log}
        """
