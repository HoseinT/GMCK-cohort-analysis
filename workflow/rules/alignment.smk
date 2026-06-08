# Step 2: Align trimmed reads with BWA mem, then merge with the unaligned BAM
# so that the RX (UMI) tag is restored in the final coordinate-sorted BAM.

REF = config["reference"]["genome"]


rule bwa_mem:
    """Align trimmed FASTQs to GRCh37 with BWA mem."""
    input:
        r1 = f"{RESULTS}/{{patient}}/trimmed/{{stype}}_R1.fastq.gz",
        r2 = f"{RESULTS}/{{patient}}/trimmed/{{stype}}_R2.fastq.gz",
    output:
        bam = temp(f"{RESULTS}/{{patient}}/aligned/{{stype}}_aligned.bam"),
    params:
        bwa      = config["tools"]["bwa"],
        samtools = config["tools"]["samtools"],
        ref      = REF,
        rg       = lambda wc: (
            f"@RG\\tID:{wc.patient}_{wc.stype}"
            f"\\tSM:{wc.patient}_{wc.stype}"
            f"\\tPL:ILLUMINA"
            f"\\tLB:{wc.patient}_{wc.stype}_lib"
            f"\\tPU:{wc.patient}"
        ),
    threads: config["threads"]["bwa"]
    log: f"{RESULTS}/{{patient}}/logs/bwa_mem_{{stype}}.log"
    shell:
        """
        {params.bwa} mem \
            -R '{params.rg}' \
            -t {threads} \
            {params.ref} \
            {input.r1} {input.r2} \
        2> {log} \
        | {params.samtools} sort \
            -@ {threads} \
            -o {output.bam} \
        2>> {log}
        {params.samtools} index {output.bam}
        """


rule merge_bam_alignment:
    """Merge aligned BAM with unaligned BAM to restore UMI (RX) tags."""
    input:
        aligned = f"{RESULTS}/{{patient}}/aligned/{{stype}}_aligned.bam",
        ubam    = f"{RESULTS}/{{patient}}/umi/{{stype}}_unaligned.bam",
    output:
        merged = temp(f"{RESULTS}/{{patient}}/merged/{{stype}}_merged.bam"),
    params:
        gatk = config["tools"]["gatk"],
        ref  = REF,
    threads: config["threads"]["gatk"]
    log: f"{RESULTS}/{{patient}}/logs/merge_bam_alignment_{{stype}}.log"
    shell:
        """
        {params.gatk} MergeBamAlignment \
            --REFERENCE_SEQUENCE {params.ref} \
            --UNMAPPED_BAM {input.ubam} \
            --ALIGNED_BAM {input.aligned} \
            --OUTPUT {output.merged} \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT \
        2> {log}
        """
