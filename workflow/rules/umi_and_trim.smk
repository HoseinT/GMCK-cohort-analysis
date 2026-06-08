# Step 1: Extract inline UMIs with fgbio, then trim adapters with Cutadapt.
# fgbio FastqToBam writes the UMI bases into the RX BAM tag so they survive
# through alignment and merge. Cutadapt then removes adapter sequences.

UMI_READ   = config["umi"]["read"]        # "R1" or "R2"
UMI_LEN    = config["umi"]["length"]
ADAPTER_R1 = config["adapters"]["r1"]
ADAPTER_R2 = config["adapters"]["r2"]
MIN_LEN    = config["params"]["cutadapt_min_length"]
QUAL       = config["params"]["cutadapt_quality"]

# fgbio ReadStructure strings: e.g. for 8-base UMI on R1: "8M+T" on R1, "+T" on R2
def _read_structure(read):
    if read == UMI_READ:
        return f"{UMI_LEN}M+T"
    return "+T"

READ_STRUCTURE_R1 = _read_structure("R1")
READ_STRUCTURE_R2 = _read_structure("R2")


rule fgbio_fastq_to_bam:
    """Convert FASTQ pair to unaligned BAM, extracting UMI into RX tag."""
    input:
        r1 = lambda wc: get_fastq(wc.patient, wc.stype, "r1"),
        r2 = lambda wc: get_fastq(wc.patient, wc.stype, "r2"),
    output:
        ubam = temp(f"{RESULTS}/{{patient}}/umi/{{stype}}_unaligned.bam"),
    params:
        rs1 = READ_STRUCTURE_R1,
        rs2 = READ_STRUCTURE_R2,
        sample = lambda wc: f"{wc.patient}_{wc.stype}",
        lib    = lambda wc: f"{wc.patient}_{wc.stype}_lib",
        fgbio  = config["tools"]["fgbio"],
    threads: config["threads"]["general"]
    log: f"{RESULTS}/{{patient}}/logs/fgbio_{{stype}}.log"
    shell:
        """
        {params.fgbio} FastqToBam \
            --input {input.r1} {input.r2} \
            --read-structures {params.rs1} {params.rs2} \
            --output {output.ubam} \
            --sample {params.sample} \
            --library {params.lib} \
            --platform ILLUMINA \
            --sort true \
        2> {log}
        """


rule ubam_to_fastq:
    """Recover adapter-trimming-ready FASTQs from the unaligned BAM."""
    input:
        ubam = f"{RESULTS}/{{patient}}/umi/{{stype}}_unaligned.bam",
    output:
        r1 = temp(f"{RESULTS}/{{patient}}/umi/{{stype}}_extracted_R1.fastq.gz"),
        r2 = temp(f"{RESULTS}/{{patient}}/umi/{{stype}}_extracted_R2.fastq.gz"),
    params:
        samtools = config["tools"]["samtools"],
    threads: config["threads"]["general"]
    log: f"{RESULTS}/{{patient}}/logs/ubam_to_fastq_{{stype}}.log"
    shell:
        """
        {params.samtools} sort -n -@ {threads} {input.ubam} \
        | {params.samtools} fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null -s /dev/null \
            -N \
        2> {log}
        """


rule cutadapt:
    """Trim adapter sequences; quality-trim 3' ends."""
    input:
        r1 = f"{RESULTS}/{{patient}}/umi/{{stype}}_extracted_R1.fastq.gz",
        r2 = f"{RESULTS}/{{patient}}/umi/{{stype}}_extracted_R2.fastq.gz",
    output:
        r1 = temp(f"{RESULTS}/{{patient}}/trimmed/{{stype}}_R1.fastq.gz"),
        r2 = temp(f"{RESULTS}/{{patient}}/trimmed/{{stype}}_R2.fastq.gz"),
        report = f"{RESULTS}/{{patient}}/trimmed/{{stype}}_cutadapt.json",
    params:
        cutadapt = config["tools"]["cutadapt"],
        adapter_r1 = ADAPTER_R1,
        adapter_r2 = ADAPTER_R2,
        min_len    = MIN_LEN,
        qual       = QUAL,
    threads: config["threads"]["cutadapt"]
    log: f"{RESULTS}/{{patient}}/logs/cutadapt_{{stype}}.log"
    shell:
        """
        {params.cutadapt} \
            -a {params.adapter_r1} \
            -A {params.adapter_r2} \
            -q {params.qual},{params.qual} \
            --minimum-length {params.min_len} \
            -j {threads} \
            --json {output.report} \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} \
        > {log} 2>&1
        """
