from os.path import join
import glob
from datetime import datetime, date, time

configfile: "./config.yaml"

# globals
now = datetime.now().strftime('%Y-%m-%d-%H-%M')
date_string = "{}".format(now)

# key step to get sample names from R1 read
SAMPLES, = glob_wildcards("projectData/{sample}_L001_R1_001.fastq.gz")
PATTERN_R1 = config['pat_r1']
PATTERN_R2 = config['pat_r2']

# overlap params
MAX_OVERLAP = config['max_overlap']
MIN_OVERLAP = config['min_overlap']
DIRS = ['flash/', 'logs/']

rule all:
    input:
        DIRS,
        expand('projectData/{sample}_L001_R1_001.fastq.gz', sample=SAMPLES),
        expand("flash/{sample}.extendedFrags.fastq", sample=SAMPLES),
        "flash_summary_results.txt"

rule clean:
    shell:
        """
        rm -f logs/*
        rm -f flash/*
        rm -f *.txt
        """

rule set_up:
    output: DIRS
    shell:
        "mkdir -p "+' '.join(DIRS)

rule flash_reads:
    input:
        r1 = join('projectData', PATTERN_R1),
        r2 = join('projectData', PATTERN_R2)
    params:
        minov=MIN_OVERLAP,
        maxov=MAX_OVERLAP
    log:
        stdout="logs/flash.{sample}.log",
        stderr="logs/flash.{sample}.stderr"
    output:
        "flash/{sample}.extendedFrags.fastq"
    shell:
        "flash --max-overlap {params.maxov} --min-overlap {params.minov} "
        "-o flash/{wildcards.sample} {input.r1} {input.r2} > {log.stdout} 2>{log.stderr}"

rule flash_results:
    input:
        expand("flash/{sample}.extendedFrags.fastq", sample=SAMPLES)
    output:
        "flash_summary_results.txt"
    shell:
        """
        for log in logs/flash.*.log; do
        comb=$(grep 'Percent' $log | cut -d':' -f2);
        name=$(echo $log | cut -d. -f2)
        echo -e "$name\t$comb" >> {output};
        done;
        find logs/ -name *.stderr -type f -empty -delete;
        """

rule quality_filter:
    input:
        expand("flash/{sample}.extendedFrags.fastq", sample=SAMPLES)
    params:
        phred=config['phred_threshold'],
        bad_run=config['bad_run'],
        fracn=config['length_fraction'],
        indicator=config['read_indicator']
    shell:
        "multiple_split_libraries_fastq.py "
        "--barcode_type not-barcoded --min_per_read_length_fraction {params.fracn} "
        "--read_indicator {params.indicator} --max_bad_run_length {params.bad_run} "
        "-i flash/ -o qiime/1_qc-fastqs/"

