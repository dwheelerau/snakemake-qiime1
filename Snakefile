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
PATTERN_CLN_R1 = config['cln_r1']
PATTERN_CLN_R2 = config['cln_r2']

# overlap params
MAX_OVERLAP = config['max_overlap']
MIN_OVERLAP = config['min_overlap']
DIRS = ['flash/', 'logs/']

rule all:
    input:
        DIRS,
        expand('projectData/{sample}_R1_001.cln.fastq.gz', sample=SAMPLES),

rule set_up:
    output: DIRS
    shell:
        "mkdir -p "+' '.join(DIRS)

rule flash_reads:
    input:
        r1 = join('projectData', PATTERN_R1),
        r2 = join('projectData', PATTERN_R2)
    output:
        "flash/"
    params:
        minov=MIN_OVERLAP,
        maxov=MAX_OVERLAP
    log:
        "logs/flash.log"
    shell:
        "flash --max-overlap {params.maxov} --min-overlap {params.minov} "
        "-o {output} {input.r1} {input.r2} > {log} 2>&1"
