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
        "flash_summary_results.txt",
        "qiime/1_qc-fastqs/seqs.fna",
        "qiime/2_usearch61_ref_chimera_checking/chimeras.txt"

rule clean:
    shell:
        """
        rm -f logs/*
        rm -f flash/*
        rm -f *.txt
        rm -rf qiime/
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
        indicator=config['read_indicator']
    output:
        "qiime/1_qc-fastqs/seqs.fna"
    shell:
        "multiple_split_libraries_fastq.py --read_indicator {params.indicator} "
        "-i flash/ -o qiime/1_qc-fastqs/ -p ./test.params"

rule identify_chimera:
    input:
        "qiime/1_qc-fastqs/seqs.fna"
    output:
        "qiime/2_usearch61_ref_chimera_checking/chimeras.txt"
    params:
        ref_db=config['ref_db']
    shell:
        "identify_chimeric_seqs.py -m usearch61 -i {input} -r {params.ref_db} "
        "-o qiime/2_usearch61_ref_chimera_checking/"

# -n is important, will remove chimeras and create a clean file
rule filter_chimera:
    input:
        fastas="qiime/1_qc-fastqs/seqs.fna",
        chimeras="qiime/2_usearch61_chimera_checking/chimeras.txt"
    output:
        "qiime/3_chimerafree-fastqs/seqs_chimeras_filtered.fna"
    shell:
        "filter_fasta.py -f {input.fastas} -o {output} "
        "-s {input.chimeras} -n"
#
#!pick_open_reference_otus.py -i ../qiime/1_qc-fastqs/seqs_chimeras_filtered.fna -o ../qiime/2_chimfree_otu-pick-uclust-97pct-16sgg -p 16S_open_params_97per.txt
#
#### below is from old script, need to integreate above files into this pipeline
## neeed to fix below to work
#rule summarize_table:
#    input:
#        "sample.otu_table_wTax.biom"
#    output:
#        "sample.otu_table_summary.txt"
#    shell:
#        """
#        biom summarize-table -i {input} -o {output}
#        cat {output}
#        echo
#        echo ' ~~~~~~ Caution! ~~~~~~ '
#        echo 'you need the above data to set the SAM_DEPTH (--sampling_depth)'
#        echo 'param in config for the core_diversity_analysis rule'
#        """
## Align sequences command
#rule align_otus:
#    input:
#        "all.otus.fasta"
#    output:
#        "pynast_aligned_seqs/all.otus_aligned.fasta"
#    shell:
#        "align_seqs.py -i {input} -o pynast_aligned_seqs/"
#
#rule filter_aln:
#    input:
#        "pynast_aligned_seqs/all.otus_aligned.fasta"
#    output:
#        "pynast_aligned_seqs/all.otus_aligned_pfiltered.fasta"
#    shell:
#        "filter_alignment.py -o pynast_aligned_seqs/ -i {input}"
#
#rule make_phylogeny:
#    input:
#        "pynast_aligned_seqs/all.otus_aligned_pfiltered.fasta"
#    output:
#        "all.otus.tre"
#    shell:
#        "make_phylogeny.py -i {input} -o {output}"
#
#rule core_div_analysis:
#    input:
#        summary="sample.otu_table_summary.txt",
#        biom="sample.otu_table_wTax.biom",
#        meta=MAP_FILE,
#        tree="all.otus.tre"
#    output:
#        "coreout/"
#    params:
#        sample_depth=SAM_DEPTH
#    shell:
#        """
#        core_diversity_analyses.py -i {input.biom} -o tmpcore/ -m {input.meta} \
#        -t {input.tree} -e {params.sample_depth}
#        mv tmpcore/* {output}
#        rmdir tmpcore/
#        """
#
#onerror:
#        print("Biom error? Have you run qiime1 from the commmand line to activate it?")
#        #shell("mail -s "an error occurred" youremail@provider.com < {log}")
#
#rule clean:
#    shell:
#        """
#        rm -f tmp/*
#        rm -f logs/*
#        rm -f readStats/*
#        rm -f *.uc
#        rm -f *.fasta
#        rm -f *.fastq
#        rm -f all*
#        rm -f *.txt
#        rm -f *.stats
#        rm -f *.biom
#        rm -rf tax*
#        rm -rf pynast_aligned_seq/
#        rm -rf coreout/
#        """
