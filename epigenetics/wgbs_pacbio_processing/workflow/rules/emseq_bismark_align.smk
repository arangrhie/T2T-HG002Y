rule emseq_bismark_align:
        input:
                trimR1 = config["short_read"]["trim_output"] + "/emseq_{emseq_sr_sample}_1.trimmed.fastq",
                trimR2 = config["short_read"]["trim_output"] + "/emseq_{emseq_sr_sample}_2.trimmed.fastq",
                bismark_indexes_dir = directory(config["reference"]["dir"] + "/bowtie2_combined")

        output:
                bam = config["short_read"]["bismark_emseq"] + "/emseq_{emseq_sr_sample}_pe.bam",
                report = config["short_read"]["bismark_emseq"] + "/emseq_{emseq_sr_sample}_PE_report.txt"
        log:
                log = "../reports/bismark/emseq_{emseq_sr_sample}.log"

        params: 
                basename = "emseq_{emseq_sr_sample}",
                outdir = config["short_read"]["bismark_emseq"],
                tempdir = "/scratch/phook",
		repo = config["software"]["dir"] + "/" + "Bismark-master" 
	
	#threads: 2 # this will use 8 threads per sample

	conda: "../envs/aligners.yaml"

        shell:
                """
                {params.repo}/bismark --bam \
                --bowtie2 \
                --genome {input.bismark_indexes_dir} \
                --basename {params.basename} \
                -1 {input.trimR1} -2 {input.trimR2} \
                --temp_dir {params.tempdir} \
                --output_dir {params.outdir} &> {log.log}
                """
