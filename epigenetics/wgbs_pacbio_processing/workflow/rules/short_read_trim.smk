rule sr_methylseq_trim:
        input:
                R1 = config["short_read"]["input"] + "/{methylseq_sr_sample}_1.fastq.gz",
                R2 = config["short_read"]["input"] + "/{methylseq_sr_sample}_2.fastq.gz"
        params:
                outdir = config["short_read"]["trim_output"],
                basename = "{methylseq_sr_sample}"
        log:
                log = "../reports/trim/methylseq_sr_{methylseq_sr_sample}_trim.log"
        output:
                trimR1 = config["short_read"]["trim_output"] + "/methylseq_{methylseq_sr_sample}_1.trimmed.fastq",
                trimR2 = config["short_read"]["trim_output"] + "/methylseq_{methylseq_sr_sample}_2.trimmed.fastq"
	threads: 16
	conda: "../envs/fastp.yaml"
        shell:
                """
                fastp --in1 {input.R1} --in2 {input.R2} \
		--out1 {output.trimR1} \
		--out2 {output.trimR2} \
		--trim_front1 10 --trim_tail1 10 --trim_front2 10 \
		-l 2 -Q -z 4 --trim_poly_g --overrepresentation_analysis --thread {threads} 2> {log.log}
                """

rule sr_emseq_trim:
	input: R1 = config["short_read"]["input"] + "/{emseq_sr_sample}_1.fastq.gz", R2 = config["short_read"]["input"] + "/{emseq_sr_sample}_2.fastq.gz"
  	params: outdir = config["short_read"]["trim_output"], basename = "{emseq_sr_sample}"
  	log: log = "../reports/trim/emseq_sr_{emseq_sr_sample}_trim.log"
  	output: trimR1 = config["short_read"]["trim_output"] + "/emseq_{emseq_sr_sample}_1.trimmed.fastq", trimR2 = config["short_read"]["trim_output"] + "/emseq_{emseq_sr_sample}_2.trimmed.fastq"
	conda: "../envs/fastp.yaml"
	threads: 16
	shell: 
		"""
		fastp --in1 {input.R1} --in2 {input.R2} \
		--out1 {output.trimR1} \
		--out2 {output.trimR2} \
		-l 2 -Q -z 4 --trim_poly_g --overrepresentation_analysis --thread {threads} 2> {log.log}
   		"""
