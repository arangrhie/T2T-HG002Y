rule emseq_mbias:
	input:
		dedup_bam = config["short_read"]["bismark_emseq"] + "/emseq_{emseq_sr_sample}_pe.deduplicated.bam",
		bismark_indexes_dir = directory(config["reference"]["dir"] + "/bowtie2_combined")
         
	params:
		outdir = config["short_read"]["bismark_emseq"] + "/mbias",
		repo = config["software"]["dir"] + "/" + "Bismark-master"
         
	log:
		log = "../reports/bismark/mbias_emseq_{emseq_sr_sample}_pe.deduplicated.log"
         
	output:
		config["short_read"]["bismark_emseq"] + "/mbias" + "/emseq_{emseq_sr_sample}_pe.deduplicated.M-bias.txt"

	conda: "../envs/aligners.yaml"

	threads: 4 
              
	shell: 
		"""
		{params.repo}/bismark_methylation_extractor -p \
		-o {params.outdir} \
		--mbias_only \
		--parallel {threads} \
		--genome_folder {input.bismark_indexes_dir} \
		{input.dedup_bam} &> {log.log}
		"""
