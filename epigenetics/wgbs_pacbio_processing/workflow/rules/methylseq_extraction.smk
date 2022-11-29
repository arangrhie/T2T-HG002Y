rule methylseq_extraction:
	input:
		dedup_bam = config["short_read"]["bismark_methylseq"] + "/methylseq_{methylseq_sr_sample}_pe.deduplicated.bam",
		bismark_indexes_dir = directory(config["reference"]["dir"] + "/bowtie2_combined")
         
	params:
		outdir = config["short_read"]["bismark_methylseq"] + "/meth_extraction",
		repo = config["software"]["dir"] + "/" + "Bismark-master"
         
	log:
		log = "../reports/bismark/extraction_methylseq_{methylseq_sr_sample}.deduplicated.log"
         
	output:
		cpg = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/CpG_context_methylseq_{methylseq_sr_sample}_pe.deduplicated.txt.gz",
		non_cpg = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/Non_CpG_context_methylseq_{methylseq_sr_sample}_pe.deduplicated.txt.gz"

	conda: "../envs/aligners.yaml"

	threads: 6 
              
	shell: 
		"""
		{params.repo}/bismark_methylation_extractor -p \
		--comprehensive \
		--mbias_off \
		--merge_non_CpG \
		-o {params.outdir} \
		--gzip \
		--parallel {threads} \
		--ignore_r2 11 \
		--genome_folder {input.bismark_indexes_dir} \
		{input.dedup_bam} &> {log.log}
		"""
