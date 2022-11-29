rule methylseq_bismark2bedgraph:
	input:
		cpg = expand(config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/CpG_context_methylseq_{methylseq_sr_sample}_pe.deduplicated.txt.gz",methylseq_sr_sample=METHYLSEQ_SR_SAMPLES)
        
	params:
		outdir = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/bedgraph",
		repo = config["software"]["dir"] + "/" + "Bismark-master",
		outname = "seqc2_methylseq_CpG",
		#inputdir = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/CpG_*"
         
	log:
		log = "../reports/bismark/bismark2bedgraph_seqc2_methylseq_CpG.log"
         
	output:
		bedGraph = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_methylseq_CpG.gz",
		cov = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_methylseq_CpG.gz.bismark.cov.gz"

	conda: "../envs/aligners.yaml"
              
	shell: 
		"""
		{params.repo}/bismark2bedGraph \
		--remove_spaces \
		--dir {params.outdir} \
		--buffer_size 20G \
		-o {params.outname} \
		{input.cpg} &> {log.log}
		"""
