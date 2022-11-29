rule concatenate_genomes:
	input:
		puc19 = config["reference"]["dir"] + "/" + "pUC19.fa",
		lamb = config["reference"]["dir"] + "/" + "lambda.fa",
		human = config["reference"]["dir"] + "/" + config["reference"]["basename"] + ".fasta"
	output:
		combined_genome = config["reference"]["dir"] + "/{aligner}_combined/methylation_reference.fa"
	shell:
		"""
    		cat {input.human} {input.puc19} {input.lamb} > {output.combined_genome}
		"""

rule index_genomes:
	input:
		combined_genome = config["reference"]["dir"] + "/{aligner}_combined/methylation_reference.fa"
	output: 
		combined_genome_index = config["reference"]["dir"] + "/{aligner}_combined/methylation_reference.fa.fai"
	conda: "../envs/samtools.yaml"
	shell:
		"""
		samtools faidx {input.combined_genome}
		"""
