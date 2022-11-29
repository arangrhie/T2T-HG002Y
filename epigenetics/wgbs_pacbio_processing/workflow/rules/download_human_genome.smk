rule download_human_genome:
	output:
		config["reference"]["dir"] + "/" + config["reference"]["basename"] + ".fasta"
	params:
		out_dir = config["reference"]["dir"],
		url = config["reference"]["url"] + "/" + config["reference"]["basename"] + ".fasta",
	shell:
		"""
		wget -nc -P {params.out_dir} {params.url}
		"""
