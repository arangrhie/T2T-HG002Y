rule download_data:
	output:
		R1 = config["data"]["dir"] + "/" + "{srr}_1.fastq.gz",
		R2 = config["data"]["dir"] + "/" + "{srr}_2.fastq.gz"
	params:
		out_dir = config["data"]["dir"],
		accession = "{srr}"

	conda: "../envs/ffq.yaml"

	shell:
		"""
		ffq --ftp {params.accession} | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs wget -P {params.out_dir}
		"""
