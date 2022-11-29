rule download_control_genomes:
	output:
		puc19 = config["reference"]["dir"] + "/" + "pUC19.fa",
		lamb = config["reference"]["dir"] + "/" + "lambda.fa"
	shell:
		"""
		wget -q -O - "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M77789.2&rettype=fasta&retmode=text" | awk 'FNR==1 {{print; next}} {{printf "%s", $0}} END {{print ""}}' > {output.puc19}
		wget -q -O - https://www.encodeproject.org/files/lambda.fa/@@download/lambda.fa.fasta.gz | gunzip -c > {output.lamb}
		"""
