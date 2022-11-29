rule bismark_download:
	output: 
		config["software"]["dir"] + "/" + "master.zip"
	params:
		out_dir = config["software"]["dir"]
	shell:
		"""
		wget -P {params.out_dir}  https://github.com/FelixKrueger/Bismark/archive/refs/heads/master.zip
		"""

rule bismark_unzip:
	input:
		config["software"]["dir"] + "/" + "master.zip"
	output: 
		dir = directory(config["software"]["dir"] + "/" + "Bismark-master"),
		command = config["software"]["dir"] + "/" + "Bismark-master" + "/" + "bismark_genome_preparation"
	params:
		out_dir = config["software"]["dir"]
	shell:
		"""
		unzip {input} -d {params.out_dir}
		"""
