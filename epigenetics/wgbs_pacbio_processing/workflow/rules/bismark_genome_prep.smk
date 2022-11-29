rule bismark_prepare_genome:
	input: fa = config["reference"]["dir"] + "/{aligner}_combined/methylation_reference.fa", genome_prep = config["software"]["dir"] + "/" + "Bismark-master" + "/" + "bismark_genome_preparation"
    	output: directory(config["reference"]["dir"] + "/{aligner}_combined/Bisulfite_Genome/")
    	log:
        	log = "../reports" + "/{aligner}.bsgenomeprep.log"
	params:
		type = "{aligner}"
	conda: "../envs/aligners.yaml" 
    	threads: 24
    	shell:
		"""
		dir=`echo "{input.fa}" | cut -d'/' -f1-7`
        	{input.genome_prep}  --{params.type} --parallel {threads} --genomic_composition --verbose ${{dir}} 2> {log.log}
		"""
