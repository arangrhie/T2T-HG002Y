rule modbam2bed_pacbio:
        input:
                modbam = config["modbam2bed_pacbio"]["input"] + "/chm13v1p1_hg002XYv2p7.HG002.bam",
		genome = config["reference"]["dir"] + "/" + config["reference"]["basename"] + ".fasta"
        output: bed = config["modbam2bed_pacbio"]["bed_output"] + "/hg002_pacbio_primrose_80.bed"
	threads: 32
	conda: "../envs/modbam2bed.yaml"
        shell:
                """
                modbam2bed -t {threads} \
		-e \
		-m 5mC \
		--cpg \
		-a 0.20 -b 0.80 \
		{input.genome} {input.modbam} > {output.bed}
                """

rule convert2bismark:
	input: 
		bed = config["modbam2bed_pacbio"]["bed_output"] + "/hg002_pacbio_primrose_80.bed"
	output:
		config["modbam2bed_pacbio"]["bed_output"] + "/hg002_pacbio_primrose_80.bismark"
	conda: "../envs/modbam2bed.yaml"
	shell:
		"""
		awk -v OFS='\t' '{{print $1,$3,$6,$13,$12,"C","CG"}}' {input.bed} > {output}
		"""
