rule modbam2bed_v20:
        input:
                modbam = config["data"]["out_v20"] + "/sorted_filtered_aligned_chm13v2.0_HG002_UL.bam",
		index = config["data"]["out_v20"] + "/sorted_filtered_aligned_chm13v2.0_HG002_UL.bam.bai",
		genome = config["reference"]["ref_v20_fa"]
        output: config["data"]["beds"] + "/chm13v2.0_HG002_UL_80_guppy.bed"
	threads: 48
	conda: "../envs/modbam2bed.yaml"
        shell:
                """
                modbam2bed -t {threads} \
		-e \
		-m 5mC \
		--cpg \
		-a 0.20 -b 0.80 \
		{input.genome} {input.modbam} > {output}
                """

rule modbam2bed_v27:
        input:
                modbam = config["data"]["out_v27"] + "/sorted_filtered_aligned_chm13v2.7XY_HG002_UL.bam",
		index = config["data"]["out_v27"] + "/sorted_filtered_aligned_chm13v2.7XY_HG002_UL.bam.bai",
                genome = config["reference"]["ref_v27"]
        output: config["data"]["beds"] + "/chm13v2.7XY_HG002_UL_80_guppy.bed"
        threads: 48
        conda: "../envs/modbam2bed.yaml"
        shell:
                """
                modbam2bed -t {threads} \
                -e \
                -m 5mC \
		--cpg \
                -a 0.20 -b 0.80 \
                {input.genome} {input.modbam} > {output}
                """

rule convert2bismark_v20:
	input: config["data"]["beds"] + "/chm13v2.0_HG002_UL_80_guppy.bed"
	output: config["data"]["beds"] + "/chm13v2.0_HG002_UL_80_guppy.bismark"
	conda: "../envs/modbam2bed.yaml"
	shell:
		"""
		awk -v OFS='\t' '{{print $1,$3,$6,$13,$12,"C","CG"}}' {input} > {output}
		"""

rule convert2bismark_v27:
        input: config["data"]["beds"] + "/chm13v2.7XY_HG002_UL_80_guppy.bed"
        output: config["data"]["beds"] + "/chm13v2.7XY_HG002_UL_80_guppy.bismark"
        conda: "../envs/modbam2bed.yaml"
        shell:
                """
                awk -v OFS='\t' '{{print $1,$3,$6,$13,$12,"C","CG"}}' {input} > {output}
                """
