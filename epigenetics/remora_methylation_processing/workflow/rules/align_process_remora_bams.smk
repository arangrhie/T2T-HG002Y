rule unaligned_modbam_fastq_conversion:
	input: lambda wc: df[df['sample'] == wc['sample']].infile
        output: config["data"]["fastq_out"] + "/" + "{sample}.fastq",
	threads: 48
	conda: "../envs/meth_alignment.yaml"
	shell:
                """
                samtools fastq -@ {threads} -T Mm,Ml {input} > {output}
                """

rule winnowmap_alignment_v20:
	input: expand(config["data"]["fastq_out"] + "/" + "{sample}.fastq",sample=df["sample"])
	output: config["data"]["out_v20"] + "/aligned_chm13v2.0_HG002_UL.sam"
	params:
                kmer_v20 = config["reference"]["kmer_v20"],
                ref_v20 = config["reference"]["ref_v20"]
	threads: 24
        conda: "../envs/meth_alignment.yaml"
	shell:
		"""
		winnowmap -t {threads} -W {params.kmer_v20} -ax map-ont -y {params.ref_v20} {input} > {output}
                """

rule winnowmap_alignment_v27:
        input: expand(config["data"]["fastq_out"] + "/" + "{sample}.fastq",sample=df["sample"])
        output: config["data"]["out_v27"] + "/aligned_chm13v2.7XY_HG002_UL.sam"
        params:
                kmer_v27 = config["reference"]["kmer_v27"],
                ref_v27 = config["reference"]["ref_v27"]
        threads: 24
        conda: "../envs/meth_alignment.yaml"
        shell:
                """
                winnowmap -t {threads} -W {params.kmer_v27} -ax map-ont -y {params.ref_v27} {input} > {output}
                """

rule sam2sortedbam_v20:
	input: config["data"]["out_v20"] + "/aligned_chm13v2.0_HG002_UL.sam"
	output: config["data"]["out_v20"] + "/sorted_filtered_aligned_chm13v2.0_HG002_UL.bam"
	threads: 24
	params: temp = config["data"]["out_v20"]
	conda: "../envs/meth_alignment.yaml"
	shell:
                """
                samtools view -@ {threads} -Sb -F 256 -F 2048 {input} | samtools sort -@ {threads} -T {params.temp} - > {output}
                """

rule sam2sortedbam_v27:
        input: config["data"]["out_v27"] + "/aligned_chm13v2.7XY_HG002_UL.sam"
        output: config["data"]["out_v27"] + "/sorted_filtered_aligned_chm13v2.7XY_HG002_UL.bam"
        threads: 24
	params: temp = config["data"]["out_v27"]
        conda: "../envs/meth_alignment.yaml"
        shell:
                """
                samtools view -@ {threads} -Sb -F 256 -F 2048 {input} | samtools sort -@ {threads} -T {params.temp} - > {output}
                """

rule index_v20:
	input: config["data"]["out_v20"] + "/sorted_filtered_aligned_chm13v2.0_HG002_UL.bam"
	output: config["data"]["out_v20"] + "/sorted_filtered_aligned_chm13v2.0_HG002_UL.bam.bai"
	threads: 48
	conda: "../envs/meth_alignment.yaml"
	shell:
		"""
		samtools index -@ {threads} {input}
		"""

rule index_v27:
        input: config["data"]["out_v27"] + "/sorted_filtered_aligned_chm13v2.7XY_HG002_UL.bam"
        output: config["data"]["out_v27"] + "/sorted_filtered_aligned_chm13v2.7XY_HG002_UL.bam.bai"
        threads: 48
        conda: "../envs/meth_alignment.yaml"
        shell:
                """
                samtools index -@ {threads} {input}
                """
