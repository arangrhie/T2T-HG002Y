rule emseq_bismark_deduplicate:
        input: 
        	bam = config["short_read"]["bismark_emseq"] + "/emseq_{emseq_sr_sample}_pe.bam"

        output:
                dedup_bam = config["short_read"]["bismark_emseq"] + "/emseq_{emseq_sr_sample}_pe.deduplicated.bam"
        log:
                log = "../reports/bismark/dedup.emseq_{emseq_sr_sample}_pe.log"

        params: 
                basename = "emseq_{emseq_sr_sample}_pe",
                outdir = config["short_read"]["bismark_emseq"],
                tempdir = "/scratch/phook",
                repo = config["software"]["dir"] + "/" + "Bismark-master" 

        conda: "../envs/samtools.yaml"

	threads: 8

        shell:
                """
                {params.repo}/deduplicate_bismark --bam \
                --paired \
                --outfile {params.basename} \
		--output_dir {params.outdir} \
                --parallel {threads} \
                {input.bam} &> {log.log}
                """
