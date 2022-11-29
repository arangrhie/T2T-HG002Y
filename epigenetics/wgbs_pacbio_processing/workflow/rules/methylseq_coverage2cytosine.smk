rule methylseq_coverage2cytosine:
        input:
                cov = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_methylseq_CpG.gz.bismark.cov.gz",
		bismark_indexes_dir = directory(config["reference"]["dir"] + "/bowtie2_combined")

        params:
                outdir = config["short_read"]["bismark_methylseq"] + "/meth_extraction" + "/bedgraph",
                repo = config["software"]["dir"] + "/" + "Bismark-master",
                outname = "seqc2_methylseq"

        log:
                log = "../reports/bismark/coverage2cytosine_methylseq_CpG.log"

        output: config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_methylseq.CpG_report.txt.gz"

        conda: "../envs/aligners.yaml"

        shell:
                """
                {params.repo}/coverage2cytosine \
                --dir {params.outdir} \
		--genome_folder {input.bismark_indexes_dir} \
		--gzip -o {params.outname} {input.cov} &> {log.log}
                """
