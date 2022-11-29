rule emseq_lab1_coverage2cytosine:
        input:
                cov = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab1_CpG.gz.bismark.cov.gz",
		bismark_indexes_dir = directory(config["reference"]["dir"] + "/bowtie2_combined")

        params:
                outdir = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph",
                repo = config["software"]["dir"] + "/" + "Bismark-master",
                outname = "seqc2_emseq_lab1"

        log:
                log = "../reports/bismark/coverage2cytosine_seqc2_emseq_lab1_CpG.log"

        output: config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab1.CpG_report.txt.gz"

        conda: "../envs/aligners.yaml"

        shell:
                """
                {params.repo}/coverage2cytosine \
                --dir {params.outdir} \
		--genome_folder {input.bismark_indexes_dir} \
		--gzip -o {params.outname} {input.cov} &> {log.log}
                """

rule emseq_lab2_coverage2cytosine:
        input:
                cov = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab2_CpG.gz.bismark.cov.gz",
                bismark_indexes_dir = directory(config["reference"]["dir"] + "/bowtie2_combined")

        params:
                outdir = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph",
                repo = config["software"]["dir"] + "/" + "Bismark-master",
                outname = "seqc2_emseq_lab2"

        log:
                log = "../reports/bismark/coverage2cytosine_seqc2_emseq_lab2_CpG.log"

        output: config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab2.CpG_report.txt.gz"

        conda: "../envs/aligners.yaml"

        shell:
                """
                {params.repo}/coverage2cytosine \
                --dir {params.outdir} \
                --genome_folder {input.bismark_indexes_dir} \
                --gzip -o {params.outname} {input.cov} &> {log.log}
                """
