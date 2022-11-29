rule emseq_lab1_bismark2bedgraph:
        input:
                cpg = expand(config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/CpG_context_emseq_{emseq_sr_lab1}_pe.deduplicated.txt.gz",emseq_sr_lab1=EMSEQ_SR_LAB1)

        params:
                outdir = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph",
                repo = config["software"]["dir"] + "/" + "Bismark-master",
                outname = "seqc2_emseq_lab1_CpG",

        log:
                log = "../reports/bismark/bismark2bedgraph_seqc2_emseq_lab1_CpG.log"

        output:
                bedGraph = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab1_CpG.gz",
                cov = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab1_CpG.gz.bismark.cov.gz"

        conda: "../envs/aligners.yaml"

        shell:
                """
                {params.repo}/bismark2bedGraph \
                --remove_spaces \
                --dir {params.outdir} \
                --buffer_size 20G \
                -o {params.outname} \
                {input.cpg} &> {log.log}
                """

rule emseq_lab2_bismark2bedgraph:
        input:
                cpg = expand(config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/CpG_context_emseq_{emseq_sr_lab2}_pe.deduplicated.txt.gz",emseq_sr_lab2=EMSEQ_SR_LAB2)

        params:
                outdir = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph",
                repo = config["software"]["dir"] + "/" + "Bismark-master",
                outname = "seqc2_emseq_lab2_CpG",

        log:
                log = "../reports/bismark/bismark2bedgraph_seqc2_emseq_lab2_CpG.log"

        output:
                bedGraph = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab2_CpG.gz",
                cov = config["short_read"]["bismark_emseq"] + "/meth_extraction" + "/bedgraph" + "/seqc2_emseq_lab2_CpG.gz.bismark.cov.gz"

        conda: "../envs/aligners.yaml"

        shell:
                """
                {params.repo}/bismark2bedGraph \
                --remove_spaces \
                --dir {params.outdir} \
                --buffer_size 20G \
                -o {params.outname} \
                {input.cpg} &> {log.log}
                """
