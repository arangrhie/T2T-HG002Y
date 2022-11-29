rule emseq_merge_sort:
        input:
                emseq = expand(config["short_read"]["bismark_emseq"] + "/emseq_{emseq_sr_samples}_pe.deduplicated.bam",emseq_sr_samples=EMSEQ_SR_SAMPLES)

        output:
                merge_bam = config["short_read"]["bismark_emseq"] + "/merge" + "/merged_HG002_emseq.bam",
		sort_bam = config["short_read"]["bismark_emseq"] + "/merge" + "/sorted_merged_HG002_emseq.bam"

        conda: "../envs/aligners.yaml"

        shell:
                """
                samtools merge -@24 -o {output.merge_bam} {input.emseq}
		samtools sort -@24 {output.merge_bam} > {output.sort_bam}
                """

rule wgbs_merge_sort:
        input:
                wgbs = expand(config["short_read"]["bismark_methylseq"] + "/methylseq_{methylseq_sr_samples}_pe.deduplicated.bam",methylseq_sr_samples=METHYLSEQ_SR_SAMPLES)

        output:
                merge_bam = config["short_read"]["bismark_methylseq"] + "/merge" + "/merged_HG002_wgbs.bam",
                sort_bam = config["short_read"]["bismark_methylseq"] + "/merge" + "/sorted_merged_HG002_wgbs.bam"

        conda: "../envs/aligners.yaml"

        shell:
                """
                samtools merge -@24 -o {output.merge_bam} {input.wgbs}
                samtools sort -@24 {output.merge_bam} > {output.sort_bam}
                """
