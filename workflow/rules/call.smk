if config["fastq"]["restrict_regions"]:
    rule BedToIntervalList:
        input:
            bed=config["fastq"]["restrict_regions"],
            dict=gatk_dict["dict"],
        output:
            "results/called/bed.interval_list"
        log:
            "logs/call/BedToIntervalList.log",
        params:
            # optional parameters
            extra="--SORT true",  # sort output interval list before writing
        # optional specification of memory usage of the JVM that snakemake will respect with global
        # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
        # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
        # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        resources:
            mem_mb=1024,
        wrapper:
            config["warpper_mirror"]+"bio/picard/bedtointervallist"

rule PreprocessIntervals:
    input:
        ref=gatk_dict["ref"],
        intervals="results/called/bed.interval_list",
    output:
        "results/called/preprocessed.interval_list"
    params:
        bin_length=config["PreprocessIntervals"]["bin_length"],
        padding=config["PreprocessIntervals"]["padding"],
        extra=config["PreprocessIntervals"]["extra"],
    log:
        "logs/call/PreprocessIntervals.log"
    conda:
        "../envs/PreprocessIntervals.yaml"
    script:
        "../scripts/PreprocessIntervals.py"

# rule RemoveChrSex:
#     input:
#         get_picard_intervals,
#     output:
#         "results/called/RemoveChrSex.interval_list"
#     shell:
#         "grep -v 'chr[XY]' {input} > {output}"

rule AnnotateIntervals:
    input:
        ref=gatk_dict["ref"],
        intervals="results/called/preprocessed.interval_list",
    output:
        "results/called/annotated.intervals.tsv"
    params:
        extra="--interval-merging-rule OVERLAPPING_ONLY"
    log:
        "logs/call/AnnotateIntervals.log"
    conda:
        "../envs/AnnotateIntervals.yaml"
    script:
        "../scripts/AnnotateIntervals.py"

rule CollectReadCounts:
    input:
        bam="results/prepared/{s}.{g}.cram",
        intervals="results/called/preprocessed.interval_list",
    output:
        counts="results/called/{s}.{g}.counts.hdf5",
    log:
        "logs/call/CollectReadCounts/{s}.{g}.log",
    params:
        extra="-R {}".format(gatk_dict["ref"]),
        mergingRule="OVERLAPPING_ONLY",  # optional
    resources:
        mem_mb=1024
    wrapper:
        config["warpper_mirror"]+"bio/gatk/collectreadcounts"

rule CreateReadCountPanelOfNormals:
    input:
        hdf5 = expand("results/called/{s}.N.counts.hdf5",s=samples.Sample.unique()),
        gc_interval = "results/called/annotated.intervals.tsv"
    output:
        "results/called/CNV.PoN.hdf5"
    threads: 32
    log:
        "logs/call/CreateReadCountPanelOfNormals.log"
    params:
        extra="--minimum-interval-median-percentile 5.0"
    conda:
        "../envs/CreateReadCountPanelOfNormals.yaml"
    script:
        "../scripts/CreateReadCountPanelOfNormals.py"

rule DenoiseReadCounts:
    input:
        hdf5=["results/called/{s}.T.counts.hdf5"],
        pon="results/called/CNV.PoN.hdf5",
        gc_interval="results/called/annotated.intervals.tsv"
    output:
        std_copy_ratio="results/called/{s}.T.standardizedCR.tsv",
        denoised_copy_ratios="results/called/{s}.T.denoisedCR.tsv",
    log:
        "logs/call/DenoiseReadCounts/{s}.T.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/denoisereadcounts"

rule CollectAllelicCounts:
    input:
        bam="results/prepared/{s}.{g}.cram",
        sam_idx="results/prepared/{s}.{g}.cram.crai",
        intervals="results/called/preprocessed.interval_list",
        ref=gatk_dict["ref"],
    output:
        counts="results/called/{s}.{g}.counts.tsv"
    log:
        "logs/call/CollectAllelicCounts/{s}.{g}.log",
    params:
        extra="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/collectalleliccounts"

rule ModelSegments:
    input:
        denoised_copy_ratios="results/called/{s}.T.denoisedCR.tsv",
        allelic_counts="results/called/{s}.T.counts.tsv",
        normal_allelic_counts="results/called/{s}.N.counts.tsv"
    output:
        "results/called/{s}.T.cr.seg",
        "results/called/{s}.T.modelFinal.seg",
        "results/called/{s}.T.hets.tsv",
    log:
        "logs/call/ModelSegments/{s}.log",
    params:
        prefix="{s}.T",
        extra="",  # optional
    conda:
        "../envs/ModelSegments.yaml"
    resources:
        mem_mb=1024
    script:
        "../scripts/ModelSegments.py"

rule CallCopyRatioSegments:
    input:
        copy_ratio_seg="results/called/{s}.T.cr.seg",
    output:
        called_copy_ratio_seg="results/called/{s}.T.called.cr.seg",
        igv_seg="results/called/{s}.T.called.igv.seg",
    log:
        "logs/call/CallCopyRatioSegments/{s}.log",
    params:
        prefix="{s}.T",
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        config["warpper_mirror"]+"bio/gatk/callcopyratiosegments"

rule PlotModeledSegments:
    input:
        allelic_counts="results/called/{s}.T.hets.tsv",
        segments="results/called/{s}.T.modelFinal.seg",
        dict=gatk_dict["dict"],
        denoised_copy_ratios="results/called/{s}.T.denoisedCR.tsv",
    output:
       report(
            "results/called/PlotModeledSegments/{s}.modeled.png",
            caption="../report/call.rst",
            category="Call",
        ), 
    log:
        "logs/call/PlotModeledSegments/{s}.log"
    params:
        output_prefix="{s}",
        extra="--minimum-contig-length 46709983"
    conda:
        "../envs/PlotModeledSegments.yaml"
    # script:
    #     "../scripts/PlotModeledSegments.py"
    shell:
        '''
        gatk PlotModeledSegments \
        --allelic-counts {input.allelic_counts} \
        --denoised-copy-ratios {input.denoised_copy_ratios} \
        --sequence-dictionary {input.dict} \
        --segments {input.segments} \
        --output-prefix {params.output_prefix} \
        {params.extra} \
        --output `dirname {output[0]}` &> {log}
        '''

rule PlotDenoisedCopyRatios:
    input:
        dict=gatk_dict["dict"],
        std_copy_ratios="results/called/{s}.T.standardizedCR.tsv",
        denoised_copy_ratios="results/called/{s}.T.denoisedCR.tsv",
    output:
        multiext("results/called/PlotDenoisedCopyRatios/{s}",
                ".denoised.png",
                ".deltaMAD.txt",
                ".denoisedMAD.txt",
                ".standardizedMAD.txt",
                ".scaledDeltaMAD.txt")
    log:
        "logs/call/PlotDenoisedCopyRatios/{s}.log"
    params:
        output_prefix="{s}",
        extra="--minimum-contig-length 46709983"
    conda:
        "../envs/PlotDenoisedCopyRatios.yaml"
    # script:
    #     "../scripts/PlotDenoisedCopyRatios.py"
    shell:
        '''
        gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios {input.std_copy_ratios} \
        --denoised-copy-ratios {input.denoised_copy_ratios} \
        --sequence-dictionary {input.dict} \
        --output-prefix {params.output_prefix} \
        {params.extra} \
        --output `dirname {output[0]}` &> {log}
        '''