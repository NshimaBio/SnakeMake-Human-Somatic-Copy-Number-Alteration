__author__ = "Victor Wang"
__copyright__ = "Copyright 2023, Victor Wang"
__email__ = "victor@bioquest.cn"
__license__ = "Apache License 2.0"

import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk --java-options '{java_opts}' PlotModeledSegments"
        " --allelic-counts {snakemake.input.allelic_counts}"
        " --denoised-copy-ratios {snakemake.input.denoised_copy_ratios}"
        " --sequence-dictionary {snakemake.input.dict}"
        " --segments {snakemake.input.segments}"
        " --output-prefix {snakemake.params.output_prefix}"
        " {extra}"
        " --tmp-dir {tmpdir}"
        " --output {snakemake.output[0]}"
        " {log}"
    )