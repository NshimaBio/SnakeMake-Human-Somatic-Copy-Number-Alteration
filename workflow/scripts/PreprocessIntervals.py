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

intervals = snakemake.input.get("intervals","")
if intervals:
    intervals = f"--intervals {intervals}"

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk --java-options '{java_opts}' PreprocessIntervals"
        " --reference {snakemake.input.ref}"
        " --bin-length {snakemake.params.bin_length}"
        " --padding {snakemake.params.padding}"
        " {extra}"
        " {intervals}"
        " --tmp-dir {tmpdir}"
        " --output {snakemake.output}"
        " {log}"
    )