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

gc_interval = snakemake.input.get("gc_interval", "")
if gc_interval:
    gc_interval = f"--annotated-intervals {gc_interval}"
    
hdf5=snakemake.input.get("hdf5",False)

if not hdf5:
    raise ValueError("hdf5 file is None")
if isinstance(hdf5,str):
    hdf5=[hdf5]
hdf5 = " ".join(["--input {}".format(x) for x in hdf5])

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk --java-options '{java_opts}' CreateReadCountPanelOfNormals"
        " {hdf5}"
        " {gc_interval}"
        " {extra}"
        " --tmp-dir {tmpdir}"
        " --output {snakemake.output}"
        " {log}"
    )