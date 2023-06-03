`GATK best practices workflow`_ Pipeline summary

=============================================
Reference
=============================================
Reference genome related files and GTAK budnle files (GATK_)

=============================================
Prepare
=============================================
1. Adapter trimming (Fastp_)
2. Aligner (`BWA mem2`_)
3. Mark duplicates (samblaster_)
4. Generates recalibration table for Base Quality Score Recalibration (BaseRecalibrator_)
5. Apply base quality score recalibration (ApplyBQSR_)
6. Merge CRAMs of every sample, repesectly (Picard_)
7. Create CRAM index (samtools_)

=============================================
Quality control report
=============================================
1. Fastp report (MultiQC_)
2. Alignment report (MultiQC_)

=============================================
Somatic copy number variants (CNA)
=============================================
1. Converts the captured regions BED file to a Picard Interval List for target seqencing (`BedToIntervalList (Picard)`_)
1. Prepares bins for coverage collection (PreprocessIntervals_)
2. AnnotateIntervals Annotates intervals with GC content, mappability, and segmental-duplication content (AnnotateIntervals_)
3. Creates a panel of normals for read-count denoising (CreateReadCountPanelOfNormals_)
4. Denoises read counts to produce denoised copy ratios (DenoiseReadCounts_)
5. Collects reference and alternate allele counts at specified sites (CollectAllelicCounts_)
6. Models segmented copy ratios from denoised copy ratios and segmented minor-allele fractions from allelic counts (ModelSegments_)
7. Calls copy-ratio segments as amplified, deleted, or copy-number neutral (CallCopyRatioSegments_)
8. Creates plots of denoised and segmented copy-ratio and minor-allele-fraction estimates (PlotModeledSegments_)
9. Creates plots of denoised copy ratios (PlotDenoisedCopyRatios_)

.. _GATK best practices workflow: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
.. _GATK: https://software.broadinstitute.org/gatk/
.. _VEP: https://www.ensembl.org/info/docs/tools/vep/index.html
.. _fastp: https://github.com/OpenGene/fastp
.. _BWA mem2: http://bio-bwa.sourceforge.net/
.. _samblaster: https://github.com/GregoryFaust/samblaster
.. _BaseRecalibrator: https://gatk.broadinstitute.org/hc/en-us/articles/13832708374939-BaseRecalibrator
.. _ApplyBQSR: https://github.com/GregoryFaust/samblaster
.. _Picard: https://broadinstitute.github.io/picard
.. _Mutect2: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-Mutect2
.. _GetPileupSummaries: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-GetPileupSummaries
.. _CalculateContamination: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-CalculateContamination
.. _LearnReadOrientationModel: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-LearnReadOrientationModel
.. _FilterMutectCalls: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-FilterMutectCalls
.. _MultiQC: https://multiqc.info
.. _samtools: http://www.htslib.org
.. _PreprocessIntervals: https://gatk.broadinstitute.org/hc/en-us/articles/13832754597915-PreprocessIntervals
.. _BedToIntervalList (Picard): https://gatk.broadinstitute.org/hc/en-us/articles/13832706340763-BedToIntervalList-Picard-
.. _AnnotateIntervals: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-AnnotateIntervals
.. _CreateReadCountPanelOfNormals: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-CreateReadCountPanelOfNormals
.. _DenoiseReadCounts: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-DenoiseReadCounts
.. _CollectAllelicCounts: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-CollectAllelicCounts
.. _ModelSegments: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-ModelSegments
.. _CallCopyRatioSegments: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-CallCopyRatioSegments
.. _PlotModeledSegments: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-PlotModeledSegments
.. _PlotDenoisedCopyRatios: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-PlotDenoisedCopyRatios
