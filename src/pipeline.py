'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='hiplexpipe')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        # Hi-Plex example: OHI031002-P02F04_S318_L001_R1_001.fastq
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq.gz'),

        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        # Hi-Plex example: OHI031002-P02F04_S318_L001_R2_001.fastq
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_R2_{lib[0]}.fastq.gz'),

        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}', '{lib[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}.clipped.bam')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.clipped.bam'),
        output='.clipped.sort.bam')

    # High quality and primary alignments
    pipeline.transform(
        task_func=stages.primary_bam,
        name='primary_bam',
        input=output_from('sort_bam_picard'),
        filter=suffix('.clipped.sort.bam'),
        output='.clipped.sort.hq.bam')

    # index bam file
    pipeline.transform(
        task_func=stages.index_sort_bam_picard,
        name='index_bam',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.clipped.sort.hq.bam.bai')

    # generate mapping metrics.
    pipeline.transform(
        task_func=stages.intersect_bed,
        name='intersect_bed',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.intersectbed.bam')

    pipeline.transform(
        task_func=stages.coverage_bed,
        name='coverage_bed',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.bedtools_hist_all.txt')

    pipeline.transform(
        task_func=stages.genome_reads,
        name='genome_reads',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.mapped_to_genome.txt')

    pipeline.transform(
        task_func=stages.target_reads,
        name='target_reads',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.mapped_to_target.txt')

    pipeline.transform(
        task_func=stages.total_reads,
        name='total_reads',
        input=output_from('align_bwa'),
        filter=suffix('.clipped.bam'),
        output='.total_raw_reads.txt')

#    pipeline.transform(
#        task_func=stages.generate_stats,
#        name='generate_stats',
#        input=output_from(['coverage_bed', 'genome_reads', 'target_reads', 'total_reads']), 
#        filter=formatter('.+/(?P<sample>.+).txt'),
#        extras=['{sample[0]}'],
#        output='all_sample.summary.txt')

    ###### GATK VARIANT CALLING ######
    # Call variants using GATK
    (pipeline.transform(
        task_func=stages.call_haplotypecaller_gatk,
        name='call_haplotypecaller_gatk',
        input=output_from('primary_bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-_]+).clipped.sort.hq.bam'),
        output='variants/gatk/{sample[0]}.g.vcf')
        .follows('index_sort_bam_picard'))

    # Genotype G.VCF files using GATK, separately genotypes replicate pairs and controls
    pipeline.merge(
        task_func=stages.genotype_gvcf_gatk_replicates,
        name='genotype_gvcf_gatk_replicates',
        input=output_from('call_haplotypecaller_gatk'),
        output='variants/gatk/replicates_controls.combined.raw.vcf')

    pipeline.transform(
         task_func=stages.variant_annotator_gatk,
         name='variant_annotator_gatk',
         input=output_from('genotype_gvcf_gatk_replicates'),
         filter=suffix('.raw.vcf'),
         output='.raw.annotate.vcf')

#### split snps and indels for filtering ####

    pipeline.transform(
        task_func=stages.select_variants_snps_gatk,
        name='select_variants_snps_gatk',
        input=output_from('variant_annotator_gatk'),
        filter=suffix('raw.annotate.vcf'),
        output='raw.annotate.snps.vcf')

    pipeline.transform(
        task_func=stages.select_variants_indels_gatk,
        name='select_variants_indels_gatk',
        input=output_from('variant_annotator_gatk'),
        filter=suffix('raw.annotate.vcf'),
        output='raw.annotate.indels.vcf')

    pipeline.transform(
        task_func=stages.apply_variant_filtration_snps_gatk,
        name='apply_variant_filtration_snps_gatk',
        input=output_from('select_variants_snps_gatk'),
        filter=suffix('raw.annotate.snps.vcf'),
        output='raw.annotate.snps.filtered.vcf')

    pipeline.transform(
        task_func=stages.apply_variant_filtration_indels_gatk,
        name='apply_variant_filtration_indels_gatk',
        input=output_from('select_variants_indels_gatk'),
        filter=suffix('raw.annotate.indels.vcf'),
        output='raw.annotate.indels.filtered.vcf')

    (pipeline.transform(
        task_func=stages.merge_filtered_vcfs_gatk,
        name='merge_filtered_vcfs_gatk',
        input=output_from('apply_variant_filtration_snps_gatk'),
        filter=suffix('.raw.annotate.snps.filtered.vcf'),
        add_inputs=add_inputs(['variants/gatk/replicates_controls.combined.raw.annotate.indels.filtered.vcf']),
        output='.raw.annotate.filtered.merged.vcf')
        .follows('apply_variant_filtration_indels_gatk'))

    return pipeline
