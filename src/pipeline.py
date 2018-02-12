'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from, regex
from stages import Stages
import re


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

    # 1. Align paired end reads in FASTQ to the reference producing a BAM file
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

    # 2. Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.clipped.bam'),
        output='.clipped.sort.bam')

    # 3. High quality and primary alignments
    pipeline.transform(
        task_func=stages.primary_bam,
        name='primary_bam',
        input=output_from('sort_bam_picard'),
        filter=suffix('.clipped.sort.bam'),
        output='.clipped.sort.hq.bam')

    # 4. index bam file
    pipeline.transform(
        task_func=stages.index_sort_bam_picard,
        name='index_bam',
        input=output_from('primary_bam'),
        filter=suffix('.clipped.sort.hq.bam'),
        output='.clipped.sort.hq.bam.bai')

    # generate mapping metrics.
    pipeline.transform(
        task_func=stages.generate_amplicon_metrics, 
        name='generate_amplicon_metrics', 
        input=output_from('primary_bam'), 
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).clipped.sort.hq.bam'),
        output='alignments/metrics/{sample[0]}.amplicon-metrics.txt',
        extras=['{sample[0]}'])

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

    pipeline.collate(
        task_func=stages.generate_stats,
        name='generate_stats',
        input=output_from('coverage_bed', 'genome_reads', 'target_reads', 'total_reads'), 
        filter=regex(r'.+/(.+BS\d\d\d\d\d\d.+S\d+)\..+\.txt'),
        output=r'all_sample.summary.\1.txt',
        extras=[r'\1', 'all_sample.summary.txt'])

    ###### GATK VARIANT CALLING ######
    # 6. Call variants using GATK
    (pipeline.transform(
        task_func=stages.call_haplotypecaller_gatk,
        name='call_haplotypecaller_gatk',
        input=output_from('primary_bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-_]+).clipped.sort.hq.bam'),
        output='variants/gatk/{sample[0]}.g.vcf')
        .follows('index_sort_bam_picard'))

    # 7. Genotype G.VCF files using GATK, separately genotypes replicate pairs
    pipeline.collate(
        task_func=stages.genotype_gvcf_gatk_replicates_collate, 
        name='genotype_gvcf_gatk_replicates_collate', 
        input=output_from('call_haplotypecaller_gatk'), 
        filter=regex(r'.+(BS\d\d\d\d\d\d).+'), 
        output=r'variants/gatk/bstp_run1-7_replicates.combined.raw.\1.vcf') 

#    # 8. Merge genotyped replicate pairs into a single vcf
#    pipeline.merge(
#        task_func=stages.merge_genotyped_replicates, 
#        name='merge_genotyped_replicates', 
#        input=output_from('genotype_gvcf_gatk_replicates_collate'), 
#        output='variants/gatk/bstp_run1-7_replicates.combined.raw.vcf')

#    # 9. Add additional annotations to the genotyped vcf
#    pipeline.transform(
#         task_func=stages.variant_annotator_gatk,
#         name='variant_annotator_gatk',
#         input=output_from('merge_genotyped_replicates'),
#         filter=suffix('.raw.vcf'),
#         output='.raw.annotate.vcf')

#### split snps and indels for filtering ####
    # 10-snps. Extract SNPs 
    pipeline.transform(
        task_func=stages.select_variants_snps_gatk,
        name='select_variants_snps_gatk',
        input=output_from('genotype_gvcf_gatk_replicates_collate'),
        filter=suffix('vcf'),
        output='snps.vcf')

    # 10-indels. Extract indels
    pipeline.transform(
        task_func=stages.select_variants_indels_gatk,
        name='select_variants_indels_gatk',
        input=output_from('genotype_gvcf_gatk_replicates_collate'),
        filter=suffix('vcf'),
        output='indels.vcf')

    # 11-snps. Filter SNPs
    pipeline.transform(
        task_func=stages.apply_variant_filtration_snps_gatk,
        name='apply_variant_filtration_snps_gatk',
        input=output_from('select_variants_snps_gatk'),
        filter=suffix('snps.vcf'),
        output='snps.filtered.vcf')

    # 11-indels. Filter indels
    pipeline.transform(
        task_func=stages.apply_variant_filtration_indels_gatk,
        name='apply_variant_filtration_indels_gatk',
        input=output_from('select_variants_indels_gatk'),
        filter=suffix('indels.vcf'),
        output='indels.filtered.vcf')

    # 12. Merge filtered SNPs and indels
    pipeline.collate(
        task_func=stages.merge_filtered_vcfs_gatk,
        name='merge_filtered_vcfs_gatk',
        input=output_from('apply_variant_filtration_snps_gatk', 'apply_variant_filtration_indels_gatk'),
        filter=regex(r'.+(BS\d\d\d\d\d\d).+'),
        output=r'variants/gatk/bstp_run1-7_replicates.combined.raw.filtered.merged.\1.vcf')

    return pipeline



#     pipeline.merge(
#         task_func=stages.genotype_gvcf_gatk_replicates,
#         name='genotype_gvcf_gatk_replicates',
#         input=output_from('call_haplotypecaller_gatk'),
#         output='variants/gatk/bstp_run1-7_replicates.combined.raw.vcf')
