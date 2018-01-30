'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os
import math

# PICARD_JAR = '$PICARD_HOME/lib/picard-1.69.jar'
PICARD_JAR = '/usr/local/picard/2.9.2/picard.jar'
SNPEFF_JAR = '/usr/local/easybuild/software/snpEff/4.1d-Java-1.7.0_80/snpEff.jar'

GATK_JAR = '/usr/local/gatk/3.7/executable/GenomeAnalysisTK.jar'

def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)

def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)

class Stages(object):
    def __init__(self, state):
        self.state = state
        self.reference = self.get_options('ref_grch37')
        self.dbsnp_hg19 = self.get_options('dbsnp_hg19')
        self.mills_hg19 = self.get_options('mills_hg19')
        self.one_k_g_snps = self.get_options('one_k_g_snps')
        self.one_k_g_indels = self.get_options('one_k_g_indels')
        self.one_k_g_highconf_snps = self.get_options('one_k_g_highconf_snps')
        self.hapmap = self.get_options('hapmap')
        self.bamclipper = self.get_options('bamclipper')
        self.coord_file = self.get_options('coord_file')
        self.interval_file = self.get_options('interval_file')
        self.primer_file = self.get_options('primer_file')
        self.primer_bedpe_file = self.get_options('primer_bedpe_file')
        self.proportionthresh = self.get_options('proportionthresh')
        self.absthresh = self.get_options('absthresh')
        self.maxvariants = self.get_options('maxvariants')
        self.other_vep = self.get_options('other_vep')

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, PICARD_JAR, mem, args)

    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, GATK_JAR, mem, args)

    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)

    def original_fastqs(self, output):
        '''Original fastq files'''
        pass

    def align_bwa(self, inputs, bam_out, sample_id, lib):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        safe_make_dir('alignments')
        read_group = '"@RG\\tID:{sample}\\tSM:{sample}\\tPU:lib1\\tPL:Illumina"' \
            .format(sample=sample_id)
        command = 'bwa mem -M -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| {bamclipper} -i -p {primer_bedpe_file} -n 1 ' \
                  '| samtools view -b -h -o {bam} -' \
                  .format(cores=cores,
                          read_group=read_group,
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          reference=self.reference,
                          bamclipper=self.bamclipper,
                          primer_bedpe_file=self.primer_bedpe_file,
                          bam=bam_out)
        run_stage(self.state, 'align_bwa', command)

    def apply_undr_rover(self, inputs, vcf_output, sample_id):
        '''Apply undr_rover to call variants from paired end fastq files'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('apply_undr_rover', 'cores')
        safe_make_dir('variants/undr_rover')
        safe_make_dir('variants/undr_rover/coverdir')
        coverdir = "variants/undr_rover/coverdir"
        coverfile = sample_id + ".coverage"

        command = 'undr_rover --primer_coords {coord_file} ' \
                  '--primer_sequences {primer_file} ' \
                  '--reference {reference} ' \
                  '--out {vcf_output} ' \
                  '--coverfile {coverdir}/{coverfile} ' \
                  '--proportionthresh {proportionthresh} ' \
                  '--absthresh {absthresh} ' \
                  '--max_variants {maxvariants} ' \
                  '{fastq_read1} {fastq_read2}'.format(
                        coord_file=self.coord_file, primer_file=self.primer_file,
                        reference=self.reference,
                        vcf_output=vcf_output,
                        coverdir=coverdir,
                        proportionthresh=self.proportionthresh,
                        absthresh=self.absthresh,
                        maxvariants=self.maxvariants,
                        coverfile=coverfile,
                        fastq_read1=fastq_read1_in,
                        fastq_read2=fastq_read2_in)
        run_stage(self.state, 'apply_undr_rover', command)

    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)

    def primary_bam(self, bam_in, sbam_out):
        '''Only keep primary alignments in the BAM file using samtools'''
        command = 'samtools view -h -q 1 -f 2 -F 4 -F 8 -F 256 -b ' \
                    '-o {sbam_out} {bam_in}'.format(
                        bam_in=bam_in, sbam_out=sbam_out)
        run_stage(self.state, 'primary_bam', command)

    # index sorted bam file
    def index_sort_bam_picard(self, bam_in, bam_index):
        '''Index sorted bam using samtools'''
        command = 'samtools index {bam_in} {bam_index}'.format(
                          bam_in=bam_in, bam_index=bam_index)
        run_stage(self.state, 'index_sort_bam_picard', command)

    ##########
    def call_haplotypecaller_gatk(self, bam_in, vcf_out):
        '''Call variants using GATK'''
        safe_make_dir('variants/gatk')
        # safe_make_dir('variants}'.format(sample=sample_id))
        gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
                    "--emitRefConfidence GVCF " \
                    "-A AlleleBalance -A AlleleBalanceBySample " \
                    "-A ChromosomeCounts -A ClippingRankSumTest " \
                    "-A Coverage -A DepthPerAlleleBySample " \
                    "-A DepthPerSampleHC -A FisherStrand " \
                    "-A GCContent -A GenotypeSummaries " \
                    "-A HardyWeinberg -A HomopolymerRun " \
                    "-A LikelihoodRankSumTest -A LowMQ " \
                    "-A MappingQualityRankSumTest -A MappingQualityZero " \
                    "-A QualByDepth " \
                    "-A RMSMappingQuality -A ReadPosRankSumTest " \
                    "-A SampleList -A SpanningDeletions " \
                    "-A StrandBiasBySample -A StrandOddsRatio " \
                    "-A TandemRepeatAnnotator -A VariantType " \
                    "-I {bam} -L {interval_list} -o {out}".format(reference=self.reference,
                                                                  bam=bam_in, interval_list=self.interval_file, out=vcf_out)
        self.run_gatk('call_haplotypecaller_gatk', gatk_args)

    def combine_gvcf_gatk(self, vcf_files_in, vcf_out):
        '''Combine G.VCF files for all samples using GATK'''
        # g_vcf_files = ' '.join(['--variant ' + vcf for vcf in vcf_files_in])
        merge_commands = []
        temp_merge_outputs = []
        for n in range(0, int(math.ceil(float(len(vcf_files_in)) / 200.0))):
            start = n * 200
            filelist = vcf_files_in[start:start + 200]
            filelist_command = ' '.join(['--variant ' + vcf for vcf in filelist])
            temp_merge_filename = vcf_out.rstrip('.vcf') + ".temp_{start}.vcf".format(start=str(start))
            gatk_args_full = "java -Xmx{mem}g -jar {jar_path} -T CombineGVCFs -R {reference} " \
                             "--disable_auto_index_creation_and_locking_when_reading_rods " \
                             "{g_vcf_files} -o {vcf_out}; ".format(reference=self.reference, 
                                                                   jar_path=GATK_JAR, 
                                                                   mem=self.state.config.get_stage_options('combine_gvcf_gatk', 'mem'), 
                                                                   g_vcf_files=filelist_command, 
                                                                   vcf_out=temp_merge_filename)
            merge_commands.append(gatk_args_full)
            temp_merge_outputs.append(temp_merge_filename)

        final_merge_vcfs = ' '.join(['--variant ' + vcf for vcf in temp_merge_outputs])
        gatk_args_full_final = "java -Xmx{mem}g -jar {jar_path} -T CombineGVCFs -R {reference} " \
                               "--disable_auto_index_creation_and_locking_when_reading_rods " \
                               "{g_vcf_files} -o {vcf_out}".format(reference=self.reference, 
                                                                   jar_path=GATK_JAR, 
                                                                   mem=self.state.config.get_stage_options('combine_gvcf_gatk', 'mem'), 
                                                                   g_vcf_files=final_merge_vcfs, 
                                                                   vcf_out=vcf_out)

        merge_commands.append(gatk_args_full_final)
        final_command = ''.join(merge_commands)
        run_stage(self.state, 'combine_gvcf_gatk', final_command)

    def genotype_gvcf_gatk(self, combined_vcf_in, vcf_out):
        '''Genotype G.VCF files using GATK'''
        cores = self.get_stage_options('genotype_gvcf_gatk', 'cores')
        gatk_args = "-T GenotypeGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--dbsnp {dbsnp} " \
                    "--num_threads {cores} --variant {combined_vcf} --out {vcf_out}" \
                    .format(reference=self.reference, dbsnp=self.dbsnp_hg19,
                            cores=cores, combined_vcf=combined_vcf_in, vcf_out=vcf_out)
        self.run_gatk('genotype_gvcf_gatk', gatk_args)

    def variant_annotator_gatk(self, vcf_in, vcf_out):
        '''Annotate G.VCF files using GATK'''
        cores = self.get_stage_options('variant_annotator_gatk', 'cores')
        gatk_args = "-T VariantAnnotator -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-A AlleleBalance -A AlleleBalanceBySample " \
                    "-A ChromosomeCounts -A ClippingRankSumTest " \
                    "-A Coverage -A DepthPerAlleleBySample " \
                    "-A DepthPerSampleHC -A FisherStrand " \
                    "-A GCContent -A GenotypeSummaries " \
                    "-A HardyWeinberg -A HomopolymerRun " \
                    "-A LikelihoodRankSumTest " \
                    "-A MappingQualityRankSumTest -A MappingQualityZero " \
                    "-A QualByDepth " \
                    "-A RMSMappingQuality -A ReadPosRankSumTest " \
                    "-A SampleList -A SpanningDeletions " \
                    "-A StrandBiasBySample -A StrandOddsRatio " \
                    "-A TandemRepeatAnnotator -A VariantType " \
                    "--num_threads {cores} --variant {vcf_in} --out {vcf_out}" \
                    .format(reference=self.reference, cores=cores, vcf_in=vcf_in, vcf_out=vcf_out)
        self.run_gatk('variant_annotator_gatk', gatk_args)

    def apply_variant_filtration_gatk(self, inputs, vcf_out):
        '''Apply Variant Filtration using gatk'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_variant_filtration_gatk', 'cores')
        gatk_args = "-T VariantFiltration --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} " \
                    "--filterExpression \"QUAL < 30.0\" --filterName \"VeryLowQual\" " \
                    "--filterExpression \"QD < 2.0\" --filterName \"LowQD\" " \
                    "--filterExpression \"DP < 10\" --filterName \"LowCoverage\" " \
                    "--filterExpression \"MQ < 30.0\" --filterName \"LowMappingQual\" " \
                    "--filterExpression \"SOR > 4.0\" --filterName \"StrandBias\" " \
                    "--filterExpression \"HRun >= 7.0\" --filterName \"HRun7\" " \
                    "--filterExpression \"MQRankSum < -12.5\" --filterName \"MQRankSum\" " \
                    "--filterExpression \"ReadPosRankSum < -8.0\" --filterName \"ReadPosRankSum\" " \
                    "--variant {vcf_in} -o {vcf_out}".format(reference=self.reference,
                                                            cores=cores, vcf_in=vcf_in, vcf_out=vcf_out)
        self.run_gatk('apply_variant_filtration_gatk', gatk_args)

    def apply_vep(self, inputs, vcf_out):
        '''Apply VEP'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_vep', 'cores')
        vep_command = "vep --cache --dir_cache {other_vep} --assembly GRCh37 --refseq --offline --fasta {reference} " \
                    "-i {vcf_in} --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs " \
                    "--format vcf -o {vcf_vep} --force_overwrite --vcf " \
                    "--fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT," \
                    "Protein_position,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,HGVSc,HGVSp,cDNA_position,CDS_position,PICK " \
                    "--fork {threads} --flag_pick".format(
                    reference=self.reference, vcf_in=vcf_in, vcf_vep=vcf_out, other_vep=self.other_vep, threads=cores)
        run_stage(self.state, 'apply_vep', vep_command)


######  stats sections

    def intersect_bed(self, bam_in, bam_out):
        '''intersect the bed file with the interval file '''
        command = "intersectBed -abam {bam_in} -b {interval_file} > {bam_out} ".format(
                     bam_in=bam_in, interval_file=self.interval_file, bam_out=bam_out)
        run_stage(self.state, 'intersect_bed', command)           


    def coverage_bed(self, bam_in, txt_out):
        ''' make coverage files '''
        command = "coverageBed -b {bam_in} -a {interval_file} -hist | grep all > {txt_out}".format(
                     bam_in=bam_in, interval_file=self.interval_file, txt_out=txt_out)
        run_stage(self.state, 'coverage_bed', command)
    
    def genome_reads(self, bam_in, txt_out):
        '''count reads that map to the genome'''
        command = 'samtools view -c -F4 {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'genome_reads', command)

    def target_reads(self, bam_in, txt_out):
        '''count reads that map to target panel'''
        command = 'samtools view -c -F4 {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'target_reads', command)

    def total_reads(self, bam_in, txt_out):
        '''count the total number of reads that we started with'''
        command = 'samtools view -c {bam_in} > {txt_out}'.format(
                        bam_in=bam_in, txt_out=txt_out)
        run_stage(self.state, 'total_reads', command)

#    def generate_stats(self, inputs, txt_out):
#        '''run R stats script'''
#        a, b, c, d, e = inputs
#        command = 'Rscript --vanilla /projects/vh83/pipelines/code/modified_summary_stat.R \
#                    {hist_in} {map_genome_in} {map_target_in} {raw_reads_in} {sample_name} \
#                    {txt_out}'.format(hist_in=a, map_genome_in=b, map_target_in=c, raw_reads_in=d , sample_name=e , txt_out=txt_out)
#        run_stage(self.state, 'generate_stats', command)

#Rscript --vanilla ~/vh83/pipelines/code/summary_stat.R \
#        metrics/${sample_run_name}.sort.bedtools_hist_all.txt \
#        metrics/${sample_run_name}.sort.mapped_to_genome.txt \
#        metrics/${sample_run_name}.sort.mapped_to_target.txt \
#        metrics/${sample_run_name}.sort.total_raw_reads.txt \
#        ${sample_run_name} \
#        ${summary_prefix}_summary_coverage.txt

    def sort_vcfs(self, vcf_in, vcf_out):
        '''sort undr_rover vcf files'''
        command = 'bcftools sort -o {vcf_out} -O z {vcf_in}'.format(vcf_out=vcf_out, vcf_in=vcf_in)
        run_stage(self.state, 'sort_vcfs', command)
    
    def index_vcfs(self, vcf_in, vcf_out):
        command = 'bcftools index -f --tbi {vcf_in}'.format(vcf_in=vcf_in)
        run_stage(self.state, 'index_vcfs', command)
    
    def concatenate_vcfs(self, vcf_files_in, vcf_out):
        join_vcf_files = ' '.join([vcf for vcf in vcf_files_in])
        command = 'bcftools concat -a -O z -o {vcf_out} {join_vcf_files} '.format(vcf_out=vcf_out, join_vcf_files=join_vcf_files)
        run_stage(self.state, 'concatenate_vcfs', command) 

    def index_final_vcf(self, vcf_in, vcf_out):
        command = 'bcftools index -f --tbi {vcf_in}'.format(vcf_in=vcf_in)
        run_stage(self.state, 'index_final_vcf', command)






