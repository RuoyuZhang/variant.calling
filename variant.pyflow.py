#!/usr/bin/env python

import os
import sys
sys.path.append("/home/fs01/rz253/.local/lib/python2.7/site-packages/pyflow/pyflow/src")
from pyflow import WorkflowRunner
import argparse
import distutils.util

def Config(config_file):
    config = dict()
    with open(config_file) as fh:
        for i,line in enumerate(fh):
            line = line.strip()
            [name,absdir] = line.split('=')
            config[name] = absdir
    
    return config


class VariantCalling(WorkflowRunner):

    def __init__(
            self,
            sample,
            fastq,
            out_dir,
            ref,
            config_file,
            stages,
    ):

        for i in ('sample','fastq','out_dir','ref','config_file','stages'):
            if i not in locals():
                raise Exception('--%s is required for %s' % (i, self.__class__.__name__))
            setattr(self, i, locals()[i])

        self.config = Config(self.config_file)  
        
        if not os.path.exists(self.out_dir + '/stat'):
            os.makedirs(self.out_dir + '/stat')
              
        self.prefix = '%s/%s' % (self.out_dir, self.sample)
            
        # trimmatic
        if len(self.fastq) == 2:
            self.trim_command = '%s PE -threads 3 %s %s %s.trimmed.R1.fastq %s.trimmed.R1.un.fastq %s.trimmed.R2.fastq %s.trimmed.R2.un.fastq ILLUMINACLIP:/home/fs01/rz253/bin/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>%s.trim.log' \
                                %(self.config['trim'], self.fastq[0], self.fastq[1], self.prefix, self.prefix, self.prefix,self.prefix,self.prefix)
        else:
            self.trim_command = '%s SE -threads 3 %s %s.trimmed.fastq ILLUMINACLIP:/home/fs01/rz253/bin/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>%s.trim.log' \
                                %(self.config['trim'], self.fastq[0], self.prefix, self.prefix)
        # bwa
        if len(self.fastq) == 2:
            self.bwa_command = '%s mem -t 6 %s %s.trimmed.R1.fastq %s.trimmed.R2.fastq > %s.bwa.sam 2>%s.bwa.log' \
                               % (self.config['bwa'], self.ref, self.prefix, self.prefix, self.prefix, self.prefix)
        else:
            self.bwa_command = '%s mem -t 6 %s %s.trimmed.fastq > %s.bwa.sam 2>%s.bwa.log' \
                               % (self.config['bwa'], self.ref, self.prefix, self.prefix, self.prefix)
    
        # sam_dedup
        self.sam_dedup_cm1 = '%s SamFormatConverter INPUT=%s.bwa.sam OUTPUT=%s.bwa.bam 2>%s.log' %(self.config['picard'],self.prefix,self.prefix,self.prefix)
        self.sam_dedup_cm2 = '%s SortSam INPUT=%s.bwa.bam OUTPUT=%s.bwa.sort.bam SORT_ORDER=coordinate 2>>%s.log'\
                        % (self.config['picard'],self.prefix,self.prefix,self.prefix)
        self.sam_dedup_cm3 = '%s AddOrReplaceReadGroups INPUT=%s.bwa.sort.bam OUTPUT=%s.sort.gp.bam SORT_ORDER=coordinate RGID=%s.gp RGLB=%s.gb RGPL=illumina RGSM=%s.gp RGPU=barcode 2>>%s.log'\
                           % (self.config['picard'],self.prefix,self.prefix,self.sample,self.sample,self.sample,self.prefix)
        self.sam_dedup_cm4 = '%s MarkDuplicates REMOVE_DUPLICATES=true INPUT=%s.sort.gp.bam OUTPUT=%s.sort.gp.rmdup.bam M=%s.duplicate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=False 2>>%s.log'\
                             % (self.config['picard'],self.prefix,self.prefix,self.prefix,self.prefix)
        self.sam_dedup_cm5 = '%s BuildBamIndex INPUT=%s.sort.gp.rmdup.bam 2>>%s.log' %(self.config['picard'],self.prefix,self.prefix) 
        
        self.sam_dedup_cm6 = '%s CollectAlignmentSummaryMetrics R=%s I=%s.sort.gp.bam O=%s/stat/alignment_metrics.txt' % (self.config['picard'],self.ref,self.prefix,self.out_dir)
        self.sam_dedup_cm7 = '%s CollectInsertSizeMetrics INPUT=%s.sort.gp.rmdup.bam OUTPUT=%s/stat/insert_metrics.txt HISTOGRAM_FILE=%s/stat/insert_size_histogram.pdf'\
                             % (self.config['picard'],self.prefix,self.out_dir,self.out_dir)
        self.sam_dedup_cm8 = '%s depth -a %s.sort.gp.rmdup.bam > %s/stat/depth_out.txt' % (self.config['samtools'],self.prefix,self.out_dir)
        
           
        # indel realigner
        self.realign_cm1 = '%s -T RealignerTargetCreator -I %s.sort.gp.rmdup.bam -R %s -o %s.indel.list 2>>%s.log' \
                           % (self.config['gatk'],self.prefix,self.ref,self.prefix,self.prefix)
        self.realign_cm2 = '%s -T IndelRealigner -l INFO -I %s.sort.gp.rmdup.bam -R %s -targetIntervals %s.indel.list -o %s.sort.gp.rmdup.realign.bam 2>>%s.log'\
                           % (self.config['gatk'],self.prefix,self.ref,self.prefix,self.prefix,self.prefix)
        
        # recalibarate
        self.recal_cm1 = '%s -T HaplotypeCaller -R %s -I %s.sort.gp.rmdup.realign.bam -o %s.raw_variants.vcf 2>>%s.log' \
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.recal_cm2 = '%s -T SelectVariants -R %s -V %s.raw_variants.vcf -selectType SNP -o %s.raw_snps.vcf 2>>%s.log'\
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.recal_cm3 = '%s -T SelectVariants -R %s -V %s.raw_variants.vcf -selectType INDEL -o %s.raw_indels.vcf 2>>%s.log'\
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.recal_cm4 = "%s -T VariantFiltration -R %s -V %s.raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName \"basic_snp_filter\" -o %s.filtered_snps.vcf 2>>%s.log" \
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.recal_cm5 = "%s -T VariantFiltration -R %s -V %s.raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName \"basic_indel_filter\" -o %s.filtered_indels.vcf 2>>%s.log" \
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.recal_cm6 = '%s -T BaseRecalibrator -R %s -I %s.sort.gp.rmdup.realign.bam -knownSites %s.filtered_snps.vcf -knownSites %s.filtered_indels.vcf -o %s.recal_data.table 2>>%s.log' \
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix)
        self.recal_cm7 = '%s -T PrintReads -R %s -I %s.sort.gp.rmdup.realign.bam -BQSR %s.recal_data.table -o %s.sort.gp.rmdup.realign.recal.bam 2>>%s.log'\
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix,self.prefix)
        
        # variants calling
        self.vcf_cm1 = '%s -T HaplotypeCaller -R %s -I %s.sort.gp.rmdup.realign.recal.bam -o %s.raw_variants_recal.vcf 2>>%s.log'\
                       % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.vcf_cm2 = '%s -T SelectVariants -R %s -V %s.raw_variants_recal.vcf -selectType SNP -o %s.raw_snps_recal.vcf 2>>%s.log'\
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.vcf_cm3 = '%s -T SelectVariants -R %s -V %s.raw_variants_recal.vcf -selectType INDEL -o %s.raw_indels_recal.vcf 2>>%s.log'\
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.vcf_cm4 = "%s -T VariantFiltration -R %s -V %s.raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName \"basic_snp_filter\" -o %s.filtered_snps_final.vcf 2>>%s.log" \
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        self.vcf_cm5 = "%s -T VariantFiltration -R %s -V %s.raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName \"basic_indel_filter\" -o %s.filtered_indels_final.vcf 2>>%s.log" \
                         % (self.config['gatk'],self.ref,self.prefix,self.prefix,self.prefix)
        
        # annotation
        
        # rm
                  
    def workflow(self):

        trim_task = None
        if 'all' in self.stages or 'trim' in self.stages:
            trim_task = self.addTask('trim', self.trim_command, isForceLocal = True, nCores=1)
        
        bwa_task = None    
        if 'all' in self.stages or 'bwa' in self.stages:
            bwa_task = self.addTask('bwa', self.bwa_command, dependencies=trim_task)
        
        sam_dedup_task = {'sam_convert':None, 'sam_sort':None, 'sam_add_RG':None, 'sam_dedup':None, 'sam_index':None}
        if 'all' in self.stages or 'sam_dedup' in self.stages:
            sam_dedup_task['sam_convert'] = self.addTask('sam_convert', self.sam_dedup_cm1, dependencies=bwa_task)
            sam_dedup_task['sam_sort'] = self.addTask('sam_sort', self.sam_dedup_cm2, dependencies=sam_dedup_task['sam_convert'])
            sam_dedup_task['sam_add_RG'] = self.addTask('sam_add_RG', self.sam_dedup_cm3, dependencies=sam_dedup_task['sam_sort'])
            sam_dedup_task['sam_dedup'] = self.addTask('sam_dedup', self.sam_dedup_cm4, dependencies=sam_dedup_task['sam_add_RG'])
            sam_dedup_task['sam_index'] = self.addTask('sam_index', self.sam_dedup_cm5, dependencies=sam_dedup_task['sam_dedup'])
            sam_dedup_task['align_stat'] = self.addTask('align_stat', self.sam_dedup_cm6, dependencies=sam_dedup_task['sam_add_RG'])
            sam_dedup_task['insert_stat'] = self.addTask('insert_stat', self.sam_dedup_cm7, dependencies=sam_dedup_task['sam_dedup'])
            sam_dedup_task['depth_stat'] = self.addTask('depth_stat', self.sam_dedup_cm8, dependencies=sam_dedup_task['sam_dedup'])
        
        realign_task = {'indel_list':None,'realign':None}
        if 'all' in self.stages or 'realign' in self.stages:
            realign_task['indel_list'] = self.addTask('indel_list', self.realign_cm1, dependencies=sam_dedup_task['sam_index'])
            realign_task['realign'] = self.addTask('realign', self.realign_cm2, dependencies=realign_task['indel_list'])
        
        recal_task = {'recal_vcf':None,'recal_take_snp':None,'recal_take_indel':None,'recal_filter_snp':None,'recal_filter_indel':None,'recal_recal':None,'print_reads':None}
        if 'all' in self.stages or 'recal' in self.stages:
            recal_task['recal_vcf'] = self.addTask('recal_vcf', self.recal_cm1, dependencies=realign_task['realign'])
            recal_task['recal_take_snp'] = self.addTask('recal_take_snp', self.recal_cm2, dependencies=recal_task['recal_vcf'])
            recal_task['recal_take_indel'] = self.addTask('recal_take_indel', self.recal_cm3, dependencies=recal_task['recal_vcf'])
            recal_task['recal_filter_snp'] = self.addTask('recal_filter_snp', self.recal_cm4, dependencies=recal_task['recal_take_snp'])
            recal_task['recal_filter_indel'] = self.addTask('recal_filter_indel', self.recal_cm5, dependencies=recal_task['recal_take_indel'])
            recal_task['recal_recal'] = self.addTask('recal_recal', self.recal_cm6, dependencies=[recal_task['recal_filter_indel'],recal_task['recal_filter_snp']])
            recal_task['print_reads'] = self.addTask('print_reads',self.recal_cm7, dependencies=recal_task['recal_recal'])
        
        vcf_task = {'vcf_vcf':None,'vcf_take_snp':None,'vcf_take_indel':None,'vcf_filter_snp':None,'vcf_filter_indel':None}
        if 'all' in self.stages or 'vcf' in self.stages:
            vcf_task['vcf_vcf'] = self.addTask('vcf_vcf', self.vcf_cm1, dependencies=recal_task['print_reads'])
            vcf_task['vcf_take_snp'] = self.addTask('vcf_take_snp', self.vcf_cm2, dependencies=vcf_task['vcf_vcf'])
            vcf_task['vcf_take_indel'] = self.addTask('vcf_take_indel', self.vcf_cm3, dependencies=vcf_task['vcf_vcf'])
            vcf_task['vcf_filter_snp'] = self.addTask('vcf_filter_snp', self.vcf_cm4, dependencies=vcf_task['vcf_take_snp'])
            vcf_task['vcf_filter_indel'] = self.addTask('vcf_filter_indel', self.vcf_cm5, dependencies=vcf_task['vcf_take_indel'])            
            

if __name__ == "__main__":
    usage = "usage: %prog [opstatistics of rawtions] \n\
            test pyflow\n\
            April 21 2017 \n\
            "

    #parser = OptionParser(usage=usage)
    parser = argparse.ArgumentParser(description = usage)
    
    parser.add_argument('--sample', type=str, required=True, help="the sample name")
    parser.add_argument('--fastq', type=str, required=True, nargs="+", help="the fastq file[s]")
    parser.add_argument('--out_dir', type=str, required=True, help="output directory")
    parser.add_argument('--ref', type=str, required=True, help="reference genome")
    parser.add_argument('--config_file', type=str, required=True, help="config")
    parser.add_argument('--stages', type=str, required=True, nargs="+", help="stages to run")
    parser.add_argument('--isContinue', type=distutils.util.strtobool, default='False')
    parser.add_argument('--isDryRun', type=distutils.util.strtobool, default='False')
            
    opts = parser.parse_args()
    
    print opts
        
    VariantCalling_wf = VariantCalling(
        opts.sample,
        opts.fastq,
        opts.out_dir,
        opts.ref,
        opts.config_file,
        opts.stages,
        
    )
    
    sys.exit(VariantCalling_wf.run(mode='local',dataDirRoot=opts.out_dir,isContinue=opts.isContinue,isDryRun=opts.isDryRun))