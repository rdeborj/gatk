from gatk.gatk import Gatk
import yaml

if __name__ == '__main__':
    gatk_obj = Gatk(output_dir='PROCESSED_BAMS')
    gatk_split = gatk_obj.split_n_cigar_reads(
        bams=['1.bam'],
        output='split.bam')
    print(gatk_split)
    knownsites=[
        '/cluster/tools/data/genomes/human/hg38/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf',
        '/cluster/tools/data/genomes/human/hg38/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf',
        '/cluster/tools/data/genomes/human/hg38/hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz']
    gatk_realign_target = gatk_obj.realigner_target_creator(
        bams=[gatk_split['output']],
        output='split.intervals',
        knownsites=knownsites)
    print(gatk_realign_target)
    gatk_indel_realign = gatk_obj.indel_realigner(
        bams=[gatk_split['output']],
        output='split.realigned.bam',
        intervals=gatk_realign_target['output'])
    print(gatk_indel_realign)
    knownsites=[
        '/cluster/tools/data/genomes/human/hg38/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf',
        '/cluster/tools/data/genomes/human/hg38/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf',
        '/cluster/tools/data/genomes/human/hg38/hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz']
    gatk_basequality_recal = gatk_obj.base_quality_recalibration(
        bam=gatk_indel_realign['output'],
        output='split.realigned.recal_data.table',
        knownsites=knownsites)
    print(gatk_basequality_recal)
    gatk_print_reads = gatk_obj.print_reads(
        bam = gatk_indel_realign['output'],
        output='split.realigned.recal.bam',
        recaldata=gatk_basequality_recal['output'])
    print(gatk_print_reads)

    # run the variant discovery suite of applications
    gatk_obj.output_dir = 'HAPLOTYPE_CALLER'
    dbsnp = '/cluster/tools/data/genomes/human/hg38/hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz'
    gatk_haplotype_caller = gatk_obj.haplotype_caller(
        bams=[gatk_print_reads['output']],
        dbsnp=dbsnp,
        output='output.raw.snps.indels.vcf')
    print(gatk_haplotype_caller)
    gatk_variant_filtration = gatk_obj.variant_filtration(
        variant=gatk_haplotype_caller['output'],
        output='output.filtered.snps.indels.vcf',
        filters={"FS":">30.0", "QD":"<2.0"})
    print(gatk_variant_filtration)

