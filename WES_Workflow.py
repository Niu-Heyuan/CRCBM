# 定义GATK可执行文件路径
GATK=~/miniforge3/bin/gatk

# 定义参考基因组和目标区域变量
GENOME=/mnt/hpc/home/liugang/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
interval=/mnt/hpc/home/liugang/references/BED/Twist_human_exome_2.0_hg38.interval_list
pon=/mnt/hpc/home/liugang/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
gnomad=/mnt/hpc/home/liugang/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz 
small_vcf=/mnt/hpc/home/liugang/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
bed=/mnt/hpc/home/liugang/references/BED/Twist_human_exome_2.0_hg38.interval_list
VEP_DIR="/mnt/hpc/home/liugang/.vep/"
# 修改样本名称
id="CRC1"
normal="CRC1N"
normalID="CRC1_N"
tumorID="CRC1"
echo "Running Sarek..."
# fq原始数据到bam，加上输出数据指控
nextflow run /mnt/hpc/home/liugang/nf-core-sarek_3.2.2/3_2_2/ --wes --input /mnt/hpc/home/liugang/DATA/fastq/${id}/${id}.csv --genome GATK.GRCh38 -profile singularity --outdir /mnt/hpc/home/liugang/DATA/result/${id} --igenomes_base /mnt/hpc/home/liugang/references/ --tools cnvkit -c /mnt/hpc/home/liugang/nf-core-sarek_3.2.2/configs/custom.config --save_mapped --save_output_as_bam

/mnt/hpc/home/liugang/DATA/result/${id}/preprocessing/recalibrated/${id}/${id}.recal.bam

echo "Running Mutect2..."

# Mutect2 变异检测
$GATK Mutect2 \
    -R ${GENOME} \
    -I /mnt/hpc/home/liugang/DATA/result/${id}/preprocessing/recalibrated/${id}/${id}.recal.bam \
    -I /mnt/hpc/home/liugang/DATA/result/${id}/preprocessing/recalibrated/${id}/${normal}.recal.bam \
    -normal ${normalID} \
    -L ${interval} \
    --panel-of-normals ${pon} \
    --germline-resource ${gnomad} \
    --genotype-germline-sites true \
    --genotype-pon-sites true \
    --f1r2-tar-gz /mnt/hpc/home/liugang/bamTOvcf/median/${id}_flr2.tar.gz \
    --native-pair-hmm-threads 10 \
    -O /mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect2.vcf.gz \
    1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1
# 检查 Mutect2 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect2.vcf.gz" ]; then
    echo "Mutect2 output file not found!"
    exit 1
fi

echo "Running LearnReadOrientationModel..."
# LearnReadOrientationModel
$GATK LearnReadOrientationModel \
    -I /mnt/hpc/home/liugang/bamTOvcf/median/${id}_flr2.tar.gz \
    -O /mnt/hpc/home/liugang/bamTOvcf/median/${id}_read-orientation-model.tar.gz \
    1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1

# 检查 LearnReadOrientationModel 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/median/${id}_read-orientation-model.tar.gz" ]; then
    echo "LearnReadOrientationModel output file not found!"
    exit 1
fi

echo "Running GetPileupSummaries..."
# GetPileupSummaries
$GATK GetPileupSummaries \
    -I /mnt/hpc/home/liugang/bamTOvcf/data/${id}.recal.bam \
    -V ${small_vcf} \
    -L ${bed} \
    -O /mnt/hpc/home/liugang/bamTOvcf/median/${id}_getpileupsummaries.table \
    1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1

# 检查 GetPileupSummaries 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/median/${id}_getpileupsummaries.table" ]; then
    echo "GetPileupSummaries output file not found!"
    exit 1
fi

echo "Running CalculateContamination..."
# CalculateContamination
$GATK CalculateContamination \
    -I /mnt/hpc/home/liugang/bamTOvcf/median/${id}_getpileupsummaries.table \
    -tumor-segmentation /mnt/hpc/home/liugang/bamTOvcf/median/${id}_segments.table \
    -O /mnt/hpc/home/liugang/bamTOvcf/median/${id}_contamination.table \
    1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1


# 检查 CalculateContamination 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/median/${id}_contamination.table" ]; then
    echo "CalculateContamination output file not found!"
    exit 1
fi

echo "Running FilterMutectCalls..."
# FilterMutectCalls
$GATK FilterMutectCalls \
    -R ${GENOME} \
    -V /mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect2.vcf.gz \
    -O /mnt/hpc/home/liugang/bamTOvcf/median/${id}_somatic.vcf.gz \
    --tumor-segmentation /mnt/hpc/home/liugang/bamTOvcf/median/${id}_segments.table \
    --contamination-table /mnt/hpc/home/liugang/bamTOvcf/median/${id}_contamination.table \
    --ob-priors /mnt/hpc/home/liugang/bamTOvcf/median/${id}_read-orientation-model.tar.gz \
    1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1

# 检查 FilterMutectCalls 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/median/${id}_somatic.vcf.gz" ]; then
    echo "FilterMutectCalls output file not found!"
    exit 1
fi

echo "Running SelectVariants..."
# SelectVariants
$GATK SelectVariants \
    -R ${GENOME} \
    -V /mnt/hpc/home/liugang/bamTOvcf/median/${id}_somatic.vcf.gz \
    --select "vc.isNotFiltered()" \
    -O /mnt/hpc/home/liugang/bamTOvcf/result/${id}_filter.vcf \
    1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1

# 检查 SelectVariants 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/result/${id}_filter.vcf" ]; then
    echo "SelectVariants output file not found!"
    exit 1
fi

echo "Running VEP annotation..."

# VEP 注释
vep --everything --offline --cache --species homo_sapiens --dir ${VEP_DIR} --fasta ${GENOME} --stats_text --stats_file /mnt/hpc/home/liugang/bamTOvcf/VEP/${id}_filter.VEP.html --vcf -i /mnt/hpc/home/liugang/bamTOvcf/result/${id}_filter.vcf  -o /mnt/hpc/home/liugang/bamTOvcf/VEP/${id}_filter.VEP.vcf --force_overwrite
         1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1

# 检查 VEP 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/VEP/${id}_filter.VEP.vcf" ]; then
    echo "VEP annotation output file not found!"
    exit 1
fi

echo "Running vcf2maf..."

#vcf转maf
vcf2maf.pl --ref-fasta ${GENOME} --vep-path /mnt/hpc/home/liugang/maf/VEP/vcf2maf/ensembl-vep-release-108.2 --ncbi-build GRCh38 --tumor-id ${tumorID} --normal-id ${normalID} --inhibit-vep --input-vcf /mnt/hpc/home/liugang/bamTOvcf/VEP/${id}_filter.VEP.vcf --output-maf /mnt/hpc/home/liugang/bamTOvcf/maf/${id}_filter.VEP.maf
1>/mnt/hpc/home/liugang/bamTOvcf/median/${id}_mutect.log 2>&1

# 检查 vcf2maf 输出文件
if [ ! -f "/mnt/hpc/home/liugang/bamTOvcf/maf/${id}_filter.VEP.maf" ]; then
    echo "vcf2maf output file not found!"
    exit 1
fi

echo "Script completed."

#OncodriveCLUST驱动突变
for g in CRCBM_BM CRCBM_P CRC ;do 
oncodriveclustl -i OncodriveCLUSTL/${g}/${g}_mut.tsv \
-r OncodriveCLUSTL/hg38_cds.tsv \
-g hg38  -c 10 --seed 123 \
-o ./OncodriveCLUSTL/${g} 1>./OncodriveCLUSTL/${g}_OncodriveCLUSTL.log 2>&1
done
##############################

echo "========================================"
echo "Process end at:"
date

# rm -rf $CURDIR/nodelist.$SLURM_JOB_ID

#OncodriveFM寻找驱动突变
for g in CRCBM_BM CRCBM_P CRC ;do 
oncodriveclustl -i OncodriveCLUSTL/${g}/${g}_mut.tsv \
-r OncodriveCLUSTL/hg38_cds.tsv \
-g hg38  -c 10 --seed 123 \
-o ./OncodriveCLUSTL/${g} 1>./OncodriveCLUSTL/${g}_OncodriveCLUSTL.log 2>&1
done
##############################

echo "========================================"
echo "Process end at:"
date

# rm -rf $CURDIR/nodelist.$SLURM_JOB_ID


#dNdScv寻找驱动突变
library(devtools)
data("dataset_simbreast", package="dndscv")
dndsout = dndscv(mutations)
head(mutations)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)
sel_cv = dndsskin$sel_cv
print(sel_cv[sel_cv$qglobal_cv<0.1,c(1:10,19)], digits = 3)

#MutSig2CV寻找驱动突变
cd ~/biosoft
mkdir MatlabMCR
cd MatlabMCR
wget -c http://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip
unzip MCR_R2016a_glnxa64_installer.zip
./install

LD_LIBRARY_PATH="${HOME}/biosoft/MATLAB/v901/runtime/glnxa64:${HOME}/biosoft/MATLAB/v901/bin/glnxa64:${HOME}/biosoft/MATLAB/v901/sys/os/glnxa64"

cd ~/biosoft
wget -c https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSigCV_1.41.zip
unzip MutSigCV_1.41.zip
cd MutSigCV_1.41

wget -c http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/gene.covariates.txt
wget -c http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/exome_full192.coverage.zip
wget -c http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/mutation_type_dictionary_file.txt
wget -c http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/chr_files_hg19.zip
wget -c https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSigCV_example_data.1.0.1.zip
wget -c https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/LUSC.MutSigCV.input.data.v1.0.zip

GENOME=/teach/database/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
mkdir chr_files_hg38
cat ${GENOME}.fai | cut -f 1  | while read id;
do 
samtools faidx ${GENOME} ${id} > chr_files_hg38/${id}.txt ; 
sed -i '1d' chr_files_hg38/${id}.txt ; 
sed -i ':a;N;$!ba;s/\n//g' chr_files_hg38/${id}.txt ;
done
