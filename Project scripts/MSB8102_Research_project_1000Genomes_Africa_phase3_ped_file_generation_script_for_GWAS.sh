#!/bin/bash

start=$SECONDS
echo $HOSTNAME

#downloading the 1000 genomes phase 3 data

wget https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1
wget https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
wget https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1

#rename the files

mv all_phase3.pgen.zst?dl=1 all_phase3.pgen.zst
mv all_phase3.pvar.zst?dl=1 all_phase3.pvar.zst
mv phase3_corrected.psam?dl=1 all_phase3.psam

#extract the .pgen files

unzstd all_phase3.pgen.zst 

extracting the African samples and converting to a bcf and then vcf

./plink2 --pfile all_phase3 vzs \
       --keep-cat-pheno SuperPop \
       --keep-cat-names AFR \
       --make-pgen \
       --out afr_phase3

mkdir 1000genomes_full_phase3_data
mv all_phase3* 1000genomes_full_phase3_data/

./plink2 --pfile afr_phase3 \
       --export bcf

#rename bcf file

mv plink2.bcf afr_phase3.bcf

#convert bcf to vcf

bcftools index afr_phase3.bcf

bcftools view -Oz -o afr_phase3.vcf.gz afr_phase3.bcf

#filtering out multi-allelic data

bcftools view -m2 -M2 -v snps afr_phase3.vcf.gz -o afr_phase3_edited.vcf.gz

#create plink files

./plink2 --vcf afr_phase3_edited.vcf.gz --make-bed --out afr_phase3

#generate ped file from plink files

#filter for MAF,LD and  relatedness

./plink --bfile afr_phase3 --maf 0.01 --make-bed --out afr_phase3_maf

./plink --bfile afr_phase3_maf --indep-pairwise 50 5 0.2

./king -b afr_phase3_maf.bed --unrelated --degree 2 --cpu 4 --prefix afr_phase3_maf_

mv afr_phase3_maf_unrelatedunrelated.txt afr_phase3_maf_unrelated.txt
mv afr_phase3_maf_unrelatedunrelated_toberemoved.txt afr_phase3_maf_unrelated_toberemoved.txt

./plink --bfile afr_phase3_maf --keep afr_phase3_maf_unrelated.txt --make-bed \
 --out afr_phase3_filtered

#split PLINK files by chromosomes for easier imputations in case of limited computational space

for chr in {1..22}; do \
./plink --bfile afr_phase3_filtered --chr $chr --make-bed --out afr_phase3_filtered${chr}; \
done

#generate ped file

for i in {1..22}; do \
./plink --bfile afr_phase3_filtered$i --recode --tab --out afr_phase3_$i \
done

end=$SECONDS
echo "duration: $((end-start)) seconds."
exit
