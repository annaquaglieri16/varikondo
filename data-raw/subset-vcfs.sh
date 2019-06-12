# write code used to produce the vcf files.
# download from GEO, align, pre-process, call, annotate.


module load bcftools

# Freebayes
bgzip -c ../inst/extdata/germline_freebayes.vcf > \
../inst/extdata/germline_freebayes.vcf.gz
tabix -p vcf ../inst/extdata/germline_freebayes.vcf.gz

bcftools filter ../inst/extdata/germline_freebayes.vcf.gz -r chr20 > \
../inst/extdata/chr20_freebayes.vcf

bgzip -c ../inst/extdata/chr20_freebayes.vcf > ../inst/extdata/chr20_freebayes.vcf.gz
tabix -p vcf ../inst/extdata/chr20_freebayes.vcf.gz


# GATK3 mutect2
bgzip -c ../inst/extdata/germline_mutect.vcf > ../inst/extdata/germline_mutect.vcf.gz
tabix -p vcf ../inst/extdata/germline_mutect.vcf.gz
bcftools filter ../inst/extdata/germline_mutect.vcf.gz -r chr20 > \
../inst/extdata/chr20_mutect.vcf

bgzip -c ../inst/extdata/chr20_mutect.vcf > ../inst/extdata/chr20_mutect.vcf.gz
tabix -p vcf ../inst/extdata/chr20_mutect.vcf.gz

# vardict
bgzip -c ../inst/extdata/germline_vardict.vcf > ../inst/extdata/germline_vardict.vcf.gz
tabix -p vcf ../inst/extdata/germline_vardict.vcf.gz
bcftools filter ../inst/extdata/germline_vardict.vcf.gz -r chr20 > \
../inst/extdata/chr20_vardict.vcf

bgzip -c ../inst/extdata/chr20_vardict.vcf > ../inst/extdata/chr20_vardict.vcf.gz
tabix -p vcf ../inst/extdata/chr20_vardict.vcf.gz


# VarScan
bgzip -c ../inst/extdata/germline_varscan.vcf > ../inst/extdata/germline_varscan.vcf.gz
tabix -p vcf ../inst/extdata/germline_varscan.vcf.gz

bcftools filter ../inst/extdata/germline_varscan.vcf.gz -r chr20 > \
../inst/extdata/chr20_varscan.vcf

bgzip -c ../inst/extdata/chr20_varscan.vcf > ../inst/extdata/chr20_varscan.vcf.gz
tabix -p vcf ../inst/extdata/chr20_varscan.vcf.gz


##
rm -rf germline_*

