#!/bin/bash
#===========================================================
# Merged script: process grouped VCF files
# Usage: bash this_script.sh <pop> <pop1> <pop2>
#   pop  : population name used for splitting by chromosome / PLINK conversion
#   pop1 : first population for XP-EHH
#   pop2 : second population for XP-EHH
#===========================================================

pop=$1
pop1=$2
pop2=$3

# ------------------------------------------------------------
# Step 1: Split VCF by chromosome (chr1..chr29)
# ------------------------------------------------------------
echo "Step 1: Splitting ${pop}.recode.vcf by chromosome..."
for i in {1..29}; do
    vcftools --vcf ${pop}.recode.vcf \
             --recode \
             --recode-INFO-all \
             --chr ${i} \
             --out ${pop}.chr${i}
done

# ------------------------------------------------------------
# Step 2: Convert per-chromosome VCF to PLINK format
# ------------------------------------------------------------
echo "Step 2: Converting per-chromosome VCF to PLINK format..."
for k in {1..29}; do
    vcftools --vcf ${pop}.chr${k}.recode.vcf \
             --plink \
             --out chr${k}.${pop}
done

echo "All files processed."

# ------------------------------------------------------------
# Step 3: Build genetic distance map files (position in Mb)
# ------------------------------------------------------------
echo "Step 3: Building genetic distance map files..."
for k in {1..29}; do
    awk 'BEGIN {OFS=" "} {print 1,".",$4/1000000,$4}' chr${k}.${pop}.map > chr${k}.${pop}.map.distance
done

# ------------------------------------------------------------
# Step 4: Run selscan XP-EHH
# ------------------------------------------------------------
echo "Step 4: Running selscan XP-EHH..."
for k in {1..29}; do
    /home/soft/selscan/src/selscan --xpehh \
                                   --vcf ${pop1}.chr${k}.phased.vcf \
                                   --vcf-ref ${pop2}.chr${k}.phased.vcf \
                                   --map chr${k}.${pop1}.map.distance \
                                   --threads 10 \
                                   --out chr${k}.${pop1}_vs_${pop2}
done

echo "All steps completed."