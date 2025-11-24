#计算grm
gcta-1.94.1 --bfile geno --autosome --maf 0.01 --make-grm --out test --thread-num 50


# 指定输入文件和输出目录
input_file="/phenotype/m841_ivn_用于计算遗传力.txt"
output_directory="hsnp/metab_ivn+covars"

# 循环处理每个表型（从第三列开始）
for ((col=3; col<=843; col++)); do
    #
    phenotype=$(awk -v col=$col 'NR==1 {print $col}' "$input_file")
    # 
    mpheno=$((col - 2))
    #
    gcta-1.94.1 --grm grm/test --pheno "$input_file" --reml --mpheno "$mpheno" --covar "phenotype/GCTA_covar.txt" --qcovar "phenotype/GCTA_qcovar.txt" --out "$output_directory/$phenotype" --thread-num 50
done








