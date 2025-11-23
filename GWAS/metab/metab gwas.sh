plink2 --bfile geno_rename --keep sampleid.txt --remove het_rmlist_for_plink.txt --maf 0.01 --hwe 1e-6 --geno 0.05 --make-bed --out geno

plink2 --bfile geno --pheno phenotype/m841_ivn_res_ivn_na.txt --no-psam-pheno --glm allow-no-covars --threads 60 --out plinkresult/m841_na/plink2

plink2 --bfile geno --pheno phenotype/ratios_input.txt --no-psam-pheno --glm allow-no-covars --threads 60 --out plinkresult/ratios_na/plink2
