# GWASFormatter
Multithreaded Java tool to generate .ped and .map files required for GWAS 

> java gwas/GWASFormatter2 --help

Usage of GWASFormatter:

java GWASFormatter2 varScan_file patient_phenotype_file -hd/-nd threshold numThreads

Input file formats:

varScan.txt has the general tab delimetted varScan layout with patient alphaNumeric (PA[0-9]{n} type): chrom	position	ref	var	normal_reads1	
normal_reads2	normal_var_freq	normal_gt	tumor_reads1	tumor_reads2	tumor_var_freq	tumor_gt	
somatic_status	variant_p_value	somatic_p_value	tumor_reads1_plus	tumor_reads1_minus	tumor_reads2_plus	
tumor_reads2_minus	normal_reads1_plus	normal_reads1_minus	normal_reads2_plus	normal_reads2_minus	patient

patient_phenotype.csv is comma separated with the general layout of: Patientfactor,PHENOTYPE,Genderfactor,Locationfactor,
Age,patientID,Gender,Location,BRAF.Alteration,Gene.Fusion

-hd to include column names in output; -nd for no header in output files

threshold specifies snp genotype and should be numeric (could be float)

numThreads indicates the number of threads requested

Output: varScanFileName.ped and varScanFileName.map files in the same directory. 

Existing files will be overwritten...

V2.7 Developed by Fadi Hariri; email: fadi.hariri@mail.mcgill.ca


