library(tidyverse)
library(data.table)
# library(xlsx)
library(parallel)

setwd("/work/sduvarcall/G35-2016-Genova/MarkDuplicates/varscan_somatic_joint_calling_default_all-samples/varscan2_maf_1000g-pon_gnomadAF_above0_0001-filtered")

read_sequenza <- function(path, ploidy_path="", manual_ploidy=F){
    cnv_files <- list.files(pattern = "segments", path = path)
    data <- map_dfr(
        .x = cnv_files,
        .f = function(file){
            sample = str_split(file, "_", simplify = T)[2]
            print(sample)
            data = read.table(paste0(path,"/",file), header = T) %>% dplyr::select(chr=chromosome, start=start.pos, stop=end.pos, cnv_tumor_total=CNt, cnv_tumor_major=A, cnv_tumor_minor=B)
            data$sample = sample
            data$sample_num = which(cnv_files==file)
            
            # Add ploidy
            if (manual_ploidy || ploidy_path==""){
                data = manually_calculate_ploidy(data)
            } else {
                ploidy_file = list.files(path = ploidy_path, pattern = paste0(sample,"_"), full.names = T)
                stats_data =  read.table(ploidy_file, header = T)
                ploidy = stats_data[2,2]
                if (is.na(ploidy)){
                    data = manually_calculate_ploidy(data)
                } else {
                    ploidy = as.numeric(ploidy)
                    data$ploidy = ploidy
                    data$ploidy_rounded = round(ploidy)
                }
                
            }
            data
        }
    )
    
    return(as.data.table(data))
}

combine_snv_cnv_data <- function(snv_data, cnv_data){
    
    snv_data = snv_data %>% dplyr::rename(sample=Tumor_Sample_Barcode, chr=Chromosome, start=Start_Position, end=End_Position)
    cnv_data = cnv_data %>% dplyr::rename(start2=start, end2=stop) #%>% as.data.table()
    
    # Set keys used for finding overlaps
    setkey(snv_data, sample, chr, start, end)
    setkey(cnv_data, sample, chr, start2, end2)
    
    combined = foverlaps(snv_data,cnv_data)
    
    return(combined)
}

merge_tumor_plasma <- function(dt1, dt2, sample2){
    dt_plasma = fread(matched_plasma)
    dt_plasma = dt_plasma[, list(chr=Chromosome, start=Start_Position, end=End_Position, #sample=Tumor_Sample_Barcode, 
                                t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count)]
    dt_plasma = dt_plasma[, ':='(t_vaf=t_alt_count/t_depth, n_vaf=n_alt_count/n_depth)]

    col_names = names(dt_plasma)[-c(1:3)]
    merged = dt1[dt_plasma, on = c("chr","start","end"), paste0(sample2,"_",col_names):=mget(paste0("i.", col_names))]
    return(merged)
}

ploidy_path="/work/sduvarcall/G35-2016-Genova/MarkDuplicates/sequenza_grouped_data/sequenza_cellularity_ploidy"
segments_path="/work/sduvarcall/G35-2016-Genova/MarkDuplicates/sequenza_grouped_data/sequenza_segments"
cnv_data = read_sequenza(segments_path, ploidy_path)




cl <- makeCluster(detectCores()-1)
clusterExport(cl, varlist = c("cnv_data","combine_snv_cnv_data","merge_tumor_plasma"))

tumors = list.files(pattern="-TUMOR-")
#parApply(cl, tumors, 1, function(tumor){
for (tumor in tumors){

    # Get sample and patient names
    sample=str_split(tumor, "_")[[1]][1]
    print(sample)
    patient=paste(str_split(sample,"-")[[1]][1:2], collapse="-")

    # Read tumor sample
    dt = fread(tumor)

    dt2 = combine_snv_cnv_data(dt, cnv_data)

    # Extract relevant columns
    dt3 = dt2[, list(sample, chr, start, end, Gene=Hugo_Symbol, Variant_Classification, Variant_Type, Ref=Reference_Allele, Alt=Allele,
                dbSNP_RS, cnv_tumor_total, cnv_tumor_major, cnv_tumor_minor, ploidy, 
                t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count)]
    # Filter variants based on read counts
    dt3 = dt3[n_alt_count==0 & n_ref_count>5 & t_alt_count>5, ]
    # Add t_vaf and n_vaf 
    dt3 = dt3[, ':='(t_vaf=t_alt_count/t_depth, n_vaf=n_alt_count/n_depth)]

    # Read matched plasma sample data
    plasma_sample = paste0(patient,"-PLASMA")
    matched_plasma = list.files(pattern=plasma_sample)
    dt_merged = merge_tumor_plasma(dt3, matched_plasma, plasma_sample)

    # Merge remaining plasma samples with data table
    remaining_plasma_samples = setdiff(list.files(pattern="PLASMA"), list.files(pattern=paste0(patient,"-PLASMA")))
    for (plasma_sample_file in remaining_plasma_samples){
        sample2=str_split(plasma_sample_file, "_")[[1]][1]
        dt_merged = merge_tumor_plasma(dt_merged, plasma_sample_file, sample2)
    }

    out.dir = "per_patient_tumor_and_all_plasma_samples/"
    dir.create(out.dir, showWarnings=F)
    # Write tsv file
    write.table(dt_merged, file=paste0(out.dir,sample,"_with_all-plasma-samples-1000g_pon-gnomadAF_above0_0001-tADmin5-nRDmin5-nADmax0.txt"), 
                row.names=F, quote=F, sep="\t")
    # Write excel file
    #write.xlsx2(dt_merged, file=paste0(out.dir,sample,"_with_all-plasma-samples-1000g_pon-gnomadAF_above0_0001-tADmin5-nRDmin5-nADmax0.xlsx"), sheetName = "Sheet1",
    #            col.names = TRUE, row.names = TRUE, append = FALSE)
#})
}