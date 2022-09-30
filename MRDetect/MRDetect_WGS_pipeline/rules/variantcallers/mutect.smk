#########################################################
####                     MuTect2                     ####
#########################################################

# Initial operations
mutect2_output_dir = output_raw+"mutect2/"
# mutect2_output_somatic = output_raw+"mutect2/somatic_variants/"
# onstart:
#     shell("mkdir -p "+mutect2_output_dir)
#     shell("mkdir -p "+mutect2_output_somatic)

# mutect2_matched_output = mutect2_output_dir
# mutect2_matched_output = main_output+"mutect2_matched/"
# onstart:
    # shell("mkdir -p "+mutect2_matched_output+"split/{{vcf,bam,f1r2}}_split")
    # shell("mkdir -p "+mutect2_matched_output+"vcf_{{files,filterFlag,filterFlag_PASS,filterFlag_norm}}")
    # shell("mkdir -p "+mutect2_matched_output+"realigned_bam/")

CHROM = ["chr"+str(i) for i in range(1,23)] + ["chrX","chrY"]

# Call Somatic Variants using Mutect2 on matched Tumor-Normal samples on per chromesome basis
rule Mutect_matched:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        vcf=mutect2_output_dir+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{chr}.vcf.gz",
        idx=mutect2_output_dir+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{chr}.vcf.gz.tbi",
        vcf_stats=mutect2_output_dir+"split/vcf_split/{tumor}_vs_{normal}_somatic_mutect2__{chr}.vcf.gz.stats",
        # bam_subfile=mutect2_output_dir+"split/bam_split/{tumor}_vs_{normal}__{chr}_mutect2.bam",
        f1r2=mutect2_output_dir+"split/f1r2_split/{tumor}_vs_{normal}_f1r2__{chr}.tar.gz"
    resources: cpus=3, mem_mb=18000
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} Mutect2 \
        -R {ref} \
        -I {input.tumor} \
        -I {input.normal} \
        -normal {wildcards.normal} \
        -pon {pon_1000g} \
        --germline-resource {gnomad} \
        --af-of-alleles-not-in-resource 0.0000025 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --native-pair-hmm-threads {resources.cpus} \
        --f1r2-tar-gz {output.f1r2} \
        -L {wildcards.chr} \
        -O {output.vcf}
        """
        # -bamout {output.bam_subfile} \ ## Debugging
        # --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
        # --af-of-alleles-not-in-resource 0.0000025 \
        # --max_alt_alleles_in_normal_count 1000000 \ # NOT IN GATK 4.1 (or 4.2)
        # --max_alt_allele_in_normal_fraction 0.10 \ # NOT IN GATK 4.1 (or 4.2)
        # --tumor-lod-to-emit  \

rule merge_somatic_vcf:
    input:
        vcf_subfile=expand(mutect2_output_dir+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{chr}.vcf.gz", chr=CHROM)
    output:
        vcf=mutect2_output_dir+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
        idx=mutect2_output_dir+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
    params:
        vcf_subfile=expand("-I "+mutect2_output_dir+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{chr}.vcf.gz", chr=CHROM)
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} GatherVcfs \
        {params.vcf_subfile} \
        -O {output.vcf}

        tabix -p vcf {output.vcf}
        """

rule merge_somatic_vcf_stats:
    input:
        vcf_subfile=expand(mutect2_output_dir+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{chr}.vcf.gz.stats", chr=CHROM)
    output:
        vcf_stats=mutect2_output_dir+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
    params:
        vcf_subfile=expand("-stats "+mutect2_output_dir+"split/vcf_split/{{tumor}}_vs_{{normal}}_somatic_mutect2__{chr}.vcf.gz.stats", chr=CHROM)
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} MergeMutectStats \
        {params.vcf_subfile} \
        -O {output.vcf_stats}
        """

# Learn Read Orientation Bias
rule learnReadOrientationModel:
    input:
        f1r2=expand(mutect2_output_dir+"split/f1r2_split/{{tumor}}_vs_{{normal}}_f1r2__{chr}.tar.gz", chr=CHROM)
    output:
        f1r2_model=mutect2_output_dir+"read_orientation_models/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
    params:
        f1r2=expand("-I "+mutect2_output_dir+"split/f1r2_split/{{tumor}}_vs_{{normal}}_f1r2__{chr}.tar.gz", chr=CHROM)
    resources:
        mem_mb=6000,
        cpus=1
    shell:
        """
        conda run --name gatk4 \
        gatk LearnReadOrientationModel \
        {params.f1r2} \
        -O {output}
        """


# Create Contamination table
rule Mutect_GetPileupSummaries_normal:
    input:
        bam=input_dir+"{normal}"+bamsuffix
    output:
        pileup=mutect2_output_dir+"contamination/{normal}_normal_pileup.table"
    resources:
        mem_mb=60000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule Mutect_GetPileupSummaries_tumor:
    input:
        bam=input_dir+"{tumor}"+bamsuffix
    output:
        pileup=mutect2_output_dir+"contamination/{tumor}_tumor_pileup.table"
    resources:
        mem_mb=150000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule Mutect_CalculateContamination:
    input:
        normal=mutect2_output_dir+"contamination/{normal}_normal_pileup.table",
        tumor=mutect2_output_dir+"contamination/{tumor}_tumor_pileup.table"
    output:
        contamination=mutect2_output_dir+"contamination/{tumor}_vs_{normal}_contamination.table"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """


# Filter Mutect2 Calls
rule FilterMutectCalls:
    input:
        vcf=mutect2_output_dir+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz",
        idx=mutect2_output_dir+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.tbi",
        stats=mutect2_output_dir+"vcf_files/{tumor}_vs_{normal}_somatic_mutect2.vcf.gz.stats",
        contamination=mutect2_output_dir+"contamination/{tumor}_vs_{normal}_contamination.table",
        read_orientation=mutect2_output_dir+"read_orientation_models/{tumor}_vs_{normal}_read-orientation-model.tar.gz"
    output:
        #vcf=mutect2_output_dir+"{tumor}_vs_{normal}_somatic_mutect2_filtered.vcf.gz",
        vcf=mutect2_output_dir+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
        #idx=mutect2_output_dir+"{tumor}_vs_{normal}_somatic_mutect2_filtered.vcf.gz.tbi"
    resources:
        mem_mb=12000,
        cpus=1
    shell:
        """
        conda run --name gatk4 \
        gatk --java-options {mem} FilterMutectCalls \
        -R {ref} \
        -V {input.vcf} \
        --contamination-table {input.contamination} \
        --orientation-bias-artifact-priors {input.read_orientation} \
        --stats {input.stats} \
        -O {output.vcf}
        """

checkpoint Mutect_copy:
    input:
        bothcombined=mutect2_output_dir+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        bothcombined=output_variants_notfiltered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf"
    params:
        folder=output_variants_notfiltered+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        mkdir -p {params.folder}
        gzip -c -d {input.bothcombined} > {output.bothcombined}
        # cp {input.bothcombined} {output.bothcombined}
        """