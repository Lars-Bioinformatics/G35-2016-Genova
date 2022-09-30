__title__ = "vcf2maf"
__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "02/06/2021"
__version__ = "1.0"

import time

#########################################################
####                       Input                     ####
#########################################################
# Matched tumour-normal samples information
configfile: "../samples-pairwise.yaml"
# configfile: "somatic_matched_samples_oneSample.yaml"

# Explicit paths for external input files
ref = "/work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta"
vep_data = "/work/sduvarcall/ensembl_vep/vep"
vep_path = "/work/miniconda3/envs/vep/bin"


# ref_build = "GRCh37"
ref_build = "GRCh38"


#########################################################
####                      Output                     ####
#########################################################
# output_mutect2 = "/work/Data/Connor/mutect2_vcf/"
output_mutect2 = "/work/Data/Connor/mutect2_joint_calling_try3/vcf_filterFlag/"
output_mutect2_pass = "/work/Data/Connor/mutect2_joint_calling_try3/vcf_filterFlag_PASS/"
output_varscan2 = "varscan_somatic_joint_calling_default/"
output_varscan2_all = "varscan_somatic_joint_calling_default_all-samples/"
# output_varscan2 = "/work/Data/Connor/varscan_somatic_joint_calling_q30/"

#########################################################
####               Sample Information                ####
#########################################################
# Mutect2 sample information
PAIR = [pair for pair in config]

# Varscan2 sample information
VARSCAN_SAMPLES = glob_wildcards(output_varscan2+"{sample}_cns_varscan2.vcf.gz")
print(VARSCAN_SAMPLES)

#########################################################
####                       Setup                     ####
#########################################################
# Timing
totim = time.time()
timeFormat = "%Y_%m_%d:%X" # year, month, day, time H:M:S

# Memory
mem = "-Xmx12g"

#########################################################
####  Define workflow start, stop and error actions  ####
#########################################################
onstart:
    # shell("mkdir -p "+output_mutect2)
    shell("mkdir -p "+output_mutect2+"maf_files/")
    shell("mkdir -p "+output_mutect2+"vcf_annotated/")
    # shell("mkdir -p "+output_mutect2+"mutect2_maf/")
    # shell("mkdir -p "+output_mutect2+"mutect2_maf_PASSonly/")
    # shell("mkdir -p "+output_varscan2+"vcf_annotated/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_1000g-pon-filtered/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomadAF_above0_001-filtered/")
    # shell("mkdir -p "+output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomadAF_above0_0001-filtered/")
    # shell("mkdir -p "+output_varscan2_all+"vcf_annotated/")
    # shell("mkdir -p "+output_varscan2_all+"varscan2_maf/")
    # shell("mkdir -p "+output_varscan2_all+"varscan2_maf_1000g-pon-filtered/")
    # shell("mkdir -p "+output_varscan2_all+"varscan2_maf_1000g-pon_gnomadAF_above0_001-filtered/")
    # shell("mkdir -p "+output_varscan2_all+"varscan2_maf_1000g-pon_gnomadAF_above0_0001-filtered/")



#########################################################
####                  Run All Rules                  ####
#########################################################
'''
Rule all
'''
rule all_pairs:
    input:
        # Multi-sample vcf
        # [expand(output_mutect2+"maf_files/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_mutect2_pass+"maf_files/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASS_norm_fromMultiSampleVcf.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_mutect2+"mutect2_maf_PASSonly/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # Varscan2 multisample per patient
        # [expand(output_varscan2+"varscan2_maf_somatic_1000g-pon-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.maf", 
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        [expand(output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomadAF_above0_001-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_001-filtered.maf", 
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR],
        [expand(output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomadAF_above0_0001-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_0001-filtered.maf", 
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_varscan2+"varscan2_maf/{tumor}_vs_{normal}_cns_varscan2.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        # Varscan2 all samples
        # [expand(output_varscan2_all+"varscan2_maf_1000g-pon-filtered/{tumor}_vs_{normal}_cns_varscan2_1000g-pon-filtered_allSamplesVcf.maf", 
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],
        [expand(output_varscan2_all+"varscan2_maf_1000g-pon_gnomadAF_above0_001-filtered/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_001-filtered_allSamplesVcf.maf", 
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR],
        [expand(output_varscan2_all+"varscan2_maf_1000g-pon_gnomadAF_above0_0001-filtered/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_0001-filtered_allSamplesVcf.maf", 
            normal=config[fam]["normal"],
            tumor=config[fam]["tumor"]) for fam in PAIR],
        # [expand(output_varscan2_all+"varscan2_maf/{tumor}_vs_{normal}_cns_varscan2_allSamplesVcf.maf",
        #     normal=config[fam]["normal"],
        #     tumor=config[fam]["tumor"]) for fam in PAIR],


#########################################################
####                  Run All Rules                  ####
#########################################################
rule vcf2maf_mutect2:
    input:
        vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz"
    output:
        maf=output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.maf",
        vep_vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vep.vcf",
        vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.vcf")
    #params:
    #    vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf")
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {wildcards.tumor} \
        --normal-id {wildcards.normal} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_mutect2_PASS:
    input:
        vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz"
    output:
        maf=output_mutect2+"maf_files/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.maf",
        vep_vcf=output_mutect2+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vep.vcf",
        vcf=temp(output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.vcf")
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {wildcards.tumor} \
        --normal-id {wildcards.normal} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_mutect2_PASS_multiSample:
    input:
        vcf=lambda wildcards: expand(output_mutect2_pass+"{sample}_somatic_mutect2_filterFlag_PASSonly.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
    output:
        maf=output_mutect2_pass+"maf_files/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASS_norm_fromMultiSampleVcf.maf",
        vep_vcf=output_mutect2_pass+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASS_norm_fromMultiSampleVcf.vep.vcf",
        vcf=temp(output_mutect2_pass+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASS_norm_fromMultiSampleVcf.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_mutect2_multiSample:
    input:
        vcf=lambda wildcards: expand(output_mutect2+"{sample}_somatic_mutect2_filterFlag.vcf.gz", sample=wildcards.normal[:5]+wildcards.normal[7:])
    output:
        maf=output_mutect2+"maf_files/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.maf",
        vep_vcf=output_mutect2+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.vep.vcf",
        vcf=temp(output_mutect2+"vcf_annotated/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_fromMultiSampleVcf.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """



rule vcf2maf_varscan2:
    input:
        vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2.vcf.gz", sample=wildcards.tumor.replace("PLASMA","TUMOR"))
    output:
        maf=output_varscan2+"varscan2_maf/{tumor}_vs_{normal}_cns_varscan2.maf",
        vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2.vep.vcf",
        vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_varscan2_somtic_pon:
    input:
        vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2_somatic_1000g-pon-filtered.vcf.gz", sample=wildcards.tumor.replace("PLASMA","TUMOR"))
    output:
        maf=output_varscan2+"varscan2_maf_somatic_1000g-pon-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.maf",
        vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.vep.vcf",
        vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon-filtered.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_varscan2_somtic_pon_gnomadAF_above0_001:
    input:
        vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_001-filtered.vcf.gz", sample=wildcards.tumor.replace("PLASMA","TUMOR"))
    output:
        maf=output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomadAF_above0_001-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_001-filtered.maf",
        vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_001-filtered.vep.vcf",
        vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_001-filtered.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_varscan2_somtic_pon_gnomadAF_above0_0001:
    input:
        vcf=lambda wildcards: expand(output_varscan2+"{sample}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_0001-filtered.vcf.gz", sample=wildcards.tumor.replace("PLASMA","TUMOR"))
    output:
        maf=output_varscan2+"varscan2_maf_somatic_1000g-pon_gnomadAF_above0_0001-filtered/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_0001-filtered.maf",
        vep_vcf=output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_0001-filtered.vep.vcf",
        vcf=temp(output_varscan2+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_somatic_1000g-pon_gnomadAF_above0_0001-filtered.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

# rule filter_PASS:
#     input
#         vcf=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag.vcf.gz",
#         maf=output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm.maf",
#     output:
#         vcf_pass=output_mutect2+"{tumor}_vs_{normal}_somatic_mutect2_filterFlag_PASSonly.vcf.gz",
#         maf_pass=output_mutect2+"mutect2_maf/{tumor}_vs_{normal}_somatic_mutect2_filterFlag_norm_PASSonly.maf",
#     shell:
#         """
#         zcat {input.vcf} | egrep '^#|PASS' | bgzip -c > {output.vcf_pass}
#         tabix -p vcf {output.vcf_pass}

#         cat {input.maf} | egrep '^#|^HUGO|PASS' | {output.maf_pass}
#         """


# ALL SAMPLES
rule vcf2maf_varscan2_all_samples:
    input:
        vcf=output_varscan2_all+"G35_all_samples_cns_varscan2.vcf.gz"
    output:
        maf=output_varscan2_all+"varscan2_maf/{tumor}_vs_{normal}_cns_varscan2_allSamplesVcf.maf",
        vep_vcf=output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2.vep.vcf",
        vcf=temp(output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_varscan2_somtic_pon_all_samples:
    input:
        vcf=output_varscan2_all+"G35_all_samples_cns_varscan2_1000g-pon-filtered.vcf.gz"
    output:
        maf=output_varscan2_all+"varscan2_maf_1000g-pon-filtered/{tumor}_vs_{normal}_cns_varscan2_1000g-pon-filtered_allSamplesVcf.maf",
        vep_vcf=output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_1000g-pon-filtered.vep.vcf",
        vcf=temp(output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_1000g-pon-filtered.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_varscan2_somtic_pon_gnomadAF_above0_001_all_samples:
    input:
        vcf=output_varscan2_all+"G35_all_samples_cns_varscan2_1000g-pon_gnomadAF_above0_001-filtered.vcf.gz"
    output:
        maf=output_varscan2_all+"varscan2_maf_1000g-pon_gnomadAF_above0_001-filtered/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_001-filtered_allSamplesVcf.maf",
        vep_vcf=output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_001-filtered.vep.vcf",
        vcf=temp(output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_001-filtered.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """

rule vcf2maf_varscan2_somtic_pon_gnomadAF_above0_0001_all_samples:
    input:
        vcf=output_varscan2_all+"G35_all_samples_cns_varscan2_1000g-pon_gnomadAF_above0_0001-filtered.vcf.gz"
    output:
        maf=output_varscan2_all+"varscan2_maf_1000g-pon_gnomadAF_above0_0001-filtered/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_0001-filtered_allSamplesVcf.maf",
        vep_vcf=output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_0001-filtered.vep.vcf",
        vcf=temp(output_varscan2_all+"vcf_annotated/{tumor}_vs_{normal}_cns_varscan2_1000g-pon_gnomadAF_above0_0001-filtered.vcf")
    params:
        new_tum_id = lambda wildcards: wildcards.tumor.split("_")[0],
        new_norm_id = lambda wildcards: str(wildcards.normal).split("_")[0]
    resources: cpus=4, mem=24000 # Since internal command of vcf2maf use fork = 4
    shell:
        """
        bgzip -d -c {input.vcf} > {output.vcf}
        vcf2maf.pl \
        --input-vcf {output.vcf} \
        --output-maf {output.maf} \
        --vep-data {vep_data} \
        --vep-path {vep_path} \
        --ncbi-build {ref_build} \
        --ref-fasta {ref} \
        --tumor-id {params.new_tum_id} \
        --normal-id {params.new_norm_id} \
        --vcf-tumor-id {wildcards.tumor} \
        --vcf-normal-id {wildcards.normal}
        """