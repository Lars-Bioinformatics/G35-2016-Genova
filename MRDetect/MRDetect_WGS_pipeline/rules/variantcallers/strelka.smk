#########################################################
####                     Strelka                     ####
#########################################################

# Initial operations
strelka_output_dir = output_raw+"strelka/"
manta_output_dir = output_raw+"manta/"
# onstart:
#     shell("mkdir -p "+strelka_output_dir)

rule Manta_config_and_run:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        wf=manta_output_dir+"{tumor}_vs_{normal}/runWorkflow.py",
        out=manta_output_dir+"{tumor}_vs_{normal}/results/variants/candidateSmallIndels.vcf.gz",
    params:
        manta_outdir=manta_output_dir+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=240000,
        cpus=threadsparallel
    threads:
        threadsparallel
    shell:
        """
        conda run --name manta \
        configManta.py \
        --normalBam {input.normal} \
        --tumorBam {input.tumor} \
        --referenceFasta {ref} \
        --callRegions {callable_regions} \
        --runDir {params.manta_outdir}

        conda run --name manta \
        {output.wf} \
        -m local \
        -j {resources.cpus}
        
        gunzip -k {output.out}
        """

rule Strelka_config_and_run:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix,
        manta=manta_output_dir+"{tumor}_vs_{normal}/results/variants/candidateSmallIndels.vcf.gz",
    output:
        wfs=strelka_output_dir+"{tumor}_vs_{normal}/runWorkflow.py",
        outsnvs=strelka_output_dir+"{tumor}_vs_{normal}/results/variants/somatic.snvs.vcf.gz",
        outindels=strelka_output_dir+"{tumor}_vs_{normal}/results/variants/somatic.indels.vcf.gz"
    params:
        outdir=strelka_output_dir+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=240000,
        cpus=threadsparallel
    threads:
        threadsparallel
    shell:
        """
        conda run --name strelka \
        configureStrelkaSomaticWorkflow.py \
        --normalBam {input.normal} \
        --tumorBam {input.tumor} \
        --referenceFasta {ref} \
        --indelCandidates {input.manta} \
        --callRegions {callable_regions} \
        --runDir {params.outdir}

        conda run --name strelka \
        {output.wfs} \
        -m local \
        -j {resources.cpus}
        
        gunzip -k {output.outsnvs} {output.outindels}
        """

checkpoint Strelka_copy:
    input:
        snvs=strelka_output_dir+"{tumor}_vs_{normal}/results/variants/somatic.snvs.vcf.gz",
        indels=strelka_output_dir+"{tumor}_vs_{normal}/results/variants/somatic.indels.vcf.gz"
    output:
        snvs=output_variants_notfiltered+"{tumor}_vs_{normal}/strelka_snvs.vcf",
        indels=output_variants_notfiltered+"{tumor}_vs_{normal}/strelka_indels.vcf"
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
        gzip -c -d {input.snvs} > {output.snvs}
        gzip -c -d {input.indels} > {output.indels}
        """