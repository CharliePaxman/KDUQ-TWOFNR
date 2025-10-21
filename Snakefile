# Snakefile

# Define the final targets
rule all:
    input:
        "FRONT",
        "outputs/PotentialsKDUQ_n_A19_Z8_E16000.pot",
        "outputs/PotentialsKDUQ_p_A19_Z8_E16000.pot"

# Rule to build TWOFNR from twofnr20.f
rule build_twofnr:
    input:
        "twofnr20.f"
    output:
        "TWOFNR"
    shell:
        "gfortran -O2 -finit-local-zero {input} -o {output}"

# Rule to build FRONT from front21.f
rule build_front:
    input:
        "front21.f"
    output:
        "FRONT"
    shell:
        "gfortran -O2 -finit-local-zero {input} -o {output}"
       
# Rule to run KDUQ
rule run_KDUQ:
    output:
        n_pot="outputs/PotentialsKDUQ_n_A{A}_Z{Z}_E{EkeV}.pot",
        p_pot="outputs/PotentialsKDUQ_p_A{A}_Z{Z}_E{EkeV}.pot"
    input:
        script="scripts/SampleKDUQ.py"
    params:
        A=lambda wildcards: wildcards.A,
        Z=lambda wildcards: wildcards.Z,
        EkeV=lambda wildcards: wildcards.EkeV,
        EMeV=lambda wildcards: float(wildcards.EkeV) / 1000.0
    shell:
        "python3 {input.script} -A {params.A} -Z {params.Z} -E {params.EMeV}"
