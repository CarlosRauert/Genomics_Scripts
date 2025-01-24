# Try test profile
cd NF_Jabba
nextflow run . -profile test 

# Does not work due to incompatible sample sheet. Ergo get some TCGA fastqs and make proper sample sheet.
cd /data/cephfs-1/home/users/rauertc_c/scratch/NF_Jabba/Bams/AB36
export GDC_TOKEN=/data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer/gdc-user-token.2025-01-22T15_32_13.891Z.txt
gdc-client download -m /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/NextFlow_Jabba_Runs/20250122_Test_Manifest_AB36.txt -t /data/cephfs-1/home/users/rauertc_c/work/GDC_Datatransfer/gdc-user-token.2025-01-22T15_32_13.891Z.txt

find "$(pwd)" -type f # for getting filepaths 
# Smaple sheet is at ~/work/Scripts_Git_Repos/Genomics_Scripts/NextFlow_Jabba_Runs/20250123_Test_samplesheet.csv

,fragcounter,dryclean,cbs,hetpileups,ascat,jabba

cp ~/work/Scripts_Git_Repos/Genomics_Scripts/NextFlow_Jabba_Runs/20250123_Test_samplesheet.csv .

nextflow run ~/work/Scripts_Git_Repos/Genomics_Scripts/NF_Jabba \
        -profile bih \
        --input './20250123_Test_samplesheet.csv' \
        --outdir './out' \
        --tools svaba,fragcounter,dryclean,cbs,hetpileups,ascat,jabba \
        --genome GRCh38 \
        --step sv_calling

# Skips most steps after GATK4. Warnings: 

WARN: Access to undefined parameter `monochromeLogs` -- Initialise it to a default value eg. `params.monochromeLogs = some_value`
WARN: Access to undefined parameter `germline_resource` -- Initialise it to a default value eg. `params.germline_resource = some_value`
WARN: Access to undefined parameter `germline_resource_tbi` -- Initialise it to a default value eg. `params.germline_resource_tbi = some_value`
WARN: Access to undefined parameter `mappability` -- Initialise it to a default value eg. `params.mappability = some_value`
WARN: Access to undefined parameter `pon` -- Initialise it to a default value eg. `params.pon = some_value`
WARN: Access to undefined parameter `pon_tbi` -- Initialise it to a default value eg. `params.pon_tbi = some_value`
WARN: Access to undefined parameter `validationSkipDuplicateCheck` -- Initialise it to a default value eg. `params.validationSkipDuplicateCheck = some_value`
WARN: Access to undefined parameter `validationS3PathCheck` -- Initialise it to a default value eg. `params.validationS3PathCheck = some_value`
WARN: Access to undefined parameter `snpeff_db` -- Initialise it to a default value eg. `params.snpeff_db = some_value`
WARN: Access to undefined parameter `vep_cache_version` -- Initialise it to a default value eg. `params.vep_cache_version = some_value`
WARN: Access to undefined parameter `vep_genome` -- Initialise it to a default value eg. `params.vep_genome = some_value`
WARN: Access to undefined parameter `vep_species` -- Initialise it to a default value eg. `params.vep_species = some_value`
WARN: Access to undefined parameter `nonintegral_jabba` -- Initialise it to a default value eg. `params.nonintegral_jabba = some_value`
WARN: Access to undefined parameter `help_jabba` -- Initialise it to a default value eg. `params.help_jabba = some_value`
WARN: There's no process matching config selector: STRELKA.*|MANTA.*
WARN: There's no process matching config selector: GATK4_MERGEVCFS
WARN: There's no process matching config selector: GRIDSS
WARN: There's no process matching config selector: SAMPLESHEET_CHECK


# Add parameter definitions in nextflow config file

# Make params file and define parameters.

-params-file /data/cephfs-1/home/users/rauertc_c/work/Scripts_Git_Repos/Genomics_Scripts/NextFlow_Jabba_Runs/20250125_Test_Params.yaml