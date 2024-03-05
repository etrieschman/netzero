from snakemake.utils import min_version
from os import path
from pathlib import Path
min_version('6.0')

# --------------------------- Workflow constants --------------------------- #
# configfile: "config.yaml"

PATH_DATA = './data/'
PATH_PROCESSED = PATH_DATA + 'processed/'
PATH_INTERIM = PATH_DATA + 'interim/'
PATH_RESOURCES = PATH_DATA + 'resources/'
PATH_RESULTS = './results/'
YEAR_START = 2006
YEAR_END = 2021

# --------------------------- Rules --------------------------- #

rule all:
    input:
        # TRANSFORM
        # PATH_PROCESSED + 'eia860_generator.parquet',
        # PATH_PROCESSED + 'eia860_plant.parquet',
        # PATH_PROCESSED + 'eia860_utility.parquet',
        # PATH_PROCESSED + 'eia860_ownership.parquet',
        # PATH_PROCESSED + 'eia923_ops.parquet',
        # PATH_PROCESSED + 'epa_facility.parquet',
        # PATH_PROCESSED + 'epa_emissions.parquet',
        # ALIGN
        outfile_generation = PATH_PROCESSED + 'df_generation.parquet',
        outfile_emissions = PATH_PROCESSED + 'df_emissions.parquet',
        # ENTITIES
        outfile_gen = PATH_PROCESSED + 'df_generators.parquet',
        outfile_plant = PATH_PROCESSED + 'df_plants.parquet',
        outfile_util = PATH_PROCESSED + 'df_utilities.parquet',
        outfile_own = PATH_PROCESSED + 'df_owners.parquet',
        # CDP FLAGS
        outfile_cdp_gen = PATH_PROCESSED + 'flag_cdp_generators.parquet',
        outfile_cdp_plant = PATH_PROCESSED + 'flag_cdp_plants.parquet',
        outfile_cdp_util = PATH_PROCESSED + 'flag_cdp_utilities.parquet',
        outfile_cdp_own = PATH_PROCESSED + 'flag_cdp_owners.parquet',

        
# ============ EXTRACT EIA ============
rule extract_eia:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END
    output:
        # directories
        eia_f860 = directory(PATH_DATA + 'raw/eia/' + f'f860/'),
        eia_f861 = directory(PATH_DATA + 'raw/eia/' + f'f861/'),
        eia_f923 = directory(PATH_DATA + 'raw/eia/' + f'f923/'),
        # flags
        flag_u = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/1___Utility_Y{YEAR_END}.xlsx',
        flat_p = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/2___Plant_Y{YEAR_END}.xlsx',
        flag_g = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/3_1_Generator_Y{YEAR_END}.xlsx',
        flag_o = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/4___Owner_Y{YEAR_END}.xlsx',
        flag_ops = PATH_DATA + 'raw/' + 'eia/' + f'f923/{YEAR_END}/EIA923_Schedules_2_3_4_5_M_12_{YEAR_END}_Final_Revision.xlsx'
    log:
        "logs/extract_eia.log"
    script:
        path.join('code', 'extract', 'eia.py')
        # Path('code/extract/eia.py')
        # 'code/extract/eia.py'

# ============ EXTRACT EPA ============
rule extract_epa:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END
    output:
        epa_emissions = directory(PATH_DATA + 'raw/epa/' + f'emissions/'),
        epa_facility = directory(PATH_DATA + 'raw/epa/' + f'facility/'),
        flag_fac = PATH_DATA + 'raw/' + 'epa/' + f'facility/{YEAR_END}/facility-{YEAR_END}.csv',
        flag_em = PATH_DATA + 'raw/' + 'epa/' + f'emissions/{YEAR_END}/emissions-daily-{YEAR_END}-wy.csv'
    log:
        "logs/extract_epa.log"
    script:
        path.join('code', 'extract', 'epa.py')

# ============ TRANSFORM EIA860 GENERATOR ============
rule transform_eia860_generator:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END,
        indir = PATH_DATA + 'raw/' + 'eia/f860/'
    input: 
        flag = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/3_1_Generator_Y{YEAR_END}.xlsx',
    output:
        intfile = PATH_INTERIM + 'eia860_generator.parquet',
        outfile = PATH_PROCESSED + 'eia860_generator.parquet'
    log:
        "logs/transform_eia860_generator.log"
    script:
        path.join('code', 'transform', 'eia860_generator.py')

# ============ TRANSFORM EIA860 PLANT ============
rule transform_eia860_plant:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END,
        indir = PATH_DATA + 'raw/' + 'eia/f860/'
    input: 
        flag = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/2___Plant_Y{YEAR_END}.xlsx',
    output:
        intfile = PATH_INTERIM + 'eia860_plant.parquet',
        outfile = PATH_PROCESSED + 'eia860_plant.parquet'
    log:
        "logs/transform_eia860_plant.log"
    script:
        path.join('code', 'transform', 'eia860_plant.py')


# ============ TRANSFORM EIA860 UTILITY ============
rule transform_eia860_utility:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END,
        indir = PATH_DATA + 'raw/' + 'eia/f860/'
    input: 
        flag = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/1___Utility_Y{YEAR_END}.xlsx',
    output:
        intfile = PATH_INTERIM + 'eia860_utility.parquet',
        outfile = PATH_PROCESSED + 'eia860_utility.parquet'
    log:
        "logs/transform_eia860_utility.log"
    script:
        path.join('code', 'transform', 'eia860_utility.py')


# ============ TRANSFORM EIA860 OWNERSHIP ============
rule transform_eia860_ownership:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END,
        indir = PATH_DATA + 'raw/' + 'eia/f860/'
    input: 
        flag = PATH_DATA + 'raw/' + 'eia/' + f'f860/{YEAR_END}/4___Owner_Y{YEAR_END}.xlsx',
    output:
        intfile = PATH_INTERIM + 'eia860_ownership.parquet',
        outfile = PATH_PROCESSED + 'eia860_ownership.parquet'
    log:
        "logs/transform_eia860_ownership.log"
    script:
        path.join('code', 'transform', 'eia860_ownership.py')


# ============ TRANSFORM EIA923 ============
rule transform_eia923_ops:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END,
        indir = PATH_DATA + 'raw/' + 'eia/f923/'
    input: 
        flag = PATH_DATA + 'raw/' + 'eia/' + f'f923/{YEAR_END}/EIA923_Schedules_2_3_4_5_M_12_{YEAR_END}_Final_Revision.xlsx'
    output:
        outfile = PATH_PROCESSED + 'eia923_ops.parquet'
    log:
        "logs/transform_eia923_ops.log"
    script:
        path.join('code', 'transform', 'eia923_ops.py')


# ============ TRANSFORM EPA ============
rule transform_epa:
    params:
        year_start = YEAR_START,
        year_end = YEAR_END,
        indir_facility = PATH_DATA + 'raw/' + 'epa/facility/',
        indir_emissions = PATH_DATA + 'raw/' + 'epa/emissions/'
    input:
        flag_fac = PATH_DATA + 'raw/' + 'epa/' + f'facility/{YEAR_END}/facility-{YEAR_END}.csv',
        flag_em = PATH_DATA + 'raw/' + 'epa/' + f'emissions/{YEAR_END}/emissions-daily-{YEAR_END}-wy.csv'        
    output:
        outfile_facility = PATH_PROCESSED + 'epa_facility.parquet',
        outfile_emissions = PATH_PROCESSED + 'epa_emissions.parquet'
    log:
        "logs/transform_epa.log"
    script:
        path.join('code', 'transform', 'epa.py')


# ============ ALIGN ENTITIES ============
rule align_entities:
    params:
        indir = PATH_PROCESSED,
        results_dir = PATH_RESULTS + 'align/entities/'
    input:
        infile_gen = PATH_PROCESSED + 'eia860_generator.parquet',
        infile_plant = PATH_PROCESSED + 'eia860_plant.parquet',
        infile_util = PATH_PROCESSED + 'eia860_utility.parquet',
        infile_own = PATH_PROCESSED + 'eia860_ownership.parquet',
    output:
        outfile_gen = PATH_PROCESSED + 'df_generators.parquet',
        outfile_plant = PATH_PROCESSED + 'df_plants.parquet',
        outfile_util = PATH_PROCESSED + 'df_utilities.parquet',
        outfile_own = PATH_PROCESSED + 'df_owners.parquet'
    log:
        'logs/align_entities.log'
    script:
        path.join('code', 'align', 'entities.py')

# ============ ALIGN GENERATION ============
rule align_generation:
    params:
        indir = PATH_PROCESSED,
        results_dir = PATH_RESULTS + 'align/generation/'
    input:
        infile_gen = PATH_PROCESSED + 'df_generators.parquet',
        infile_ops = PATH_PROCESSED + 'eia923_ops.parquet'
    output:
        outfile = PATH_PROCESSED + 'df_generation.parquet'
    log:
        'logs/align_generation.log'
    script:
        path.join('code', 'align', 'generation.py')


# ============ ALIGN EMISSIONS ============
rule align_emissions:
    params:
        indir = PATH_DATA,
        results_dir = PATH_RESULTS + 'align/emissions/'
    input:
        infile_gen = PATH_PROCESSED + 'df_generators.parquet',
        infile_ops = PATH_PROCESSED + 'df_generation.parquet',
        infile_epa_emissions = PATH_PROCESSED + 'epa_emissions.parquet',
        crosswalk = PATH_DATA + 'resources/epa_eia_crosswalk.csv',
        emission_factors = PATH_DATA + 'resources/se_emissions_factors.csv'
    output:
        outfile = PATH_PROCESSED + 'df_emissions.parquet'
    log:
        'logs/align_emissions.log'
    script:
        path.join('code', 'align', 'emissions.py')

# ============ FLAG CDP ============
rule cdp_flag:
    params:
        threshold = 90,
        indir = PATH_DATA,
        results_dir = PATH_RESULTS + 'align/cdp_flag/'
    input:
        infile_gen = PATH_PROCESSED + 'df_generators.parquet',
        infile_plant = PATH_PROCESSED + 'df_plants.parquet',
        infile_util = PATH_PROCESSED + 'df_utilities.parquet',
        infile_own = PATH_PROCESSED + 'df_owners.parquet',
        infile_csids = PATH_RESOURCES + 'out_parent_subsidiary_mapping.csv'
    output:
        intfile_util = PATH_INTERIM + 'map_utility_cdp.csv',
        intfile_own = PATH_INTERIM + 'map_owner_cdp.csv',
        outfile_util = PATH_PROCESSED + 'flag_cdp_utilities.parquet',
        outfile_own = PATH_PROCESSED + 'flag_cdp_owners.parquet',
        outfile_plant = PATH_PROCESSED + 'flag_cdp_plants.parquet',
        outfile_gen = PATH_PROCESSED + 'flag_cdp_generators.parquet'
    log:
        'logs/cdp_flag.log'
    script:
        path.join('code', 'align', 'cdp_flag.py')

# ============ FLAG CDP LLm ============
rule cdp_flag_llm:
    params:
        threshold = 90,
        indir = PATH_DATA,
        results_dir = PATH_RESULTS + 'align/cdp_flag/'
    input:
        infile_gen = PATH_PROCESSED + 'df_generators.parquet',
        infile_plant = PATH_PROCESSED + 'df_plants.parquet',
        infile_util = PATH_PROCESSED + 'df_utilities.parquet',
        infile_own = PATH_PROCESSED + 'df_owners.parquet',
        infile_csids = PATH_RESOURCES + 'out_parent_subsidiary_mapping.csv'
    output:
        intfile_util = PATH_INTERIM + 'map_utility_cdp_llm.csv',
        intfile_own = PATH_INTERIM + 'map_owner_cdp_llm.csv',
        outfile_util = PATH_PROCESSED + 'flag_cdp_utilities_llm.parquet',
        outfile_own = PATH_PROCESSED + 'flag_cdp_owners_llm.parquet',
        outfile_plant = PATH_PROCESSED + 'flag_cdp_plants_llm.parquet',
        outfile_gen = PATH_PROCESSED + 'flag_cdp_generators_llm.parquet'
    log:
        'logs/cdp_flag_llm.log'
    script:
        path.join('code', 'align', 'cdp_flag_llm.py')


