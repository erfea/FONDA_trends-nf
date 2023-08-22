
process TSA_SMA_RBF_parameters {

    // inputs:
    // AOI, Endmembers' AOI, endmember+variant, endmember's no, RMSE{TRUE,FALSE}, RBF sigma1, sigma2, sigma3, OUTPUT_TSI{TRUE,FALSE}, OUTPUT_SPL{TRUE,FALSE}, OUTPUT_LSP{TRUE,FALSE}
    input:
    tuple val(aoi), val(aoipar), val(name), val(endmember), val(rmse), val(sigma1), val(sigma2), val(sigma3), val(output_tsi), val(spl), val(lsp)

    output:
    tuple val(name), path( "${name}_config.prm" ), val(aoi), val(aoipar) 

    """

    touch ${name}_config.prm

    echo "++PARAM_TSA_START++" >> ${name}_config.prm

    echo "# INPUT/OUTPUT DIRECTORIES" >> ${name}_config.prm
    echo "DIR_LOWER = level2_norm" >> ${name}_config.prm
    echo "DIR_HIGHER = output/${aoi}_${aoipar}/$name" >> ${name}_config.prm
    echo "DIR_PROVENANCE = prov" >> ${name}_config.prm

    echo "# MASKING" >> ${name}_config.prm
    echo "DIR_MASK =  cubed" >> ${name}_config.prm
    echo "BASE_MASK = GRA_2018_10m_conv_noData.tif" >> ${name}_config.prm

    echo "# OUTPUT OPTIONS" >> ${name}_config.prm
    echo "OUTPUT_FORMAT = GTiff" >> ${name}_config.prm
    echo "FILE_OUTPUT_OPTIONS = NULL" >> ${name}_config.prm
    echo "OUTPUT_EXPLODE = FALSE" >> ${name}_config.prm

    echo "# PARALLEL PROCESSING" >> ${name}_config.prm
    echo "NTHREAD_READ = 15" >> ${name}_config.prm
    echo "NTHREAD_COMPUTE = 12" >> ${name}_config.prm
    echo "NTHREAD_WRITE = 4" >> ${name}_config.prm

    echo "# PROCESSING EXTENT AND RESOLUTION" >> ${name}_config.prm
    echo "X_TILE_RANGE = 15 108" >> ${name}_config.prm
    echo "Y_TILE_RANGE = 22 103" >> ${name}_config.prm
    echo "FILE_TILE = vectors/${aoi}_tiles.txt" >> ${name}_config.prm
    echo "BLOCK_SIZE = 0" >> ${name}_config.prm
    echo "RESOLUTION = 30" >> ${name}_config.prm
    echo "REDUCE_PSF = TRUE" >> ${name}_config.prm
    echo "USE_L2_IMPROPHE = FALSE" >> ${name}_config.prm

    echo "# SENSOR ALLOW-LIST" >> ${name}_config.prm
    echo "SENSORS = LND04 LND05 LND07 LND08 LND09 SEN2A SEN2B" >> ${name}_config.prm
    echo "SPECTRAL_ADJUST = FALSE" >> ${name}_config.prm

    echo "# QAI SCREENING" >> ${name}_config.prm
    echo "SCREEN_QAI = NODATA CLOUD_OPAQUE CLOUD_BUFFER CLOUD_CIRRUS CLOUD_SHADOW SNOW SUBZERO SATURATION" >> ${name}_config.prm
    echo "ABOVE_NOISE = 3" >> ${name}_config.prm
    echo "BELOW_NOISE = 1" >> ${name}_config.prm

    echo "# PROCESSING TIMEFRAME" >> ${name}_config.prm
    echo "DATE_RANGE = 1984-01-01 2022-12-31" >> ${name}_config.prm
    echo "DOY_RANGE = 1 365" >> ${name}_config.prm

    echo "# SPECTRAL INDEX" >> ${name}_config.prm
    echo "INDEX = SMA" >> ${name}_config.prm
    echo "STANDARDIZE_TSS = NONE" >> ${name}_config.prm
    echo "OUTPUT_TSS = TRUE" >> ${name}_config.prm

    echo "# SPECTRAL MIXTURE ANALYSIS" >> ${name}_config.prm
    echo "FILE_ENDMEM  = endm/${aoipar}_endmembers.txt" >> ${name}_config.prm
    echo "SMA_SUM_TO_ONE = TRUE" >> ${name}_config.prm
    echo "SMA_NON_NEG = TRUE" >> ${name}_config.prm
    echo "SMA_SHD_NORM = FALSE" >> ${name}_config.prm
    echo "SMA_ENDMEMBER = $endmember" >> ${name}_config.prm
    echo "OUTPUT_RMS = $rmse" >> ${name}_config.prm

    echo "# INTERPOLATION PARAMETERS" >> ${name}_config.prm
    echo "INTERPOLATE = RBF" >> ${name}_config.prm
    echo "MOVING_MAX = 16" >> ${name}_config.prm
    echo "RBF_SIGMA = $sigma1 $sigma2 $sigma3" >> ${name}_config.prm
    echo "RBF_CUTOFF = 0.95" >> ${name}_config.prm
    echo "HARMONIC_MODES = 3" >> ${name}_config.prm
    echo "HARMONIC_FIT_RANGE = 1984-01-01 2022-12-31" >> ${name}_config.prm
    echo "OUTPUT_NRT = FALSE" >> ${name}_config.prm
    echo "INT_DAY = 16" >> ${name}_config.prm
    echo "STANDARDIZE_TSI = NONE" >> ${name}_config.prm
    echo "OUTPUT_TSI = $output_tsi" >> ${name}_config.prm

    echo "# PYTHON UDF PARAMETERS" >> ${name}_config.prm
    echo "FILE_PYTHON = NULL" >> ${name}_config.prm
    echo "PYTHON_TYPE = PIXEL" >> ${name}_config.prm
    echo "OUTPUT_PYP = FALSE" >> ${name}_config.prm

    echo "# SPECTRAL TEMPORAL METRICS" >> ${name}_config.prm
    echo "OUTPUT_STM = FALSE" >> ${name}_config.prm
    echo "STM = Q25 Q50 Q75 AVG STD" >> ${name}_config.prm

    echo "# FOLDING PARAMETERS" >> ${name}_config.prm
    echo "FOLD_TYPE = AVG" >> ${name}_config.prm
    echo "STANDARDIZE_FOLD = NONE" >> ${name}_config.prm
    echo "OUTPUT_FBY = FALSE" >> ${name}_config.prm
    echo "OUTPUT_FBQ = FALSE" >> ${name}_config.prm
    echo "OUTPUT_FBM = FALSE" >> ${name}_config.prm
    echo "OUTPUT_FBW = FALSE" >> ${name}_config.prm
    echo "OUTPUT_FBD = FALSE" >> ${name}_config.prm
    echo "OUTPUT_TRY = FALSE" >> ${name}_config.prm
    echo "OUTPUT_TRQ = FALSE" >> ${name}_config.prm
    echo "OUTPUT_TRM = FALSE" >> ${name}_config.prm
    echo "OUTPUT_TRW = FALSE" >> ${name}_config.prm
    echo "OUTPUT_TRD = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAY = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAQ = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAM = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAW = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAD = FALSE" >> ${name}_config.prm

    echo "# LAND SURFACE PHENOLOGY PARAMETERS - SPLITS-BASED" >> ${name}_config.prm
    echo "LSP_DOY_PREV_YEAR = 273" >> ${name}_config.prm
    echo "LSP_DOY_NEXT_YEAR = 91" >> ${name}_config.prm
    echo "LSP_HEMISPHERE = NORTH" >> ${name}_config.prm
    echo "LSP_N_SEGMENT = 4" >> ${name}_config.prm
    echo "LSP_AMP_THRESHOLD = 0.2" >> ${name}_config.prm
    echo "LSP_MIN_VALUE = 500" >> ${name}_config.prm
    echo "LSP_MIN_AMPLITUDE = 500" >> ${name}_config.prm
    echo "LSP = DSS DES" >> ${name}_config.prm
    echo "STANDARDIZE_LSP = NONE" >> ${name}_config.prm
    echo "OUTPUT_SPL = $spl" >> ${name}_config.prm
    echo "OUTPUT_LSP = $lsp" >> ${name}_config.prm
    echo "OUTPUT_TRP = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAP = FALSE" >> ${name}_config.prm

    echo "# LAND SURFACE PHENOLOGY PARAMETERS - POLAR-BASED" >> ${name}_config.prm
    echo "POL_START_THRESHOLD = 0.2" >> ${name}_config.prm
    echo "POL_MID_THRESHOLD = 0.5" >> ${name}_config.prm
    echo "POL_END_THRESHOLD = 0.8" >> ${name}_config.prm
    echo "POL_ADAPTIVE = TRUE" >> ${name}_config.prm
    echo "POL = VSS VPS VES VSA RMR IGS" >> ${name}_config.prm
    echo "STANDARDIZE_POL = NONE" >> ${name}_config.prm
    echo "OUTPUT_PCT = FALSE" >> ${name}_config.prm
    echo "OUTPUT_POL = FALSE" >> ${name}_config.prm
    echo "OUTPUT_TRO = FALSE" >> ${name}_config.prm
    echo "OUTPUT_CAO = FALSE" >> ${name}_config.prm

    echo "# TREND PARAMETERS" >> ${name}_config.prm
    echo "TREND_TAIL = TWO" >> ${name}_config.prm
    echo "TREND_CONF = 0.95" >> ${name}_config.prm
    echo "CHANGE_PENALTY = FALSE" >> ${name}_config.prm

    echo "++PARAM_TSA_END++" >> ${name}_config.prm

    """


}
