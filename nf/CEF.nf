process CEF {
  // debug true

  tag { "${aoi_aoipar}" }
  // publishDir "${params.outdata}/output/cef/${aoi_aoipar}", mode:'copy'

  conda "/data/Jakku/users/lewinska/miniconda3/envs/envpy39/"

  input:
  tuple val(sraData), val(aoi), val(aoiprm), val(aoi_aoipar), path(in_path)
  path(cefCode)

  output:
  tuple val(sraData), val(aoi), val(aoiprm), val(aoi_aoipar), path("${in_path}/cef")

  """
  mkdir -p ${in_path}/cef
  mkdir -p ${params.outdata}/output/cef
  # conda activate envpy39gdal
  python ${cefCode} ${in_path}

  """
}
