process FNF {
  // debug true

  tag { "${aoi_aoipar}" }
  publishDir "${params.outdata}/output/${aoi_aoipar}", mode:'copy'

  conda "/data/Jakku/users/lewinska/miniconda3/envs/envpy39gdal/"

  input:
  tuple val(sraData), val(aoi), val(aoiprm), val(aoi_aoipar), path(in_path)
  // path(path)
  path(fnfCode)

  output:
  tuple val(sraData), val("${aoi[0]}"), val("${aoiprm[0]}"), val("${aoi_aoipar}"), path("${sraData[0]}")

  // path("${fnfpath[0]}/**")// , val("sraData[0]"), val(tile)
  // fnfpath[0] because we only want a single path
  // tuple val(fnfpath), path("${fnfpath[0]}/")

  """
  # conda activate envpy39gdal
  python ${fnfCode} "."

  """
}
