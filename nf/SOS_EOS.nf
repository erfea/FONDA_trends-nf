process SOS_EOS {
  // debug true

  tag { "(${aoi}-${aoipar})" }
  publishDir "${params.outdata}/", mode:'copy'

  // tag { "$tile ${sraData[1].baseName}" }
  conda "/data/Jakku/users/lewinska/miniconda3/envs/envpy39/"

  input:
  // tuple path(sosParent), path(sospath)
  tuple val(sraData), val(aoi), val(aoipar), val(aoi_aoipar), path(sospath)
  path(soseosCode)

  output:
  tuple val("${aoi}"), val("${aoipar}"), val("${aoi_aoipar}"),
  // path( "output/${aoi}_${aoipar}/${sraData}/$tile/*" )
  path( "${sospath}/")

  """
   python ${soseosCode} ${sospath}

  """

  //  # mkdir -p output/${aoi}_${aoipar}/${sraData}/$tile/*
  // # conda activate envpy39

}
