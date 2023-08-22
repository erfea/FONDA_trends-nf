process GLS {
  debug true

  tag { "${aoi_aoipar}" }
  publishDir "${params.outdata}", mode:'copy'
  container "lewinska/rdoc:rparts"

  cpus 20
  // memory '40.GB'

  input:
  tuple val(aoi), val(aoipar), val(aoi_aoipar), path(in_path)
  path( glsCode )

  output:
  tuple val("${aoi}"), val("${aoipar}"), val("${aoi_aoipar}"), path("GLS/${aoi_aoipar}")

  """
  mkdir -p GLS/${aoi_aoipar}
  Rscript ${glsCode} ${in_path} ${aoi} GLS/${aoi_aoipar}
  """
  // docker run -v /data:/data -v /mnt:/mnt -u $(id -u):$(id -g) -it --rm --cpuset-cpus="0-19" --name rPARTS lewinska/rdoc:rparts Rscript ${arCode} ${in_path

}
