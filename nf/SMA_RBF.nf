process SMA_RBF {

    // debug true

    // tag { "$tile ${sraData[1].baseName} (${sraData[2]}-${sraData[3]})" }
    tag {"${sraData[0]} ${sraData[1].baseName} (${sraData[2]}-${sraData[3]})" }
    publishDir "${params.outdata}/", mode:'copy'
    container 'davidfrantz/force'
    // , saveAs: {"${sraData[2]}"_"${sraData[3]}/${sraData[0]}/${tile}/$it"},

    cpus 10
    memory '40.GB'

    input:
    // tuple val(tile), path("level2_norm/$tile/*"), path("mask/$tile/*")
    // tuple val(tile), path("level2_norm/*/*"), path("mask/*/*")
    path( "*" )
    path( "*" ) 
    each sraData
    path( "*" )
    path( "*" )
    // path( "level2_norm/datacube-definition.prj" )

    output:
    // tuple val(tile), val("$sraData[0]"), path( "output/${sraData[2]}_${sraData[3]}/${sraData[0]}/$tile/*" ), val("$sraData[2]"), val("$sraData[3]")
    tuple val("${sraData[0]}"), val("${sraData[2]}"), val("${sraData[3]}"), val("${sraData[2]}_${sraData[3]}"), path( "output/${sraData[2]}_${sraData[3]}/${sraData[0]}/*/*" )
    """
    mkdir -p output/${sraData[2]}_${sraData[3]}/${sraData[0]}/ prov 
    force-higher-level ${sraData[1]}
    """

// # for tl in ${tile}; do mkdir -p output/${sraData[2]}_${sraData[3]}/${sraData[0]}/$tl prov ; done
//     # mkdir -p output/${sraData[2]}_${sraData[3]}/${sraData[0]}/prov

}
