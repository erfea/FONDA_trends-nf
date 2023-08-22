
process SMA_RBF {


    tag { "$tile ${sraData[1].baseName}" }
    container 'davidfrantz/force'

    cpus 2
    memory '20.GB'

    input:
    tuple val(tile), path("level2_norm/$tile/*"), path("mask/$tile/*")
    each sraData
    path( "vectors/*" )
    path( "endm/*" )
    path( "level2_norm/datacube-definition.prj" )

    output:
    tuple val(tile), val("sraData[0]"), path( "output/${sraData[0]}/$tile/*" )

    """
    mkdir -p output/${sraData[0]} prov
    force-higher-level ${sraData[1]}
    """

}

process WriteSMA {

    input:
    val( name )
    val( numberOfEndmembers )

    output:
    tuple val(name), path( "config.prm" )

    """
    touch config.prm
    echo "++PARAM_TRAIN_START++" >> config.prm
    echo "DIR_HIGHER = output/$name"
    echo "SMA_ENDMEMBER = $numberOfEndmembers"
    """


}



workflow {

    //normalizedLevel2 = Channel.fromPath( "${params.level2norm}/X????_Y????/*" ).map{ [it.parent.name, it] }.groupTuple().map{ [it[0],it[1],file("${params.masks}/${it[0]}/${params.maskName}")] }
    normalizedLevel2 = Channel.fromPath( "${params.level2norm}/X0058_Y0087/*" ).map{ [it.parent.name, it] }.groupTuple().map{ [it[0],it[1],file("${params.masks}/${it[0]}/${params.maskName}")] }
    SRAParamFiles = Channel.of( 
        ["gv",file(params.sraGV)]
    )
    //SRAParamFiles = WriteSMA( 
    //    Channel.of(
    //        ["gv",1]
    //    )
    //).out
    SATiles = file( params.SATiles )
    endmember = file( params.endmembers )
    datacube = file( params.datacube )
    SRAParamFiles.view()
    SMA_RBF( normalizedLevel2, SRAParamFiles, SATiles, endmember, datacube )
    
    timeSeries = SMA_RBF.out.filter{ it[1] in ["gv","npv"] }
    rbf = SMA_RBF.out.filter{ it[1] in ["gv","npv","soil","shade"] }


    WriteSMA( 
        Channel.of(
            ["gv",1]
        )
    )


}