
include { TSA_SMA_RBF_parameters } from './TSA_SMA_RBF_parameters'
include { SMA_RBF } from './SMA_RBF'
include { SOS_EOS } from './SOS_EOS'
include { FNF } from './FNF'
include { CEF } from './CEF'
include { AR } from './AR'
include { GLS } from './GLS'

workflow {

    // normalizedLevel2 = Channel.fromPath( "${params.level2norm}/X????_Y????/*" ).map{ [it.parent.name, it] }.groupTuple().map{ [it[0],it[1],file("${params.masks}/${it[0]}/${params.maskName}")] }
    // normalizedLevel2 = Channel.fromPath( "${params.level2norm}/X005?_Y008?/*" ).map{ [it.parent.name, it] }.groupTuple().map{ [it[0],it[1],file("${params.masks}/${it[0]}/${params.maskName}")] }
    normalizedLevel2 = Channel.fromPath( "${params.level2norm}" )//.map{ [it.parent.name, it] }.groupTuple().map{ [it[0],it[1]]}//.view()
    masksChannel = Channel.fromPath( "${params.masks}")//.map{ [it.parent.name, it] }//.view()
    // ${param.maskname}").map{ [it.name, it] }.groupTuple().map{ [it[0],it[1],file("${params.masks}/*/${params.maskName}")] }
    // cCh = normalizedLevel2.map{[it[0]]}.view()
    
    aoiChannel = Channel.of( "SA" )
    aoiparChannes = Channel.of("SA")
    parametersChannel = Channel.of(
          // AOI, Endmembers' AOI, endmember+variang, endmember's no, RMSE{TRUE,FALSE}, RBF: sigma1,
          // sigma2, sigma3, OUTPUT_TSI{TRUE,FALSE}, OUTPUT_SPL{TRUE,FALSE}, OUTPUT_LSP{TRUE,FALSE}
            [ "gv",         1,"TRUE",  8,  16, 32, "TRUE", "TRUE",  "TRUE" ],
            [ "npv",        2,"FALSE", 8,  16, 32, "TRUE", "FALSE", "FALSE"],
            [ "soil",       3,"FALSE", 8,  16, 32, "TRUE", "FALSE", "FALSE"],
            [ "shade",      4,"FALSE", 8,  16, 32, "TRUE", "FALSE", "FALSE"],
            [ "gv_wide",    1,"FALSE", 16, 48, 96, "TRUE", "FALSE", "FALSE"],
            [ "npv_wide",   2,"FALSE", 16, 48, 96, "TRUE", "FALSE", "FALSE"],
            [ "soil_wide",  3,"FALSE", 16, 48, 96, "TRUE", "FALSE", "FALSE"],
            [ "shade_wide", 4,"FALSE", 16, 48, 96, "TRUE", "FALSE", "FALSE"]
        )

    aoiCombinations = aoiChannel.combine(aoiparChannes)

    TSA_SMA_RBF_parameters_input = aoiCombinations.combine(parametersChannel)


    SRAParamFiles = TSA_SMA_RBF_parameters( TSA_SMA_RBF_parameters_input )//.view()

    // // SATiles = file( params.SATiles )
    // // endmembers = file( params.endmembers )

    Tiles = Channel.fromPath( "${params.DirTiles}" )//.view()
    endmembers = Channel.fromPath( "${params.DirEndm}")

    datacube = file( params.datacube )

    // // Tiles.view{ "Tile: $it" }
    // SRAParamFiles.view()

    SMA_RBF( normalizedLevel2, masksChannel, SRAParamFiles, Tiles, endmembers)//.view() //, datacube ).view()
    // timeSeries = SMA_RBF.out.filter{ it[1] in ["gv","npv"] }
    // rbf = SMA_RBF.out.filter{ it[0] in ["gv","npv","soil","shade","gv_wide","npv_wide","soil_wide","shade_wide"]}.view()
    //       .map(it.parent.parent.unique()).unique().collect().view()


    // soseos = SMA_RBF.out
    // .map{  [ it[0], it[1],  ( it[2] .findAll{it =~"S-LSP.tif"} ).parent.parent.unique() ] }
    // .findAll{ !it[2].isEmpty() } .view() 

    smaDataGV = SMA_RBF.out.filter{it[0] == 'gv'}.map{ [ it[0], it[1], it[2], it[3],
                                  (it[4].findAll{it =~"S-LSP.tif"}).parent.parent.unique()]}//.view()

    // smaData = SMA_RBF.out.map{ [it[0], it[1], it[2], it[3],
    //                               (it[4].findAll{it =~"S-LSP.tif"}).parent.parent.unique()]}//.view()
                                //   .branch{data: !it[5].isEmpty()
                                //           nodata: it[5].isEmpty()}
                                //   .set{phenChannel}
    // phenChannel.data.view()
    SOS_EOS(smaDataGV, file(params.soeosCode)).view()
    // SOS_EOS(smaData.data , file(params.soeosCode)).view()

    // // // // soseos = SMA_RBF.out.map{[ it[0], it[1], it[3], it[4],
    // // // //                         ( it[2].findAll{it =~"S-LSP.tif"}).parent.parent.unique() ] }
    // // // //                             .branch{data: !it[2].isEmpty()
    // // // //                                     nodata: it[2].isEmpty()}
    // // // //                                     .set{soseosin}
    // // // //                             soseosin.data.view()
    // // // //

    rbfChannel = SMA_RBF.out.map{[it[0], it[1], it[2], it[3], 
                                 (it[4].parent.parent.unique()).collect()]}
                                 .groupTuple(by: 3)
                                 .map{[it[0], it[1], it[2], it[3], it[4].flatten()]}
                                //  .view()


    FNF(rbfChannel, file(params.fnfCode))//.view()

    cefChannel = FNF.out.map{[it[0], it[1], it[2], it[3],
                             (it[4].parent)]}.view()

    CEF(cefChannel, file(params.cefCode))//.view()

    arChannel = CEF.out.map{[it[1], it[2], it[3], it[4]]}//.view()

    AR(arChannel, file(params.arCode))//.view()
    
    glsChannel = AR.out.map{[it[0], it[1], it[2], it[3]]}.view()

    GLS(glsChannel, params.glsCode).view()




// ## ARCHIVE ## \\

    // soseos = SMA_RBF.out.map{[it[0], it[1], ((it[2].findAll{it =~"S-LSP.tif"}).parent.parent.unique()).findAll{!it[2].isEmpty()}] }.view()
    // soseos = SMA_RBF.out.map{[it[0], it[1], ((it[2].findAll{it =~"S-LSP.tif"}).parent.parent.unique())] }.view()
    // soseos = SMA_RBF.out
    // .map{
    //    [ it[0], it[1], ( it[2]
    //                                       .findAll{it =~"S-LSP.tif"}
    //                                 ).parent.parent.unique()
    //    ] }
    //    . view{ it.size }
    // .view()

    // params.fnfCode.view{"fnfCode: $it" }


    // // ! does not work as negation O_o
    // // cefin1 = SMA_RBF.out.map{it[2].parent.parent.unique().findAll{!it.endsWith("gv")}}.collect().view()
    //
    // CEFin = (SMA_RBF.out.map{it[2].parent.parent.unique().findAll{!it.endsWith("gv")}})
    //         .join(SOS_EOS.out.map{it[0]}, remainder: true).collect().view()
    //         // .map{[!it[0].isEmpty()]}.view()
    //
    //       // merge(SMA_RBF.out.map{it[2].parent.parent.unique().findAll{it.endsWith("_npv")}}).view()
    //
    //                     // SMA_RBF.out.map{it[2].parent.parent.unique().findAll{it.endsWith("_soil")}},
    //                     // SMA_RBF.out.map{it[2].parent.parent.unique().findAll{it.endsWith("_shade")}},
    //                     //
    //
    // //
    // // // (CEFin.map{it[2]}).join(rbf).view()
    // // //cef = FnF.out.map{it[1].findAll{it =~"FnF"}}.view()
    // //
    // CEF(CEFin, file(params.cefCode)).view()
    //


}
