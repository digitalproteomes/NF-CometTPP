nextflow.enable.dsl=2

include {tpp_main;
	 tpp_summaries;
	 tpp_exports;
	 run_proteinprophet;
	 tpp_peter;
	 apply_progenesis_patch	 
} from './nfmodule_tpp/tpp_workflows.nf'

workflow {
    main:
    log.info("++++++++++========================================")
    log.info("Executing CometTPP workflow")
    log.info("")
    log.info("Parameters:")
    log.info(" DDA folder:\t $params.dda_folder")
    log.info(" Protein DB:\t $params.protein_db")
    log.info(" Decoy tag:\t $params.decoy")
    log.info(" TPP params:\t $params.tpp")
    if(params.no_pool) {
	log.info("  All searches analyzed separately in TPP")
    }
    else {
	log.info("  All searches merged for TPP analysis")
    }
    log.info(" Comet threads:\t $params.comet_threads")
    log.info(" xinteract threads:\t $params.xinteract_threads")
    log.info("++++++++++========================================")
    
    tpp_main(params.dda_folder,
	     params.comet_params,
	     params.protein_db,
	     params.comet_threads,
	     params.tpp,
	     params.decoy,
	     params.no_pool,
	     params.xinteract_threads)
    
    tpp_summaries(tpp_main.out.pepxmls,
		  params.comet_params,
		  params.protein_db,
		  params.decoy)

    if(params.no_pool) {
	// StPeter only makes sense if did not pool pep.xmls
	//
	// Channels in nextflow are FIFO and tpp_main.out.pepxmls and
	// tpp_main.out.protxmls are loaded with the corresponding pep
	// and prot xml files in the corresponding positions

	tpp_peter(tpp_main.out.pepxmls,
		  tpp_main.out.protxmls,
		  channel.fromPath("$params.dda_folder/*").collect())
    }

    if(params.export_to_tsv) {
	tpp_exports(tpp_main.out.pepxmls,
		    tpp_main.out.pepxmlmodels,
		    tpp_main.out.protxmls,
		    tpp_main.out.protxmlmodels)
    }

    if(params.patch_for_progenesis) {
	apply_progenesis_patch(tpp_main.out.pepxmls)
    }
}
