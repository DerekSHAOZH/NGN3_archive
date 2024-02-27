#!/usr/bin/env bash

bash PrepareData/prepareData.sh

bash IsoformMS/getIsoform.sh

bash RBPs/run.sh

mkdir -p Result/GO

if [ ! -f Result/GO/Figure_EV2B.pdf ] || [ ! -f Result/GO/Figure_EV2C.pdf ] || [ ! -f Result/GO/Figure_3G.pdf ];then
	Rscript GO/plotGO.r Result/GO Data/GOfiles
fi
