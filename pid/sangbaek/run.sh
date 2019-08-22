#!/usr/bin/bash

export pdir=`pwd`
export groovy=$COATJAVA"/bin/run-groovy"
data_path="/volatile/clas12/rg-a/production/recon/pass0/v15/"
dir=`ls $data_path`
#for run in $dir
#do
#       run=${run:2}
#       $groovy electron_pID.groovy $run `find $data_path"00"$run -name "*.hipo"`
#done
#JYPATH=$PWD $groovy electron_pID.groovy 5038 10 `find /volatile/clas12/rg-a/production/recon/pass0/v15/005038/ -name "*.hipo"`
#JYPATH=$PWD $groovy electron_pID.groovy 5038 100 `find /volatile/clas12/rg-a/production/recon/pass0/v15/005038/ -name "*.hipo"`
#JYPATH=$PWD $groovy electron_pID.groovy 5038 500 `find /volatile/clas12/rg-a/production/recon/pass0/v15/005038/ -name "*.hipo"`

JYPATH=$PWD $groovy electron_pID.groovy 5037 10 /work/clas12/sangbaek/generators/recon_dvcs.hipo
JYPATH=$PWD $groovy electron_pID.groovy 5037 100 /work/clas12/sangbaek/generators/recon_dvcs.hipo
JYPATH=$PWD $groovy electron_pID.groovy 5037 500 /work/clas12/sangbaek/generators/recon_dvcs.hipo


