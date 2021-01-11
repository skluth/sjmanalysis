#!/bin/bash

LD_LIBRARY_PATH=$PWD:/home/iwsatlas1/skluth/Downloads/root/root_v6.12.04/lib:/home/iwsatlas1/skluth/qcd/hepmc/hepmc2.06.09/install/lib
GFORTRAN_UNBUFFERED_PRECONNECTED=y
EXPORTLIST="LD_LIBRARY_PATH GFORTRAN_UNBUFFERED_PRECONNECTED"

( export $EXPORTLIST ; \
    ./runjob sjmconfig_91_all.cfg >& runsjm91_all.log ; ) &
( export $EXPORTLIST ; \
    ./runjob sjmconfig_130.cfg >& runsjm130.log ; \
    ./runjob sjmconfig_136.cfg >& runsjm136.log ; \
    ./runjob sjmconfig_161.cfg >& runsjm161.log ; \
    ./runjob sjmconfig_172.cfg >& runsjm172.log ; \
    ./runjob sjmconfig_183.cfg >& runsjm183.log ; ) &
( export $EXPORTLIST ; \
    ./runjob sjmconfig_189.cfg >& runsjm189.log ; \
    ./runjob sjmconfig_192.cfg >& runsjm192.log ; \
    ./runjob sjmconfig_196.cfg >& runsjm196.log ; \
    ./runjob sjmconfig_200.cfg >& runsjm200.log ; ) &
( export $EXPORTLIST ; \
    ./runjob sjmconfig_202.cfg >& runsjm202.log ; \
    ./runjob sjmconfig_205.cfg >& runsjm205.log ; \
    ./runjob sjmconfig_207.cfg >& runsjm207.log ; ) &

exit
