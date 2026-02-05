#!/bin/bash

# Load shared libs from root and fastjet
LD_LIBRARY_PATH=$PWD:/opt/root/lib:/opt/fastjet/fastjet-3.5.1/lib/
# Unbuffered Fortran IO for readability
GFORTRAN_UNBUFFERED_PRECONNECTED=y
EXPORTLIST="LD_LIBRARY_PATH GFORTRAN_UNBUFFERED_PRECONNECTED"

#( export $EXPORTLIST ; echo "Start LEP1 job" ; \
#    ./runjob sjmconfig_91_all.cfg >& runsjm91_all.log ; ) &
( export $EXPORTLIST ; echo "Start LEP1 96-98 job" ; \
  ./runjob sjmconfig_91_96-98.cfg >& runsjm91_96-98.log ; ) &
( export $EXPORTLIST ; echo "Start LEP1 99-2k job" ; \
  ./runjob sjmconfig_91_99-2k.cfg >& runsjm91_99-2k.log ; ) &
( export $EXPORTLIST ; \
    echo "Start 130 job" ; ./runjob sjmconfig_130.cfg >& runsjm130.log ; \
    echo "Start 136 job" ; ./runjob sjmconfig_136.cfg >& runsjm136.log ; \
    echo "Start 161 job" ; ./runjob sjmconfig_161.cfg >& runsjm161.log ; \
    echo "Start 172 job" ; ./runjob sjmconfig_172.cfg >& runsjm172.log ; \
    echo "Start 183 job" ; ./runjob sjmconfig_183.cfg >& runsjm183.log ; ) &
( export $EXPORTLIST ; echo "Start 189, 192, 196 200 jobs" ; \
    ./runjob sjmconfig_189.cfg >& runsjm189.log ; \
    ./runjob sjmconfig_192.cfg >& runsjm192.log ; \
    ./runjob sjmconfig_196.cfg >& runsjm196.log ; \
    ./runjob sjmconfig_200.cfg >& runsjm200.log ; ) &
( export $EXPORTLIST ; echo "Start 202, 205, 207 jobs" ; \
    ./runjob sjmconfig_202.cfg >& runsjm202.log ; \
    ./runjob sjmconfig_205.cfg >& runsjm205.log ; \
    ./runjob sjmconfig_207.cfg >& runsjm207.log ; ) &

exit
