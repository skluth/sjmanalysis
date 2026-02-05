sjmanalysis
===========

Analysis of event shapes, jet rates and moments with e+e- data

The commandline to run h2root on the hbook ntuples, da* for data, mc* for the MC ntuples

ls da*.histo | xargs -I{} sh -c 'echo h2root {} `echo {} | sed -e "s/histo/root/g"` ; h2root {} `echo {}|sed -e "s/histo/root/g"`' > h2rootlog.txt

