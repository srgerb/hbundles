for i in {1..144}; do echo "/software/rosetta/latest/bin/score_jd2.hdf5.linuxgccrelease -symmetry_definition /home/srgerb/hBundles/SymFiles/C4_Z.sym -in:file:silent_struct_type binary -out:file:scorefile "$i"_relax.sc -in:file:silent " | tr '\r\n' ' '; for j in sg180508_"$i"_*; do echo $j |tr '\r\n' ' ' ;done ; echo " "; done > array_tasks.list