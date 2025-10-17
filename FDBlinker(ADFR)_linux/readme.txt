#Installation guidance
install ADFR  
conda install openbabel

"ADFR_HOME" ADFR installation path

ftlk    $ADFR_HOME/bin/    running ftlk.py and convert parameters into ftlk.py

autogrid4 path   $ADFR_HOME/bin/autogrid4

ftlk.py   $ADFR_HOME/CCSBpckgs/ADFR/bin      line46 set "ADFR_HOME   = os.environ.get("ADFR_HOME", "your_path/ADFRsuite_x86_64Linux_1.0")"

pdbqt_tools.py  $ADFR_HOME/CCSBpckgs/ADFR/bin   ftlk.py  call its function

partial_minimize_ob3.py $ADFR_HOME/CCSBpckgs/ADFR/bin   not minimize the cycle of fragment 

MakeGrids.py   $ADFR_HOME/CCSBpckgs/ADFR/util  line81  set "   known_path = "your_path/ADFRsuite_x86_64Linux_1.0/bin/autogrid4"  "   #autogrid4 path   $ADFR_HOME/bin/autogrid4

runADFR.py   $ADFR_HOME/CCSBpckgs/ADFR/utils/   

#run commandline 
./ftlk --complex /input_pdb_path/0.pdb --ligand /input_frag_path/Toluene.pdbqt --atomsA A1 A2 A3 --atomsB B1 --jobName job --outputPath /output_path --maxEvals 400000 --nbRuns 8

