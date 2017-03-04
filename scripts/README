## Submitting jobs in batch from csv files in pion_data/

Run 'python ScanMomentaSubmit.py'. It will grab the momenta values from the csv files in pion_data and prepare scripts to run piscat for those momenta with variations of the FSI parameters

The flag to submit or not is in ScanMomentaSubmit.py

The details of the grid of FSI parameters to be scanned is in ParameterVariations.py (or ParameterVariationsExtra.py)

## Submitting jobs indvidually

Running 'ParameterVariations.build_pbs_files("single","c",mom,211,601,SubmitFlag)' will run the particle guns for the momenta specified as mom (211 means it is piP)

## Dealing with '-nan' results

Sometimes a piscat process will exit without error code but somehow without any events recorded. Here is how to re-run that piscat run, extract the results and fix the summary tables.

1) Use grep to find the lines with 'nan' as a result. Dump then in a text file. This can take a couple minutes.

   grep -r 'nan' c/summary/c_*.txt fe/summary/fe_pi* cu/summary/cu_pi* pb/summary/pb_pi* > nan_scan_072516.dat

2) Run ReSubmitNan.py on this file. It will read the nuclei, polarity and FSI pars info from each line and re-run piscat. It will also find replace the incorrect line on the summary file

   python ReRunNan.py nan_scan_072516.dat

## To merge the results from all runs use something like this command

   cat fe/summary/*.dat cu/summary/*.dat pb/summary/*.dat > scan_all_fe_cu_pb.dat

## Read into root file for FSIFitter
To read the summary .dat file from the previous step into a root file use the simple macro read_summary_file.C