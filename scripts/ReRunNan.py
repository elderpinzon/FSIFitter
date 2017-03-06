import os
import sys
import re
import fileinput
import elementinfo

## Argument should be the nan_scan file
if len(sys.argv) < 2:
    print "Usage: python ReRunNan.py [scan_file_name]        (See README for details)"
    sys.exit(0)

os.system('echo -n "Started at: "\n')
os.system('date\n')

os.system('source /home/s/sbhadra/elder/installed_NEUT/v0/setup_defaultGEANT.sh')

## Folder for storage (this path)
GENERAL_SCRATCH = os.path.dirname(os.path.realpath(__file__)) + "/"

nan_file = open(sys.argv[1], 'r')

for line in nan_file:

    words = re.split('/|:| ',line)

    ## Grab element inifo
    element = words[0]
    element_Z = words[3]

    ## Grab file name
    filename = words[2]
    
    ## Grab momentum
    momenta = words[4]

    ## Grab pion polarity
    pid = words[5]

    if pid == "211":
        pion = "pip"
    if pid == "-211":
        pion = "pim"
        
    ## Grab fsi parameters
    fsipars = words[6:13]
    ## Remove list formatting
    fsipars_under = '_'.join(fsipars)
    fsipars = ' '.join(fsipars)

    ## Basename for outputs
    basename = element + "_" + momenta + "_" + fsipars_under + "_" + pion

    ## Check that this hasn't been run already
    if os.path.isfile("piscat_" + basename +".txt"):
        print "Already ran this... skip"
        sys.exit(0)
    
    ## Re-build command piece by piece

    preface = "python " + GENERAL_SCRATCH + "BuildCardFileArgs.py "
    command = preface + element + " " + momenta + " " + pion + " " + fsipars
    run_command = " && " + GENERAL_SCRATCH + "/./run_single_process.sh "
    command += run_command + basename + ".c " + pid
    xsec_command = " 500 && " + GENERAL_SCRATCH + "/./extract_xsec.exe piscat_"
    command += xsec_command + basename + ".root " + element_Z
    command += " " + momenta + " " + pid + " " +  str(elementinfo.eftarget[element]) + " " + fsipars + ";"

    ## And then run it
    print command + "\n"
    os.system(command)

    ## Grab fixed summary line (result from running piscat)
    fixed_line = "oops"
    result_summary_file = open("piscat_" + basename +".txt", 'r')
    for fline in result_summary_file:
        fixed_line = fline

    print "Fixed line: " + fixed_line

    ## Create backup copy of summary filw with wrong line
    nan_summary = element + "/summary/" + filename
    os.system('cp %s %s.bak' % (nan_summary,nan_summary))
    nan_summary_file = open(nan_summary + ".bak", 'r')

    ##Open fixed file
    fixed_summary_file = open(nan_summary, 'w')
    for fline in nan_summary_file:
        if fsipars in fline:
            fixed_summary_file.write(fixed_line)
        else:
            fixed_summary_file.write(fline)
            

    result_summary_file.close()
    nan_summary_file.close()
    fixed_summary_file.close()

os.system('echo -n "Finished at: "\n')
os.system('date\n')
