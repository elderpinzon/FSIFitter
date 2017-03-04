import re,os
from itertools import product
#import BuildCardFile
import elementinfo

## Folder to store everything
GENERAL_SCRATCH = os.path.dirname(os.path.realpath(__file__)) + "/"

## Nominal parameter values
fefqe_nom = 0.9
fefabs_nom = 1.25
fefcx_nom = 0.8
fefinel_nom = 1.0
fefqeh_nom = 1.8
fefcxh_nom = 1.8
fefall_nom = 1.0

## Grid parameters
grid_step_LE = 0.1
grid_step_HE = 0.2
#grid_size = 7 ##Sould be an odd number!!

## List of parameters to do
## This are the ones that I hadn't done already
## (Remove the if conditions to do the full grid)

## (0.1 ~ 0.5) and (1.3 ~1.7)
lfefqe = [x*grid_step_LE for x in range(1,18)]
## (0.35 ~ 0.85) and (1.65 ~ 1.95)
lfefabs = [(x*grid_step_LE)+0.05 for x in range(3,20)]
## (0.1 ~ 0.4) and (1.2 ~ 1.6)
lfefcx = [x*grid_step_LE for x in range(1,17)]
## (0.2) and (1.8 ~ 2.6)
lfefinel = [x*grid_step_HE for x in range(1,14)]
## (0.8 ~ 1.0) and (2.6 ~ 2.8)
lfefqeh = [x*grid_step_HE for x in range(4,15)]

## Build file name from parameters
def file_name(element,mom,pid_name,fsipars):

    return "%s_%.3f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%s" % (element,mom,fsipars[0],fsipars[1],fsipars[2],fsipars[3],fsipars[4],fsipars[5],fsipars[6],pid_name)


## Prints the command that will be added to the parallel stack inside the script
def PrintActualCommand(counter,Element,fmomentum,PID_NAME,fsipars,pid):

    ## Time-out waiting period in seconds
    TOSEC = "500"

    ## Name for file, etc.
    PREFIX_NAME = "piscat_" +file_name(Element,fmomentum,PID_NAME,fsipars) + ".root"

    ## Build card file on the fly
    command = "python " + GENERAL_SCRATCH + "BuildCardFileArgs.py %s %s %s %.2f %.2f %.2f %.2f %.2f %.2f %.2f && " % (Element,fmomentum,PID_NAME,fsipars[0],fsipars[1],fsipars[2],fsipars[3],fsipars[4],fsipars[5],fsipars[6])

    ## Run bash script with actual NEUT particle gun instructions
    command += GENERAL_SCRATCH + "./run_single_process.sh " + file_name(Element,fmomentum,PID_NAME,fsipars)+ ".c " + str(pid) + " " + TOSEC +" && "

    ## Run app to get cross section values
    command += GENERAL_SCRATCH + "./extract_xsec.exe " + PREFIX_NAME + " %d %.3f %d %.3f %.2f %.2f %.2f %.2f %.2f %.2f %.2f;\n" % (elementinfo.NUMBNDP[Element],fmomentum,pid,elementinfo.eftarget[Element],fsipars[0],fsipars[1],fsipars[2],fsipars[3],fsipars[4],fsipars[5],fsipars[6])

    return command

## Builds and (if asked to) submits jobs for given parameters
def build_pbs_files(option,Element,momentum,pid,processes_per_file,submit=False):

    ## Check that option is supported
    if(option == "all_nominal" or option == "single"):
        print "\nBuilding pbs submission scripts. Option %s. Processes per file: %d" % (option,processes_per_file)
    else:
        print "\nOption: %s not supported... exiting" % option
        return

    ## Create folders
    SCRATCH = GENERAL_SCRATCH + Element.lower()
    os.system("mkdir -p " + SCRATCH + " " + SCRATCH + "/logs " + SCRATCH + "/files " + SCRATCH +  "/summary " + SCRATCH +  "/card " + SCRATCH + "/pbs")
    
    if pid == 211:
        PID_NAME = "pip"
    elif pid == -211:
        PID_NAME = "pim"
    else:
        print "Invalid pid :%s" % pid

    ## Total number of permutations
    if(option == "all_nominal"):
        total = len(momentum)
        print "Building pbs for A:%s. All momenta. Total: %d. Divide into %d pieces" % (Element,total,total/processes_per_file+1)

    elif(option == "single"):
        total = len(lfefqe) * len(lfefabs) * len(lfefcx) * len(lfefinel) * len(lfefqeh)
        if(momentum>500):
            total = len(lfefinel) * len(lfefqeh)
            print total
        if(momentum < 400):
            total = len(lfefqe) * len(lfefabs) * len(lfefcx)
        print "Building pbs for A:%s p:%.3f. Total %d. Divide into %d pieces" % (Element,momentum,total,total/processes_per_file+1)

    ## List to store file names
    job_file_name = []
    
    for job in range(total/processes_per_file + 1):
        
        ## Open job template file
        default_file = open("job.template",'r')

        ## Create new pbs script 
        if(option == "all_nominal"):
            job_file_name.append("%s_%s_indiv_%d" % (Element,PID_NAME,job))
        elif(option == "single"):
            job_file_name.append("%s_%s_%.3f_%d" % (Element,PID_NAME,momentum,job))

        write_file = open("%s/pbs/%s.pbs" % (SCRATCH,job_file_name[job]), 'w')
        print "\nWorking on job: %s/pbs/%s.pbs " % (SCRATCH,job_file_name[job])

        ## Loop over template and make mods
        for line in default_file:

            if line.find("#PBS -N piscat-aho") != -1:
                if(option == "all_nominal"):
                    line = re.sub("aho", "%s-%s-indiv-%d" % (Element,PID_NAME,job), line)
                elif(option == "single"):
                    line = re.sub("piscat-aho", "%s-%s-%.3f-%d" % (Element,PID_NAME,momentum,job), line)

            elif line.find("insert_parallel") != -1:
                line = "parallel -j 8 <<EOF\n"
                counter = 0

                ## In this case loop over momenta:
                if(option == "all_nominal"):
                
                    for imom in range(len(momentum)):
                        
                        fmomentum = momentum[imom]
                        fefqe = fefqe_nom
                        fefabs = fefabs_nom
                        fefcx = fefcx_nom
                        fefinel = fefinel_nom
                        fefqeh = fefqeh_nom
                        fefcxh = fefcxh_nom
                        fefall = fefall_nom

                        fsipars = [fefqe,fefabs,fefcx,fefinel,fefqeh,fefcxh,fefall]
                        if job == counter/processes_per_file:
                            line += PrintActualCommand(counter,Element,fmomentum,PID_NAME,fsipars,pid)
                            
                        counter += 1
                
                ## In this case loop over LE parameter combinations
                elif(option == "single" and momentum < 400):

                    for iqe,iabs,icx in product(range(len(lfefqe)),range(len(lfefabs)),range(len(lfefcx))):
                    
                        fmomentum = momentum
                        fefqe  = lfefqe[iqe]
                        fefabs = lfefabs[iabs]
                        fefcx  = lfefcx[icx]
                        fefinel = fefinel_nom
                        fefqeh = fefqeh_nom
                        fefcxh = fefcxh_nom
                        fefall = fefall_nom
                        
                        fsipars = [fefqe,fefabs,fefcx,fefinel,fefqeh,fefcxh,fefall]
                        if job == counter/processes_per_file:
                            line += PrintActualCommand(counter,Element,fmomentum,PID_NAME,fsipars,pid)
                                                        
                        counter += 1

                ## In this case loop over HE parameter combinations
                elif(option == "single" and momentum > 500):

                    for iinel,iqeh in product(range(len(lfefinel)),range(len(lfefqeh))):
                    
                        fmomentum = momentum
                        fefqe = fefqe_nom
                        fefabs = fefabs_nom
                        fefcx = fefcx_nom
                        fefinel  = lfefinel[iinel]
                        fefqeh = lfefqeh[iqeh]
                        fefcxh = fefcxh_nom
                        fefall = fefall_nom
                                                
                        fsipars = [fefqe,fefabs,fefcx,fefinel,fefqeh,fefcxh,fefall]
                        if job == counter/processes_per_file:
                            line += PrintActualCommand(counter,Element,fmomentum,PID_NAME,fsipars,pid)
                            
                        counter += 1

                ## In this case loop over LE,HE parameter combinations
                elif(option == "single" and momentum > 400 and momentum < 500):

                    for iqe,iabs,icx,iinel,iqeh in product(range(len(lfefqe)),range(len(lfefabs)),range(len(lfefcx)),range(len(lfefinel)),range(len(lfefqeh))):
                        
                        fmomentum = momentum
                        fefqe  = lfefqe[iqe]
                        fefabs = lfefabs[iabs]
                        fefcx  = lfefcx[icx]
                        fefinel  = lfefinel[iinel]
                        fefqeh = lfefqeh[iqeh]
                        fefcxh = fefcxh_nom
                        fefall = fefall_nom
                                            
                        fsipars = [fefqe,fefabs,fefcx,fefinel,fefqeh,fefcxh,fefall]
                        if job == counter/processes_per_file:
                            line += PrintActualCommand(counter,Element,fmomentum,PID_NAME,fsipars,pid)
                            
                        counter += 1

                line += "EOF\n"

            elif line.find("insert_summary") != -1:
                line = "    cat *.txt > %s/summary/%s.dat\n" % (SCRATCH,job_file_name[job])
                #line += '    for f in *.log; do (echo "Log-File: $f"; cat "${f}"; echo) >> %s/logs/logs_%s.logs; done\n' % (SCRATCH,job_file_name[job])
                #line += "    tar cf %s/files/files_%s.tar *.root\n" % (SCRATCH,job_file_name[job])
        
            write_file.write(line)
           
        print "Closing default file"
        default_file.close()
        write_file.close()

    ## Now that pbs files are build, go back and add them to the queue
    for job in range(total/processes_per_file+1):
        print "Submitting job %s.pbs to the queue" % job_file_name[job]
        #submit_command = "cd " + SCRATCH + "/pbs; qsub -q debug " + job_file_name[job] + ".pbs"
        submit_command = "cd " + SCRATCH + "/pbs; qsub " + job_file_name[job] + ".pbs"
        ## If last job send email when it is done
        if job == total/processes_per_file:
            submit_command = "cd " + SCRATCH + "/pbs; qsub -m ae -M elderpinzon@gmail.com " + job_file_name[job] + ".pbs"
        print submit_command
        if submit == True:
            os.system(submit_command)
