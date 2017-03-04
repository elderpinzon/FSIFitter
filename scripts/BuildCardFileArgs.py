#!/usr/bin/python

import re,sys,os,commands
import elementinfo

import sys
import locale

GENERAL_SCRATCH="/scratch/s/sbhadra/elder/NIWG/FSIFitter/scripts/"

#reac = sys.argv[1]
element = sys.argv[1]
mom = locale.atof(sys.argv[2])
pid_name = sys.argv[3]
fsipars = sys.argv[4:11]
fsipars = [float(x) for x in fsipars]

def file_name(element,mom,pid_name,fsipars):

    return "%s_%.3f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%s" % (element,mom,fsipars[0],fsipars[1],fsipars[2],fsipars[3],fsipars[4],fsipars[5],fsipars[6],pid_name)

def neut_card_file(element,mom,pid_name,fsipars):
    
    # Generate neut card file from template
    default_file = open(GENERAL_SCRATCH + "neut.template", 'r')
    write_file = open(file_name(element,mom,pid_name,fsipars)+".c", 'w')

    for line in default_file:
        # if line.find("EVCT-NEVT  100000") != -1 and reac == 1:
        #     line = re.sub("100000", "25000", line) ## Make less MC for reactive
        if line.find("EVCT-PV aho") != -1:
            line = re.sub("aho", "%.3f" % mom, line)
        if line.find("NEUT-NUMBNDN aho") != -1:
            line = re.sub("aho", str(elementinfo.NUMATOM[element]-elementinfo.NUMBNDP[element]), line)
        elif line.find("NEUT-NUMBNDP aho") != -1:
            line = re.sub("aho", str(elementinfo.NUMBNDP[element]), line)
        elif line.find("NEUT-NUMATOM aho") != -1:
            line = re.sub("aho", str(elementinfo.NUMATOM[element]), line)
        elif line.find("NEUT-FEFQE aho") != -1:
            line = re.sub("aho", str(fsipars[0]), line)
        elif line.find("NEUT-FEFABS aho") != -1:
            line = re.sub("aho", str(fsipars[1]), line)
        elif line.find("NEUT-FEFCX aho") != -1:
            line = re.sub("aho", str(fsipars[2]), line)
        elif line.find("NEUT-FEFINEL aho") != -1:
            line = re.sub("aho", str(fsipars[3]), line)
        elif line.find("NEUT-FEFQEH aho") != -1:
            line = re.sub("aho", str(fsipars[4]), line)
        elif line.find("NEUT-FEFCXH aho") != -1:
            line = re.sub("aho", str(fsipars[5]), line)  
        elif line.find("NEUT-FEFALL aho") != -1:
            line = re.sub("aho", str(fsipars[6]), line)  
        write_file.write(line)

    default_file.close()
    write_file.close()

## Build Card File on the fly
neut_card_file(element,mom,pid_name,fsipars)
