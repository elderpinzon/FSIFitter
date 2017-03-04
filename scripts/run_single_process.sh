#!/bin/bash

PISCAT=/home/s/sbhadra/elder/installed_NEUT/v0/NEUT/v5r3p5/Linux-x86_64/neut_5.3.5/src/pionsmpl/./piscat
CARD_NAME=$1
BASE_NAME=`basename ${CARD_NAME} .c`

PID=$2
if [[ $PID == 211 ]]; then
    PID_NAME=pip
elif [[ $PID == -211 ]]; then
    PID_NAME=pim
else
    echo `printf "Wrong PID: %s" ${PID}`
    exit
fi

PISCAT_NAME=`printf piscat_%s.hbk ${BASE_NAME}`
RUN_SINGLE_PISCAT=`printf "%s %s %s %s" ${PISCAT} ${CARD_NAME} ${PISCAT_NAME} ${PID}`

echo ${RUN_SINGLE_PISCAT}

## Sometimes the piscat process gets hung reading the card file
## Check after given # of secs and if so restart the process 
timeout $3 $RUN_SINGLE_PISCAT > log_${BASE_NAME}.log 2>&1
if [ $? -eq 124 ]; then
    # Timeout occurred
    echo "Process hung... Try again"
    rm ${PISCAT_NAME}
    $RUN_SINGLE_PISCAT > log_${BASE_NAME}.log 2>&1
    echo "... recovered."
else
    # Did not hang
    echo "Process finished successfully!"
fi

## Convert hbk file to ROOT file
h2root ${PISCAT_NAME} > /dev/null 2>&1

## Remove hbk file and used temporary card file
#rm ${PISCAT_NAME} ${CARD_NEW} ${CARD_NAME}
rm ${PISCAT_NAME} ${CARD_NAME}

echo `printf "run_single_process.sh for %s done..." ${BASE_NAME}`
