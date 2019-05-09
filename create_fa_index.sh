#!/bin/bash

VERSION="1.0"
AUTHOR="Yibin Wang"
PROG="$(basename $0)"
function info {
    echo
    echo "To create fasta index, to variants calling."
    echo "                                Copyright:$AUTHOR"
    echo
}

function usage {
    echo -e "Usage : $PROG -f fasta [-t]"
    echo -e "Use option -h for more information."
}

function help {
    info;
    usage;
    echo
    echo "$PROG $VERSION"
    echo
    echo " -f : input fasta file"
    echo " [-t] : num_thread [default:1]"
    echo " [-h] : help"
    exit;
}

#Get options
while getopts "f:t:vh" opt; do
    case $opt in
        f)fasta=$OPTARG;;
        t)nct=$OPTARG;;
        v)version ;;
        h)help;;
        ?)echo "Invalid option : -$OPTARG"
            help;
            exit;;
        :)echo "Option -$OPTARG requires an argument."
            help;
            exit;;
    esac
done

if [[ -z $fasta ]]; then
    help;
    exit
fi

if [[ -z $nct ]]; then
    nct=1
fi

#Determine the fasta file if exist
if [ ! -e $fasta ]; then
    echo "[ERROR] No such file of $fasta"
    exit 1
fi

#Determine the fasta index fai dict whether if exist
if [ ! -e ${fasta}.fai ]; then
    echo "[INFO] Create ${fasta}.fai ..."
    samtools faidx $fasta
    echo "Done"
fi

fa_prefix=$(echo $fasta |sed "s/.fa.*//")
if [ ! -e ${fa_prefix}.dict ]; then
    echo -e "[INFO] Create ${fa_prefix}.dict..."
    java -jar /public1/home/stu_wangyibin/software/picard.jar CreateSequenceDictionary R=$fasta O=${fa_prefix}.dict
    echo "Done"
fi
:<<!
if [ ! -e ${fasta}.bwt ]; then
    echo "Create ${fasta}.bwt ..."
    bwa index -a bwtsw ${fasta}
    echo "Done"
fi
!
echo -e "[INFO] $fasta index was created."

