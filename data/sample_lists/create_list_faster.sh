
#!/bin/bash

# Uses direct SAM URL lookup instead of samweb locate-file and pnfsToXRootD
# Runs in parallel
# Keeps optional file limit

source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode v10_10_03_01 -q e26:prof
htgettoken -a htvaultprod.fnal.gov -i sbnd

SAMDEF=${1}                 # Required: SAM definition name 
OUTDIR=${2}                 # Required: Output directory 
EXP=${3}                    # Required: Experiment (sbnd/icarus) 
LIMIT=${4:-""}              # Optional: Number of files
USER_PARALLEL=${5:-""}      # Optional: Number of jobs in parallel 

if [ -z "${SAMDEF}" ] || [ -z "${OUTDIR}" ] || [ -z "${EXP}" ]; then
    echo "Usage:"
    echo "  $0 <SAM definition> <output dir> <experiment> [file limit] [parallel jobs]"
    exit 1
fi

mkdir -p ${OUTDIR}

LIST_CMD="samweb list-definition-files -e ${EXP} ${SAMDEF}"
if [ -n "${LIMIT}" ]; then
    LIST_CMD="${LIST_CMD} | head -n ${LIMIT}"
fi

if [ -n "${LIMIT}" ]; then
    NFILES=${LIMIT}
else
    echo "Counting files in SAM definition..."
    NFILES=$(samweb list-definition-files -e ${EXP} ${SAMDEF} | wc -l)
fi

NCPU=$(nproc)

if [ -n "${USER_PARALLEL}" ]; then
    PARALLEL=${USER_PARALLEL}
else
    MAX_CAP=10   # don't hog shared login nodes
    PARALLEL=$(( NCPU / 2 ))

    if [ ${PARALLEL} -gt ${MAX_CAP} ]; then
        PARALLEL=${MAX_CAP}
    fi
    if [ ${PARALLEL} -gt ${NFILES} ]; then
        PARALLEL=${NFILES}
    fi
    if [ ${PARALLEL} -lt 1 ]; then
        PARALLEL=1
    fi
fi

eval ${LIST_CMD} \
| xargs -P ${PARALLEL} -I {} samweb get-file-access-url --schema=root -e ${EXP} {} \
2>&1 | tee ${OUTDIR}/${SAMDEF}.txt

echo "Done."

