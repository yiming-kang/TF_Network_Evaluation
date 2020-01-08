#!/bin/bash
#SBATCH --mem=5G
#SBATCH --ntasks-per-node=1

## Example usage:
## sbatch --mail-type=FAIL,END --mail-user=yiming.kang@wustl.edu evaluate_network_w_nonbinary_benchmark.sh <network_file> <regulators_file> <genes_file> 32 20 <output_directory>
##

NETWORK=$1
REGULATORS=$2
GENES=$3
MAX_RANK=$4 ## unit in 1000, i.e. MAX_RANK*1000
NUM_BINS=$5
OUT_DIR=$6
CHIP_NET=$7
PWM_NET=$8

if [ -z "${CHIP_NET}" ]; then
	CHIP_NET=/scratch/mblab/yiming.kang/proj_netprophet_2.0/resources/yeast_benchmark_network/chip_net.txt
fi
if [ -z "${PWM_NET}" ]; then
	PWM_NET=/scratch/mblab/yiming.kang/proj_netprophet_2.0/resources/yeast_benchmark_network/motif_net.txt
fi

ml R/3.2.1

network_base=${NETWORK##*/}
network_name=${network_base%.adjmtr}
rand_str=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c 10`
adjlst_network=/tmp/${network_name}_${rand_str}.txt
min_rank=$(bc <<< "scale=3; $MAX_RANK/$NUM_BINS")

sed -i 's/ /\t/g' ${NETWORK}

echo -e "REGULATOR\tTARGET\tCONFIDENCE" > ${adjlst_network}
./adjmtr2interactions.rb -a ${NETWORK} -r ${REGULATORS} -c ${GENES} >> ${adjlst_network}

Rscript analyze_binding_overlap.r ${adjlst_network} ${CHIP_NET} ${PWM_NET} ${OUT_DIR}/analysis.${NUM_BINS}bins.top${min_rank}to${MAX_RANK}k.${network_name}.txt ${MAX_RANK} ${NUM_BINS}

rm ${adjlst_network}
