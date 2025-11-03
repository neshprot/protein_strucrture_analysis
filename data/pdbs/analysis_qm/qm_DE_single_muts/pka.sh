#!/bin/bash

P=$1
Project=${P::-4}
~/script_library/new_rhodopsin_optimization_mm/genpsf.sh ${Project}.pdb
rm ${Project}_new.psf
~/script_library/new_rhodopsin_optimization_mm/pkacalc.sh ${Project}_new.pdb
pka88=$(grep " 88 A" D88E_D118L_qm_new.pka | awk '{print $1, $2, $4}' | tail -1)
pka215=$(grep " 215 A" D88E_D118L_qm_new.pka | awk '{print $1, $2, $4}' | tail -1)
pka197=$(grep " 197 A" D88E_D118L_qm_new.pka | awk '{print $1, $2, $4}' | tail -1)
pka207=$(grep " 207 A" D88E_D118L_qm_new.pka | awk '{print $1, $2, $4}' | tail -1)
echo $Project $pka88 $pka215 $pka197 $pka207 >> pkas
rm A.pdb ${Project}_new* chains.txt PROTLIST
