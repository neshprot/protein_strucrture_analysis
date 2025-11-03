#!/bin/bash

propkadir=/home/ibsbnicigt/nikolaev_d/bin

Help()
{
   # Display Help
   echo ""
   echo ""
   echo "Function for calculation of pKa of titratable residues"
   echo ""
   echo "./pkacalc.sh Project.pdb" 
   echo ""
   echo "where Project is the name of your protein pdb file." 
   echo "Output: Project.pka with pKa values, PROTLIST with ASP, GLU, HIS pKas"
   echo "h     Print this Help."
   echo ""
   echo ""
}

while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;   
   esac           
done              
                  
P=$1
Project=${P::-4}

sed  "s/HID/HIS/g; s/HIE/HIS/g; s/HIP/HIS/g; s/GLH/GLU/; s/ASH/ASP/; s/CLA/CL /g; s/BRA/BR /" ${Project}.pdb > ${Project}_0.pdb
$propkadir/propka31 ${Project}_0.pdb
mv ${Project}_0.pka ${Project}.pka 
rm ${Project}_0.pdb 

awk '$0 == "       Group      pKa  model-pKa   ligand atom-type" {i=1;next};i && i++ <= 30' ${Project}.pka > PROTLIST
grep "ASP" PROTLIST > PASP
grep "GLU" PROTLIST > PGLU
grep "HIS" PROTLIST > PHIS
cat PGLU >> PASP
cat PHIS >> PASP
mv PASP PROTLIST
rm PGLU PHIS
grep "\*" ${Project}.pka >> PROTLIST
