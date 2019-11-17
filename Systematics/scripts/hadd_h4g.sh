#!/bin/bash

#for i in 50 45 40 35 30 25 20 15 10;
for i in 40;
do
    hadd signal_m_${i}.root /eos/user/t/twamorka/Nov11/output_SUSYGluGluToHToAA_AToGG_M-${i}_*.root
    mv signal_m_${i}.root /eos/user/t/twamorka/Nov11_wScalesandSmearings/Signal/
    hadd_workspaces signal_m_${i}_ws.root /eos/user/t/twamorka/Nov11/output_SUSYGluGluToHToAA_AToGG_M-${i}_*.root
    mv signal_m_${i}_ws.root /eos/user/t/twamorka/Nov11_wScalesandSmearings/Signal/
done
