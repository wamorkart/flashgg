#!/bin/bash

# for m in 60 45 35 25 15;
# do
#    echo ${m}
#    for n in 2016 2017 2018;
#    do
#      echo ${n}
#      for o in CatMVA_PhoMVA_Only CatMVA_PhoMVA_KinVars CatMVA_PhoMVA_KinVars_Many
#      do
#        echo ${o}
#        python ApplyMVATransformation.py -id /eos/user/t/twamorka/h4g_fullRun2/noSystematics/${n}/hadd/ReducedSignal_CatMVA/${o}/ -iF signal_m_${m}.root -m ${m} -y ${n}
#      done
#    done
# done


for n in 2016 2017 2018;
do
  echo ${n}
  for o in CatMVA_PhoMVA_Only CatMVA_PhoMVA_KinVars CatMVA_PhoMVA_KinVars_Many
  do
    echo ${o}
    python ApplyMVATransformation.py -id /eos/user/t/twamorka/h4g_fullRun2/noSystematics/${n}/hadd/ReducedSignal_CatMVA/${o}/ -iF data_${n}.root  -y ${n}
  done
done


for n in 2016 2017 2018;
do
  echo ${n}
  for o in CatMVA_PhoMVA_Only CatMVA_PhoMVA_KinVars CatMVA_PhoMVA_KinVars_Many
  do
    echo ${o}
    python ApplyMVATransformation.py -id /eos/user/t/twamorka/h4g_fullRun2/noSystematics/${n}/hadd/ReducedSignal_CatMVA/${o}/ -iF data_mix.root  -y ${n}
  done
done
