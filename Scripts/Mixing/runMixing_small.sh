#!/bin/bash
for i in 0 4;
do
  echo $i $(( $i + 1 )) $(( $i + 2 )) $(( $i + 3 ))
  python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 $i -i2 $(( $i + 1 )) -i3 $(( $i + 2 )) -i4 $(( $i + 3 )) -y "2017" -o /eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/mixing/data_mix_${i}$(( $i + 1 ))$(( $i + 2 ))$(( $i + 3 )).root
  # echo $(( $i + 1 ))
done
