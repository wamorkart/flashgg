#!/bin/bash
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 0 -i2 1 -i3 2 -i4 3 -y "2017" -o data_mix_0123.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 4 -i2 5 -i3 6 -i4 7 -y "2017" -o data_mix_4567.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 8 -i2 9 -i3 10 -i4 11 -y "2017" -o data_mix_891011.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 12 -i2 13 -i3 14 -i4 15 -y "2017" -o data_mix_12131415.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 16 -i2 17 -i3 18 -i4 19 -y "2017" -o data_mix_16171819.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 20 -i2 21 -i3 22 -i4 23 -y "2017" -o data_mix_20212223.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 24 -i2 25 -i3 26 -i4 27 -y "2017" -o data_mix_24252627.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 28 -i2 29 -i3 30 -i4 31 -y "2017" -o data_mix_28293031.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 32 -i2 33 -i3 34 -i4 35 -y "2017" -o data_mix_32333435.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 36 -i2 37 -i3 38 -i4 39 -y "2017" -o data_mix_36373839.root
python H4GTreeMixing.py  -i /eos/user/t/twamorka/h4g_fullRun2/2017/hadd/data_2017_skim.root -t 'Data_13TeV_H4GTag_0' -i1 40 -i2 41 -i3 42 -i4 43 -y "2017" -o data_mix_40414243.root
#hadd /eos/user/t/twamorka/h4g_fullRun2/2016/hadd/data_mix.root data_mix_0123.root data_mix_4567.root data_mix_891011.root data_mix_12131415.root data_mix_16171819.root data_mix_20212223.root data_mix_24252627.root data_mix_28293031.root data_mix_32333435.root data_mix_36373839.root data_mix_40414243.root
