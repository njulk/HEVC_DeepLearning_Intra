#!/bin/bash
qp=(22 22 27 32 37)
num=${#qp[@]}
sequence=$1
ResultDir=/njulk/HEVC/data/propose/
OrigiDir=/njulk/HEVC/data/original/
cfgDir=/njulk/HEVC/allModels/sequences/
for((i=2;i<num;i++));
do
	 sed -i "s/${qp[i]}/32/" ${cfgDir}encoder_intra_main.cfg
done
ps aux | grep "TAppEncoderStatic" |grep -v grep| cut -c 9-15 | xargs kill -9
/njulk/HEVC/original/HM-16.12+SCM-8.3/bin/TAppEncoderStatic -c  ${cfgDir}encoder_intra_main.cfg -c ${cfgDir}${sequence}.cfg  



