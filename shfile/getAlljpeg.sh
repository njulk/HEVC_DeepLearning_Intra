#! /bin/bash
qp=(22 22 27 32 37)
num=${#qp[@]}
BinPath=/njulk/HEVC/HM-16.12+SCM-8.3/bin/
SavePath=/njulk/HEVC/data/jpeg_save/
CfgDir=/njulk/HEVC/allModels/sequences/
SequenceName=$1
rm -rf ${SavePath}${SequenceName}/*
mkdir ${SavePath}${SequenceName}
for((i=1;i<5;i++));do
	sed -i "s/${qp[i]}/32/" ${CfgDir}encoder_intra_main.cfg
done
${BinPath}TAppEncoderStatic -c ${CfgDir}encoder_intra_main.cfg -c ${CfgDir}${SequenceName}.cfg ${SavePath}${SequenceName}/
 sed -i "s/32/22/" ${CfgDir}encoder_intra_main.cfg
