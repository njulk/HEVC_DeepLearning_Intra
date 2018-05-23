#! /bin/bash
qp=(22 22 27 32 37)
num=${#qp[@]}
sequence=$1
ResultDir=/njulk/HEVC/data/propose/
OrigiDir=/njulk/HEVC/data/original/
cfgDir=/njulk/HEVC/allModels/sequences/
for((i=2;i<num;i++));
do
	 sed -i "s/${qp[i]}/22/" ${cfgDir}encoder_intra_main.cfg
done
mkdir $ResultDir$sequence
#for((i=1;i<num;i++));
#do	
#	mkdir $ResultDir$sequence/${qp[i]}
#	sed -i "s/${qp[i-1]}/${qp[i]}/" ${cfgDir}encoder_intra_main.cfg 
#	ps aux | grep "TAppEncoderStatic" |grep -v grep| cut -c 9-15 | xargs kill -9
#	rm -rf $ResultDir$sequence/${qp[i]}/*_Result.txt
#	rm -rf $ResultDir$sequence/${qp[i]}/result.txt
#	/njulk/HEVC/HM-16.12+SCM-8.3/bin/TAppEncoderStatic -c  ${cfgDir}encoder_intra_main.cfg -c ${cfgDir}${sequence}.cfg $ResultDir$sequence/${qp[i]}/ $sequence | tee $ResultDir$sequence/${qp[i]}/result.txt	
#	rm -rf rec.yuv
#	rm -rf str.bin
#done

#sed -i "s/37/22/" ${cfgDir}encoder_intra_main.cfg
mkdir $OrigiDir$sequence
for((i=1;i<num;i++));
do
	mkdir $OrigiDir$sequence/${qp[i]}
        sed -i "s/${qp[i-1]}/${qp[i]}/" ${cfgDir}encoder_intra_main.cfg
        ps aux | grep "TAppEncoderStatic" |grep -v grep| cut -c 9-15 | xargs kill -9
        rm -rf $OrigiDir$sequence/${qp[i]}/result.txt
        /njulk/HEVC/original/HM-16.12+SCM-8.3/bin/TAppEncoderStatic -c  ${cfgDir}encoder_intra_main.cfg -c ${cfgDir}${sequence}.cfg | tee $OrigiDir$sequence/${qp[i]}/result.txt  
        rm -rf rec.yuv
        rm -rf str.bin
done



