#! /bin/bash
qp=(22 22 27 32 37)
num=${#qp[@]}
#./TAppEncoderStatic -c  encoder_intra_main.cfg -c /njulk/HEVC/HM-16.12+SCM-8.3/cfg/per-sequence/SlideEditing.cfg /njulk/HEVC/data/32/
for((i=1;i<num;i++));
do	
	sed -i "s/${qp[i-1]}/${qp[i]}/" encoder_intra_main.cfg 
	ps aux | grep "/TAppEncoderStatic" |grep -v grep| cut -c 9-15 | xargs kill -9
	rm -rf /njulk/HEVC/data/propose/SlideEditing/${qp[i]}/SlideEditing_1280x720_30.yuv_Result.txt
	rm -rf /njulk/HEVC/data/propose/SlideEditing/${qp[i]}/result.txt
	./TAppEncoderStatic -c  encoder_intra_main.cfg -c SlideEditing.cfg /njulk/HEVC/data/propose/SlideEditing/${qp[i]}/SlideEditing_1280x720_30.yuv_Result.txt | tee /njulk/HEVC/data/propose/SlideEditing/${qp[i]}/result.txt	
done
