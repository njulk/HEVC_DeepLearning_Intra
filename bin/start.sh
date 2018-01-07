#! /bin/bash
 ps aux | grep "/TAppEncoderStatic" |grep -v grep| cut -c 9-15 | xargs kill -9
./TAppEncoderStatic -c  encoder_intra_main.cfg -c /njulk/HEVC/HM-16.12+SCM-8.3/cfg/per-sequence/SlideEditing.cfg /njulk/HEVC/data/32/

