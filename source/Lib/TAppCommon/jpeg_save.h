#ifndef JPEG_SAVE
#define JPEG_SAVE
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sys/types.h>
#include<sys/stat.h>
extern "C" {
#include <jpeglib.h>
}
	void GeneJpegFile(const char* jpegFileName, unsigned char* inputData, int nWidth, int nHeight, int nChannel, int nQuality);
#endif
