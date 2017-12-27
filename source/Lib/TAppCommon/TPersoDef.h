#ifndef DEEP_LEARNING
#define DEEP_LEARNING
#include<iostream>
#include<stdio.h>
#include <sys/types.h>  
#include <sys/stat.h>
#include<map>
#include<errno.h>
extern char*** INTRAMODE_DATA;
extern long long **MODECOUNT;
extern std::map<unsigned int,int> CuMap;
extern FILE** LabelFile;
extern FILE* ResultLog;
void Mkdirs(const char* prefixFilepath);
extern const char* prefixPath;
void freeData();
#endif // DEEP_LEARNING
