#ifndef DEEP_LEARNING
#define DEEP_LEARNING
#include<iostream>
#include<stdio.h>
#include <sys/types.h>  
#include <sys/stat.h>
#include<map>
#include<errno.h>
#include<string>
#include"TAppCommon/classification.h"
using namespace std;
extern char**** INTRAMODE_DATA;
extern long long ***MODECOUNT;
extern std::map<unsigned int,int> CuMap;
extern FILE*** LabelFile;
extern FILE* ResultLog;
void Mkdirs(const char* prefixFilepath);
extern const char* prefixPath;
extern int IndexCurFrame;
extern int NumSeperation;
void freeData();
#endif // DEEP_LEARNING
#undef DEEP_CLASSIFY
#ifdef DEEP_CLASSIFY
extern Classifier* classifier8;
extern Classifier* classifier16;
extern Classifier* classifier32;
extern Classifier* classifier64;
void buildClassifier(const char* sequence);
#endif // DEEP_CLASSIFY

