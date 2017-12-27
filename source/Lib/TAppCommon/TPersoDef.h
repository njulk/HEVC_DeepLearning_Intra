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
void Mkdirs();
void freeData();
#endif // DEEP_LEARNING
