#include"TPersoDef.h"
#include<string.h>
using namespace std;
const int num=5;
char*** INTRAMODE_DATA = new char**[num];
long long **MODECOUNT = new long long*[num];
map<unsigned int, int> CuMap;
FILE** LabelFile=new FILE*[num];
FILE* ResultLog;
const char* prefixPath=NULL;
void init(const char* prefixFilepath){
	for(int i=0;i<num;i++){
		INTRAMODE_DATA[i]=new char*[35];
		for(int j=0;j<35;j++){
			INTRAMODE_DATA[i][j]=new char[100];
			memset(INTRAMODE_DATA[i][j],0,100);
		}
		MODECOUNT[i]=new long long[35];
	}
	CuMap[64]=0;
	CuMap[32]=1;
	CuMap[16]=2;
	CuMap[8]=3;
	CuMap[4]=4;
	const char* tmp[num]={"/label_64x64.txt","/label_32x32.txt","/label_16x16.txt","/label_8x8.txt","/label_4x4.txt"};
	for(int i=0;i<num;i++){
		string tmpPath=string(prefixFilepath)+string(tmp[i]);
		LabelFile[i]=fopen(tmpPath.c_str(),"w+");
	}

}
void Mkdirs(const char* prefixFilepath) {
	prefixPath = prefixFilepath;
	init(prefixPath);
	int cuSize[num]={64,32,16,8,4};
	for(int j=0;j<num;j++){
		char tmp[100];
		memset(tmp,0,100);
		sprintf(tmp,"%s/%d",prefixPath,cuSize[j]);
		int rval=mkdir(tmp,S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
		if(rval!=0){
			perror("mkdir error");
		}
		for(int i = 0; i < 35; i++) {
			char str[100];
			memset(str, 0, 100);
			sprintf(str,"%s%s%d",tmp,"/",i);
			mkdir(str,S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
			memcpy(INTRAMODE_DATA[j][i],str,100);
			MODECOUNT[j][i]=0;
		}
	}
}

void freeData(){
	if(MODECOUNT != NULL){
		for(int i=0;i<num;i++){
			delete MODECOUNT[i];
			fclose(LabelFile[i]);
		}
	}
	if(INTRAMODE_DATA != NULL){
		for(int j=0;j<num;j++){
			for(int i=0; i<35; i++){
				delete INTRAMODE_DATA[j][i];
			}
			delete INTRAMODE_DATA[j];
		}
	}
	delete LabelFile;
	fclose(ResultLog);
}
