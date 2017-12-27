#include"TPersoDef.h"
#include<string.h>
using namespace std;
char*** INTRAMODE_DATA = new char**[4];
long long **MODECOUNT = new long long*[4];
map<unsigned int, int> CuMap;
void init(){
	for(int i=0;i<4;i++){
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
	
}
void Mkdirs() {
	const char* prefix = "";
	init();
	int cuSize[4]={64,32,16,8};
	for(int j=0;j<4;j++){
		char tmp[100];
		memset(tmp,0,100);
		sprintf(tmp,"%s%d",prefix,cuSize[j]);
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
		for(int i=0;i<4;i++)
			delete MODECOUNT[i];
	}
	if(INTRAMODE_DATA != NULL){
		for(int j=0;j<4;j++){
			for(int i=0; i<35; i++){
				delete INTRAMODE_DATA[j][i];
			}
		}
	}
}
