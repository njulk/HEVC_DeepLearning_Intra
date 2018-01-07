#include"TPersoDef.h"

#include"unistd.h"
using namespace std;
int IndexCurFrame = 0;
const int num = 4;
const int numTestTrain = 2;
int NumSeperation = 30;
char**** INTRAMODE_DATA = new char***[numTestTrain];
long long ***MODECOUNT = new long long**[numTestTrain];
map<unsigned int, int> CuMap;
FILE*** LabelFile = new FILE**[numTestTrain];
FILE* ResultLog;
const char* prefixPath = NULL;
int createDirectory(const char* directoryPath, mode_t mode);
void init(const char* prefixFilepath) {
	for (int x = 0; x < numTestTrain; x++) {
		INTRAMODE_DATA[x] = new char**[num];
		MODECOUNT[x] = new long long*[num];
		for (int i = 0; i < num; i++) {
			INTRAMODE_DATA[x][i] = new char*[35];
			for (int j = 0; j < 35; j++) {
				INTRAMODE_DATA[x][i][j] = new char[100];
				memset(INTRAMODE_DATA[x][i][j], 0, 100);
			}
			MODECOUNT[x][i] = new long long[35];
		}
		LabelFile[x]=new FILE*[num];
	}
	CuMap[64] = 0;
	CuMap[32] = 1;
	CuMap[16] = 2;
	CuMap[8] = 3;
	const char* tmp_train[num] = { "/train/train_label_64x64.txt","/train/train_label_32x32.txt","/train/train_label_16x16.txt","/train/train_label_8x8.txt" };
	const char* tmp_test[num] = { "/test/test_label_64x64.txt","/test/test_label_32x32.txt","/test/test_label_16x16.txt","/test/test_label_8x8.txt" };
	for (int i = 0; i < num; i++) {
		string tmpPath = string(prefixFilepath) + string(tmp_train[i]);
		createDirectory(tmpPath.c_str(), 0666);
		LabelFile[0][i] = fopen(tmpPath.c_str(), "w+");
		if(LabelFile[0][i]==NULL){
			perror(tmpPath.c_str());
		}
		tmpPath = string(prefixFilepath) + string(tmp_test[i]);
		createDirectory(tmpPath.c_str(), 0666);
		LabelFile[1][i] = fopen(tmpPath.c_str(), "w+");
		if(LabelFile[1][i]==NULL){
                        perror(tmpPath.c_str());
                }
	}

}


//避免多级目录的失败
int createDirectory(const char* directoryPath, mode_t mode)
{
	int dirPathLen = strlen(directoryPath);
	char tmpDirPath[256] = { 0 };
	for (int i = 0; i < dirPathLen; ++i)
	{
		tmpDirPath[i] = directoryPath[i];
		if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
		{
			if (access(tmpDirPath, 0) != 0)
			{
				int ret = mkdir(tmpDirPath, mode);
				if (ret != 0)
				{
					return ret;
				}
			}
		}
	}
	return 0;
}

void Mkdirs(const char* prefixFilepath) {
	prefixPath = prefixFilepath;
	init(prefixPath);
	int cuSize[num] = { 64,32,16,8 };
	const char* testTrain[2] = { "train","test" };
	for (int x = 0; x < 2; x++) {
		for (int j = 0; j < num; j++) {
			char tmp[100];
			memset(tmp, 0, 100);
			sprintf(tmp, "%s/%s/%d", prefixPath, testTrain[x], cuSize[j]);
			int rval = mkdir(tmp, 0666);
			if (rval != 0) {
				perror(tmp);
			}
			for (int i = 0; i < 35; i++) {
				char str[100];
				memset(str, 0, 100);
				sprintf(str, "%s%s%d", tmp, "/", i);
				mkdir(str, 0666);
				memcpy(INTRAMODE_DATA[x][j][i], str, 100);
				MODECOUNT[x][j][i] = 0;
			}
		}
	}
}

void freeData() {
	if (MODECOUNT != NULL) {
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < num; i++) {
				delete MODECOUNT[j][i];
				fclose(LabelFile[j][i]);
			}
			delete MODECOUNT[j];
			delete LabelFile[j];
		}		
	}
	if (INTRAMODE_DATA != NULL) {
		for (int x = 0; x < 2; x++) {
			for (int j = 0; j < num; j++) {
				for (int i = 0; i < 35; i++) {
					delete INTRAMODE_DATA[x][j][i];
				}
				delete INTRAMODE_DATA[x][j];
			}
			delete INTRAMODE_DATA[x];
		}		
	}
	delete LabelFile;
	fclose(ResultLog);
}


#ifdef DEEP_CLASSIFY
Classifier* classifier8=new Classifier("/njulk/HEVC/sequences/prototxt/SlideEditing8_deploy.prototxt","/njulk/HEVC/sequences/prototxt/8_iter_250000.caffemodel","/njulk/HEVC/sequences/lmdb/SlideEditing_8x8_mean_train.binaryproto","/njulk/HEVC/sequences/lmdb/label.txt");
Classifier* classifier16= new Classifier("/njulk/HEVC/sequences/prototxt/SlideEditing16_deploy.prototxt","/njulk/HEVC/sequences/prototxt/16_iter_10000.caffemodel","/njulk/HEVC/sequences/lmdb/SlideEditing_16x16_mean_train.binaryproto","/njulk/HEVC/sequences/lmdb/label.txt");
Classifier* classifier32 = new Classifier("/njulk/HEVC/sequences/prototxt/SlideEditing32_deploy.prototxt","/njulk/HEVC/sequences/prototxt/32_iter_5000.caffemodel","/njulk/HEVC/sequences/lmdb/SlideEditing_32x32_mean_train.binaryproto","/njulk/HEVC/sequences/lmdb/label.txt");
Classifier* classifier64 = new Classifier("/njulk/HEVC/sequences/prototxt/SlideEditing64_deploy.prototxt","/njulk/HEVC/sequences/prototxt/64_iter_5000.caffemodel","/njulk/HEVC/sequences/lmdb/SlideEditing_64x64_mean_train.binaryproto","/njulk/HEVC/sequences/lmdb/label.txt");
#endif
