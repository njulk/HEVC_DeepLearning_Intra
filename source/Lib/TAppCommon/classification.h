#ifndef DEEP_CLASSIFY
#define DEEP_CLASSIFY
#include <caffe/caffe.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <algorithm>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>
#include<memory.h>
using std::string;
//using namespace std;
using namespace caffe;
typedef std::pair<string, float> Prediction;

class Classifier {
public:
	Classifier(const string& model_file,
		const string& trained_file,
		const string& mean_file,
		const string& label_file);

	std::vector<Prediction> Classify(const cv::Mat& img, int N = 5);

private:
	void SetMean(const string& mean_file);

	std::vector<float> Predict(const cv::Mat& img);

	void WrapInputLayer(std::vector<cv::Mat>* input_channels);

	void Preprocess(const cv::Mat& img,
		std::vector<cv::Mat>* input_channels);

private:
        caffe::shared_ptr<Net<float>> net_;
	cv::Size input_geometry_;
	int num_channels_;
	cv::Mat mean_;
	std::vector<string> labels_;
};


Prediction getPrediction(Classifier* classifier, string file);

#endif // !DEEP_CLASSIFY
