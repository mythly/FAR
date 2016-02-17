#ifndef far_h
#define far_h

#include "cv_face.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
using namespace std;

#include "Eigen/Dense"
using namespace Eigen;

const float PI_2 = acos(-1.0f) * 0.5f;

const int cell_min = 2;
const int fast_n = 25;
const int fast_step = 2;
const int fine_n = 160;
const int cell_n = 40;
const int fine_steps[] = { 9, 3, 1 };
const int max_iteration = 8;
const float padding = 1.6f;
const float sigmoid_factor = 7.141f;
const float sigmoid_bias = 0.482f;
const float translate_eps = 0.005f;
const float threshold_error = 0.4f;

typedef Matrix<float, 8, 1> Vector8f;
typedef Matrix<float, 32, 1> Vector32f;

class Surf
{
public:
	Surf(int width, int height);

	Vector4f kernel(float angle);
	void process(const unsigned char *gray, float angle);
	void set_cell(float cell);
	void set_step(int step);

	float* cell_hist(int x, int y);
	float cell_norm(int x, int y);
	void descriptor(float x, float y, float *f);
	void gradient(float x, float y, float *f, float *dx, float *dy);
	void descriptor4(float x, float y, float *f);
	void gradient4(float x, float y, float *f, float *dx, float *dy);

public:
	float angle, tx[4], ty[4];
	Vector4f kx, ky;
	int W, H, C, step;

private:
	MatrixXf sum, hist, norm;
	MatrixXi flag;
	Vector8f zero;
};

class Warp
{
public:
	Warp(int width, int height);

	void setr(Vector3f rotate);
	void sett(Vector3f translate);

	Vector2f project(Vector3f p);
	Vector3f transform(Vector3f p);
	Vector2f transform2(Vector3f p);

	Matrix<float, 2, 6> gradient(Vector3f p);
	void steepest(Matrix<float, 6, 1> parameters);
	void euler(float &roll, float &yaw, float &pitch);

public:
	Vector2f c;
	float f;
	Vector3f r;
	Vector3f t;

private:
	Matrix3f R, Dx, Dy, Dz;
};

class FART
{
public:
	FART(const unsigned char *gray, int width, int height, cv_rect_t rect, ostream *os = NULL);
	cv_rect_t track(const unsigned char *gray);
    void restart(const unsigned char *gray, cv_rect_t rect);

private:
	Vector3f locate(cv_rect_t rect);
	cv_rect_t window(Vector3f translate);

	void fast_train(Warp warp);
	void fine_train(Warp warp);
	Vector3f fast_test(Warp warp);
	Warp fine_test(Warp warp);

	float sigmoid(float x);
	Warp Lucas_Kanade(Warp warp);
	float evaluate(Warp warp);

public:
	ostream *log;
    int image_width, image_height;
    float window_width, window_height;
	Surf feature;
    Warp warp;
	vector<Vector2i> fast_samples;
	vector<Vector3f> fine_samples;
	vector<float> fine_errors;
	MatrixXf fast_model, fine_model;
    float error, roll, yaw, pitch;
	int N;
};

ostream& operator<<(ostream& cout, const cv_rect_t &rect);

#endif /* far_h */

