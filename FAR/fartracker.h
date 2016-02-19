#ifndef fartracker_h
#define fartracker_h

// C interface

typedef struct far_rect_t {
	float x;
	float y;
	float width;
	float height;
} far_rect_t;

typedef void* far_tracker_t;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

	far_tracker_t far_init(const unsigned char *gray, int width, int height, far_rect_t rect);
	far_rect_t far_track(far_tracker_t tracker, const unsigned char *gray);
	far_rect_t far_restart(far_tracker_t tracker, const unsigned char *gray, far_rect_t rect);
	bool far_check(far_tracker_t tracker);
	void far_release(far_tracker_t tracker);

#ifdef __cplusplus
}
#endif /* __cplusplus */

// C++ implementation

#ifdef __cplusplus

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cfloat>
using namespace std;

#include "Eigen/Dense"
using namespace Eigen;

const float PI_2 = acos(-1.0f) * 0.5f;

const int cell_min = 2;
const int fast_n = 25;
const int fast_step = 4;
const int fine_n = 160;
const int cell_n = 40;
const int fine_steps[] = { 27, 9, 3, 1 };
const int max_iteration = 9;
const float padding = 1.5f;
const float sigmoid_factor = 7.141f;
const float sigmoid_bias = 0.482f;
const float translate_eps = 0.005f;
const float threshold_error = 0.4f;

typedef Matrix<float, 8, 1> Vector8f;
typedef Matrix<float, 32, 1> Vector32f;

float rectArea(far_rect_t rect);
float rectOverlap(far_rect_t a, far_rect_t b);
far_rect_t rectBound(far_rect_t rect);
ostream& operator<<(ostream& cout, const far_rect_t &rect);

class Surf
{
public:
	Surf(int width, int height);

	Vector4f kernel(float angle);
	void process(const unsigned char *gray, float angle);
	void set_cell(float cell);
	void set_step(int step);

	Ref<Vector8f> cell_hist(int x, int y);
	void descriptor(float x, float y, Ref<Vector8f> f);
	void gradient(float x, float y, Ref<Vector8f> f, Ref<Vector8f> dx, Ref<Vector8f> dy);
	void descriptor4(float x, float y, Ref<Vector32f> f);
	void gradient4(float x, float y, Ref<Vector32f> f, Ref<Vector32f> dx, Ref<Vector32f> dy);

public:
	float angle, tx[4], ty[4];
	Vector4f kx, ky;
	int W, H, C, step;

private:
	MatrixXf sum, hist;
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

class FARTracker
{
public:
	FARTracker(const unsigned char *gray, int width, int height, far_rect_t rect, ostream *os = NULL);
	far_rect_t track(const unsigned char *gray);
	far_rect_t restart(const unsigned char *gray, far_rect_t rect);
	bool check();

private:
	Vector3f locate(far_rect_t rect);
	far_rect_t window(Vector3f translate);

	void fast_train(Warp warp);
	void fine_train(Warp warp);
	Vector3f fast_test(Warp warp);
	Warp fine_test(Warp warp);

	float sigmoid(float x);
	Warp Lucas_Kanade(Warp warp);
	float evaluate(Warp warp);

public:	
	int image_width, image_height;
	float window_width, window_height;
	Warp warp;
	float roll, yaw, pitch;
	vector<float> fine_errors;

private:
	Surf feature;
	vector<Vector2i> fast_samples;
	vector<Vector3f> fine_samples;		
	MatrixXf fast_model, fine_model;
	ostream *log;
	int N;
};

#endif /* __cplusplus */

#endif /* fartracker_h */
