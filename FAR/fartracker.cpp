#include "fartracker.h"

// C interface

far_tracker_t far_init(const unsigned char *gray, int width, int height, far_rect_t rect)
{
	return new FARTracker(gray, width, height, rect, &cerr);
}

far_rect_t far_track(far_tracker_t tracker, const unsigned char *gray)
{
	assert(tracker != NULL);
	FARTracker *t = static_cast<FARTracker*>(tracker);
	return t->track(gray);
}

far_rect_t far_retrack(far_tracker_t tracker, const unsigned char *gray, const far_rect_t rects[], int n_rects)
{
	assert(tracker != NULL);
	FARTracker *t = static_cast<FARTracker*>(tracker);	
	return t->retrack(gray, vector<far_rect_t>(rects, rects + n_rects));
}

void far_transform(far_tracker_t tracker, far_rect_t start_rect, float *x, float *y)
{
	assert(tracker != NULL);
	FARTracker *t = static_cast<FARTracker*>(tracker);
	assert(x != NULL && y != NULL);
	Vector3f p3(*x - (start_rect.x + start_rect.width * 0.5f), *y - (start_rect.y + start_rect.height * 0.5f), 0.0f);
	Vector2f p2 = t->warp.transform2(p3);
	*x = p2.x();
	*y = p2.y();
}

void far_info(far_tracker_t tracker, float *error, float *roll, float *yaw, float *pitch)
{
    assert(tracker != NULL);
    FARTracker *t = static_cast<FARTracker*>(tracker);
    assert(error != NULL && roll != NULL && yaw != NULL && pitch != NULL);
    *error = t->fine_errors.back();
    t->warp.euler(*roll, *yaw, *pitch);
    *roll *= 90.0f / PI_2;
    *yaw *= 90.0f / PI_2;
    *pitch *= 90.0f / PI_2;
}

bool far_check(far_tracker_t tracker)
{
	assert(tracker != NULL);
	FARTracker *t = static_cast<FARTracker*>(tracker);
	return t->check();
}

void far_release(far_tracker_t tracker)
{
	if (tracker != NULL) {
		FARTracker *t = static_cast<FARTracker*>(tracker);
		delete t;
	}
}

// C++ implementation

float rectArea(far_rect_t rect)
{
	if (rect.width <= 0.0f || rect.height <= 0.0f)
		return 0.0f;
	return rect.width * rect.height;
}

float rectOverlap(far_rect_t a, far_rect_t b)
{
	float sa = rectArea(a), sb = rectArea(b);
	if (sa <= 0.0f || sb <= 0.0f)
		return 0.0f;
	float lx = max(a.x, b.x);
	float ly = max(a.y, b.y);
	float rx = min(a.x + a.width, b.x + b.width);
	float ry = min(a.y + a.height, b.y + b.height);
	return max(rx - lx, 0.0f) * max(ry - ly, 0.0f);
}

ostream& operator<<(ostream& cout, const far_rect_t&rect)
{
	cout << "[" << int(rect.width) << " x " << int(rect.height);
	cout << " from (" << int(rect.x) << ", " << int(rect.y) << ")]";
	return cout;
}

Surf::Surf(int width, int height):
W(width), H(height), C(0), step(0),
img(H + 2, W + 2),
flag(H, W),
sum(H + 1, W + 1),
hist(H, W),
zero(1, 1)
{
	sum.set(0.0f);
	zero.set(0.0f);	
}

void Surf::rotate(float angle, float kernel[4])
{
	float c = cos(angle), s = sin(angle);
	float wx = c * (1.0f - abs(s));
	float wy = s * (1.0f - abs(c));
	float wu = 0.0f, wv = 0.0f;
	if (c >= 0)
		(s >= 0 ? wu : wv) = c * s;
	else
		(s >= 0 ? wv : wu) = -c * s;
	kernel[0] = wx;
	kernel[1] = wy;
	kernel[2] = wu;
	kernel[3] = wv;	
}

void Surf::process(const unsigned char *gray, float angle)
{		
	for (int y = 0; y < img.rows; ++y) {
		int yy = y - 1;
		if (y == 0)
			yy = 0;
		if (y == img.rows - 1)
			yy = H - 1;
		float *f = img.ptr(y);
		const unsigned char *g = gray + yy * W;		
		f[0] = float(g[0]);
		f[img.cols - 1] = float(g[W - 1]);
		for (int x = 1; x <= W; ++x)
			f[x] = float(g[x - 1]);
	}	
	A = angle;
	float kx[4], ky[4];
	rotate(A, kx);
	rotate(A + PI_2, ky);	
	for (int y = 1; y <= H; ++y) {
		float *f0 = img.ptr(y - 1);
		float *f = img.ptr(y);
		float *f1 = img.ptr(y + 1);
		float *s0 = sum.ptr(y - 1);
		float *s1 = sum.ptr(y);
		float s[8] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		for (int x = 1; x <= W; ++x) {
			++f0;
			++f;
			++f1;
			s0 += 8;
			s1 += 8;
			float gx = f[1] - f[-1];
			float gy = f1[0] - f0[0];
			float gu = f1[1] - f0[-1];
			float gv = f1[-1] - f0[1];
			float dx = kx[0] * gx + kx[1] * gy + kx[2] * gu + kx[3] * gv;
			float dy = ky[0] * gx + ky[1] * gy + ky[2] * gu + ky[3] * gv;
			if (dx > 0.0f) {
				if (dy > 0.0f) {
					s[1] += dx;
					s[3] += dx;
					s[5] += dy;
					s[7] += dy;
				}
				else {
					s[0] += dx;
					s[2] += dx;
					s[5] += dy;
					s[7] += -dy;
				}
			}
			else {
				if (dy > 0.0f) {
					s[1] += dx;
					s[3] += -dx;
					s[4] += dy;
					s[6] += dy;
				}
				else {
					s[0] += dx;
					s[2] += -dx;
					s[4] += dy;
					s[6] += -dy;
				}
			}
			s1[0] = s0[0] + s[0];
			s1[1] = s0[1] + s[1];
			s1[2] = s0[2] + s[2];
			s1[3] = s0[3] + s[3];
			s1[4] = s0[4] + s[4];
			s1[5] = s0[5] + s[5];
			s1[6] = s0[6] + s[6];
			s1[7] = s0[7] + s[7];
		}
	}

	C = 0;
	step = 0;
	flag.set(0);	
}

void Surf::set_cell(float cell)
{
	cell = cell * 0.5f;
	X[0] = -cell; Y[0] = -cell;
	X[1] = cell; Y[1] = -cell;
	X[2] = cell; Y[2] = cell;
	X[3] = -cell; Y[3] = cell;
	for (int i = 0; i < 4; ++i) {
		float x = cos(A) * X[i] - sin(A) * Y[i];
		float y = sin(A) * X[i] + cos(A) * Y[i];
		X[i] = x;
		Y[i] = y;
	}
	C = max(int(floor(cell)), cell_min);
}

void Surf::set_step(int step)
{
	this->step = step;
}

float* Surf::cell_hist(int x, int y)
{	
	if (x < 0 || x >= flag.cols || y < 0 || y >= flag.rows)
		return zero.data();
	if (flag(y, x) != C) {
		int x0 = max(x - C, 0);
		int x1 = min(x + C + 1, W);
		int y0 = max(y - C, 0);
		int y1 = min(y + C + 1, H);
		float *s00 = sum.ptr(y0, x0);
		float *s01 = sum.ptr(y0, x1);
		float *s10 = sum.ptr(y1, x0);
		float *s11 = sum.ptr(y1, x1);
		float *h = hist.ptr(y, x);
		for (int i = 0; i < 8; ++i)
			h[i] = s11[i] + s00[i] - s01[i] - s10[i];
		flag(y, x) = C;
	}
	return hist.ptr(y, x);
}

void Surf::descriptor(float x, float y, float *f)
{
	x = x / step;
	y = y / step;
	int ixp = (int)floor(x);
	int iyp = (int)floor(y);
	float wx1 = x - ixp, wx0 = 1.0f - wx1;
	float wy1 = y - iyp, wy0 = 1.0f - wy1;
	float w00 = wx0 * wy0;
	float w01 = wx0 * wy1;
	float w10 = wx1 * wy0;
	float w11 = wx1 * wy1;
	ixp *= step;
	iyp *= step;
	float *f00 = cell_hist(ixp, iyp);
	float *f01 = cell_hist(ixp, iyp + step);
	float *f10 = cell_hist(ixp + step, iyp);
	float *f11 = cell_hist(ixp + step, iyp + step);
	for (int i = 0; i < 8; ++i)
		f[i] = f00[i] * w00 + f01[i] * w01 + f10[i] * w10 + f11[i] * w11;	
}

void Surf::gradient(float x, float y, float *f, float *dx, float *dy)
{
	x = x / step;
	y = y / step;
	int ixp = (int)floor(x);
	int iyp = (int)floor(y);
	float wx1 = x - ixp, wx0 = 1.0f - wx1;
	float wy1 = y - iyp, wy0 = 1.0f - wy1;
	float w00 = wx0 * wy0;
	float w01 = wx0 * wy1;
	float w10 = wx1 * wy0;
	float w11 = wx1 * wy1;
	wx0 /= step;
	wy0 /= step;
	wx1 /= step;
	wy1 /= step;
	ixp *= step;
	iyp *= step;
	float *f00 = cell_hist(ixp, iyp);
	float *f01 = cell_hist(ixp, iyp + step);
	float *f10 = cell_hist(ixp + step, iyp);
	float *f11 = cell_hist(ixp + step, iyp + step);
	for (int i = 0; i < 8; ++i) {
		f[i] = f00[i] * w00 + f01[i] * w01 + f10[i] * w10 + f11[i] * w11;
		dx[i] = (f10[i] - f00[i]) * wy0 + (f11[i] - f01[i]) * wy1;
		dy[i] = (f01[i] - f00[i]) * wx0 + (f11[i] - f10[i]) * wx1;
	}
}

void Surf::descriptor4(float x, float y, float *f)
{
	for (int i = 0; i < 4; ++i)
		descriptor(x + X[i], y + Y[i], f + i * 8);
	float S = 0.0f;
	for (int i = 0; i < 32; ++i)
		S += f[i] * f[i];
	float iS = S < 1.0f ? 0.0f : 1.0f / sqrt(S);
	for (int i = 0; i < 32; ++i)
		f[i] *= iS;
}

void Surf::gradient4(float x, float y, float *f, float *dx, float *dy)
{
	for (int i = 0; i < 4; ++i)
		gradient(x + X[i], y + Y[i], f + i * 8, dx + i * 8, dy + i * 8);
	float S = 0.0f, Sx = 0.0f, Sy = 0.0f;
	for (int i = 0; i < 32; ++i) {
		S += f[i] * f[i];
		Sx += f[i] * dx[i];
		Sy += f[i] * dy[i];
	}
	float iS = S < 1.0f ? 0.0f : 1.0f / sqrt(S);
	float iSx = Sx * iS * iS * iS;
	float iSy = Sy * iS * iS * iS;
	for (int i = 0; i < 32; ++i) {
		dx[i] = dx[i] * iS - f[i] * iSx;
		dy[i] = dy[i] * iS - f[i] * iSy;
		f[i] *= iS;
	}
}

Warp::Warp(int width, int height) :
c(width * 0.5f, height * 0.5f),
f(float(max(width, height)))
{
	setr(Vector3f(0.0f, 0.0f, 0.0f));
	sett(Vector3f(0.0f, 0.0f, 0.0f));
}

void Warp::setr(Vector3f rotate)
{
	r = rotate;
	double rx = r(0), ry = r(1), rz = r(2);
	double theta = sqrt(rx * rx + ry * ry + rz * rz);
	double R[9], J[27];
	if (theta < DBL_EPSILON)
	{
		memset(R, 0, sizeof(R));
		R[0] = R[4] = R[8] = 1;
		memset(J, 0, sizeof(J));
		J[5] = J[15] = J[19] = -1;
		J[7] = J[11] = J[21] = 1;
	}
	else
	{
		const double I[] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

		double c = cos(theta);
		double s = sin(theta);
		double c1 = 1. - c;
		double itheta = theta ? 1. / theta : 0.;

		rx *= itheta; ry *= itheta; rz *= itheta;

		double rrt[] = { rx*rx, rx*ry, rx*rz, rx*ry, ry*ry, ry*rz, rx*rz, ry*rz, rz*rz };
		double _r_x_[] = { 0, -rz, ry, rz, 0, -rx, -ry, rx, 0 };
		for (int k = 0; k < 9; k++)
			R[k] = c*I[k] + c1*rrt[k] + s*_r_x_[k];
		double drrt[] = { rx + rx, ry, rz, ry, 0, 0, rz, 0, 0,
			0, rx, 0, rx, ry + ry, rz, 0, rz, 0,
			0, 0, rx, 0, 0, ry, rx, ry, rz + rz };
		double d_r_x_[] = { 0, 0, 0, 0, 0, -1, 0, 1, 0,
			0, 0, 1, 0, 0, 0, -1, 0, 0,
			0, -1, 0, 1, 0, 0, 0, 0, 0 };
		for (int i = 0; i < 3; i++)
		{
			double ri = i == 0 ? rx : i == 1 ? ry : rz;
			double a0 = -s*ri, a1 = (s - 2 * c1*itheta)*ri, a2 = c1*itheta;
			double a3 = (c - s*itheta)*ri, a4 = s*itheta;
			for (int k = 0; k < 9; k++)
				J[i * 9 + k] = a0*I[k] + a1*rrt[k] + a2*drrt[i * 9 + k] +
				a3*_r_x_[k] + a4*d_r_x_[i * 9 + k];
		}
	}

	float _R[9], _J[3][9];
	for (int i = 0; i < 9; ++i)
		_R[i] = float(R[i]);
	for (int i = 0; i < 27; ++i)
		_J[i / 9][i % 9] = float(J[i]);
	this->R << _R[0], _R[1], _R[2], _R[3], _R[4], _R[5], _R[6], _R[7], _R[8];
	Dx << _J[0][0], _J[1][0], _J[2][0],
		_J[0][3], _J[1][3], _J[2][3],
		_J[0][6], _J[1][6], _J[2][6];
	Dy << _J[0][1], _J[1][1], _J[2][1],
		_J[0][4], _J[1][4], _J[2][4],
		_J[0][7], _J[1][7], _J[2][7];
	Dz << _J[0][2], _J[1][2], _J[2][2],
		_J[0][5], _J[1][5], _J[2][5],
		_J[0][8], _J[1][8], _J[2][8];
}

void Warp::sett(Vector3f translate)
{
	t = translate;
}

Vector2f Warp::project(Vector3f p)
{
	return Vector2f(p.x() / p.z() * f, p.y() / p.z() * f) + c;
}

Vector3f Warp::transform(Vector3f p)
{
	return R * p + t;
}

Vector2f Warp::transform2(Vector3f p)
{
	return project(transform(p));
}

Vector2f Warp::gradient(Vector3f p, Matrix<float, 2, 6> &dW)
{
	Matrix3f D1 = p.x() * Dx + p.y() * Dy + p.z() * Dz;	
	Vector3f tp = transform(p);		
	float fz = f / tp.z(), fzz = f / (tp.z() * tp.z());
	dW(0, 3) = fz;
	dW(0, 4) = 0.0f;
	dW(0, 5) = -tp.x() * fzz;
	dW(1, 3) = 0.0f;
	dW(1, 4) = fz;
	dW(1, 5) = -tp.y() * fzz;
	dW.leftCols(3).noalias() = dW.rightCols(3) * D1;
	return project(tp);
}

void Warp::steepest(Matrix<float, 6, 1> parameters)
{
	float rx = r.x() + parameters(0);
	float ry = r.y() + parameters(1);
	float rz = r.z() + parameters(2);
	float tx = t.x() + parameters(3);
	float ty = t.y() + parameters(4);
	float tz = t.z() + parameters(5);
	setr(Vector3f(rx, ry, rz));
	sett(Vector3f(tx, ty, tz));
}

void Warp::euler(float &roll, float &yaw, float &pitch)
{
	if (abs(1 - abs(R(2, 1))) > 1.0e-7f) {
		roll = atan2(-R(0, 1), R(1, 1));
		yaw = atan2(-R(2, 0), R(2, 2));
		pitch = asin(R(2, 1));
	}
	else {
		roll = atan2(R(1, 0), R(0, 0));
		yaw = 0.0f;
		pitch = R(2, 1) > 0 ? PI_2 : -PI_2;
	}
}

FARTracker::FARTracker(const unsigned char *gray, int width, int height, far_rect_t rect, ostream *os) :
image_width(width), image_height(height),
window_width(rect.width), window_height(rect.height),
warp(width, height),
feature(width, height),
log(os)
{
	warp.sett(locate(rect));
	float fine_stride = sqrt(window_width * window_height / fine_n);
	int W = int(floor(window_width / (2.0f * fine_stride)));
	int H = int(floor(window_height / (2.0f * fine_stride)));
	for (int y = 0; y <= 2 * H; ++y)
	for (int x = 0; x <= 2 * W; ++x)
		fine_samples.push_back(Vector3f((x - W) * fine_stride, (y - H) * fine_stride, 0.0f));

	feature.process(gray, 0.0f); 
	N = 0;
	fine_train(warp);
	fast_train(warp);
	fine_errors.push_back(0.0f);
	roll = yaw = pitch = 0.0f;
}

far_rect_t FARTracker::track(const unsigned char *gray)
{
	if (log != NULL) {
		(*log) << "roll = " << roll * 90.0f / PI_2 << endl;
		(*log) << "yaw = " << yaw * 90.0f / PI_2 << endl;
		(*log) << "pitch = " << pitch * 90.0f / PI_2 << endl;
	}	
	feature.process(gray, roll);

	Warp w = warp;
	if (log != NULL)
		(*log) << "track at " << w.t.transpose() << " " << window(w.t) << endl;
	w = fine_test(w);
	float e = evaluate(w);
	if (e > threshold_error) {
		Warp w2 = warp;		
		w2.sett(fast_test(warp));
		if (log != NULL)
			(*log) << "search at " << w2.t.transpose() << " " << window(w2.t) << endl;
		w2 = fine_test(w2);
		float e2 = evaluate(w2);
		if (e2 < e) {
			w = w2;
			e = e2;
		}
	}
	update(w, e);

	return window(w.t);
}

far_rect_t FARTracker::retrack(const unsigned char *gray, const vector<far_rect_t> &detections)
{
	if (log != NULL) {
		(*log) << "roll = " << roll * 90.0f / PI_2 << endl;
		(*log) << "yaw = " << yaw * 90.0f / PI_2 << endl;
		(*log) << "pitch = " << pitch * 90.0f / PI_2 << endl;
	}
	feature.process(gray, 0.0f);

	Warp w = warp;
	w.setr(Vector3f(0.0f, 0.0f, 0.0f));
	if (log != NULL)
		(*log) << "track at " << w.t.transpose() << " " << window(w.t) << endl;
	w = fine_test(w);
	float e = evaluate(w);
	Warp w2 = warp;
	w2.setr(Vector3f(0.0f, 0.0f, 0.0f));
	w2.sett(fast_test(warp));
	if (log != NULL)
		(*log) << "search at " << w2.t.transpose() << " " << window(w2.t) << endl;
	w2 = fine_test(w2);
	float e2 = evaluate(w2);
	if (e2 < e) {
		w = w2;
		e = e2;
	}
	for (auto d : detections) {
		Warp w3 = warp;
		w3.setr(Vector3f(0.0f, 0.0f, 0.0f));
		w3.sett(locate(d));
		if (log != NULL)
			(*log) << "detect at " << w3.t.transpose() << " " << window(w3.t) << endl;
		w3 = fine_test(w3);
		float e3 = evaluate(w3);
		if (e3 < e) {
			w = w3;
			e = e3;
		}
	}	 	
	update(w, e);

	return window(w.t);
}

bool FARTracker::check()
{
	auto iter = fine_errors.rbegin();
	int n = 0;
	float s = 0.0f;
	while (n < 10 && iter != fine_errors.rend()) {
		s += *iter;
		++iter;
		++n;
	}
	return s / n < threshold_error;
}


Vector3f FARTracker::locate(far_rect_t rect)
{
	float scale = sqrt(window_width * window_height / rectArea(rect));
	float x = rect.x + rect.width * 0.5f - warp.c.x();
	float y = rect.y + rect.height * 0.5f - warp.c.y();
	return scale * Vector3f(x, y, warp.f);
}

far_rect_t FARTracker::window(Vector3f translate)
{
	Vector2f center = warp.project(translate);
	float scale = warp.f / translate.z();
	far_rect_t ret;
	ret.width = window_width * scale;
	ret.height = window_height * scale;
	ret.x = center.x() - ret.width * 0.5f;
	ret.y = center.y() - ret.height * 0.5f;
	return ret;
}

void FARTracker::update(Warp w, float e)
{
	if (log != NULL) {
		(*log) << "final translation = " << w.t.transpose() << " " << window(w.t) << endl;
		(*log) << "final rotation = " << w.r.transpose() << endl;
		(*log) << "final error = " << e << endl;
	}
	fine_errors.push_back(e);
	if (e < threshold_error) {
		warp = w;
		warp.euler(roll, yaw, pitch);
		fine_train(warp);
	}
	else
		warp.t = w.t * (warp.t.z() / w.t.z());
	fast_train(warp);
}

void FARTracker::fast_train(Warp warp)
{
	far_rect_t rect = window(warp.t);
	float fast_stride = sqrt(rectArea(rect) / fast_n);
	feature.set_cell(fast_stride);
	int W = int(rect.width * 0.5f / fast_stride);
	int H = int(rect.height * 0.5f / fast_stride);
	int ox = int(rect.width * 0.5f + 0.5f);
	int oy = int(rect.height * 0.5f + 0.5f);
	int stride = int(fast_stride + 0.5f);
	fast_samples.clear();
	for (int y = 0; y <= 2 * H; ++y)
	for (int x = 0; x <= 2 * W; ++x)
		fast_samples.push_back(Vector2i(ox + (x - W) * stride, oy + (y - H) * stride));

	fast_model = MatrixXf::Zero(8, fast_samples.size());
	int x = int(rect.x + 0.5f);
	int y = int(rect.y + 0.5f);
	for (int i = 0; i < fast_samples.size(); ++i) {
		int tx = x + fast_samples[i].x();
		int ty = y + fast_samples[i].y();
		memcpy(fast_model.col(i).data(), feature.cell_hist(tx, ty), 8 * sizeof(float));		
	}
}

void FARTracker::fine_train(Warp warp)
{
	far_rect_t rect = window(warp.t);
	float fine_cell = sqrt(rectArea(rect) / cell_n);
	feature.set_cell(fine_cell);
	feature.set_step(1);

	MatrixXf model(32, fine_samples.size());
	for (int i = 0; i < fine_samples.size(); ++i) {
		Vector2f p = warp.transform2(fine_samples[i]);
		feature.descriptor4(p.x(), p.y(), model.col(i).data());
	}
	if (N == 0) {
		N = 1;
		fine_model = model;
	}
	else {
		++N;
		fine_model = (float(N - 1) / N) * fine_model + (1.0f / N) * model;
	}
}

Vector3f FARTracker::fast_test(Warp warp)
{
	far_rect_t rect = window(warp.t);
	float fast_stride = sqrt(rectArea(rect) / fast_n);
	feature.set_cell(fast_stride);
	far_rect_t region = window(warp.t / (1.0f + padding));
	float minminx = -rect.width * 0.5f;
	float minminy = -rect.height * 0.5f;
	float maxmaxx = image_width + rect.width * 0.5f;
	float maxmaxy = image_height + rect.height * 0.5f;
	int minx = int(max(region.x, minminx) + 0.5f);
	int miny = int(max(region.y, minminy) + 0.5f);
	int maxx = int(min(region.x + region.width, maxmaxx) - rect.width + 0.5f);
	int maxy = int(min(region.y + region.height, maxmaxy) - rect.height + 0.5f);

	float best_score = 0.0f;
	Vector3f best_translate = warp.t;
	MatrixXf model(8, fast_samples.size());
	for (int y = miny; y <= maxy; y += fast_step)
	for (int x = minx; x <= maxx; x += fast_step) {		
		for (int i = 0; i < fast_samples.size(); ++i) {
			int tx = x + fast_samples[i].x();
			int ty = y + fast_samples[i].y();
			memcpy(model.col(i).data(), feature.cell_hist(tx, ty), 8 * sizeof(float));
		}
		float S = model.squaredNorm();
		float score = S < 1.0f ? 0.0f : model.cwiseProduct(fast_model).sum() / sqrt(S);		
		if (score > best_score) {
			far_rect_t best_rect = rect;
			best_rect.x = float(x);
			best_rect.y = float(y);
			best_translate = locate(best_rect);
			best_score = score;
		}
	}
	return best_translate;
}

Warp FARTracker::fine_test(Warp warp)
{
	far_rect_t rect = window(warp.t);
	float fine_cell = sqrt(rectArea(rect) / cell_n);
	feature.set_cell(fine_cell);
	for (auto fine_step : fine_steps) {
		if (fine_step > 2.0f * fine_cell)
			continue;
		feature.set_step(fine_step);
		if (log != NULL)
			(*log) << "\tcell = " << fine_cell << " step = " << fine_step << " : ";
		warp = Lucas_Kanade(warp);
	}
	return warp;
}

float FARTracker::sigmoid(float x)
{
	return 1.0f / (1.0f + exp(-sigmoid_factor * (x - sigmoid_bias)));
}

void FARTracker::hessian(Matrix<float, 6, 6> &H, float w, const Matrix<float, 2, 6> &dW, const Matrix<float, 32, 2> &dF)
{
	float F00 = w * dF.col(0).squaredNorm();
	float F11 = w * dF.col(1).squaredNorm();
	float F01 = w * dF.col(0).dot(dF.col(1));
	float x0 = dW(0, 0) * F00 + dW(1, 0) * F01;
	float y0 = dW(0, 0) * F01 + dW(1, 0) * F11;
	H(0, 0) += x0 * dW(0, 0) + y0 * dW(1, 0);
	H(0, 1) += x0 * dW(0, 1) + y0 * dW(1, 1);
	H(0, 2) += x0 * dW(0, 2) + y0 * dW(1, 2);
	H(0, 3) += x0 * dW(0, 3) + y0 * dW(1, 3);
	H(0, 4) += x0 * dW(0, 4) + y0 * dW(1, 4);
	H(0, 5) += x0 * dW(0, 5) + y0 * dW(1, 5);
	float x1 = dW(0, 1) * F00 + dW(1, 1) * F01;
	float y1 = dW(0, 1) * F01 + dW(1, 1) * F11;
	H(1, 1) += x1 * dW(0, 1) + y1 * dW(1, 1);
	H(1, 2) += x1 * dW(0, 2) + y1 * dW(1, 2);
	H(1, 3) += x1 * dW(0, 3) + y1 * dW(1, 3);
	H(1, 4) += x1 * dW(0, 4) + y1 * dW(1, 4);
	H(1, 5) += x1 * dW(0, 5) + y1 * dW(1, 5);
	float x2 = dW(0, 2) * F00 + dW(1, 2) * F01;
	float y2 = dW(0, 2) * F01 + dW(1, 2) * F11;
	H(2, 2) += x2 * dW(0, 2) + y2 * dW(1, 2);
	H(2, 3) += x2 * dW(0, 3) + y2 * dW(1, 3);
	H(2, 4) += x2 * dW(0, 4) + y2 * dW(1, 4);
	H(2, 5) += x2 * dW(0, 5) + y2 * dW(1, 5);
	float x3 = dW(0, 3) * F00 + dW(1, 3) * F01;
	float y3 = dW(0, 3) * F01 + dW(1, 3) * F11;
	H(3, 3) += x3 * dW(0, 3) + y3 * dW(1, 3);
	H(3, 4) += x3 * dW(0, 4) + y3 * dW(1, 4);
	H(3, 5) += x3 * dW(0, 5) + y3 * dW(1, 5);
	float x4 = dW(0, 4) * F00 + dW(1, 4) * F01;
	float y4 = dW(0, 4) * F01 + dW(1, 4) * F11;
	H(4, 4) += x4 * dW(0, 4) + y4 * dW(1, 4);
	H(4, 5) += x4 * dW(0, 5) + y4 * dW(1, 5);
	float x5 = dW(0, 5) * F00 + dW(1, 5) * F01;
	float y5 = dW(0, 5) * F01 + dW(1, 5) * F11;
	H(5, 5) += x5 * dW(0, 5) + y5 * dW(1, 5);
}

Warp FARTracker::Lucas_Kanade(Warp warp)
{
	float last_E = 1.0f;
	for (int iter = 0; iter < max_iteration; ++iter) {
		Matrix<float, 6, 1> G = Matrix<float, 6, 1>::Constant(0.0f);
		Matrix<float, 6, 6> H = Matrix<float, 6, 6>::Constant(0.0f);
		float E = 0.0f;
		for (int i = 0; i < fine_samples.size(); ++i) {
			Matrix<float, 2, 6> dW;
			Vector2f p = warp.gradient(fine_samples[i], dW);			 
			Vector32f F;
			Matrix<float, 32, 2> dF;
			feature.gradient4(p.x(), p.y(), F.data(), dF.col(0).data(), dF.col(1).data());
			F -= fine_model.col(i);
			float e = sigmoid(F.squaredNorm());
			float w = sigmoid_factor * e * (1.0f - e);			
			G += w * (dW.transpose() * (dF.transpose() * -F));
			//H.triangularView<Upper> += w * (dW.transpose() * (dF.transpose() * dF) * dW);
			hessian(H, w, dW, dF);
			E += e;
		}
		E = E / fine_samples.size();
		H.triangularView<Lower>() = H.transpose();
		Matrix<float, 6, 1> D = H.fullPivHouseholderQr().solve(G);
		warp.steepest(D);
		if (log != NULL)
			(*log) << E << " ";
		if (iter > 1 && D.segment<3>(3).squaredNorm() < translate_eps && last_E - E < error_eps)
			break;
		last_E = E;
	}
	if (log != NULL)
		(*log) << endl;
	return warp;
}

float FARTracker::evaluate(Warp warp)
{
	far_rect_t rect = window(warp.t);
	float fine_cell = sqrt(rectArea(rect) / cell_n);
	feature.set_cell(fine_cell);
	feature.set_step(1);

	float E = 0.0f;
	for (int i = 0; i < fine_samples.size(); ++i) {		
		Vector2f p = warp.transform2(fine_samples[i]);
		Vector32f F;
		feature.descriptor4(p.x(), p.y(), F.data());		
		E = E + sigmoid((F - fine_model.col(i)).squaredNorm());
	}
	return E / fine_samples.size();
}
