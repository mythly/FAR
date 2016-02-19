#include "fartracker.h"

// C interface

far_tracker_t far_init(const unsigned char *gray, int width, int height, far_rect_t rect)
{
	return new FARTracker(gray, width, height, rect, NULL);
}

far_rect_t far_track(far_tracker_t tracker, const unsigned char *gray)
{
	assert(tracker != NULL);
	FARTracker *t = static_cast<FARTracker*>(tracker);
	return t->track(gray);
}

far_rect_t far_restart(far_tracker_t tracker, const unsigned char *gray, far_rect_t rect)
{
	assert(tracker != NULL);
	FARTracker *t = static_cast<FARTracker*>(tracker);
	return t->restart(gray, rect);
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

far_rect_t rectBound(far_rect_t rect)
{
	rect.width = ceil(rect.x + rect.width) - floor(rect.x);
	rect.x = floor(rect.x);
	rect.height = ceil(rect.y + rect.height) - floor(rect.y);
	rect.y = floor(rect.y);
	return rect;
}

ostream& operator<<(ostream& cout, const far_rect_t&rect)
{
	cout << "[" << int(rect.width) << " x " << int(rect.height);
	cout << " from (" << int(rect.x) << ", " << int(rect.y) << ")]";	
	return cout;
}

Surf::Surf(int width, int height)
{
	W = width;
	H = height;
	C = 0;
	step = 0;	
	sum = MatrixXf((W + 1) * 8, H + 1);	
	hist = MatrixXf(W * 8, H);	
	flag = MatrixXi(W, H);
	zero = Vector8f::Zero();
}

Vector4f Surf::kernel(float angle)
{
	float c = cos(angle), s = sin(angle);
	float wx = c * (1.0f - abs(s));
	float wy = s * (1.0f - abs(c));
	float wu = 0.0f, wv = 0.0f;
	if (c >= 0)
		(s >= 0 ? wu : wv) = c * s;
	else
		(s >= 0 ? wv : wu) = -c * s;
	return Vector4f(wx, wy, wu, wv);
}

void Surf::process(const unsigned char *gray, float angle)
{
	this->angle = angle;
	Vector4f kx = kernel(angle), ky = kernel(angle + PI_2);
	sum.setConstant(0.0f);		
	for (int y = 0; y < H; ++y) {
		int y0 = y > 0 ? y - 1 : y;
		int y1 = y < H - 1 ? y + 1 : y;
		const unsigned char *ptr_y0 = gray + y0 * W;
		const unsigned char *ptr_y = gray + y * W;
		const unsigned char *ptr_y1 = gray + y1 * W;
		for (int x = 0; x < W; ++x) {
			int x0 = x > 0 ? x - 1 : x;
			int x1 = x < W - 1 ? x + 1 : x;
			float gx = float(ptr_y[x1]) - float(ptr_y[x0]);
			float gy = float(ptr_y1[x]) - float(ptr_y0[x]);
			float gu = float(ptr_y1[x1]) - float(ptr_y0[x0]);
			float gv = float(ptr_y1[x0]) - float(ptr_y0[x1]);
			Vector4f g(gx, gy, gu, gv);
			float dx = kx.dot(g), dy = ky.dot(g);
			float *f = sum.data() + ((y + 1) * (W + 1) + x + 1) * 8;
			switch ((dx > 0.0f ? 1 : 0) + (dy > 0.0f ? 2 : 0)) {
			case 0:
				f[0] = dx;
				f[2] = -dx;
				f[4] = dy;
				f[6] = -dy;
				break;
			case 1:
				f[0] = dx;
				f[2] = dx;
				f[5] = dy;
				f[7] = -dy;
				break;
			case 2:
				f[1] = dx;
				f[3] = -dx;
				f[4] = dy;
				f[6] = dy;
				break;
			case 3:
				f[1] = dx;
				f[3] = dx;
				f[5] = dy;
				f[7] = dy;
				break;
			}
		}
	}
	for (int y = 1; y <= H; ++y) {		
		for (int x = 1; x <= W; ++x)
			sum.col(y).segment<8>(x * 8) += sum.col(y).segment<8>((x - 1) * 8);
		sum.col(y) += sum.col(y - 1);
	}	
	C = 0;
	step = 0;
	flag.setConstant(0);
	zero.setConstant(0.0f);
}

void Surf::set_cell(float cell)
{
	cell = cell * 0.5f;
	tx[0] = -cell; ty[0] = -cell;
	tx[1] = cell; ty[1] = -cell;
	tx[2] = cell; ty[2] = cell;
	tx[3] = -cell; ty[3] = cell;
	for (int i = 0; i < 4; ++i) {
		float x = cos(angle) * tx[i] - sin(angle) * ty[i];
		float y = sin(angle) * tx[i] + cos(angle) * ty[i];
		tx[i] = x;
		ty[i] = y;
	}
	C = max(int(floor(cell)), cell_min);
}

void Surf::set_step(int step)
{
	this->step = step;
}

Ref<Vector8f> Surf::cell_hist(int x, int y)
{
	if (x < 0 || x >= W || y < 0 || y >= H)
		return zero;
	if (flag(x, y) != C) {
		int x0 = max(x - C, 0);
		int x1 = min(x + C + 1, W);
		int y0 = max(y - C, 0);
		int y1 = min(y + C + 1, H);	
		Ref<Vector8f> s00 = sum.col(y0).segment<8>(x0 * 8);
		Ref<Vector8f> s01 = sum.col(y0).segment<8>(x1 * 8);
		Ref<Vector8f> s10 = sum.col(y1).segment<8>(x0 * 8);
		Ref<Vector8f> s11 = sum.col(y1).segment<8>(x1 * 8);
		hist.col(y).segment<8>(x * 8) = s11 + s00 - s01 - s10;
		flag(x, y) = C;
	}
	return hist.col(y).segment<8>(x * 8);
}

void Surf::descriptor(float x, float y, Ref<Vector8f> f)
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
	Ref<Vector8f> f00 = cell_hist(ixp * step, iyp * step);
	Ref<Vector8f> f01 = cell_hist(ixp * step, (iyp + 1) * step);
	Ref<Vector8f> f10 = cell_hist((ixp + 1) * step, iyp * step);
	Ref<Vector8f> f11 = cell_hist((ixp + 1) * step, (iyp + 1) * step);
	f = w00 * f00 + w01 * f01 + w10 * f10 + w11 * f11;
}

void Surf::gradient(float x, float y, Ref<Vector8f> f, Ref<Vector8f> dx, Ref<Vector8f> dy)
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
	Ref<Vector8f> f00 = cell_hist(ixp * step, iyp * step);
	Ref<Vector8f> f01 = cell_hist(ixp * step, (iyp + 1) * step);
	Ref<Vector8f> f10 = cell_hist((ixp + 1) * step, iyp * step);
	Ref<Vector8f> f11 = cell_hist((ixp + 1) * step, (iyp + 1) * step);
	f = w00 * f00 + w01 * f01 + w10 * f10 + w11 * f11;
	dx = wy0 * (f10 - f00) + wy1 * (f11 - f01);
	dy = wx0 * (f01 - f00) + wx1 * (f11 - f10);	
}

void Surf::descriptor4(float x, float y, Ref<Vector32f> f)
{
	for (int i = 0; i < 4; ++i)
		descriptor(x + tx[i], y + ty[i], f.segment<8>(i * 8));
	float S = f.squaredNorm();	
	float iS = S < 1.0f ? 0.0f : 1.0f / sqrt(S);
	f *= iS;	
}

void Surf::gradient4(float x, float y, Ref<Vector32f> f, Ref<Vector32f> dx, Ref<Vector32f> dy)
{
	for (int i = 0; i < 4; ++i)
		gradient(x + tx[i], y + ty[i], f.segment<8>(i * 8), dx.segment<8>(i * 8), dy.segment<8>(i * 8));
	float S = f.squaredNorm(), Sx = f.dot(dx), Sy = f.dot(dy);	
	float iS = S < 1.0f ? 0.0f : 1.0f / sqrt(S);
	float iSx = Sx * iS * iS * iS;
	float iSy = Sy * iS * iS * iS;
	dx = iS * dx - iSx * f;	
	dy = iS * dy - iSy * f;	
	f *= iS;
}

Warp::Warp(int width, int height):
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

Matrix<float, 2, 6> Warp::gradient(Vector3f p)
{
	Matrix3f D1 = p.x() * Dx + p.y() * Dy + p.z() * Dz;
	Vector3f tp = transform(p);
	Matrix<float, 2, 3> D2;
	D2 <<  f / tp.z(), 0.0f, -f * tp.x() / (tp.z() * tp.z()), 
		0.0f, f / tp.z(), -f * tp.y() / (tp.z() * tp.z());
	Matrix<float, 2, 3> D3 = D2 * D1;
	Matrix<float, 2, 6> G;

	G << D3(0, 0), D3(0, 1), D3(0, 2), D2(0, 0), D2(0, 1), D2(0, 2),
		D3(1, 0), D3(1, 1), D3(1, 2), D2(1, 0), D2(1, 1), D2(1, 2);
	return G;
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
	vector<Vector3f> candidates;
	candidates.push_back(warp.t);
	candidates.push_back(fast_test(warp));
	
	Warp best_warp(image_width, image_height);
	float best_error = 1.0f;
	for (auto& t : candidates) {
		if (log != NULL)
			(*log) << "track at " << t.transpose() << " " << window(t) << endl;
		Warp w(image_width, image_height);		
		w.setr(warp.r);
		w.sett(t);
		w = fine_test(w);
		float e = evaluate(w);
		if (log != NULL) {
			(*log) << "final translation = " << w.t.transpose() << " " << window(w.t) << endl;
			(*log) << "final rotation = " << w.r.transpose() << endl;
			(*log) << "final error = " << e << endl;
		}
		if (e < best_error) {
			best_warp = w;
			best_error = e;
		}
	}	

	fine_errors.push_back(best_error);
	if (best_error < threshold_error) {		
		warp = best_warp;
		warp.euler(roll, yaw, pitch);
		fine_train(warp);
	}
	else 		
		warp.t = best_warp.t * (warp.t.z() / best_warp.t.z());	
	fast_train(warp);

	return window(best_warp.t);
}

far_rect_t FARTracker::restart(const unsigned char *gray, far_rect_t rect)
{	
	roll = yaw = pitch = 0.0f;
	feature.process(gray, roll);

	if (log != NULL)
		(*log) << "restart at " << locate(rect).transpose() << " " << rect << endl;
	Warp w(image_width, image_height);
	w.sett(locate(rect));
	w = fine_test(w);
	float e = evaluate(w);	
	if (log != NULL) {
		(*log) << "final translation = " << w.t.transpose() << " " << window(w.t) << endl;
		(*log) << "final rotation = " << w.r.transpose() << endl;
		(*log) << "final error = " << e << endl;
	}
	Warp best_warp = w;
	float best_error = e;

	fine_errors.push_back(best_error);
	if (best_error < threshold_error) {
		warp = best_warp;
		warp.euler(roll, yaw, pitch);
		fine_train(warp);
	}
	else
		warp.t = best_warp.t * (warp.t.z() / best_warp.t.z());
	fast_train(warp);

	return window(best_warp.t);
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
		fast_model.col(i) = feature.cell_hist(tx, ty);		
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
		feature.descriptor4(p.x(), p.y(), model.col(i));
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
	for (int y = miny; y <= maxy; y += fast_step)
	for (int x = minx; x <= maxx; x += fast_step) {
		float S = 0.0f, score = 0.0f;
		for (int i = 0; i < fast_samples.size(); ++i) {
			int tx = x + fast_samples[i].x();
			int ty = y + fast_samples[i].y();
			Ref<Vector8f> f = fast_model.col(i);
			Ref<Vector8f> g = feature.cell_hist(tx, ty);
			S += g.squaredNorm();
			score += f.dot(g);
		}
		score *= S < 1.0f ? 0.0f : 1.0f / sqrt(S);
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

Warp FARTracker::Lucas_Kanade(Warp warp)
{	
	for (int iter = 0; iter < max_iteration; ++iter) {
		Matrix<float, 6, 1> G;
		Matrix<float, 6, 6> H;
		G.setConstant(0.0f);
		H.setConstant(0.0f);		
		float E = 0.0f;
		for (int i = 0; i < fine_samples.size(); ++i) {
			Vector32f T(fine_model.col(i)), F;
			Matrix<float, 32, 2> dF;
			Matrix<float, 2, 6> dW = warp.gradient(fine_samples[i]);
			Vector2f p = warp.transform2(fine_samples[i]);
			feature.gradient4(p.x(), p.y(), F, dF.col(0), dF.col(1));	
			T -= F;
			float e = sigmoid(T.dot(T));
			E += e;
			float w = sigmoid_factor * e * (1.0f - e);
			G += w * ((T.transpose() * dF) * dW).transpose();			
			H += w * (dW.transpose() * (dF.transpose() * dF) * dW);
		}
		E = E / fine_samples.size();		
		Matrix<float, 6, 1> D = H.colPivHouseholderQr().solve(G);
		warp.steepest(D);
		if (log != NULL)
			(*log) << E << " ";
		if (iter > 1 && D(3) * D(3) + D(4) * D(4) + D(5) * D(5) < translate_eps) 			
			break;
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
		Vector32f T(fine_model.col(i)), I;
		Vector2f p = warp.transform2(fine_samples[i]);
		feature.descriptor4(p.x(), p.y(), I);
		T -= I;
		E = E + sigmoid(T.dot(T));
	}
	return E / fine_samples.size();
}
