#include "far.h"

Surf::Surf(int width, int height)
{
	W = width;
	H = height;
	C = 0;
	step = 0;	
	sum = MatrixXf((W + 1) * 8, H + 1);	
	hist = MatrixXf(W * 8, H);
	norm = MatrixXf(W, H);
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
	float *f = sum.data();
	for (int y = 0; y < H; ++y) {
		int y0 = y > 0 ? y - 1 : y;
		int y1 = y < H - 1 ? y + 1 : y;
		const unsigned char *ptr_y0 = gray + W * y0;
		const unsigned char *ptr_y = gray + W * y;
		const unsigned char *ptr_y1 = gray + W * y1;
		for (int x = 0; x < W; ++x) {
			int x0 = x > 0 ? x - 1 : x;
			int x1 = x < W - 1 ? x + 1 : x;
			float gx = float(ptr_y[x1]) - float(ptr_y[x0]);
			float gy = float(ptr_y1[x]) - float(ptr_y0[x]);
			float gu = float(ptr_y1[x1]) - float(ptr_y0[x0]);
			float gv = float(ptr_y1[x0]) - float(ptr_y0[x1]);
			Vector4f g(gx, gy, gu, gv);
			float dx = kx.dot(g), dy = ky.dot(g);			
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
			f += 8;
		}
	}
	for (int y = 1; y <= H; ++y) {
		Vector8f s;
		s.setConstant(0.0f);
		for (int x = 1; x <= W; ++x) {			
			s += sum.col(y).segment<8>(x * 8);
			sum.col(y).segment<8>(x * 8) = sum.col(y - 1).segment<8>(x * 8) + s;
		}
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

float* Surf::cell_hist(int x, int y)
{
	if (x < 0 || x >= W || y < 0 || y >= H)
		return zero.data();
	if (flag(x, y) != C) {
		int x0 = max(x - C, 0);
		int x1 = min(x + C + 1, W);
		int y0 = max(y - C, 0);
		int y1 = min(y + C + 1, H);
		float *s00 = &sum(x0 * 8, y0);
		float *s01 = &sum(x0 * 8, y1);
		float *s10 = &sum(x1 * 8, y0);
		float *s11 = &sum(x1 * 8, y1);
		float *f = &hist(x * 8, y);
		float S = 0.0f;
		for (int i = 0; i < 8; ++i) {
			f[i] = s11[i] + s00[i] - s01[i] - s10[i];
			S += f[i] * f[i];
		}
		norm(x, y) = S;
		flag(x, y) = C;
	}
	return &hist(x * 8, y);
}

float Surf::cell_norm(int x, int y)
{
	if (x < 0 || x >= W || y < 0 || y >= H)
		return 0.0f;
	else
		return norm(x, y);
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
	float *f00 = cell_hist(ixp * step, iyp * step);
	float *f01 = cell_hist(ixp * step, (iyp + 1) * step);
	float *f10 = cell_hist((ixp + 1) * step, iyp * step);
	float *f11 = cell_hist((ixp + 1) * step, (iyp + 1) * step);
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
	float *f00 = cell_hist(ixp * step, iyp * step);
	float *f01 = cell_hist(ixp * step, (iyp + 1) * step);
	float *f10 = cell_hist((ixp + 1) * step, iyp * step);
	float *f11 = cell_hist((ixp + 1) * step, (iyp + 1) * step);
	wx0 /= step;
	wy0 /= step;
	wx1 /= step;
	wy1 /= step;
	for (int i = 0; i < 8; ++i) {
		f[i] = f00[i] * w00 + f01[i] * w01 + f10[i] * w10 + f11[i] * w11;
		dx[i] = (f10[i] - f00[i]) * wy0 + (f11[i] - f01[i]) * wy1;
		dy[i] = (f01[i] - f00[i]) * wx0 + (f11[i] - f10[i]) * wx1;
	}
}

void Surf::descriptor4(float x, float y, float *f)
{
	for (int i = 0; i < 4; ++i)
		descriptor(x + tx[i], y + ty[i], f + i * 8);
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
		gradient(x + tx[i], y + ty[i], f + i * 8, dx + i * 8, dy + i * 8);
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
	double theta = sqrt(rx*rx + ry*ry + rz*rz);
	double R[9], J[27];	
	if (theta < 1e-9)
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
		// R = cos(theta)*I + (1 - cos(theta))*r*rT + sin(theta)*[r_x]
		// where [r_x] is [0 -rz ry; rz 0 -rx; -ry rx 0]
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
	this->R << R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8];
	double _J[3][9];
	for (int i = 0; i < 27; ++i)
		_J[i / 9][i % 9] = J[i];	
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
	return Vector2f(p(0) / p(2) * f, p(1) / p(2) * f) + c;
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
	Matrix3f D1 = p(0) * Dx + p(1) * Dy + p(2) * Dz;
	Vector3f tp = transform(p);
	Matrix<float, 2, 3> D2;
	D2 <<  f / tp(2), 0.0f, -f * tp(0) / (tp(2) * tp(2)), 0.0f, f / tp(2), -f * tp(1) / (tp(2) * tp(2));
	Matrix<float, 2, 3> D3 = D2 * D1;
	Matrix<float, 2, 6> G;	
	G << D3(0, 0), D3(0, 1), D3(0, 2), D2(0, 0), D2(0, 1), D2(0, 2),
		D3(1, 0), D3(1, 1), D3(1, 2), D2(1, 0), D2(1, 1), D2(1, 2);
	return G;
}

void Warp::steepest(Matrix<float, 6, 1> parameters)
{
	float rx = r(0) + parameters(0);
	float ry = r(1) + parameters(1);
	float rz = r(2) + parameters(2);
	float tx = t(0) + parameters(3);
	float ty = t(1) + parameters(4);
	float tz = t(2) + parameters(5);
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

FART::FART(const unsigned char *gray, int width, int height, cv_rect_t rect, ostream *os) :
log(os),
image_width(width), image_height(height), 
window_width(rect.width), window_height(rect.height),
feature(width, height),
warp(width, height)
{
	warp.sett(locate(rect));
	float fine_stride = sqrt(window_width * window_height / fine_n);
	int W = int(floor(window_width / (2.0f * fine_stride)));
	int H = int(floor(window_height / (2.0f * fine_stride)));
	for (int y = 0; y <= 2 * H; ++y)
	for (int x = 0; x <= 2 * W; ++x)
		fine_samples.push_back(Vector3f((x - W) * fine_stride, (y - H) * fine_stride, 0.0f));

	feature.process(gray, 0.0f);
	fine_train(warp);
	fast_train(warp);
	error = 0.0f;
	roll = yaw = pitch = 0.0f;
	N = 0;
}

cv_rect_t FART::track(const unsigned char *gray)
{
	if (log != NULL) {
		(*log) << "roll = " << roll * 90.0f / PI_2 << endl;
		(*log) << "yaw = " << yaw * 90.0f / PI_2 << endl;
		(*log) << "pitch = " << pitch * 90.0f / PI_2 << endl;
	}
	feature.process(gray, roll);
	vector<Vector3f> candidates;
	//candidates.push_back(warp.t);
	candidates.push_back(fast_test(warp));

	Warp best_warp(image_width, image_height);
	float best_error = 1.0f;
	for (auto& t : candidates) {
		if (log != NULL)
			(*log) << "candidate " << t.transpose() << " " << window(t) << endl;
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
	candidates.clear();

	error = best_error;
	if (error < threshold_error) {
		warp = best_warp;
		warp.euler(roll, yaw, pitch);
		fine_train(warp);
	}
	else
		warp.t = best_warp.t * (warp.t(2) / best_warp.t(2));
	fast_train(warp);

	return window(best_warp.t);
}

void FART::restart(const unsigned char *gray, cv_rect_t rect)
{
    //to do
}

Vector3f FART::locate(cv_rect_t rect)
{
	float scale = sqrt((window_width * window_height) / (rect.width * rect.height));
	float x = rect.x + rect.width * 0.5f - warp.c(0);
	float y = rect.y + rect.height * 0.5f - warp.c(1);
	return Vector3f(x, y, warp.f) * scale;
}

cv_rect_t FART::window(Vector3f translate)
{
	Vector2f center = warp.project(translate);
	float scale = warp.f / translate(2);
    cv_rect_t ret;
	ret.width = window_width * scale;
	ret.height = window_height * scale;
    ret.x = center(0) - ret.width * 0.5f;
    ret.y = center(1) - ret.height * 0.5f;
    return ret;
}

void FART::fast_train(Warp warp)
{
	cv_rect_t rect = window(warp.t);
	float fast_stride = sqrt(rect.width * rect.height / fast_n);
	feature.set_cell(fast_stride);
	int W = int(rect.width * 0.5f / fast_stride);
	int H = int(rect.height * 0.5f / fast_stride);
	int ox = int(rect.width * 0.5f + 0.5f);
	int oy = int(rect.height * 0.5f + 0.5f);
	int stride = int(round(fast_stride));
	fast_samples.clear();
	for (int y = 0; y <= 2 * H; ++y)
	for (int x = 0; x <= 2 * W; ++x)
		fast_samples.push_back(Vector2i(ox + (x - W) * stride, oy + (y - H) * stride));

	fast_model = MatrixXf::Zero(8, fast_samples.size());
	int x = int(round(rect.x));
	int y = int(round(rect.y));
	for (int i = 0; i < fast_samples.size(); ++i) {
		int tx = x + fast_samples[i](0);
		int ty = y + fast_samples[i](1);
		fast_model.col(i) << Map<Vector8f>(feature.cell_hist(tx, ty));
	}
}

void FART::fine_train(Warp warp)
{
	cv_rect_t rect = window(warp.t);
	float fine_cell = sqrt(rect.width * rect.height / cell_n);
	feature.set_cell(fine_cell);
	feature.set_step(1);

	MatrixXf model(32, fine_samples.size());	
	for (int i = 0; i < fine_samples.size(); ++i) {
		Vector2f p = warp.transform2(fine_samples[i]);
		feature.descriptor4(p(0), p(1), &model(0, i));
	}
	if (fine_model.size() == 0) {
		N = 1;
		fine_model = model;
	}
	else {
		++N;
		fine_model = (float(N - 1) / N) * fine_model + (1.0f / N) * model;
	}
}

Vector3f FART::fast_test(Warp warp)
{
	cv_rect_t rect = window(warp.t);
	float fast_stride = sqrt(rect.width * rect.height / fast_n);
	feature.set_cell(fast_stride);
	cv_rect_t region = window(warp.t / (1.0f + padding));
	int minx = int(max(region.x, -rect.width * 0.5f) + 0.5f);
	int miny = int(max(region.y, -rect.height * 0.5f) + 0.5f);
	int maxx = int(min(region.x + region.width, image_width + rect.width * 0.5f) - rect.width + 0.5f);
	int maxy = int(min(region.y + region.height, image_height + rect.height * 0.5f) - rect.height + 0.5f);

	float best_score = 0.0f;
	Vector3f best_translate = warp.t;
	for (int y = miny; y <= maxy; y += fast_step)
	for (int x = minx; x <= maxx; x += fast_step) {
		float S = 0.0f, score = 0.0f;
		for (int i = 0; i < fast_samples.size(); ++i) {
			int tx = x + fast_samples[i](0);
			int ty = y + fast_samples[i](1);
			Vector8f f = fast_model.col(i);
			Map<Vector8f> g(feature.cell_hist(tx, ty));
			S += feature.cell_norm(tx, ty);
			score += f.dot(g);
		}
		score *= S < 1.0f ? 0.0f : 1.0f / sqrt(S);
		if (score > best_score) {
            cv_rect_t best_rect;
            best_rect.x = x;
            best_rect.y = y;
            best_rect.width = rect.width;
            best_rect.height = rect.height;
            best_translate = locate(best_rect);
			best_score = score;
		}
	}
	return best_translate;
}

Warp FART::fine_test(Warp warp)
{
	cv_rect_t rect = window(warp.t);
	float fine_cell = sqrt(rect.width * rect.height / cell_n);
	feature.set_cell(fine_cell);
	for (auto fine_step : fine_steps) {
		if (fine_step > 2.0f * fine_cell)
			continue;
		feature.set_step(fine_step);
		if (log != NULL)
			(*log) << "\tcell = " << fine_cell << " step = " << fine_step << endl;
		warp = Lucas_Kanade(warp);
	}
	return warp;
}

float FART::sigmoid(float x)
{
	return 1.0f / (1.0f + exp(-sigmoid_factor * (x - sigmoid_bias)));
}

Warp FART::Lucas_Kanade(Warp warp)
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
			Matrix<float, 6, 2> dW = warp.gradient(fine_samples[i]).transpose();
			Vector2f p = warp.transform2(fine_samples[i]);
			feature.gradient4(p(0), p(1), F.data(), dF.col(0).data(), dF.col(1).data());	
			T -= F;
			float e = sigmoid(T.dot(T));
			E += e;
			float w = sigmoid_factor * e * (1.0f - e);
			G += w * (dW * (dF.transpose() * T));
			H += w * (dW * (dF.transpose() * dF) * dW.transpose());
		}
		E = E / fine_samples.size();
		Matrix<float, 6, 1> D = H.colPivHouseholderQr().solve(G);
		warp.steepest(D);
		if (iter > 1 && D(3) * D(3) + D(4) * D(4) + D(5) * D(5) < translate_eps) {
			if (log != NULL)
				(*log) << "\terror in iteration " << iter << " = " << E << endl;
			break;
		}
	}
	return warp;
}

float FART::evaluate(Warp warp)
{
	cv_rect_t rect = window(warp.t);
	float fine_cell = sqrt(rect.width * rect.height / cell_n);
	feature.set_cell(fine_cell);
	feature.set_step(1);

	float E = 0.0f;
	for (int i = 0; i < fine_samples.size(); ++i) {
		Matrix<float, 32, 1> T(fine_model.col(i)), I;
		Vector2f p = warp.transform2(fine_samples[i]);
		feature.descriptor4(p(0), p(1), I.data());
		T -= I;
		E = E + sigmoid(T.dot(T));
	}
	return E / fine_samples.size();
}

ostream& operator<<(ostream& cout, const cv_rect_t &rect)
{
    cout << "[(" << rect.x << ", " << rect.y << ") " << rect.width << "x" << rect.height << "]";
	return cout;
}
