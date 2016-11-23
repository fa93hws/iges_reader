#include "shape.h"
#include <math.h>		//sin,cos,atan2
#include <algorithm>    // std::reverse
#include <typeinfo>


#define PI 3.14159265359
NSI_BEG
void Shape::MatrixTransform(const std::vector<std::vector<double>> &R, const std::vector<double> &T, NS::Point3D &pt) const {
	double x, y, z;
	x = R[0][0] * pt.x() + R[0][1] * pt.y() + R[0][2] * pt.z() + T[0];
	y = R[1][0] * pt.x() + R[1][1] * pt.y() + R[1][2] * pt.z() + T[1];
	z = R[2][0] * pt.x() + R[2][1] * pt.y() + R[2][2] * pt.z() + T[2];
	pt.set(x, y, z);
}

double Shape::linear_interpolation(double seg0, double seg1, double pt) const {
	// parameter space from -1 to 1
	// find parameter for pt1 on [seg0,seg1]
	return (pt * 2 + 2) / (seg1 - seg0) - 1;
}

/***************************************** Shape *******************************************/
void Shape::set_tran_matrix(std::vector<std::vector<double>> &_R, std::vector<double> &_T){
	std::vector<std::vector<double>> R_res; R_res.resize(3);
	for (int i = 0; i < 3; i++)	R_res[i].resize(3);
	std::vector<double> T_res;	T_res.resize(3);
	for (int i = 0; i < 3; i++) {
		T_res[i] = _R[i][0] * T[0] + _R[i][1] * T[1] + _R[i][2] * T[2] +_T[i];
		for (int j = 0; j < 3; j++) {
			R_res[i][j] = _R[i][0] * R[0][j] + _R[i][1] * R[1][j] + _R[i][2] * R[2][j];
		}
	}
	R = R_res;
	T = T_res;
}


/***************************************** Circlular arc *******************************************/
CircularArc::CircularArc(NS::Point3D &pt1, NS::Point3D &pt2, NS::Point3D &pt3) {
	cnt = pt1;
	radius = pt1.distanceTo(pt2);
	theta.resize(2);
	theta[0] = atan2(pt2.y() - pt1.y(), pt2.x() - pt1.x());
	theta[1] = atan2(pt3.y() - pt1.y(), pt3.x() - pt1.x());
	if (pt2.distanceTo(pt3) < 1e-8) {
		_isClosed = true;
		theta[0] = -PI;
		theta[1] = PI;
	}
	if (abs(theta[0] - PI) < 1e-8) 
		theta[0] = -PI;
}
CircularArc::CircularArc(const NS::Point3D &_cnt, const double _radius,const double theta0, const double theta1) {
	cnt = _cnt; 
	theta.push_back(theta0); 
	theta.push_back(theta1); 
	radius = _radius;
}
void CircularArc::get_pt_from_para(const double t, NS::Point3D &pt) const {
	double x, y, z;
	double t1=theta[0], t2=theta[1];
	if (abs(theta[0] - theta[1]) < 1e-10 || theta[0] > theta[1]) {
		t2 += 2 * PI;
	}
	double _theta = (t2-t1)*(t + 1) / 2 +t1;
	x = radius * cos(_theta) + cnt.x();
	y = radius * sin(_theta) + cnt.y();
	pt = NS::Point3D(x, y, cnt[2]);
	MatrixTransform(R, T, pt);
}
void CircularArc::discrete(const int num_pt, std::vector<NS::Point3D> &pts)  {
	pts.clear();
	double para;
	for (int i = 0; i < num_pt; i++) {
		para = 2.0 / (num_pt - 1) * i - 1;
		NS::Point3D pt_tmp;
		get_pt_from_para(para, pt_tmp);
		pts.push_back(pt_tmp);
	}
	discretePts = pts;
	if (_direction == -1)
		std::reverse(pts.begin(), pts.end());
	if (_isClosed) pts.front() = pts.back();
}	
void CircularArc::discrete(const double TOL, std::vector<NS::Point3D> &pts)  {
	if (discretePts.size() > 0) {
		pts = discretePts;
		if (_direction == -1)
			std::reverse(pts.begin(), pts.end());
		return;
	}
	double a = 1.0 - TOL;
	// 2(1-cosx)-a^2x^2 = 0, x not equal to 0
	// cosx = 1 - x^2/2 + x^4/24 - o(x^6)
	// f(x)  = -x^2/12  + 1 - a^2 = 0
	// x = sqrt(12*(1-a^2))
	double dt2 = 1.0 / 12.0 - sqrt(1.0 / 144.0 - (1 - a*a) / 90.0);
	dt2 *= 180;
	double dt = sqrt(dt2);
	double dtheta = theta[1] - theta[0];
	if (dtheta < 0) dtheta += 2 * PI;
	int n = (int)(dtheta / dt + 1);
	assert(n > 0);
	if (_isClosed) n = (int)(2 * PI / dt);
	n = n < 2 ? 2 : n;
	discrete(n, pts);		
}
void CircularArc::set_tran_matrix(std::vector<std::vector<double>> &_R, std::vector<double> &_T) {
	std::vector<NS::Point3D> pts(4);
	pts[0] = cnt;
	get_pt_from_para(-1, pts[1]);
	get_pt_from_para(0, pts[2]);
	get_pt_from_para(1, pts[3]);
	double t00 = theta[0], t01 = theta[1];
	for (NS::Point3D &pt:pts)
		MatrixTransform(_R, _T, pt);
	cnt = pts[0];
	radius = pts[1].distanceTo(pts[0]);

	theta[0] = atan2(pts[1].y() - pts[0].y(), pts[1].x() - pts[0].x());
	theta[1] = atan2(pts[3].y() - pts[0].y(), pts[3].x() - pts[0].x());
	if (_isClosed) return;
	NS::Vector3D v0(pts[1], pts[2]), v1(pts[2], pts[3]);
	NS::Vector3D v = v0.crossProduct(v1);
	if (v.z() < 0) { //clock-wise
		std::swap(theta[0], theta[1]);
		_direction = -_direction;
	}
	if (abs(theta[0] - PI) < 1e-8)
		theta[0] = -PI;
}

bool CircularArc::isConincideWithArc(const CircularArc &arc,const double TOL) const {
	if (cnt.distanceTo(arc.getCnt()) > TOL) return false;
	if (abs(radius) - arc.getRadius() > TOL) return false;
	std::vector<double> tmp; arc.getTheta(tmp);
	double _t[2][2] = { {theta[0],theta[1]}, {tmp[0], tmp[1]} };
	if (_isClosed) { _t[0][0] = -PI; _t[0][1] = PI; }
	if (arc.isClosed()) { _t[1][0] = -PI; _t[1][1] = PI; }
	std::vector<std::vector<double>> t(2, std::vector<double>(2));
	for (int i = 0; i < 2; i++) {
		t[i][0] = _t[i][0]; t[i][1] = _t[i][1];
		if (_t[i][1] < _t[i][0] && _t[i][0] < PI) {
			t[i][1] = -PI; t[i].push_back(PI), t[i].push_back(_t[i][1]);
		}
	}
	for (int i = 0; i < 2; i++) {
		if (t[i].size() == 2) continue;
		for (int j = 0; j < 2; j++) {
			if (abs(t[i][2 * j] - t[i][2 * j + 1]) < TOL) {
				std::vector<double> tmp(2);
				tmp[0] = t[i][2 - 2 * j];
				tmp[1] = t[i][3 - 2 * j];
				t[i] = tmp;
			}
		}
	}
	for (size_t i = 0; i < t[0].size() / 2; i++) {
		for (size_t j = 0; j < t[1].size() / 2; j++) {
			if (t[1][2 * j + 1]>t[0][2 * i] + TOL && t[0][2 * i] - TOL > t[1][2 * j]) return true;
			if (t[1][2 * j + 1]>t[0][2 * i + 1] + TOL && t[0][2 * i + 1] - TOL > t[1][2 * j]) return true;
			if (t[0][2 * i + 1]>t[1][2 * j] + TOL && t[1][2 * j] - TOL > t[0][2 * i]) return true;
			if (t[0][2 * i + 1]>t[1][2 * j + 1] + TOL && t[1][2 * j + 1] - TOL > t[0][2 * i]) return true;
			if (abs(t[0][2 * i] - t[1][2 * j]) < TOL && abs(t[0][2 * i + 1] - t[1][2 * j + 1]) < TOL) return true;
		}
	}
	return false;
} // end of isConincideWithArc
/***************************************** Line *******************************************/	
void Line::get_pt_from_para(const double t, NS::Point3D &pt)  const {
	pt = (pt1 - pt0) / 2 * (t + 1) + pt0;
	MatrixTransform(R, T, pt);
}
void Line::discrete(const int num_pt, std::vector<NS::Point3D> &pts)  {
	pts.clear();
	double para;
	for (int i = 0; i < num_pt; i++) {
		para = 2.0 / (num_pt - 1) * i - 1;
		NS::Point3D pt_tmp;
		get_pt_from_para(para, pt_tmp);
		pts.push_back(pt_tmp);
	}
	if (_direction == -1)
		std::reverse(pts.begin(), pts.end());
}
bool Line::isConincideWithLine(const Line &line,const double TOL) const {
	// is parallel
	std::vector<NS::Point3D> pts(4);
	pts[0] = pt0, pts[1] = pt1; pts[2] = line.pt0; pts[3] = line.pt1;
	NS::Vector3D v0(pts[0], pts[1]);
	NS::Vector3D v1(pts[2], pts[3]);
	NS::Vector3D v3 = v0.crossProduct(v1);
	if (abs(v3.z()) > TOL) return false;
	double d0 = NS::Geometry::PointDistanceToSeg(pts[2], pts[0], pts[1]);
	double d1 = NS::Geometry::PointDistanceToSeg(pts[3], pts[0], pts[1]);
	double d3 = NS::Geometry::PointDistanceToSeg(pts[0], pts[2], pts[3]);
	double d4 = NS::Geometry::PointDistanceToSeg(pts[1], pts[2], pts[3]);
	if (d0 > TOL && d1 > TOL && d3>TOL && d4 > TOL) return false;
	std::vector<std::vector<double>> l(4, std::vector<double>(4));
	for (int m = 0; m<4; m++)
		for (int k = 0; k<4; k++)
			l[m][k] = pts[m].distanceTo(pts[k]);
	if (l[0][2] < TOL && l[1][3] < TOL) return true;
	if (l[1][2] < TOL && l[0][3] < TOL) return true;
	if (l[0][2] < TOL && abs(l[1][3] - l[1][0] - l[3][0]) < TOL) return false;
	if (l[0][3] < TOL && abs(l[1][2] - l[1][0] - l[2][0]) < TOL) return false;
	if (l[1][2] < TOL && abs(l[0][3] - l[1][0] - l[1][3]) < TOL) return false;
	if (l[1][3] < TOL && abs(l[0][2] - l[1][0] - l[1][2]) < TOL) return false;
	return true;
}
void Line::set_tran_matrix(std::vector<std::vector<double>> &_R, std::vector<double> &_T) {
	MatrixTransform(_R, _T, pt0);
	MatrixTransform(_R, _T, pt1);
}

/******************************************* ShapeGroup **********************************************/	
void ShapeGroup::add_line(const Line &line) {
	lines.push_back(line);
	_idx.push_back(0);
}
void ShapeGroup::add_arc(const CircularArc &arc) {
	arcs.push_back(arc);
	_idx.push_back(1);
}
void ShapeGroup::add_nurbs_curves(const IGES_NURBSCurve &curve) {
	curves.push_back(curve);
	_idx.push_back(2);
}
void ShapeGroup::discrete(const int num_pt, std::vector<NS::Point3D> &pts)  {
	std::cout << "not supported for discrete GroupShapes or Surface2D into std::vector<NS::Point3D>" << std::endl;
}
void ShapeGroup::discrete(const double TOL, std::vector<NS::Point3D> &pts)  {
	std::cout << "not supported for discrete GroupShapes or Surface2D into std::vector<NS::Point3D>" << std::endl;
}
void ShapeGroup::discrete(const double TOL, std::vector<std::list<int>> &sequence)  {
	std::vector<std::vector<NS::Point3D>> allShapesPts;
	int i_line = 0, i_arc = 0, i_nurbs_curve = 0;
	std::vector<NS::Point3D> member_pts;
	for (size_t  i = 0; i < _idx.size(); i++) {
		switch (_idx[i]) {
		case 0: {	// lines
			lines[i_line].discrete(2, member_pts);
			i_line++;
			break;
		}
		case 1: {	// arcs
			arcs[i_arc].discrete(TOL, member_pts);
			i_arc++;
			break;
		}
		case 2: {	// nurbs_curve
					//curves[i_nurbs_curve].generate_pts(num_pt + 1, member_pts);
			curves[i_nurbs_curve].discrete(TOL, member_pts);
			i_nurbs_curve++;
			break;
		}
		}
		allShapesPts.push_back(member_pts);
	}
	// connect pts
	std::vector<bool> connected(allShapesPts.size(), false); connected[0] = true;
	sequence.resize(1); sequence[0].push_back(0);
	std::vector<bool> isClosed; isClosed.push_back(false);
	do {
		bool newCurveFlag = true;
		for (size_t  i = 1; i < allShapesPts.size() && newCurveFlag; i++) {
			if (connected[i]) continue;
			for (size_t  j = 0; j < sequence.size() && newCurveFlag; j++) {
				if (isClosed[j]) continue;
				NS::Point3D firstPt = allShapesPts[sequence[j].front()][0];
				NS::Point3D lastPt = allShapesPts[sequence[j].back()].back();
				if (newCurveFlag && allShapesPts[i].back().distanceTo(firstPt) < 1e-8) {
					connected[i] = true;
					sequence[j].push_front(i);
					newCurveFlag = false;
				}
				//if (newCurveFlag && allShapesPts[i].front().distanceTo(firstPt) < 1e-8) {
				//	connected[i] = true;
				//	sequence[j].push_front(i);
				//	newCurveFlag = false;
				//	std::reverse(allShapesPts[i].begin(), allShapesPts[i].end());
				//}
				if (newCurveFlag && (*allShapesPts[i].begin()).distanceTo(lastPt) < 1e-8) {
					connected[i] = true;
					sequence[j].push_back(i);
					newCurveFlag = false;
				}
				//if (newCurveFlag && allShapesPts[i].back().distanceTo(lastPt) < 1e-8) {
				//	connected[i] = true;
				//	sequence[j].push_back(i);
				//	newCurveFlag = false;
				//	std::reverse(allShapesPts[i].begin(), allShapesPts[i].end());
				//}
				if (!newCurveFlag) {
					NS::Point3D pt1 = allShapesPts[sequence[j].front()][0];
					NS::Point3D pt2 = allShapesPts[sequence[j].back()].back();
					if (pt1.distanceTo(pt2) < 1e-8)
						isClosed[j] = true;
					break;
				}
			}
		}
		if (newCurveFlag) {
			bool exitFlag = true;
			for (size_t  i = 0; i < connected.size(); i++) {
				if (!connected[i]) {
					std::list<int> temp; temp.push_back(i);
					sequence.push_back(temp);
					isClosed.push_back(false);
					exitFlag = false;
					connected[i] = true;
					break;
				}
			}
			if (exitFlag)
				break;
		}
	} while (true);
}
void ShapeGroup::discrete(const double TOL, std::vector<std::vector<NS::Point3D>> &pts)  {
	pts.clear();
	std::vector<std::vector<NS::Point3D>> allShapesPts;
	int i_line = 0, i_arc = 0, i_nurbs_curve = 0;
	std::vector<NS::Point3D> member_pts;
	for (size_t  i = 0; i < _idx.size(); i++) {
		switch (_idx[i]) {
			case 0: {	// lines
				lines[i_line].discrete(2, member_pts);
				i_line++;
				break;
			}
			case 1: {	// arcs
				arcs[i_arc].discrete(TOL, member_pts);
				i_arc++;
				break;
			}
			case 2: {	// nurbs_curve
						//curves[i_nurbs_curve].generate_pts(num_pt + 1, member_pts);
				curves[i_nurbs_curve].discrete(TOL, member_pts);
				i_nurbs_curve++;
				break;
			}
		}
		allShapesPts.push_back(member_pts);
	}
	// connect pts
	std::vector<std::list<int>> sequence;
	discrete(TOL, sequence);
	for (size_t  i = 0; i < sequence.size(); i++) {
		std::vector<NS::Point3D> pts_row;
		for (auto k = sequence[i].begin(); k != sequence[i].end(); k++) {
			for (size_t  j = 0; j < allShapesPts[*k].size()-1; j++)
				pts_row.push_back(allShapesPts[*k][j]);
			pts_row.push_back(allShapesPts[*k].back());				
		}
		pts.push_back(pts_row);
	}
}
void ShapeGroup::get_pt_from_para(const double t, NS::Point3D &pt) const {

}
void ShapeGroup::find_dir(const NS::Point3D &pt ) {
	int direction = -1;
	return;
}
void ShapeGroup::discreteByIdx(const int idx, const double TOL,  std::vector<NS::Point3D> &pts, std::string &layer, std::vector<double> &color)  {
	pts.clear();
	ASSERT((size_t)idx <= _idx.size());
	int typ = _idx[idx];
	int seq = 0;
	for (int i = 0; i != idx; i++) 
		if (_idx[idx] == _idx[i])
			seq++;
	switch (typ) {
		case 0: {	// lines
			lines[seq].discrete(2, pts);
			layer = lines[seq].get_layer();
			lines[seq].get_color(color);
			return;
		}
		case 1: {	// arcs
			arcs[seq].discrete(TOL, pts);
			layer = arcs[seq].get_layer();
			arcs[seq].get_color(color);
			return;
		}
		case 2: {	// nurbs_curve
					//curves[i_nurbs_curve].generate_pts(num_pt + 1, member_pts);
			curves[seq].discrete(TOL, pts);
			layer = curves[seq].get_layer();
			curves[seq].get_color(color);
			return;
		}
		default:
			ASSERT(0);
	}
}
/******************************************* 2D-surface-143 **********************************************/
void Surface2D::discrete(const double TOL, std::vector<std::vector<NS::Point3D>> &pts) {
	pts.clear();
	int i_line = 0, i_arc = 0, i_nurbs_curve = 0;
	std::vector<NS::Point3D> member_pts;
	for (size_t i = 0; i < _idx.size(); i++) {
		switch (_idx[i]) {
			case 0: {	// lines
				lines[i_line].discrete(2, member_pts);
				i_line++;
				break;
			}
			case 1: {	// arcs
				arcs[i_arc].discrete(TOL, member_pts);
				i_arc++;
				break;
			}
			case 2: {	// nurbs_curve
						//curves[i_nurbs_curve].generate_pts(num_pt + 1, member_pts);
				curves[i_nurbs_curve].discrete(TOL, member_pts);
				i_nurbs_curve++;
				break;
			}
		}
		pts.push_back(member_pts);
	}
}
void Surface2D::add_sg(const ShapeGroup &sg) {
	_sg.push_back(sg);
	_idx.push_back(3);
}
/***************************************** IGES_NURBSCurve *******************************************/
//bool _sort(int i, int j) { return i < j; };

void IGES_NURBSCurve::discrete(const int num_pt, std::vector<NS::Point3D> &pts)  {
	pts.clear();
	double para;
	for (int i = 0; i < num_pt; i++) {
		para = 2.0 / (num_pt - 1) * i - 1;
		NS::Point3D pt_tmp;
		get_pt_from_para(para, pt_tmp);
		pts.push_back(pt_tmp);
	}
	discretePts = pts;
	if (_direction == -1)
		std::reverse(pts.begin(), pts.end());
}
void IGES_NURBSCurve::get_pt_from_para(const double t, NS::Point3D &pt)  const {
	c.getCoord(t, pt);
}

double IGES_NURBSCurve::polyline_length(const NS::Point3D &pt0, const NS::Point3D &pt1, const NS::Point3D &pt2) const {
	return pt0.distanceTo(pt1) + pt1.distanceTo(pt2);
}
double IGES_NURBSCurve::arc_len_pt_loc_convex(const double u0, const double u1) const {
	double TOL = 0.025;
	double un = (u0 + u1) / 2, uo;
	std::vector<NS::Point3D> pts(3);
	c.getCoord(u0, pts[1]);
	c.getCoord(u1, pts[2]);
	//if (pts[2].distanceTo(NS::Point3D(-104.4457, -56.6572,0.0)) < 1e-3)
	//	auto sss = 1;
	double dif;
	do {
		uo = un;
		double dp = 0, dpp = 0;
		std::vector<NS::Vector3D> temp;
		c.ratCurveDerives(2, uo, temp);
		pts[0] = NS::Point3D(temp[0].x(), temp[0].y(), temp[0].z());
		NS::Vector3D fp = temp[1];
		NS::Vector3D fpp = temp[2];
		for (int i = 0; i < 2; i++) {
			NS::Vector3D f(pts[i + 1], pts[0]);
			double sum = f.dotProduct(fp);
			dp += sum / f.length();
			dpp -= sum * sum * pow(f.lengthSqrd() , -1.5);
			dpp += (fp.lengthSqrd() + f.dotProduct(fpp)) / f.length();
		}
		dif = dp / dpp;
		un = uo - dif;
		if (un < u0 || un > u1)
			return uo;
	} while (dif / (u1 - u0) > TOL);
	return un;
}
double IGES_NURBSCurve::arc_len_pt_loc_concave1(const double u0, const double u1,const double ui) const {
	assert(c.basis.order > 2);
	double TOL = 0.025;
	double un = ui, uo;
	std::vector<NS::Point3D> pts(3);
	c.getCoord(u0, pts[1]);
	c.getCoord(u1, pts[2]);
	int ntime = 0;
	double dif;
	do {
		uo = un;
		double dpp = 0, dppp = 0;
		std::vector<NS::Vector3D> temp;
		c.ratCurveDerives(3, uo, temp);		
		pts[0] = NS::Point3D(temp[0].x(), temp[0].y(), temp[0].z());
		NS::Vector3D fp = temp[1];
		NS::Vector3D fpp = temp[2];
		NS::Vector3D fppp = temp[3];
		for (int i = 0; i < 2; i++) {
			pts[0] = NS::Point3D(temp[0].x(), temp[0].y(), temp[0].z());
			NS::Vector3D f(pts[i + 1], pts[0]);
			double fdfp  = f.dotProduct(fp);
			double fdfpp = f.dotProduct(fpp);
			dpp  -= pow(fdfp,2) * pow(f.lengthSqrd(), -1.5);
			dpp  += (fp.lengthSqrd() + f.dotProduct(fpp)) / f.length();

			dppp += 3 * pow(fdfp, 3) / pow(f.lengthSqrd(), 2.5);
			dppp -= 3 * fdfp * (fdfpp + fp.lengthSqrd()) / pow(f.lengthSqrd(), 1.5);
			dppp += (3 * fp.dotProduct(fpp) + f.dotProduct(fppp)) / f.length();
		}
		dif = dpp / dppp;
		un = uo - dif;
		if (un > u1) un = (u1*10+u0)/11;
		if (un<u0) un = (u0*10+u1)/11;
		ntime++;
		if (ntime > 20) return (u0 + u1) / 2;
	} while (abs(dif) / (u1 - u0) > TOL);
	//debug
	NS::Point3D debugPt;
	c.getCoord(un, debugPt);
	//debug
	return un;
}
double IGES_NURBSCurve::arc_len_pt_loc_general(const double u0, const double u1,const double arc_len) const {
	double TOL = 0.05;
	double u_up = u1;
	double u_down = u0;
	double u,u_mu,u_md,dist_mu,dist_md,dist;
	NS::Point3D pt_up, pt_down,pt_mu,pt_md,pt,pt0,pt1;
	c.getCoord(u_up, pt0);
	c.getCoord(u_down, pt1);
	do {
		u = (u_up + u_down) / 2;
		u_mu = (u * 10 + u_up) / 11.0;
		u_md = (u * 10 + u_down) / 11.0;
		c.getCoord(u_up,pt_up);
		c.getCoord(u_down,pt_down);
		c.getCoord(u_mu, pt_mu);
		c.getCoord(u_md, pt_md);
		c.getCoord(u, pt);
		dist_mu = polyline_length(pt0, pt_mu, pt1);
		dist_md = polyline_length(pt0, pt_md, pt1);
		dist    = polyline_length(pt0, pt, pt1);
		if (dist_mu > dist)
			u_down = u;
		else if (dist_md > dist)
			u_up = u;
		else
			break;
	} while (u_up-u_down > TOL*2*(u1-u0));
	return u;
}
double IGES_NURBSCurve::arc_len_pt_loc(const double u0, const double u1, const double arc_len) const {
	assert(false);
	return 0.0;
	//int mode = 0;
	//int idx = -1;
	//for (int i = 0; i < c.controlPoints.size() - 2 && mode < 2; i++) {
	//	std::vector<NS::Point3D> pts(c.controlPoints.begin() + i, c.controlPoints.end() + i + 4);
	//	std::vector<NS::Vector3D> vecs(3);
	//	for (int j = 0; j < 3; j++) vecs[j] = NS::Vector3D(pts[j], pts[j + 1]);
	//	NS::Vector3D v0 = vecs[0].crossProduct(vecs[1]);
	//	NS::Vector3D v1 = vecs[1].crossProduct(vecs[2]);
	//	if (v0*v1 < 0) {
	//		if (idx == -1) idx = i;
	//		else mode = 2;
	//	}
	//}
	//switch (mode) {
	//case 0:
	//	return arc_len_pt_loc_convex(u0, u1);
	//case 1: {
	//	double ui = 0;
	//	return arc_len_pt_loc_concave1(u0, u1,ui);
	//}
	//case 2:
	//	return arc_len_pt_loc_general(u0, u1, arc_len);
	//}
}
	
	
void IGES_NURBSCurve::discrete(const double TOL, std::vector<NS::Point3D> &pts)  {
	if (discretePts.size() > 0) {
		pts = discretePts;
		if (_direction == -1)
			std::reverse(pts.begin(), pts.end());
		return;
	}
	pts.clear();
	assert(c.basis.order < 4);
	int mode = 0;
	int idx = -1;
	
	for (size_t i = 0;mode < 2; i++) {
		if (c.controlPoints.size() <= 3) break;
		std::vector<NS::Point3D> _pts(c.controlPoints.begin() + i, c.controlPoints.begin() + i + 4);
		std::vector<NS::Vector3D> vecs(3);
		for (int j = 0; j < 3; j++) vecs[j] = NS::Vector3D(_pts[j], _pts[j + 1]);
		NS::Vector3D v0 = vecs[0].crossProduct(vecs[1]);
		NS::Vector3D v1 = vecs[1].crossProduct(vecs[2]);
		if (v0.z()*v1.z() < 0) {
			if (idx == -1) { idx = i, mode = 1; }
			else mode = 2;
		}
	}
	assert(mode < 2);
	std::vector<std::vector<double>> k;
	if (mode == 1) {
		double uc = arc_len_pt_loc_concave1(-1.0, 1.0, 0);
		k.push_back(std::vector<double>(2));
		k.push_back(std::vector<double>(2));
		k[0][0] = -1.0; k[0][1] = uc; k[1][0] = uc; k[1][1] = 1.0;
	}
	else if (mode == 0) {
		k.push_back(std::vector<double>(2)); 
		k[0][0] = -1.0; k[0][1] = 1.0;
	}
	int maxDiv = 5;
	//std::vector<double> k_int;
	//for (size_t  i = 0; i < c.basis.knotVector.size(); i++) {
	//	if (abs(c.basis.knotVector[i] + 1) < 1e-6)
	//		continue;
	//	if (abs(c.basis.knotVector[i] - 1) < 1e-6)
	//		break;
	//	bool multiFlag = false;
	//	for (size_t  j = 0; j < k_int.size(); j++)
	//		if (abs(c.basis.knotVector[i] - k_int[j]) < 1e-6) {
	//			multiFlag = true;
	//			break;
	//		}
	//	if (multiFlag) continue;
	//	k_int.push_back(c.basis.knotVector[i]);
	//}
	//std::vector<std::vector<double>> k(k_int.size()+1);		
	//for (size_t  i = 0; i < k.size(); i++) {
	//	k[i].resize(2);
	//	k[i][0] = i == 0           ? -1 : k_int[i - 1];
	//	k[i][1] = i == k.size() - 1 ? 1 : k_int[i];
	//}
	std::vector<bool> pass(k.size(), false);
	bool need_refine = true;
	int div = 0;
	do {
		need_refine = false;
		int k_siz = k.size();
		for (int i = 0; i < k_siz; i++) {
			if (pass[i]) continue;
			NS::Point3D pt0, pt1;
			c.getCoord(k[i][0], pt0);
			c.getCoord(k[i][1], pt1);
			double dist			= pt0.distanceTo(pt1);
			double arc_length	= c.arc_length(k[i][0], k[i][1]);
			if (abs(1 - dist / arc_length) < TOL)
				pass[i] = true;
			else {
				need_refine = true;
				double u = arc_len_pt_loc_convex(k[i][0], k[i][1]);
				std::vector<double> temp; temp.push_back(u); temp.push_back(k[i][1]);
				k.push_back(temp);
				k[i][1] = u;
				pass.push_back(false);
			}
		}
		div++;
	} while (need_refine && div < maxDiv);

	std::vector<double> k_tmp;
	for (std::vector<double> tmp : k)
		k_tmp.push_back(tmp[0]);
	k_tmp.push_back(1.0);
	std::sort(k_tmp.begin(), k_tmp.end());
	for (double u : k_tmp) {
		NS::Point3D tmp; c.getCoord(u, tmp);
		pts.push_back(tmp);
	}
	discretePts = pts;
	if (_direction == -1)
		std::reverse(discretePts.begin(), discretePts.end());
	return;
}
int IGES_NURBSCurve::isConincideWith(const IGES_NURBSCurve &_cc, const double TOL) const{
	NURBS::Curve _c;
	_cc.get_curve(_c);
	bool conincide = true;
	if (_c.basis.knotVector.size() != c.basis.knotVector.size()) return 0;
	for (size_t i = 0; i < _c.basis.knotVector.size() && conincide; i++)
		if (c.basis.knotVector[i] != _c.basis.knotVector[i])
			conincide = false;
	for (size_t i = 0; i < c.controlPoints.size() && conincide; i++) {
		if (c.controlPoints[i].distanceTo(_c.controlPoints[i]) > TOL)
			conincide = false;
		if (abs(c.weightVector[i] - _c.weightVector[i])>1e-8)
			conincide = false;
	}
	if (conincide) return 1;

	conincide = true;
	_c.reverse();
	for (size_t i = 0; i < _c.basis.knotVector.size() && conincide; i++)
		if (c.basis.knotVector[i] != _c.basis.knotVector[i])
			conincide = false;
	for (size_t i = 0; i < c.controlPoints.size() && conincide; i++) {
		double d = c.controlPoints[i].distanceTo(_c.controlPoints[i]);
		if ( d> TOL)
			conincide = false;
		if (abs(c.weightVector[i] - _c.weightVector[i])>1e-8)
			conincide = false;
	}
	if (conincide) return -1;
	else return 0;
}
/***************************************** Shape-Manger2DSurfs *******************************************/
void ShapeManager2DSurfs::addLine(const int sid, const Line &l,const int dir) {
	_surfs[sid].push_back(Member(_lines.size(), 0, dir));
	_lines.push_back(l);		
}
void ShapeManager2DSurfs::addArc(const int sid, CircularArc &c, const int dir) {
	if (dir == c.get_dir())
		_surfs[sid].push_back(Member(_arcs.size(), 1, 1));
	else
		_surfs[sid].push_back(Member(_arcs.size(), 1, -1));
	c.set_dir(1);
	_arcs.push_back(c);
}
void ShapeManager2DSurfs::addNURBSCurve(const int sid, const IGES_NURBSCurve &c, const int dir) {
	NURBS::Curve temp;
	c.get_curve(temp);
	NURBS::ParentCurve curve(temp);
	if (dir == -1)	curve.reverse();
	curve.subDivision();
	std::vector<NURBS::Curve> children;
	curve.getChildren(children);
	for (NURBS::Curve child : children) {
		_surfs[sid].push_back(Member(_curves.size(), 2, 1));
		_curves.push_back(IGES_NURBSCurve(child));
	}
	//NURBS::Curve curve;
	//c.get_curve(curve);
	//double dist = curve.arc_length(-1, 1);
	//std::cout << dist << std::endl;
}
void ShapeManager2DSurfs::discrete(const int i, const int j, const double TOL, std::vector<NS::Point3D> &pts)  {
	Member m = _surfs[i][j];
	switch (m._typ) {
	case 0:
		_lines[m._idx].discrete(TOL, pts);
		break;
	case 1:
		_arcs[m._idx].discrete(TOL, pts);
		break;
	case 2:
		_curves[m._idx].discrete(TOL, pts);
		break;
	default:
		std::cout << "typ not supported in ShapeManager2DSurfs";
	}
}
void ShapeManager2DSurfs::discrete(const int i, const double TOL, std::vector<std::vector<NS::Point3D>> &pts)  {
	pts.clear();
	for (size_t j = 0; j < _surfs[i].size(); j++) {
		std::vector<NS::Point3D> temp;
		discrete(i, j, TOL, temp);
		pts.push_back(temp);
	}
}
		/**********Remove Conincide Objects**********/
void ShapeManager2DSurfs::removeConincideLines() { // O(n^2)
	if (_lines.size() == 0) return;
	double TOL = 1e-8;
	// Uniformlize
	std::vector<NS::Point3D> tmp;
	_lines[0].discrete(2, tmp);
	double xmin = tmp[0].x(), xmax = tmp[0].x();
	double ymin = tmp[0].y(), ymax = tmp[0].y();
	for (size_t i = 0; i < _lines.size(); i++) {
		_lines[0].discrete(2, tmp);
		for (NS::Point3D pt : tmp) {
			if (pt.x() > xmax) xmax = pt.x();
			if (pt.x() < xmin) xmin = pt.x();
			if (pt.y() > ymax) ymax = pt.y();
			if (pt.y() < ymin) ymin = pt.y();
		}
	}
	double scale = xmax - xmin > ymax - ymin ? (xmax - xmin) / 10.0 : (ymax - ymin) / 10.0;
	TOL *= scale;
	//
	std::vector<int> idx2group(_lines.size(),-1);
	std::vector<std::vector<int>> group;
	for (size_t i = 0; i < _lines.size() - 1; i++) { // O(n^2)
		for (size_t j = i + 1; j < _lines.size(); j++) {
			if (!_lines[i].isConincideWithLine(_lines[j],TOL))
				continue;
			if (idx2group[i] == -1 && idx2group[j] == -1) {
				idx2group[i] = group.size();
				idx2group[j] = group.size();
				group.push_back(std::vector<int>(2));
				group.back().front() = i; group.back().back() = j;
			}
			else if (idx2group[i] == -1) {
				group[idx2group[j]].push_back(i);
				idx2group[i] = idx2group[j];
			}
			else if (idx2group[j] == -1) {
				group[idx2group[i]].push_back(j);
				idx2group[j] = idx2group[i];
			}
		}
	}
	if (group.size() == 0) return;
	// split lines in the linegroup
	std::vector<std::vector<int>> lineModify(_lines.size());
	std::vector<int> dir(_lines.size(), 1);
	std::vector<bool> removeidx(_lines.size());
	int lineSize = _lines.size();
	for (auto i = group.cbegin(); i != group.cend(); i++) { // O(n^2)
		std::vector<NS::Point3D> keyPts;
		std::vector<std::vector<int>> linePtIdx((*i).size(),std::vector<int>(2,-1));
		// group all points
		for (size_t j = 0; j < (*i).size(); j++) {
			assert(lineModify[(*i)[j]].size() == 0);
			_lines[(*i)[j]].discrete(2, tmp);
			for (int k = 0; k < 2; k++) {
				int ptidx = -1;
				for (size_t l = 0; l < keyPts.size(); l++) {
					if (keyPts[l].distanceTo(tmp[k]) < TOL) {
						ptidx = l;
						break;
					}
				}
				if (ptidx > -1)
					linePtIdx[j][k] = ptidx;
				else {
					linePtIdx[j][k] = keyPts.size();
					keyPts.push_back(tmp[k]);
				}
			}
		}
		for (auto j = (*i).cbegin(); j != (*i).cend(); j++) removeidx[*j] = true;
		//std::vector<int> pointsGroup;
		//for (size_t j = 0; j < keyPts.size(); j++)
		//	pointsGroup.push_back(j);
		// find two furthest points
		double maxD = 0;
		std::vector<int> idx(keyPts.size(), -1);
		for (auto j = keyPts.cbegin(); j != keyPts.cend() - 1; j++) {
			for (auto k = j + 1; k != keyPts.cend(); k++) {
				double dist = (*j).distanceTo((*k));
				if (dist > maxD) {
					maxD = dist;
					idx.front() = j - keyPts.cbegin();
					idx.back() = k - keyPts.cbegin();
				}
			}
		}
		// sort points based on distance
		std::vector<double> intDist;
		std::vector<int> intidx;
		for (size_t j = 0; j < keyPts.size(); j++) intidx.push_back(j);
		for (auto j = keyPts.cbegin(); j != keyPts.cend(); j++)
			intDist.push_back((*j).distanceTo(keyPts[idx[0]] ));
		std::sort(intidx.begin(), intidx.end(), [&intDist](int i, int j) {return intDist[i] < intDist[j]; });
		for (size_t j = 1; j < keyPts.size() - 1; j++) idx[j] = intidx[j];
		// build lineModify and reverse
		std::vector<int> lineAddedTo(keyPts.size() - 1, -1);
		for (size_t j = 0; j < (*i).size(); j++) {
			int pb = std::find(idx.begin(), idx.end(), linePtIdx[j][0]) - idx.begin();
			int pe = std::find(idx.begin(), idx.end(), linePtIdx[j][1]) - idx.begin();
			assert(pb != pe);
			if (pb > pe) { dir[(*i)[j]] = -1; std::swap(pb, pe); }
			for (int k = pb; k < pe; k++) {
				lineModify[(*i)[j]].push_back(_lines.size() + k);
			}
		}
		for (auto j = idx.cbegin(); j != idx.cend() - 1; j++) 
			_lines.push_back( Line(keyPts[*j],keyPts[*(j+1)]) );
	}
	// finish find lineModify and reverse
	// modify _surface and _reversePolyLines
	for (size_t i = 0; i < _surfs.size(); i++) {
		for (size_t j = 0; j < _surfs[i].size(); j++) {
			if (_surfs[i][j]._typ != 0) continue;
			int lineid = _surfs[i][j]._idx;
			if (lineid >= lineSize) continue;
			if (removeidx[_surfs[i][j]._idx]) {
				bool isReversed = _surfs[i][j]._dir != dir[lineid];
				_surfs[i][j]._dir = isReversed? -1 : 1;
				_surfs[i][j]._idx = isReversed ? lineModify[lineid].back() : lineModify[lineid].front();
				for (size_t k = 1; k < lineModify[lineid].size(); k++) {
					if (isReversed)
						_surfs[i].insert(_surfs[i].begin() + j + k, Member(*(lineModify[lineid].end() - k - 1), 0, _surfs[i][j]._dir));
					else
						_surfs[i].insert(_surfs[i].begin() + j + k, Member(lineModify[lineid][k], 0, _surfs[i][j]._dir));
				}
			}
		}
	}
	// shrink _polylines and _surfaces to fit
	int count = 0;
	std::vector<int> shrinkPolyIdx(_lines.size(), -1);
	std::vector<Line> newLines;
	for (size_t i = 0; i < _lines.size(); i++)
		if (i >= removeidx.size() || !removeidx[i]) { shrinkPolyIdx[i] = count++; newLines.push_back(_lines[i]); }
	_lines = newLines;
	// surface
	for (auto i = _surfs.begin(); i != _surfs.end(); i++)
		for (auto j = (*i).begin(); j != (*i).end(); j++) {
			if ((*j)._typ != 0) continue;
			(*j)._idx = shrinkPolyIdx[(*j)._idx];
		}
}//removeConincideLines_end
void ShapeManager2DSurfs::removeConincideArcs() {
	if (_arcs.size() == 0) return;
	double TOL = 1e-8;
	// Uniformlize
	double xmin = _arcs[0].getCnt().x() , xmax = xmin;
	double ymin = _arcs[0].getCnt().y(), ymax = ymin;
	for (size_t i = 0; i < _arcs.size(); i++) {
		if (_arcs[i].getCnt().x() + _arcs[i].getRadius() > xmax) xmax = _arcs[i].getCnt().x() + _arcs[i].getRadius();
		if (_arcs[i].getCnt().x() - _arcs[i].getRadius() < xmin) xmin = _arcs[i].getCnt().x() - _arcs[i].getRadius();
		if (_arcs[i].getCnt().y() + _arcs[i].getRadius() > ymax) ymax = _arcs[i].getCnt().y() + _arcs[i].getRadius();
		if (_arcs[i].getCnt().y() - _arcs[i].getRadius() > ymin) ymin = _arcs[i].getCnt().y() - _arcs[i].getRadius();
	}
	double scale = xmax - xmin > ymax - ymin ? (xmax - xmin) / 10.0 : (ymax - ymin) / 10.0;
	TOL *= scale;
	//
	std::vector<int> idx2group(_arcs.size(), -1);
	std::vector<std::vector<int>> group;
	for (size_t i = 0; i < _arcs.size() - 1; i++) { // O(n^2)
		for (size_t j = i + 1; j < _arcs.size(); j++) {
			if (!_arcs[i].isConincideWithArc(_arcs[j],TOL)) continue;
			if (idx2group[i] == -1 && idx2group[j] == -1) {
				idx2group[i] = group.size();
				idx2group[j] = group.size();
				group.push_back(std::vector<int>(2));
				group.back().front() = i; group.back().back() = j;
			}
			else if (idx2group[i] == -1) {
				group[idx2group[j]].push_back(i);
				idx2group[i] = idx2group[j];
			}
			else if (idx2group[j] == -1) {
				group[idx2group[i]].push_back(j);
				idx2group[j] = idx2group[i];
			}
		}
	}
	if (group.size() == 0) return;
	std::vector<std::vector<int>> arcModify(_arcs.size());
	std::vector<bool> removeidx(_arcs.size());
	int arcSize = _arcs.size();
	for (auto i = group.cbegin(); i != group.cend(); i++) { // O(n^2)
		NS::Point3D cnt = _arcs[(*i)[0]].getCnt();
		double radius = _arcs[(*i)[0]].getRadius();
		std::vector<double> thetas;
		//thetas.push_back(-PI); thetas.push_back(PI);
		std::vector<std::vector<int>> arcAngleIdx((*i).size(), std::vector<int>(2, -1));
		// group all angles
		for (size_t j = 0; j < (*i).size(); j++) {
			std::vector<double> angles;
			_arcs[(*i)[j]].getTheta(angles);
			for (int k = 0; k < 2; k++) {
				int ptidx = -1;
				for (size_t l = 0; l < thetas.size(); l++) {
					if (abs(thetas[l] - angles[k]) < 1e-8) {
						ptidx = l;
						break;
					}
				}
				if (ptidx > -1)
					arcAngleIdx[j][k] = ptidx;
				else {
					arcAngleIdx[j][k] = thetas.size();
					thetas.push_back(angles[k]);
				}
			}
		}
		for (auto j = (*i).cbegin(); j != (*i).cend(); j++) {
			assert(removeidx[*j] == false);
			removeidx[*j] = true; 
		}

		int pidx = -1, ngpidx = -1;
		bool isPiExist = false;
		bool isAssistedPiUsed = false;
		for (size_t j = 0; j < thetas.size(); j++) {
			if (pidx == -1 && abs(thetas[j] - PI) < 1e-8) { pidx = j; isPiExist = true; }
			if (ngpidx == -1 && abs(thetas[j] + PI) < 1e-8) { ngpidx = j; isPiExist = true;}
		}
		if (ngpidx == -1) { ngpidx = thetas.size(); thetas.push_back(-PI); }
		if (pidx == -1) { pidx = thetas.size(); thetas.push_back(PI); }			
		for (size_t j = 0; j < arcAngleIdx.size(); j++) {
			if (thetas[arcAngleIdx[j][1]] < thetas[arcAngleIdx[j][0]] && thetas[arcAngleIdx[j][0]] < PI) {
				arcAngleIdx[j].push_back(pidx);
				arcAngleIdx[j].push_back(arcAngleIdx[j][1]);
				arcAngleIdx[j][1] = ngpidx;
			}
		}
		// find two furthest points
		std::vector<int> idx(thetas.size(), -1);
		idx[0] = ngpidx; idx.back() = pidx;
		std::vector<double> intDist;
		std::vector<int> intidx;
		for (size_t j = 0; j < thetas.size(); j++) intidx.push_back(j);
		for (auto j = thetas.cbegin(); j != thetas.cend(); j++)
			intDist.push_back((*j) - thetas[idx[0]]);
		std::sort(intidx.begin(), intidx.end(), [&intDist](int i, int j) {return intDist[i] < intDist[j]; });
		for (size_t j = 1; j < thetas.size() - 1; j++) idx[j] = intidx[j];
		// build lineModify and reverse
		std::vector<int> arcAddedTo(thetas.size() - 1, -1);
		for (size_t j = 0; j < (*i).size(); j++) {
			for (size_t k = 0; k < arcAngleIdx[j].size() / 2; k++) {
				if (k > 0) isAssistedPiUsed = true;
				int pb = std::find(idx.begin(), idx.end(), arcAngleIdx[j][k*2]) - idx.begin();
				int pe = std::find(idx.begin(), idx.end(), arcAngleIdx[j][k*2+1]) - idx.begin();
				assert(pb != pe);
				for (int k = pb; k < pe; k++) {
					if (isPiExist && k == idx.size() - 1) continue;
					arcModify[(*i)[j]].push_back(_arcs.size() + k);
				}
			}
		}
		bool used0 = false, used1 = false;
		for (size_t j = 0; j < (*i).size() && !used0; j++)
			for (int tmp : arcModify[(*i)[j]])
				if (tmp == _arcs.size())
					used0 = true;
		if (!used0) {
			for (auto j = idx.begin(); j != idx.end() - 1; j++)	*j = *(j + 1);
			idx.pop_back();
			for (size_t j = 0; j < (*i).size() && !used0; j++)
				for (int &tmp : arcModify[(*i)[j]])
					tmp--;
		}
		for (size_t j = 0; j < (*i).size() && !used1; j++)
			for (int tmp : arcModify[(*i)[j]])
				if (tmp == _arcs.size() + idx.size()-2)
					used1 = true;
		if (!used1) idx.pop_back();
		if (!isPiExist && isAssistedPiUsed) {
			assert(used0);
			assert(used1);
			idx[0] = *(idx.end() - 2); idx.pop_back();
		}
		for (auto j = idx.cbegin(); j != idx.cend() - 1; j++) 
			_arcs.push_back(CircularArc(cnt,radius, thetas[*j], thetas[*(j + 1)]));
	}
	// finish find lineModify and reverse
	// modify _surface and _reversePolyLines
	for (size_t i = 0; i < _surfs.size(); i++) {
		for (size_t j = 0; j < _surfs[i].size(); j++) {
			if (_surfs[i][j]._typ != 1) continue;
			int arcid = _surfs[i][j]._idx;
			if (arcid >= arcSize) continue;
			if (removeidx[_surfs[i][j]._idx]) {
				bool isReversed = _surfs[i][j]._dir == -1;
				_surfs[i][j]._idx = isReversed ? arcModify[arcid].back() : arcModify[arcid].front();
				for (size_t k = 1; k < arcModify[arcid].size(); k++) {
					if (isReversed)
						_surfs[i].insert(_surfs[i].begin() + j + k, Member(*(arcModify[arcid].end() - k - 1), 1, _surfs[i][j]._dir));
					else
						_surfs[i].insert(_surfs[i].begin() + j + k, Member(arcModify[arcid][k], 1, _surfs[i][j]._dir));
				}
			}
		}
	}
	// shrink _polylines and _surfaces to fit
	int count = 0;
	std::vector<int> shrinkPolyIdx(_arcs.size(), -1);
	std::vector<CircularArc> newArcs;
	for (size_t i = 0; i < _arcs.size(); i++)
		if (i >= removeidx.size() || !removeidx[i]) { shrinkPolyIdx[i] = count++; newArcs.push_back(_arcs[i]); }
	_arcs = newArcs;
	// surface
	for (auto i = _surfs.begin(); i != _surfs.end(); i++)
		for (auto j = (*i).begin(); j != (*i).end(); j++) {
			if ((*j)._typ != 1) continue;
			(*j)._idx = shrinkPolyIdx[(*j)._idx];
		}
} //removeConincideArcs
void ShapeManager2DSurfs::removeConcincideNURBSCrrves() {
	if (_curves.size() == 0) return;
	double TOL = 1e-8;
	// Uniformlize
	NURBS::Curve c;
	_curves[0].get_curve(c);
	NS::Point3D pt; c.getCoord(-1, pt);
	double xmin = pt.x(), xmax = xmin;
	double ymin = pt.y(), ymax = ymin;
	for (size_t i = 0; i < _curves.size(); i++) {
		double _k[2] = { -1.0,1.0 };
		for (int k = 0; k<2; k++) {
			_curves[i].get_curve(c);
			c.getCoord(_k[k], pt);
			if (pt.x() > xmax) xmax = pt.x();
			if (pt.x() < xmin) xmin = pt.x();
			if (pt.y() > ymax) ymax = pt.y();
			if (pt.y() > ymin) ymin = pt.y();
		}
	}
	double scale = xmax - xmin > ymax - ymin ? (xmax - xmin) / 10.0 : (ymax - ymin) / 10.0;
	TOL *= scale;
	//
	std::vector<int> idx2group(_curves.size(), -1);
	std::vector<std::vector<int>> group;
	std::vector<std::vector<int>> dir;
	for (size_t i = 0; i < _curves.size() - 1; i++) { // O(n^2)
		for (size_t j = i + 1; j < _curves.size(); j++) {
			int conDir = _curves[i].isConincideWith(_curves[j], TOL);
			if (conDir ==0) continue;
			if (idx2group[i] == -1 && idx2group[j] == -1) {
				idx2group[i] = group.size();
				idx2group[j] = group.size();
				group.push_back(std::vector<int>(2));
				group.back().front() = i; group.back().back() = j;
				dir.push_back(std::vector<int>(2,1));
				if (conDir == -1) dir.back().back() = -1;
			}
			else if (idx2group[i] == -1) {
				group[idx2group[j]].push_back(i);
				idx2group[i] = idx2group[j];
				int idx = std::find(group[idx2group[j]].cbegin(), group[idx2group[j]].cend(), i) - group[idx2group[j]].cbegin();
				dir[idx2group[j]].push_back(conDir * dir[idx2group[j]][idx]);
			}
			else if (idx2group[j] == -1) {
				group[idx2group[i]].push_back(j);
				idx2group[j] = idx2group[i];
				int idx = std::find(group[idx2group[i]].cbegin(), group[idx2group[i]].cend(), j) - group[idx2group[i]].cbegin();
				dir[idx2group[i]].push_back(conDir * dir[idx2group[i]][idx]);
			}
		}
	}
	if (group.size() == 0) return;
	// update surf
	std::vector<int> splineModify(_curves.size());
	std::vector<int> splDir(_curves.size(), 0);
	std::vector<bool> modIdx(_curves.size(), false);
	for (size_t i = 0; i < group.size(); i++) 
		for (size_t j = 1; j < group[i].size(); j++) {
			modIdx[group[i][j]] = true;
			splineModify[group[i][j]] = group[i][0];
			splDir[group[i][j]] = dir[i][j];
		}
	for (auto i = _surfs.begin(); i != _surfs.end(); i++) {
		for (auto j = (*i).begin(); j != (*i).end(); j++) {
			if ((*j)._typ != 2) continue;
			if (modIdx[(*j)._idx]) {
				if (splDir[(*j)._idx] == -1)
					(*j)._dir = -(*j)._dir;
				(*j)._idx = splineModify[(*j)._idx];
			}
		}
	}
	std::vector<int> newCurvesIdx;
	std::vector<IGES_NURBSCurve> newCurves;
	int count = 0;
	for (size_t i = 0; i < _curves.size(); i++) {
		if (modIdx[i]) newCurvesIdx.push_back(-1);
		else {
			newCurvesIdx.push_back(count++);
			newCurves.push_back(_curves[i]);
		}
	}
	_curves = newCurves;
	for (auto i = _surfs.begin(); i != _surfs.end(); i++) 
		for (auto j = (*i).begin(); j != (*i).end(); j++) {
			if ((*j)._typ != 2) continue;
			(*j)._idx = newCurvesIdx[(*j)._idx];
		}
} // removeConcincideNURBSCrrves
NS_END