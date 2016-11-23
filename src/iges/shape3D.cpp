#include "shape3D.h"
#include <typeinfo>

NSI_BEG


/**************************************** Shape3D **************************************************/
void Shape3D::matrixTran( NS::Point3D &pt) const {
	arma::vec coor(3, 1);
	coor[0] = pt.x();
	coor[1] = pt.y();
	coor[2] = pt.z();
	coor = _R*coor + _T;
	pt.set(coor[0], coor[1], coor[2]);
}
Shape3D::Shape3D() {
	_R = arma::mat(3, 3, arma::fill::eye);
	_T = arma::vec(3, 1, arma::fill::zeros);
}
void Shape3D::setMatTran(const arma::mat &R, const arma::vec T) {
	return;
}
/**************************************** NurbsSurf 128 ************************************************/
void Surf128::setSurf(const NURBS::ParentSurface &surf) { 
	_surf = surf;
	_subClass = 1; 
}
void Curve126::setCurve(const NURBS::ParentCurve &cc) {
	_curve = cc; 
	_subClass = 2;
}
/**************************************** Shell 514 **************************************************/
void Shell514::addFace(const Face510 &face, const int dir) { 
	_faces.push_back(face);
	_direction.push_back(dir);
}
/**************************************** Loop 508 **************************************************/
void Loop508::addEdge(Edge504 const &edge, int const dir) {
	_edge.push_back(edge);
	_direction.push_back(dir);
}
/**************************************** Group **************************************************/
template class ShapesGroup<NS::Point3D, std::pair<int, int>>;
template class ShapesGroup<std::shared_ptr<Surface3D>, int>;
template class ShapesGroup<std::shared_ptr<Curve3D>, int>;
template class ShapesGroup<Solid186, int>;
template<class COMP, class DIR>
int ShapesGroup<COMP, DIR>::count() const {
	return _comp.size() - 1; 
}

template<class COMP, class DIR>
int ShapesGroup<COMP, DIR>::findDirIdx(const DIR & dir) {
	auto it = std::find(_directoryIdx.begin(), _directoryIdx.end(), dir);
	if (it != _directoryIdx.end())
		return it - _directoryIdx.begin();
	else
		return -1;
}

/**************************************** Surfaces Group **************************************************/
void SurfacesGroup::addComp(const std::shared_ptr<Surface3D> &comp, const int &idx) {
	switch (comp->getSurfaceTyp()) {
	case 1: {
		auto temp = std::static_pointer_cast<Surf128>(comp);
		NURBS::ParentSurface ss;
		temp->getSurf(ss);
		if (ss.is_plane)






			//something
			auto s = 1;
		else
			_comp.push_back(comp);
		break;
	}
	default:
		assert(0);
		break;
	}	
	_directoryIdx.push_back(idx);
}
NS_END