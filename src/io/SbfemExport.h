#include <vector>
#include "../geo/point3d.h"
#include "../iges/shape.h"
#include "../iges/IGESReader.h"
#include "../nurbs/curve.h"
#include <assert.h> 

namespace IGES {
	class SBFEMExport{
	private:
		std::vector<NS::Point3D> _keyPoints;
		std::vector<std::vector<int>> _polyLines;
		struct Member {
		public:
			int _idx;
			bool _reversed;
			Member(const int idx, const bool rev) { _idx = idx; _reversed = rev; }
		};
		std::vector<std::vector<Member>> _surfaces;
		std::vector<bool> _isHardpt;
		std::vector<int> _level;
		std::vector<std::string> _layerName;
		std::vector<std::vector<double>> _color;
		std::vector<int> _matColor;
		double calculateAngle3pt(const NS::Point3D &pt0, const NS::Point3D &pt1, const NS::Point3D &pt2);

		void parseSurf2D(ShapeManager2DSurfs &s2d, const double TOL);					// how many points? how to locate?
		void findHardPt(const double angle);

		//void parseShapeGrouptoPolylines(const ShapeGroup &gs, const double TOL);
		int buildPolyLines(const std::vector<NS::Point3D> &member_pts, int pt_idx);
	public:
		void getLevel(std::vector<int> &l)const { l = _level; }
		void getLayerName(std::vector<std::string> &l) const { l = _layerName; }
		void getMatColor(std::vector<int> &c) const { c = _matColor; }
		void setMatColor(const std::vector<int> &c) { _matColor = c; }
		void getKPT(std::vector<NS::Point3D> &keyPoints) const { keyPoints = _keyPoints; }
		void getPolyLines(std::vector<std::vector<int>> &poly) const { poly = _polyLines; }
		void getSurfaces(std::vector<std::vector<int>> &surf) const;
		void getReversed(std::vector<std::vector<bool>> &rev) const;
		void getIsHardPt(std::vector<bool> &isHardPt) const { isHardPt = _isHardpt; }

		SBFEMExport(ShapeManager2DSurfs &in,const double TOL, const double angle);
		void export2vtk(const std::string &pathOut, const std::string &name) const;
		void export2MFile2D(const std::string &pathOut, const std::string &name) const;
		void filterByLayer(const std::vector<std::string> &layer, const char mode);
		void setLayerColor(const std::string &layer, const int color);
		void setLayerColor(const std::vector<int> &level, const std::vector<int> &color);
	};
}