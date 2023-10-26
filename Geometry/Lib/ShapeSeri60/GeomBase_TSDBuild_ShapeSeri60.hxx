#pragma once
#include<TopoDS_Shape.hxx>
#include<TopoDS_Edge.hxx>
#include<vector>
#include<gp_Pnt.hxx>
#include<Geom_BSplineCurve.hxx>
#include<Geom_Curve.hxx>
#include<Geom_Surface.hxx>
#include <Geom_Geometry.hxx>
#include<TColgp_Array1OfPnt.hxx>


class Offset
{

private:

	int theDegree_;
	double theCB_;
	double theLCB_;

public:


	Offset(int deg, double cb, double lcb);


	void PlotCurve(std::vector<Handle(Geom_Curve)> curve, std::string filename, int n);



	std::vector<double> CalculateKnots
	(
		const std::vector<gp_Pnt>& ctrlPts,
		int degree,
		bool isPeriodic
	);

	std::vector<int> CalculateMultiplicities
	(
		const std::vector<gp_Pnt>& ctrlPts,
		int degree,
		bool isPeriodic
	);




	std::vector<double> CalcBasisFunctionVector
	(
		const std::vector<double>& knots,
		int i,
		int p,
		double u
	) const;

	double CalcBasisFunction
	(
		const std::vector<double>& knots,
		int i,
		int p,
		double u
	) const;

	int FindSpan
	(
		int n,
		int p,
		double u,
		const std::vector<double>& knots
	);

	std::vector<double> CalcKnotBarEquallySpaced
	(
		const std::vector<gp_Pnt>& Q
	);

	std::vector<double> CalcKnotBarChordLength
	(
		const std::vector<gp_Pnt>& Q
	);

	std::vector<double> CalcKnotBarCentripetal
	(
		const std::vector<gp_Pnt>& Q
	);

	std::vector<double> CalcKnotVector
	(
		const std::vector<gp_Pnt>& Q,
		const std::vector<double>& knotsBar,
		int degree
	);

	std::vector<gp_Pnt> Interpolate
	(
		const std::vector<gp_Pnt>& Q,
		int degree
	);

	Handle(Geom_BSplineCurve) CreateBSPLineCurve
	(
		const std::vector<gp_Pnt>& ctrlPts,
		const std::vector<double>& theKnots,
		const std::vector<int>& theMultiplicities,
		int degree,
		bool isPerodic
	);


	Handle(Geom_BSplineCurve) SetData
	(
		const std::vector<gp_Pnt>& ctrlPts,
		int degree
	);


	std::vector<std::vector<gp_Pnt>> ReadOffsetOne();




	bool areEqual(const std::vector<gp_Pnt>& v1, const std::vector<gp_Pnt>& v2);
	std::vector<std::vector<gp_Pnt>> ReadOffsetOneFindMidUp()const;
	std::vector<std::vector<gp_Pnt>> ReadOffsetOneFindMidDown()const;
	std::vector<double> FindeMidCb()const;


	void SetDegree(int deg) { theDegree_ = deg; }
	void SetCB(double cbval);
	void SetLCB(double lcb) { theLCB_ = lcb; }


	int GetDegree()const { return theDegree_; }
	double GetCB()const { return theCB_; }
	double GetLCB()const { return theLCB_; }

	std::string doubleToString(double value)const;

	std::vector<Handle(Geom_Curve)> CreatCurveOffsetByMySpline();
	std::vector<Handle(Geom_Curve)> CreatCurveOffsetOneByInterpolateNewMethod();
	std::vector<Handle(Geom_Curve)> CreatCurveOffsetOneByInterpolate();
	std::vector<Handle(Geom_Curve)> CreatCurvByOCC();

	TopoDS_Shape SurfaceOffsetOneByOccBsplineCurve();
	TopoDS_Shape shapebyinterpolate1();
	TopoDS_Shape ShapeOfGeomTrimCurve();

	bool HandeValueOfLCB();

	void creatIGSfile(TopoDS_Shape shape);

};




