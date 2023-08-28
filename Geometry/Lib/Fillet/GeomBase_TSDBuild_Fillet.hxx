#pragma once
#include<TopoDS_Shape.hxx>
#include<TopoDS_Edge.hxx>
#include<vector>
#include<GeomBaseGlobal_TSDBuild.hxx>



class TonbSoftDev_GeomTools_Fillet
{
private:

	TopoDS_Shape  theShape_;
	double  theRadius_;
	std::vector<int> theIndex_;

public:

	Geometry_EXPORT TonbSoftDev_GeomTools_Fillet(const TopoDS_Shape& shape, double radius, const std::vector<int>& in);
	Geometry_EXPORT TopoDS_Shape ApplyFillet();
	Geometry_EXPORT void PlotFillet(std::string filename);
	
	//TonbSoftDev_Fillet_EXPORT TopoDS_Shape ApplyFilletTwo();
	//TonbSoftDev_Fillet_EXPORT void PlotFilletTwo(std::string filename);

private:

	Geometry_EXPORT void SetShape( TopoDS_Shape sh);
	Geometry_EXPORT void SetRadius( double rad);
	Geometry_EXPORT void SetVectorOfEdge(std::vector<int> vector);

	Geometry_EXPORT TopoDS_Shape GetShape()const;
	Geometry_EXPORT double GetRadius()const;
	Geometry_EXPORT std::vector<int> GetVectorOfEdge()const;
	
	Geometry_EXPORT std::vector<TopoDS_Edge> FindEdgesFromNumber(std::vector<int> index);
	Geometry_EXPORT int NumberOfEdge();
	Geometry_EXPORT std::vector<TopoDS_Edge> FindEdgesFromShape();

};