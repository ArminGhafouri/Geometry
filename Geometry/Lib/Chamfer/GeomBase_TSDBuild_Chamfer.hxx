#pragma once
#include<TopoDS_Shape.hxx>
#include<TopoDS_Edge.hxx>
#include<vector>
#include<GeomBaseGlobal_TSDBuild.hxx>



class TonbSoftDev_GeomTools_Chamfer
{
private:

	TopoDS_Shape  theShape_;
	double  theRadius_;
	std::vector<int> theIndex_;

public:

	Geometry_EXPORT TonbSoftDev_GeomTools_Chamfer(const TopoDS_Shape& shape, double radius, const std::vector<int>& in);
	Geometry_EXPORT TopoDS_Shape ApplyChamfer();
	Geometry_EXPORT void PlotChamfer(std::string filename);
	
	

private:

	Geometry_EXPORT void SetShape(TopoDS_Shape sh);
	Geometry_EXPORT void SetRadius(double rad);
	Geometry_EXPORT void SetVectorOfEdge(std::vector<int> vector);

	Geometry_EXPORT TopoDS_Shape GetShape()const;
	Geometry_EXPORT double GetRadius()const;
	Geometry_EXPORT std::vector<int> GetVectorOfEdge()const;

	Geometry_EXPORT std::vector<TopoDS_Edge> FindEdgesFromNumber(std::vector<int> index);
	Geometry_EXPORT int NumberOfEdge();
	Geometry_EXPORT std::vector<TopoDS_Edge> FindEdgesFromShape();

};