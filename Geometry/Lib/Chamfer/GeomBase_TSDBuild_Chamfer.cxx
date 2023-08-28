#include"GeomBase_TSDBuild_Chamfer.hxx"
#include<GeomBase_TSDBuild_Tools.hxx>

#include<BRepFilletAPI_MakeChamfer.hxx>
#include<TopoDS.hxx>
#include<TopoDS_Face.hxx>
#include<TopExp_Explorer.hxx>
#include<TopExp.hxx>
#include<BRepBuilderAPI_MakeSolid.hxx>
#include<ShapeFix_Shape.hxx>
#include<BRepBuilderAPI_MakeWire.hxx>
#include<BRepBuilderAPI_MakeEdge.hxx>
#include<BRepBuilderAPI_MakeFace.hxx>




TonbSoftDev_GeomTools_Chamfer::TonbSoftDev_GeomTools_Chamfer(const TopoDS_Shape& shape, double radius,
	const std::vector<int>& in)
{
	SetShape(shape);
	SetRadius(radius);
	SetVectorOfEdge(in);

}


void TonbSoftDev_GeomTools_Chamfer::SetShape(const TopoDS_Shape shape)
{

	if (shape.IsNull())
	{
		throw std::exception("The Shape is null!");

	}
	else
	{
		theShape_ = shape;
	}

}


void TonbSoftDev_GeomTools_Chamfer::SetRadius(const double rad)
{

	if (rad <= 0.0)
	{
		throw std::exception("Radius is not valid!");

	}

	theRadius_ = rad;

}

void TonbSoftDev_GeomTools_Chamfer::SetVectorOfEdge(std::vector<int> vector)
{
	if (vector.size() == 0.0)
	{
		throw std::exception("number of Edge is null!");

	}

	for (int i = 0; i < vector.size(); i++)
	{
		if (vector[i] > NumberOfEdge())
		{
			throw std::exception("number of Edge is not valid!");
		}

	}

	theIndex_ = vector;
}



TopoDS_Shape TonbSoftDev_GeomTools_Chamfer::GetShape()const
{
	//if (theShape_.IsNull())
	//{ 
	//	throw std::exception("The Shape is null!");
	//}	//throw;
	return theShape_;
}

double TonbSoftDev_GeomTools_Chamfer::GetRadius()const
{
	//if (theRadius_ < 0.0)
	//{
	//	throw std::exception("Radis is not valid!");
	//	/*throw;*/
	//}
	//
	return theRadius_;
}

std::vector<int> TonbSoftDev_GeomTools_Chamfer::GetVectorOfEdge() const
{
	//if (theIndex_.size() == 0.0)
	//{
	//	throw std::exception("number of Edge null!");
	//	
	//}
	return theIndex_;
}





int TonbSoftDev_GeomTools_Chamfer::NumberOfEdge(/*const TopoDS_Shape shape*/)
{

	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
	int num = mymap.Size();

	return num;
}


TopoDS_Shape TonbSoftDev_GeomTools_Chamfer::ApplyChamfer()
{

	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
	BRepFilletAPI_MakeChamfer mkChamfer(theShape_);

	for (int i = 0; i < theIndex_.size(); i++)
	{

		mkChamfer.Add(theRadius_, TopoDS::Edge(mymap(theIndex_[i])));
		//mkFillet.Build();

	}

	TopoDS_Shape shape;
	try
	{
		shape = mkChamfer.Shape();
	}

	catch (const Standard_Failure& ex)
	{
		throw std::exception(ex.GetMessageString());
	}

	return shape;

}


void TonbSoftDev_GeomTools_Chamfer::PlotChamfer(std::string filename)
{

	Tools::PlotShapeTwo(ApplyChamfer(), filename);

}


std::vector<TopoDS_Edge> TonbSoftDev_GeomTools_Chamfer::FindEdgesFromNumber(std::vector<int> index)
{

	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
	std::vector<TopoDS_Edge> edgs;


	for (int i = 0; i < index.size(); i++)
	{
		TopoDS_Edge edge = TopoDS::Edge(mymap(index[i]));
		edgs.push_back(edge);
	}


	return edgs;
}


std::vector<TopoDS_Edge> TonbSoftDev_GeomTools_Chamfer::FindEdgesFromShape()
{
	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
	std::vector<TopoDS_Edge> edgs;


	for (int i = 1; i < mymap.Size() + 1; i++)
	{
		TopoDS_Edge edge = TopoDS::Edge(mymap(i));
		edgs.push_back(edge);
	}

	//TopExp_Explorer Ex;
	//std::vector<TopoDS_Edge> edgs;
	//for (Ex.Init(GetShape(), TopAbs_EDGE); Ex.More(); Ex.Next())
	//{
	//	const TopoDS_Edge& edge = TopoDS::Edge(Ex.Current());
	//	edgs.push_back(edge);
	//}

	std::cout << "num  of edge " << edgs.size();
	return  edgs;
}