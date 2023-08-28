#include"GeomBase_TSDBuild_Fillet.hxx"
#include<GeomBase_TSDBuild_Tools.hxx>

#include<BRepFilletAPI_MakeFillet.hxx>
#include<TopoDS.hxx>
#include<TopoDS_Face.hxx>
#include<TopExp_Explorer.hxx>
#include<TopExp.hxx>
#include<BRepBuilderAPI_MakeSolid.hxx>
#include<ShapeFix_Shape.hxx>
#include<BRepBuilderAPI_MakeWire.hxx>
#include<BRepBuilderAPI_MakeEdge.hxx>
#include<BRepBuilderAPI_MakeFace.hxx>




TonbSoftDev_GeomTools_Fillet::TonbSoftDev_GeomTools_Fillet(const TopoDS_Shape& shape, double radius, 
	const std::vector<int>& in)
{
	SetShape(shape);
	SetRadius(radius);
	SetVectorOfEdge(in);
	
}


void TonbSoftDev_GeomTools_Fillet::SetShape(const TopoDS_Shape shape)
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


void TonbSoftDev_GeomTools_Fillet::SetRadius(const double rad)
{
	
	if (rad <= 0.0)
	{
		throw std::exception("Radius is not valid!");

	}

      theRadius_ = rad;
	
}

 void TonbSoftDev_GeomTools_Fillet::SetVectorOfEdge(std::vector<int> vector)
{
	if (vector.size() == 0.0 )
	{
		throw std::exception("number of Edge is null!");

	}
	
	for (int i = 0 ; i < vector.size() ; i++)
	{
		if (vector[i] > NumberOfEdge())
		{
			throw std::exception("number of Edge is not valid!");
		}

	}

	theIndex_ = vector;
}



TopoDS_Shape TonbSoftDev_GeomTools_Fillet::GetShape()const
{
	//if (theShape_.IsNull())
	//{ 
	//	throw std::exception("The Shape is null!");
	//}	//throw;
		return theShape_;
}

double TonbSoftDev_GeomTools_Fillet::GetRadius()const
{
	//if (theRadius_ < 0.0)
	//{
	//	throw std::exception("Radis is not valid!");
	//	/*throw;*/
	//}
	//
		return theRadius_;
}

std::vector<int> TonbSoftDev_GeomTools_Fillet::GetVectorOfEdge() const
{
	//if (theIndex_.size() == 0.0)
	//{
	//	throw std::exception("number of Edge null!");
	//	
	//}
	return theIndex_;
}





int TonbSoftDev_GeomTools_Fillet::NumberOfEdge(/*const TopoDS_Shape shape*/)
{

	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
	int num = mymap.Size();

	return num;
}


TopoDS_Shape TonbSoftDev_GeomTools_Fillet::ApplyFillet()
{

		TopTools_IndexedMapOfShape mymap;
		TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
		BRepFilletAPI_MakeFillet mkFillet(theShape_);

		for (int i = 0; i < theIndex_.size(); i++)
		{

			mkFillet.Add(theRadius_, TopoDS::Edge(mymap(theIndex_[i])));
			//mkFillet.Build();

		}
		
		TopoDS_Shape shape;
			try
			{
				 shape = mkFillet.Shape();
			}
	
			catch (const Standard_Failure& ex)
			{
				throw std::exception(ex.GetMessageString());
			}

	return shape;

}


void TonbSoftDev_GeomTools_Fillet::PlotFillet(std::string filename)
{

	Tools::PlotShapeTwo(ApplyFillet(), filename);

}


std::vector<TopoDS_Edge> TonbSoftDev_GeomTools_Fillet::FindEdgesFromNumber(std::vector<int> index)
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


//TopoDS_Shape TonbSoftDev_GeomTools_Fillet::ApplyFilletTwo()
//{
//	/*****/
//	std::vector<TopoDS_Edge> AllEdgs = FindEdgesFromShape();
//	BRepFilletAPI_MakeFillet mkFillet(GetShape());
//	
//
//		//BRepFilletAPI_MakeFillet2d mkFillet();
//	for (int i = 0; i < theIndex_.size(); i++)
//	{
//		//mkFillet().ModifyFillet(AllEdgs[theIndex_[i]], theRadius_);
//	    mkFillet.Add(theRadius_, AllEdgs[theIndex_[i]]);
//	    //mkFillet.Build();
//	}
//
//	TopoDS_Shape shape = mkFillet.Shape();
//	return shape;
//
//}




//void TonbSoftDev_GeomTools_Fillet::PlotFilletTwo(std::string filename)
//{
//	Tools::PlotShapeTwo(ApplyFilletTwo(), filename );
//
//}





std::vector<TopoDS_Edge> TonbSoftDev_GeomTools_Fillet::FindEdgesFromShape()
{
	TopTools_IndexedMapOfShape mymap;
	TopExp::MapShapes(GetShape(), TopAbs_EDGE, mymap);
	std::vector<TopoDS_Edge> edgs;
	

	for (int i = 1 ; i < mymap.Size()+1 ; i++)
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

	std::cout << "num  of edge " <<  edgs.size();
	return  edgs;
}