#include<iostream>
#include<GeomBase_TSDBuild_Fillet.hxx>
#include<GeomBase_TSDBuild_Chamfer.hxx>
#include<BRepPrimAPI_MakeBox.hxx>
#include<TopTools_IndexedMapOfShape.hxx>
#include<Bnd_Box.hxx>


using namespace std;

int main()
{
	try
	{


		{

			gp_Pnt corner1(0.0, 0.0, 0.0);
			gp_Pnt corner2(5.0, 5.0, 5.0);
			BRepPrimAPI_MakeBox make(corner1, corner2);
			TopoDS_Shape shape0 = make.Shape();


			TonbSoftDev_GeomTools_Fillet fil(shape0,0.5 , {1 , 9, 5} );
			TonbSoftDev_GeomTools_Chamfer cham(shape0,0.5 , {1 , 9, 3} );
			

			fil.PlotFillet("Fillet Edge");
			cham.PlotChamfer("Chamfer Edge");

		}

	}// <-- this bracket for try



	catch (const Standard_Failure& ex)
	{

		std::cout << ex.GetMessageString() << std::endl;

	}



	return 0;

}
