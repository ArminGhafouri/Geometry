#include"GeomBase_TSDBuild_ShapeSeri60.hxx"
#include<BRep_Tool.hxx>
#include<GeomBase_TSDBuild_Tools.hxx>
#include<GeomBase_TSDBuild_Bspline.hxx>



#include <algorithm>
#include <vector>
#include<iostream>
#include <fstream>  
#include <string> 

#include <GeomFill_Coons.hxx> 
#include <GeomFill_Curved.hxx > 
#include <GeomFill_Stretch.hxx> 
#include <GeomConvert_CompBezierSurfacesToBSplineSurface.hxx> 
#include<Geom_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include<Geom2dAPI_InterCurveCurve.hxx>
#include <GeomAPI.hxx>
#include <gp_Pln.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomConvert.hxx>
#include <Geom_Plane.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <gp_Pnt2d.hxx>
#include <BRepBuilderAPI_MakeShape.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <NCollection_Array1.hxx>
#include <GeomFill_Sweep.hxx>
#include <GeomFill_SweepSectionGenerator.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <GeomFill_AppSweep.hxx>



#include<BrepOffsetAPI_MakeOffset.hxx>
#include <STEPControl_Writer.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_BezierCurve.hxx>
#include <Geom2d_Geometry.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepAlgoAPI_Fuse.hxx>


#include <GeomFill_BezierCurves.hxx>
#include <GeomFill_BSplineCurves.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepFill_CurveConstraint.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shape.hxx>
#include <Geom_BoundedSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Circle.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomFill_Pipe.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom2dAPI_PointsToBSpline.hxx>
#include <TColGeom_Array1OfCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <Geom_Line.hxx>
#include <gp_Lin.hxx>
#include <GeomFill.hxx>
#include <TopoDS.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <STEPControl_Writer.hxx>
#include <Interface_Static.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <BRepOffsetAPI_Sewing.hxx>
#include <BRepPrimAPI_MakeSweep.hxx>
#include <GeomFill_SectionGenerator.hxx>
#include <GeomFill_AppSurf.hxx>
#include <GeomFill_Line.hxx>


#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <STEPControl_Controller.hxx>
#include <STEPControl_Writer.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <IGESControl_Controller.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <Interface_Static.hxx>
#include <TopoDS.hxx>

#include <Geom_Curve.hxx>
#include <Geom_BezierCurve.hxx>
#include <GeomAPI.hxx>
#include <Geom_Curve.hxx>
#include <GeomFill_Pipe.hxx>
#include <Geom_Surface.hxx>
#include <Geom_BezierCurve.hxx>
#include <GeomAPI.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomFill_Sweep.hxx>
#include <TopoDS_Shell.hxx>

#include <TopoDS_Shape.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <Poly_Triangulation.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <Poly_Array1OfTriangle.hxx>
#include <Poly_CoherentTriangulation.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <GeomFill_NSections.hxx>

#include <FrgBase_Vcr.hxx>
#include <FrgBase_FullMtx.hxx>
#include <FrgBase_AbsMtx.hxx>

#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <TopoDS_Compound.hxx>



Offset::Offset(int deg, double cb, double lcb)
{
	theDegree_ = deg;        /*    degree for creat curve     */
	SetCB(cb);               /* 0.6 <= cb <= 0.8 */
	theLCB_ = lcb;

}


void Offset::SetCB(double cbval)
{
	if (cbval > 0.800 || cbval < 0.600)
	{

		throw std::exception("CB not Valid ");
	}

	theCB_ = cbval;
}



std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetOneByInterpolate()
{

	std::vector<Handle(Geom_Curve)> curvetotal0;

	Handle(TColgp_HArray1OfPnt) Pts0 =
		new TColgp_HArray1OfPnt;

	std::vector <std::vector < gp_Pnt >> poinOff = ReadOffsetOne();

	for (int i = 0; i < poinOff.size(); i++)
	{

		Pts0->Resize(1, poinOff[i].size(), false);


		for (int j = 0; j < poinOff[i].size(); j++)
		{

			/*double x = poinOff[i][j].X();
			double y = poinOff[i][j].Y();
			double z = poinOff[i][j].Z();
			gp_Pnt pointfinall(x, y, z);*/

			Pts0->SetValue(j + 1, poinOff[i][j].XYZ());


		}


		GeomAPI_Interpolate inter0(Pts0, false, 1e-2);
		inter0.Perform();
		Handle(Geom_Curve) cur0 = inter0.Curve();

		curvetotal0.push_back(cur0);

	}

	if (curvetotal0.size() < 2)
	{
		throw std::exception(" Vector of Geom_Curve Size is Small");
	}

	return curvetotal0;

}



TopoDS_Shape Offset::shapebyinterpolate1()
{

	NCollection_List < Handle(Geom_Curve)> col;

	for (int i = 0; i < CreatCurveOffsetByMySpline().size(); i++)
	{

		col.Append(CreatCurveOffsetByMySpline()[i]);
	}

	std::cout << "Ncollection_list size: " << col.Size();

	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PConfusion());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());


	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt2 = 0;
	Standard_Real aTol3d = 1e-4, aTol2d = Precision::Parametric(aTol3d);


	GeomFill_AppSurf anAlgo2(1, 2, aTol3d, aTol2d, aNbIt2, false);


	anAlgo2.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;


	aRes = new Geom_BSplineSurface(anAlgo2.SurfPoles(), anAlgo2.SurfWeights(),
		anAlgo2.SurfUKnots(), anAlgo2.SurfVKnots(), anAlgo2.SurfUMults(), anAlgo2.SurfVMults(),
		anAlgo2.UDegree(), anAlgo2.VDegree());



	BRepBuilderAPI_MakeFace face(aRes, 1, 10, 1, 10, 10);
	TopoDS_Shape shape = face.Shape();

	if(shape.IsNull())
	{
		throw std::exception("The Shape is Null");
	}
	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	//	gp_Pnt po0(ReadOffsetOne()[0][0].X(), ReadOffsetOne()[0][0].Y(), ReadOffsetOne()[0][0].Z());
	gp_Pnt po0(ReadOffsetOne()[0][0].XYZ());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, false);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	return ship;
}





std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetByMySpline()
{

	std::vector <std::vector <gp_Pnt>> poinOff = ReadOffsetOne();

	std::vector<Handle(Geom_Curve)> curvetotal;

	for (int i = 0; i < poinOff.size(); i++)
	{
		std::vector<gp_Pnt> Poin;

		for (int j = 0; j < poinOff[i].size(); j++)
		{

			gp_Pnt pnt = poinOff[i][j];
			Poin.push_back(pnt);

		}

		GeomBase_TSDBuild_Bspline bs(theDegree_, Poin);
		Handle(Geom_Curve) cur = bs.CreateOCCurve();

		curvetotal.push_back(cur);

	}

	if (curvetotal.size() < 2)
	{
		throw std::exception(" Vector of Geom_Curve Size is Small");
	}

	return  curvetotal;

}








std::vector<std::vector<gp_Pnt>> Offset::ReadOffsetOne()
{

	double length, beam, draft/*, valuecb, valuecm, valuelcb*/;

	std::cout << "Enter Length Of Ship: " << std::endl;
	std::cin >> length;
	std::cout << "Enter Beam Of Ship: " << std::endl;
	std::cin >> beam;
	std::cout << "Enter ShipDepth : " << std::endl;
	std::cin >> draft;


	std::vector<double> cbvalue{ 0.600,0.610,0.620,0.630,0.640,0.650,0.650,0.660,0.670,0.670
	,0.680,0.690,0.700,0.710,0.720,0.730,0.740,0.750,0.760,0.770,0.780,0.790,0.800 };



	for (int ij = 0; ij < cbvalue.size(); ij++)
	{

		while (theCB_ == cbvalue[ij])
		{

			if (!HandeValueOfLCB())
			{

				throw std::exception("LCB not valid!");

				std::cout << "   LCB NOT VALID   " << std::endl;

			}


			std::string str1 = doubleToString(theCB_);
			//std::string str2 = doubleToString(theCM_);
			std::string str3 = doubleToString(theLCB_);

			std::string space = " ";
			std::string theFileName_ = "cb=" + str1 + space + "lcb=" + str3 + ".txt";
			std::cout << "File_Name:" << theFileName_ << std::endl;

			std::fstream fin(theFileName_, std::ios::in);

			if (!fin.is_open())
			{
				std::exception e(" file not open (LCB  amount is not available) ");
				throw e;
			}

			std::string a;
			int numberofWaterline, numberofeStation;
			double  CB, CM, LCB;

			fin >> numberofWaterline;
			getline(fin, a);
			fin >> numberofeStation;

			/*std::cout << "numberofWaterline = "<< numberofWaterline << "     :    " <<"numberofeStation = "
			<< numberofeStation << std::endl;*/

			for (int i = 0; i < 6; i++)
			{

				getline(fin, a);

			}

			fin >> CB;
			getline(fin, a);

			fin >> CM;
			getline(fin, a);

			fin >> LCB;
			getline(fin, a);

			std::cout << "CB=" << CB << " , " << "CM=" << CM << " , " << "LCB=" << LCB << std::endl;

			std::vector<long double> Z_Core;
			std::vector<long double>X_Core;
			std::vector< std::vector<long float>>Y_Coretotal;

			for (int i = 0; i < numberofWaterline; i++)
			{
				double z;
				fin >> z;
				double zPrim = z * ((draft / 1.5) / 4);

				Z_Core.push_back(zPrim);
				//std::cout << Z_Core[i] << std::endl;

			}
			// std::cout << Z_Core[15] << std::endl;



			for (int i = 0; i < 2; i++)
			{

				getline(fin, a);

			}
			//std::cout << a << std::endl;



			for (int i = 0; i < numberofeStation; i++)
			{

				double x;
				fin >> x;
				double xPrim = x * (length / 20);
				X_Core.push_back(xPrim);

				//std::cout << x << std::endl;

				std::vector< double>Y_Core;

				for (int j = 0; j < numberofWaterline; j++)
				{
					double y;
					fin >> y;
					double yPrim = y * beam;

					Y_Core.push_back(yPrim);

				}

				Y_Coretotal.push_back(Y_Core);

				//std::cout << Y_Coretotal[i][1] << std::endl;
				getline(fin, a);

			}



			std::vector< std::vector<gp_Pnt> > pointss;


			for (int i = 0; i < numberofeStation; i++)
			{


				std::vector<gp_Pnt> points;
				double xdummy = X_Core[i];


				for (int j = 0; j < numberofWaterline; j++)
				{

					double Ysec = Y_Coretotal[i][j];
					double zdummy = Z_Core[j];

					/*if (j == 0)
					{
						Ysec = 0.0;
					}*/

					gp_Pnt point1(xdummy, Ysec, zdummy);

					points.push_back(point1);


				}


				pointss.push_back(points);

			}


			std::vector< std::vector<gp_Pnt> > pointss2;

			for (int i = 0; i < pointss.size(); i++)
			{


				std::vector<gp_Pnt> points0;
				double xdummy = X_Core[i];


				for (int j = 0; j < pointss[i].size() /*-1*/; j++)
				{

					double y1 = Y_Coretotal[i][j];


					if (std::abs(y1) != 0.0 && j != 0)
					{


						double Ysec000 = Y_Coretotal[i][j - 1];
						double zdummy = Z_Core[j - 1];

						gp_Pnt point00(xdummy, Ysec000, zdummy);
						points0.push_back(point00);

					}


				}

				int lastpoint = pointss[i].size() - 1;
				double Ysec00 = Y_Coretotal[i][lastpoint];
				double zdummy = Z_Core[lastpoint];
				gp_Pnt point00(xdummy, Ysec00, zdummy);
				points0.push_back(point00);
				pointss2.push_back(points0);

			}


			std::vector< std::vector<gp_Pnt> > pointssFinall;


			for (int i = 0; i < pointss2.size(); i++)
			{


				std::vector<gp_Pnt> points;
				long double xdummy = X_Core[i];


				for (int j = 0; j < pointss2[i].size(); j++)
				{

					long float Ysec0 = pointss2[i][j].Y();
					long double zdummy0 = pointss2[i][j].Z();

					if (j == 0)
					{
						Ysec0 = 0.0;
					}

					gp_Pnt point1(xdummy, Ysec0, zdummy0);

					points.push_back(point1);


				}


				pointssFinall.push_back(points);

			}


			if (pointssFinall.size() == 0)
			{
				throw std::exception("The size of vector Of Points is Zero ");
			}
			return pointssFinall;

		}


	}


	if (theCB_ != 0.600 || 0.610 || 0.620 || 0.630 || 0.640
		|| 0.650 || 0.650 || 0.660 || 0.670 || 0.670
		|| 0.680 || 0.690 || 0.700 || 0.710 || 0.720
		|| 0.730 || 0.740 || 0.750 || 0.760 || 0.770
		|| 0.780 || 0.790 || 0.800)
	{
		std::vector<std::vector<gp_Pnt>> upvec = ReadOffsetOneFindMidUp();
		std::vector<std::vector<gp_Pnt>> downvec = ReadOffsetOneFindMidDown();


		std::vector<std::vector<gp_Pnt>> pointssfin;

		double cbup = FindeMidCb()[1];
		double cbdown = FindeMidCb()[0];

		for (int i = 0; i < upvec.size(); i++)
		{
			std::vector<gp_Pnt> pointssedit;
			for (int j = 0; j < upvec[i].size(); j++)
			{

				double x2 = upvec[i][j].X();
				double x1 = downvec[i][j].X();
				double y2 = upvec[i][j].Y();
				double y1 = downvec[i][j].Y();
				double z2 = upvec[i][j].Z();
				double z1 = downvec[i][j].Z();

				double xprim = x1 + (((x2 - x1) / (cbup - cbdown)) * (cbup - cbdown));
				double yprim = y1 + (((y2 - y1) / (cbup - cbdown)) * (cbup - cbdown));
				double zprim = z1 + (((z2 - z1) / (cbup - cbdown)) * (cbup - cbdown));

				double scalex = xprim * (length / 20);
				double scaley = yprim * beam;
				double scalez = zprim * ((draft / 1.5) / 4);

				gp_Pnt edit(scalex, scaley, scalez);

				pointssedit.push_back(edit);

			}
			pointssfin.push_back(pointssedit);

		}
		if (pointssfin.size() == 0)
		{
			throw std::exception("The size of vector Of Points is Zero ");
		}
		return pointssfin;
	}
	//}
}

//for (int i = 0; i < pointss2.size(); i++)
//{
//	std::cout << "*************** NEW SEC *************" << std::endl;
//	std::cout << "*********************" << std::endl;

//	for (int j = 0; j < pointssFinall[i].size(); j++)
//	{
//		//std::cout << " size two (" << i << ") = " << pointss2[i].size() << std::endl;
//		std::cout << " Section [" << i << "] = " << pointssFinall[i][j].X()
//			<< "  " << pointssFinall[i][j].Y() << "  " << pointssFinall[i][j].Z() << std::endl;
//	}
//	

//}









Handle(Geom_BSplineCurve) Offset::CreateBSPLineCurve(const std::vector<gp_Pnt>& ctrlPts,
	const std::vector<double>& theKnots,
	const std::vector<int>& theMultiplicities,
	int degree,
	bool isPerodic)
{
	/*if (theCurve_)
	{
		theCurve_->Delete();
		theCurve_ = nullptr;
	}*/

	TColStd_Array1OfReal knots;
	TColStd_Array1OfInteger mult;

	knots.Resize(1, theKnots.size(), false);
	mult.Resize(1, theMultiplicities.size(), false);

	for (int i = 0; i < theKnots.size(); i++)
	{
		knots.SetValue(i + 1, theKnots[i]);
		mult.SetValue(i + 1, theMultiplicities[i]);
	}

	/*if constexpr (Dim == 2)
	{
		TColgp_Array1OfPnt2d Pts(1, ctrlPts.size());
		for (int i = 0; i < ctrlPts.size(); i++)
			Pts.SetValue(i + 1, gp_Pnt2d(ctrlPts[i].X(), ctrlPts[i].Y()));
		theCurve_ = opencascade::handle<Geom2d_BSplineCurve>(new Geom2d_BSplineCurve(Pts, knots, mult, degree, isPerodic));
	}*/
	/*else if constexpr (Dim == 3)
	{*/
	TColgp_Array1OfPnt Pts(1, ctrlPts.size());
	for (int i = 0; i < ctrlPts.size(); i++)
		Pts.SetValue(i + 1, gp_Pnt(ctrlPts[i].X(), ctrlPts[i].Y(), ctrlPts[i].Z()));

	Handle(Geom_BSplineCurve)theCurve_  /*opencascade::handle<Geom_BSplineCurve>*/(new Geom_BSplineCurve(Pts, knots, mult, degree, isPerodic));
	/*}*/

		/*Handle(Geom_Curve) curve = theCurve_;*/
	return theCurve_;

}

Handle(Geom_BSplineCurve) Offset::SetData(const std::vector<gp_Pnt>& ctrlPts, int degree)
{

	//theDegree_ = degree;

	/*if (ctrlPts.empty())
		return ;*/

		//if (ctrlPts.size() == 1 || ctrlPts.size() == 2)
		//{
		//	FrgVisual_PolylineActor<Dim>::SetData(ctrlPts);
		//	return;
		//}

	int newDegree = degree;
	if (degree >= ctrlPts.size())
	{
		newDegree = ctrlPts.size() - 1;
	}

	bool isPeriodic = false;
	if (ctrlPts[0].IsEqual(ctrlPts[ctrlPts.size() - 1], 1e-4))
	{
		isPeriodic = true;
	}

	Handle(Geom_BSplineCurve) curve = CreateBSPLineCurve(ctrlPts, CalculateKnots(ctrlPts, newDegree, isPeriodic), CalculateMultiplicities(ctrlPts, newDegree, isPeriodic), newDegree, isPeriodic);

	/*UpdateActor();*/
	return curve;
}




void Offset::PlotCurve(std::vector<Handle(Geom_Curve)> curve, std::string filename, int n)
{

	std::fstream My_File(filename + ".plt", std::ios::out);
	if (!My_File.is_open())
	{
		std::cout << "file not open";
		return;
	}


	My_File << "VARIABLES = X Y Z" << std::endl;
	// My_File << "ZONE T = Curve" << std::endl;

	auto curve0 = curve;


	for (int j = 0; j < curve0.size(); j++)
	{

		My_File << "ZONE T = Curve" << std::endl;
		double du = (curve0[j]->LastParameter() - curve0[j]->FirstParameter()) / n;
		/* std::cout << curve0[j]->FirstParameter() << std::endl;
		 std::cout << curve0[j]->LastParameter();*/

		for (int i = 0; i <= n; i++)
		{
			double u = (i * du) + curve0[j]->FirstParameter();
			gp_Pnt p = curve0[j]->Value(u);


			My_File << p.X() << "  " << p.Y() << "  " << p.Z() << std::endl;

		}


	}

	My_File.close();


}



TopoDS_Shape Offset::SurfaceOffsetOneByOccBsplineCurve()
{
	std::vector<TColgp_Array1OfPnt> points;
	std::vector<Handle(Geom_BSplineCurve)> BScurve;
	std::vector<std::vector<gp_Pnt>> all = ReadOffsetOne();

	for (int i = 0; i < all.size(); i++)
	{

		TColgp_Array1OfPnt po /*= TColgp_Array1OfPnt(1 , vectorOfvector()[n].size() )*/;
		po.Resize(1, all[i].size(), false);



		for (int j = 0; j < all[i].size(); j++)
		{
			/*std::cout << "siZe vecofvec[I] =  " << vectorOfvector()[num].size() << std::endl;
			std::cout << "siZe vecofvec =  " << vectorOfvector().size() << std::endl;*/

			po.SetValue(j + 1, all[i][j]/*.XYZ()*/);


		}

		std::cout << " // ********************* //  " << std::endl;
		std::cout << " //   *****************   //  " << std::endl;

		points.push_back(po);

	}



	for (int num = 0; num < points.size(); num++)
	{

		/*GeomAbs_C0,
		  GeomAbs_G1,
		  GeomAbs_C1,
		  GeomAbs_G2,
		  GeomAbs_C2,
		  GeomAbs_C3,
		  GeomAbs_CN*/

		GeomAPI_PointsToBSpline Bs;

		//int DegMax = theDegree_;
		//Bs = GeomAPI_PointsToBSpline(points[num],  1 , 2 , GeomAbs_C0, 0.00001);
		//Bs = GeomAPI_PointsToBSpline(points[num],  1 , 2, GeomAbs_G1, 0.00001);
		Bs = GeomAPI_PointsToBSpline(points[num], 1, 3, GeomAbs_G2, 0.1);
		//Bs = GeomAPI_PointsToBSpline(points[num],  1 , 2 , GeomAbs_C1, 0.1);
		//Bs = GeomAPI_PointsToBSpline(points[num],  1 , 2 , GeomAbs_C2, 1.0e-1);
		//Bs = GeomAPI_PointsToBSpline(points[num],  1 , 2 , GeomAbs_C3, 1.0e-3);
		//Bs = GeomAPI_PointsToBSpline(points[num], 1, 2 , GeomAbs_CN, 0.000001);

		const Handle(Geom_BSplineCurve)& cur = Bs.Curve();

		BScurve.push_back(cur);

	}


	NCollection_List < Handle(Geom_BSplineCurve)> col;
	//NCollection_TListIterator < Handle(Geom_BSplineCurve)> col0;

	for (int i = 0; i < BScurve.size(); i++)
	{

		col.Append(BScurve[i]);

	}

	/*Handle(Geom_Surface) ACISAlgo::MakeSkinSurface(const NCollection_List<Handle(Geom_Curve)>&theSections)
	{*/
	std::cout << "Ncollection_list size: " << col.Size() << std::endl;
	//populate section generator


	GeomFill_SectionGenerator aSecGenerator;

	for (NCollection_List<Handle(Geom_Curve)>::Iterator anIt(col); anIt.More(); anIt.Next())
	{
		const Handle(Geom_Curve)& aCurve = anIt.Value();
		aSecGenerator.AddCurve(aCurve);
	}
	aSecGenerator.Perform(Precision::PIntersection());
	Handle(GeomFill_Line) aLine = new GeomFill_Line(col.Size());

	//parameters
	const Standard_Integer aMinDeg = 1, aMaxDeg = BSplCLib::MaxDegree(), aNbIt = 0;
	Standard_Real aTol3d = 0.1, aTol2d = Precision::Parametric(aTol3d);

	//algorithm
	//GeomFill_AppSweep anAlgo(aMinDeg, aMaxDeg, aTol3d, aTol2d, aNbIt, false);
	GeomFill_AppSurf anAlgo(1, 2, 1e-6, 1e-7, aNbIt, false);

	//anAlgo.Init(1 , 3, aTol3d, aTol2d, aNbIt, false);
	//anAlgo.Perform(aLine, aSecGenerator);

	anAlgo.Perform(aLine, aSecGenerator);
	Handle(Geom_Surface) aRes;

	/*if (!anAlgo.IsDone())
	{
		return aRes;
	}*/

	aRes = new Geom_BSplineSurface(anAlgo.SurfPoles(), anAlgo.SurfWeights(),
		anAlgo.SurfUKnots(), anAlgo.SurfVKnots(), anAlgo.SurfUMults(), anAlgo.SurfVMults(),
		anAlgo.UDegree(), anAlgo.VDegree());


	BRepBuilderAPI_MakeFace face(aRes, 1, 10, 1, 10, 10);
	TopoDS_Shape shape = face.Shape();

	if (shape.IsNull())
	{
		throw std::exception("The Shape is null!");

	}
	
	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, shape);


	gp_Pnt po0(all[0][0].XYZ());
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, false);
	TopoDS_Shape halfshape = build.Shape();

	builder.Add(ship, halfshape);


	return ship;
}



void Offset::creatIGSfile(TopoDS_Shape shape )
{

	//TopoDS_Shape shapeforconvert = shape;



	IGESControl_Controller::Init();
	BRep_Builder builder;

	/*if (!BRepTools::Read(shape, "input.brep", builder))
	{
		std::cout << "Failed to read the input shape." << std::endl;
		return  ;
	}*/


	IGESControl_Writer writer(Interface_Static::CVal("XSTEP.iges.unit"),
		Interface_Static::IVal("XSTEP.iges.writebrep.mode"));

	writer.AddShape(shape);

	writer.ComputeModel();
	writer.Write("ship.igs");

	if (!writer.Write("ship.igs"))
	{
		throw std::exception("can't Write file");
	}
	//std::cout << "IGES file successfully created." << std::endl;

}


std::vector<gp_Pnt> Offset::Interpolate(const std::vector<gp_Pnt>& Q, int degree)
{
	auto knotsBar = CalcKnotBarCentripetal(Q);

	int n = Q.size() - 1;
	FrgBase_FullMtx<double> A(n + 1, n + 1);
	FrgBase_Vcr<double> rhsX(n + 1);
	FrgBase_Vcr<double> rhsY(n + 1);
	FrgBase_Vcr<double> rhsZ(n + 1);

	for (int i = 0; i < Q.size(); i++)
	{
		rhsX[i] = Q[i].X();
		rhsY[i] = Q[i].Y();


		rhsZ[i] = Q[i].Z();
	}

	auto knots = CalcKnotVector(Q, knotsBar, degree);

	for (int i = 0; i <= n; i++)
	{
		int span = FindSpan(n, degree, knotsBar[i], knots);
		auto N = CalcBasisFunctionVector(knots, span, degree, knotsBar[i]);
		for (int j = 0; j <= degree; j++)
			A[i][j + span - degree] = N[j];
	}

	A[0][0] = 1.0;
	A[n][n] = 1.0;

	A.GaussElim(rhsX);
	A.GaussElim(rhsY);
	A.GaussElim(rhsZ);

	std::vector<gp_Pnt> ctrlPts;
	for (int i = 0; i < rhsX.size(); i++)
	{

		ctrlPts.emplace_back(rhsX[i], rhsY[i], rhsZ[i]);
	}

	return std::move(ctrlPts);

}



std::vector<double> Offset::CalcKnotBarChordLength(const std::vector<gp_Pnt>& Q)
{
	std::vector<double> knots;
	knots.push_back(0.0);

	int n = Q.size() - 1;
	double d = 0.0;
	for (int k = 1; k <= n; k++)
		d += Q[k].Distance(Q[k - 1]);

	for (int k = 1; k <= n - 1; k++)
		knots.push_back(knots[k - 1] + (Q[k].Distance(Q[k - 1]) / d));

	knots.push_back(1.0);

	return std::move(knots);

}

std::vector<double> Offset::CalcKnotBarCentripetal(const std::vector<gp_Pnt>& Q)
{
	std::vector<double> knots;
	knots.push_back(0.0);

	int n = Q.size() - 1;
	double d = 0.0;
	for (int k = 1; k <= n; k++)
	{


		d += Q[k].Distance(Q[k - 1]);
	}

	for (int k = 1; k <= n - 1; k++)
		knots.push_back(knots[k - 1] + (Q[k].Distance(Q[k - 1]) / d));

	knots.push_back(1.0);

	return std::move(knots);
}

std::vector<double> Offset::CalcKnotBarEquallySpaced(const std::vector<gp_Pnt>& Q)
{
	std::vector<double> knots;
	knots.push_back(0.0);

	int n = Q.size() - 1;
	double d = 0.0;
	for (int k = 1; k <= n; k++)
	{
		d += (Q[k].Distance(Q[k - 1]));
	}
	for (int k = 1; k <= n - 1; k++)
		knots.push_back(knots[k - 1] + (Q[k].Distance(Q[k - 1]) / d));

	knots.push_back(1.0);

	return std::move(knots);

}

std::vector<double> Offset::CalcKnotVector
(const std::vector<gp_Pnt>& Q, const std::vector<double>& knotsBar, int degree)
{
	std::vector<double> knots;

	for (int i = 0; i <= degree; i++)
		knots.push_back(0.0);

	int n = Q.size() - 1;
	for (int j = 1; j <= n - degree; j++)
	{
		double A = 0.0;
		for (int i = j; i <= j + degree - 1; i++)
			A += knotsBar[i];
		knots.push_back(1.0 / (static_cast<double>(degree)) * A);
	}

	for (int i = 0; i <= degree; i++)
		knots.push_back(1.0);

	return std::move(knots);


}

int Offset::FindSpan(int n, int p, double u, const std::vector<double>& knots)
{
	if (u >= knots[n + 1])
		return n;

	int low = p;
	int high = n + 1;
	int mid = (low + high) / 2;
	while (u < knots[mid] || u >= knots[mid + 1])
	{
		if (u < knots[mid])
			high = mid;
		else
			low = mid;

		mid = (low + high) / 2.0;
	}
	return mid;


}

std::vector<double> Offset::CalcBasisFunctionVector
(const std::vector<double>& U, int i, int p, double u) const
{
	std::vector<double> N(p + 1);
	std::vector<double> left(p + 1);
	std::vector<double> right(p + 1);

	N[0] = 1.0;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		double saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			double temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}

	return N;


}

std::vector<int> Offset::CalculateMultiplicities
(const std::vector<gp_Pnt>& ctrlPts, int degree, bool isPeriodic)
{
	std::vector<int> mult;

	if (!isPeriodic)
	{
		int sizeOfOther = ctrlPts.size() - degree - 1;
		mult.resize(2 + sizeOfOther);

		for (int i = 1; i < mult.size() - 1; i++)
			mult[i] = 1;

		mult[0] = degree + 1;
		mult[mult.size() - 1] = degree + 1;
	}
	else
	{
		mult.push_back(degree);
		int sizeOfOther = ctrlPts.size() - mult[0];
		for (int i = 0; i < sizeOfOther; i++)
			mult.push_back(1);
		mult.push_back(degree);
	}

	return std::move(mult);

}

double Offset::CalcBasisFunction(const std::vector<double>& U, int i, int p, double u) const
{
	int m = U.size() - 1;

	if ((i == 0 && u == U[0]) || (i == m - p - 1 && u == U[m]))
	{
		return 1;
	}

	if (u < U[i] || u >= U[i + p + 1])
	{
		return 0;
	}

	std::vector<double> N(p + 1);

	for (int j = 0; j <= p; j++)
		if (u >= U[i + j] && u < U[i + j + 1]) N[j] = 1;
		else
			N[j] = 0;

	double saved, temp;
	double Uleft, Uright;

	for (int k = 1; k <= p; k++)
	{
		if (N[0] == 0.) saved = 0;
		else
		{
			saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
		}

		for (int j = 0; j < p - k + 1; j++)
		{
			Uleft = U[i + j + 1];
			Uright = U[i + j + k + 1];

			if (N[j + 1] == 0)
			{
				N[j] = saved;
				saved = 0;
			}
			else
			{
				temp = N[j + 1] / (Uright - Uleft);
				N[j] = saved + (Uright - u) * temp;
				saved = (u - Uleft) * temp;
			}
		}
	}
	return N[0];

	/*int m = knots.size() - 1;

	if (p == 0)
	{
		if (u >= knots[i] && u < knots[i + 1])
			return 1.0;
		else
			return 0.0;
	}

	double A = (u - knots[i]) / (knots[i + p] - knots[i]);
	double B = (knots[i + p + 1] - u) / (knots[i + p + 1] - knots[i + 1]);

	std::cout <<"u = " << u << ", A = " << A << ", B = " << B << std::endl;
	system("pause");

	return A * CalcBasisFunction(knots, i, p - 1, u) + B * CalcBasisFunction(knots, i + 1, p - 1, u);*/
}



std::vector<double> Offset::CalculateKnots(const std::vector<gp_Pnt>& ctrlPts, int degree, bool isPeriodic)
{
	std::vector<double> knots;

	if (!isPeriodic)
	{
		knots.resize(2 + (ctrlPts.size() - degree - 1));
		knots[0] = 0.0;
		knots[knots.size() - 1] = 1.0;

		double dK = 1.0 / ((double)knots.size() - 1.0);

		for (int i = 1; i < knots.size() - 1; i++)
			knots[i] = i * dK;
	}
	else
	{
		knots.push_back(0.0);
		int sizeOfOther = ctrlPts.size() - degree;
		double dK = 1.0 / ((double)sizeOfOther + 1.0);
		for (int i = 0; i < sizeOfOther; i++)
			knots.push_back((i + 1.0) * dK);
		knots.push_back(1.0);
	}

	return std::move(knots);

}

std::string Offset::doubleToString(double value)const
{
	std::vector<char> Allchar;
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(3) << value;
	std::string str = oss.str();
	for (int i = 0; i < str.size(); i++)
	{
		char ch = str[i];
		Allchar.push_back(ch);

	}
	std::vector<char> charVector = Allchar;
	std::string strfinall(charVector.begin(), charVector.end());

	return strfinall;
}


std::vector<Handle(Geom_Curve)> Offset::CreatCurveOffsetOneByInterpolateNewMethod()
{
	std::vector<Handle(Geom_Curve)> curvetotal0;
	std::vector<std::vector< gp_Pnt>> poinOff = ReadOffsetOne();

	//std::vector<gp_Pnt> points;


	for (int i = 0; i < poinOff.size(); i++)
	{

		/*for (int j = 0; j < poinOff[i].size(); j++)
		{
			gp_Pnt poi = (poinOff[i][j]);
			points.push_back(poi);

		}*/

		Handle(Geom_Curve) cur = SetData(Interpolate(poinOff[i], theDegree_), theDegree_);
		curvetotal0.push_back(cur);
		//points.clear();
	}

	if (curvetotal0.size() < 2)
	{
		throw std::exception(" Vector of Geom_Curve Size is Small");
	}

	return curvetotal0;

}





bool Offset::areEqual(const std::vector<gp_Pnt>& v1, const std::vector<gp_Pnt>& v2)
{
	if (v1.size() != v2.size())
		return false;

	for (int i = 0; i < v1.size(); ++i)
	{
		double x = v1[i].X();
		double y = v1[i].Y();
		double z = v1[i].Z();
		double xp = v2[i].X();
		double yp = v2[i].Y();
		double zp = v2[i].Z();
		// Assuming you have defined the equality operator for gp_pnt
		if (x != xp || y != yp || z != zp)
			return false;
	}

	return true;

}


TopoDS_Shape Offset::ShapeOfGeomTrimCurve()
{
	std::vector<Handle(Geom_Curve)> curveVector = CreatCurveOffsetOneByInterpolateNewMethod();
	TColGeom_SequenceOfCurve curveSeq;
	//Create a composite curve using the curves from the vector
	for (int i = 0; i < curveVector.size(); i++)
	{

		Standard_Real u1 = curveVector[i]->FirstParameter();
		Standard_Real u2 = curveVector[i]->LastParameter();

		//Create a Geom_TrimmedCurve using the composite curve and the parameter range
		Handle(Geom_TrimmedCurve) trimmedCurve = new Geom_TrimmedCurve(curveVector[i], u1, u2);


		curveSeq.Append(trimmedCurve);
	}
	// Geom_TrimmedCurve object that represents the trimmed composite curve.
	GeomFill_NSections fillOp(curveSeq);
	fillOp.ComputeSurface();

	auto aSurf = fillOp.BSplineSurface();

	// convert it into a TopoDS_Shape
	BRepBuilderAPI_MakeFace mkFace(aSurf, 1e-4);
	//if (!mkFace.IsDone())
	//{
	//	// Handle error
	//	std::cerr << "Error in face creation." << std::endl;
	//	return;
	//}

	TopoDS_Shape aShape = mkFace.Shape();

	if (aShape.IsNull())
	{
		throw std::exception("The Shape is null!");

	}

	BRep_Builder builder;
	TopoDS_Compound ship;
	builder.MakeCompound(ship);
	builder.Add(ship, aShape);


	gp_Pnt po0(curveVector[0]->Value(0).XYZ()*(-1));
	gp_Dir dir0(0.0, 1.0, 0.0);
	gp_Dir dir1(0.0, 0.0, 1.0);
	gp_Ax2 ax0(po0, dir0, dir1);
	gp_Trsf trsf2;
	trsf2.SetMirror(ax0);

	BRepBuilderAPI_Transform build(ship, trsf2, false);
	TopoDS_Shape halfshape = build.Shape();
	builder.Add(ship, halfshape);

	return ship;
}


std::vector<Handle(Geom_Curve)> Offset::CreatCurvByOCC()
{
	std::vector<TColgp_Array1OfPnt> points;
	std::vector<Handle(Geom_Curve)> BScurve;
	std::vector<std::vector<gp_Pnt>> all = ReadOffsetOne();


	for (int i = 0; i < all.size(); i++)
	{
		TColgp_Array1OfPnt po;

		po.Resize(1, all[i].size(), false);

		for (int j = 0; j < all[i].size(); j++)
		{

			po.SetValue(j + 1, all[i][j].XYZ());

		}

		points.push_back(po);
	}

	for (int num = 0; num < points.size(); num++)
	{

		GeomAPI_PointsToBSpline Bs;

		Bs = GeomAPI_PointsToBSpline(points[num],  1 , 1 , GeomAbs_C0, 1e-6);
		//Bs = GeomAPI_PointsToBSpline(points[num],  1 ,1, GeomAbs_G1, 1e-6);
		//Bs = GeomAPI_PointsToBSpline(points[num], 1, 2, GeomAbs_G2, 1e-4);
		//Bs = GeomAPI_PointsToBSpline(points[num], 1 , 2 , GeomAbs_C1, 1e-4);
	    //Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_C2, 1e-4);
		//Bs = GeomAPI_PointsToBSpline(points[num],  2 , 5 , GeomAbs_C3, 1e-4);
		//Bs = GeomAPI_PointsToBSpline(points[num], 1,2, GeomAbs_CN, 1e-2);

		const Handle(Geom_BSplineCurve)& cur = Bs.Curve();
		BScurve.push_back(cur);

	}

	if (BScurve.size() < 2)
	{
		throw std::exception(" Vector of Geom_Curve Size is Small");
	}

	return BScurve;
}


std::vector<std::vector<gp_Pnt>> Offset::ReadOffsetOneFindMidUp()const
{
	std::vector<double> value = FindeMidCb();

	//std::string str1 = doubleToString(value[0]);
	std::string str2 = doubleToString(value[1]);

	std::string theFileName_ = "cb=" + str2 + ".txt";
	std::cout << "file1:" << theFileName_ << std::endl;

	std::fstream fin(theFileName_, std::ios::in);

	if (!fin.is_open())
	{
		std::exception e("file not open");
		throw e;
	}

	std::string a;
	int numberofWaterline, numberofeStation;
	double  CB, CM, LCB;

	fin >> numberofWaterline;
	getline(fin, a);
	fin >> numberofeStation;

	for (int i = 0; i < 6; i++)
	{

		getline(fin, a);

	}


	fin >> CB;
	getline(fin, a);

	fin >> CM;
	getline(fin, a);

	fin >> LCB;
	getline(fin, a);

	std::cout << "CB up=" << CB << " , " << "CM up=" << CM << " , " << "LCB up=" << LCB << std::endl;

	std::vector<long double> Z_Core;
	std::vector<long double>X_Core;
	std::vector< std::vector<long float>>Y_Coretotal;

	for (int i = 0; i < numberofWaterline; i++)
	{
		double z;
		fin >> z;
		//double zPrim = z * ((draft / 1.5) / 4);

		Z_Core.push_back(z/*Prim*/);

	}


	for (int i = 0; i < 2; i++)
	{

		getline(fin, a);

	}



	for (int i = 0; i < numberofeStation; i++)
	{

		double x;
		fin >> x;
		//double xPrim = x * (length / 20);
		X_Core.push_back(x/*Prim*/);

		//std::cout << x << std::endl;

		std::vector< double>Y_Core;

		for (int j = 0; j < numberofWaterline; j++)
		{
			double y;
			fin >> y;
			// yPrim = y * beam;

			Y_Core.push_back(y/*Prim*/);

		}

		Y_Coretotal.push_back(Y_Core);
		getline(fin, a);


	}

	//fin.close();

	std::vector< std::vector<gp_Pnt> > pointss;


	for (int i = 0; i < numberofeStation; i++)
	{


		std::vector<gp_Pnt> points;
		double xdummy = X_Core[i];


		for (int j = 0; j < numberofWaterline; j++)
		{

			double Ysec = Y_Coretotal[i][j];
			double zdummy = Z_Core[j];

			/*if (j == 0)
			{
				Ysec = 0.0;
			}*/

			gp_Pnt point1(xdummy, Ysec, zdummy);

			points.push_back(point1);


		}


		pointss.push_back(points);

	}


	std::vector< std::vector<gp_Pnt> > pointss2;

	for (int i = 0; i < pointss.size(); i++)
	{

		std::vector<gp_Pnt> points0;
		double xdummy = X_Core[i];


		for (int j = 0; j < pointss[i].size() /*-1*/; j++)
		{

			double y1 = Y_Coretotal[i][j];


			if (std::abs(y1) != 0.0 && j != 0)
			{


				double Ysec000 = Y_Coretotal[i][j - 1];
				double zdummy = Z_Core[j - 1];

				gp_Pnt point00(xdummy, Ysec000, zdummy);
				points0.push_back(point00);

			}


		}

		int lastpoint = pointss[i].size() - 1;
		double Ysec00 = Y_Coretotal[i][lastpoint];
		double zdummy = Z_Core[lastpoint];
		gp_Pnt point00(xdummy, Ysec00, zdummy);
		points0.push_back(point00);
		pointss2.push_back(points0);

	}


	std::vector< std::vector<gp_Pnt> > pointssFinall;


	for (int i = 0; i < pointss2.size(); i++)
	{


		std::vector<gp_Pnt> points;
		long double xdummy = X_Core[i];


		for (int j = 0; j < pointss2[i].size(); j++)
		{

			long float Ysec0 = pointss2[i][j].Y();
			long double zdummy0 = pointss2[i][j].Z();

			if (j == 0)
			{
				Ysec0 = 0.0;
			}

			gp_Pnt point1(xdummy, Ysec0, zdummy0);

			points.push_back(point1);


		}

		pointssFinall.push_back(points);

	}
	std::cout << "sizeUP:" << pointssFinall.size() << std::endl;
	if (pointssFinall.size() == 0)
	{
		throw std::exception(" The size vector Of Points for up bound is Zero ");
	}
	return pointssFinall;

}


std::vector<std::vector<gp_Pnt>> Offset::ReadOffsetOneFindMidDown() const
{

	std::vector<double> value = FindeMidCb();

	std::string str1 = doubleToString(value[0]);
	//std::string str2 = doubleToString(value[1]);

	std::string theFileName_ = "cb=" + str1 + ".txt";
	std::cout << "file2:" << theFileName_ << std::endl;

	std::fstream fin(theFileName_, std::ios::in);

	if (!fin.is_open())
	{
		std::exception e("file not open");
		throw e;
	}

	std::string a;
	int numberofWaterline, numberofeStation;
	double  CB, CM, LCB;

	fin >> numberofWaterline;
	getline(fin, a);
	fin >> numberofeStation;

	for (int i = 0; i < 6; i++)
	{

		getline(fin, a);

	}


	fin >> CB;
	getline(fin, a);

	fin >> CM;
	getline(fin, a);

	fin >> LCB;
	getline(fin, a);

	std::cout << "CB down=" << CB << " , " << "CM down=" << CM << " , " << "LCB down=" << LCB << std::endl;

	std::vector<long double> Z_Core;
	std::vector<long double>X_Core;
	std::vector< std::vector<long float>>Y_Coretotal;

	for (int i = 0; i < numberofWaterline; i++)
	{
		double z;
		fin >> z;
		//double zPrim = z * ((draft / 1.5) / 4);

		Z_Core.push_back(z/*Prim*/);

	}


	for (int i = 0; i < 2; i++)
	{

		getline(fin, a);

	}



	for (int i = 0; i < numberofeStation; i++)
	{

		double x;
		fin >> x;
		//double xPrim = x * (length / 20);
		X_Core.push_back(x/*Prim*/);

		//std::cout << x << std::endl;

		std::vector< double>Y_Core;

		for (int j = 0; j < numberofWaterline; j++)
		{
			double y;
			fin >> y;
			// yPrim = y * beam;

			Y_Core.push_back(y/*Prim*/);

		}

		Y_Coretotal.push_back(Y_Core);
		getline(fin, a);


	}

	//fin.close();

	std::vector< std::vector<gp_Pnt> > pointss;


	for (int i = 0; i < numberofeStation; i++)
	{


		std::vector<gp_Pnt> points;
		double xdummy = X_Core[i];


		for (int j = 0; j < numberofWaterline; j++)
		{

			double Ysec = Y_Coretotal[i][j];
			double zdummy = Z_Core[j];

			/*if (j == 0)
			{
				Ysec = 0.0;
			}*/

			gp_Pnt point1(xdummy, Ysec, zdummy);

			points.push_back(point1);


		}


		pointss.push_back(points);

	}


	std::vector< std::vector<gp_Pnt> > pointss2;

	for (int i = 0; i < pointss.size(); i++)
	{


		std::vector<gp_Pnt> points0;
		double xdummy = X_Core[i];


		for (int j = 0; j < pointss[i].size() /*-1*/; j++)
		{

			double y1 = Y_Coretotal[i][j];


			if (std::abs(y1) != 0.0 && j != 0)
			{


				double Ysec000 = Y_Coretotal[i][j - 1];
				double zdummy = Z_Core[j - 1];

				gp_Pnt point00(xdummy, Ysec000, zdummy);
				points0.push_back(point00);

			}


		}

		int lastpoint = pointss[i].size() - 1;
		double Ysec00 = Y_Coretotal[i][lastpoint];
		double zdummy = Z_Core[lastpoint];
		gp_Pnt point00(xdummy, Ysec00, zdummy);
		points0.push_back(point00);
		pointss2.push_back(points0);

	}


	std::vector< std::vector<gp_Pnt> > pointssFinall;


	for (int i = 0; i < pointss2.size(); i++)
	{


		std::vector<gp_Pnt> points;
		long double xdummy = X_Core[i];


		for (int j = 0; j < pointss2[i].size(); j++)
		{

			long float Ysec0 = pointss2[i][j].Y();
			long double zdummy0 = pointss2[i][j].Z();

			if (j == 0)
			{
				Ysec0 = 0.0;
			}

			gp_Pnt point1(xdummy, Ysec0, zdummy0);

			points.push_back(point1);


		}


		pointssFinall.push_back(points);
		//.close();
	}
	std::cout << "sizeDOWN:" << pointssFinall.size() << std::endl;
	if (pointssFinall.size() == 0)
	{
		throw std::exception(" The size vector Of Points for down bound is Zero ");
	}
	return pointssFinall;
}






std::vector<double> Offset::FindeMidCb(/*double value*/) const
{
	std::vector<double> cbvalue{ 0.600,0.610,0.620,0.630,0.640,0.650,0.650,0.660,0.670,0.670
	,0.680,0.690,0.700,0.710,0.720,0.730,0.740,0.750,0.760,0.770,0.780,0.790,0.800 };

	std::vector<double> results;

	for (int i = 0.0; i < cbvalue.size(); i++)
	{
		double per = theCB_ - cbvalue[i];
		if (per < 0.0)
		{

			results.push_back(cbvalue[i - 1]);
			results.push_back(cbvalue[i]);
			break;
		}

	}


	return results;
}


bool Offset::HandeValueOfLCB()
{

	double valuelcb = theLCB_;


	if (GetCB() == 0.600 || 0.610)
	{

		if (valuelcb == -1.000 || -1.500 || -2.000)
		{
			return true;

		}

		//return false;
	}

	return false;
	if (GetCB() == 0.620 || 0.630 || 0.640)
	{

		if (valuelcb == -0.500 || -1.500 || -1.000)
		{

			return true;
		}
		//return false;
	}
	return false;
	if (GetCB() == 0.650 || 0.660 || 0.670 || 0.680)
	{
		if (valuelcb == 0.000 || -0.500 || -1.000)
		{

			return true;
		}
		//return false;
	}
	return false;
	if (GetCB() == 0.690 || 0.700)
	{
		if (valuelcb == -0.500 || 0.000 || 0.500)
		{

			return true;
		}

		//return false;
	}
	return false;
	if (GetCB() == 0.710 || 0.720)
	{
		if (valuelcb == 0.000 || 0.500 || 1.000)
		{

			return true;
		}

		//return false;
	}
	return false;
	if (GetCB() == 0.730 || 0.740 || 0.750)
	{

		if (valuelcb == 0.500 || 1.500 || 1.000)
		{

			return true;
		}

		//return false;
	}
	return false;

	if (GetCB() == 0.760 || 0.770 || 0.780 || 0.790)
	{

		if (valuelcb == 1.000 || 1.500 || 2.000)
		{

			return true;
		}

		//return false;
	}
	return false;

	if (GetCB() == 0.800)
	{

		if (valuelcb == 2.000 || 1.500 || 2.500)
		{

			return true;
		}

		//return false;
	}

	return false;


}
