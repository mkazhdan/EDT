/*
Copyright (c) 2021, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "Misha/CmdLineParser.h"
#include "Misha/Miscellany.h"
#include "Misha/RegularGrid.h"
#include "Misha/IsoSurface3D.h"
#include "Misha/MarchingCubes.h"
#include "Misha/IsoSurface3D.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineReadable FullCaseTable( "full" ) ,  Polygons( "polygons" ) , NonManifold( "nonManifold" );
Misha::CmdLineParameter< float > Value( "value" , 0.f );
Misha::CmdLineParameter< int > InterpolationType( "fit" , IsoSurface3D< float >::INTERPOLATE_LINEAR );


Misha::CmdLineReadable* params[] =
{
	&In , &Out , &Value , &FullCaseTable , &InterpolationType , &Polygons , &NonManifold , NULL
};

void ShowUsage( char* ex )
{
	std::cout << "Usage: " << std::string(ex) << std::endl;
	std::cout << "\t --" << In.name  << " <input image list>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output mesh>]" << std::endl;
	std::cout << "\t[--" << Value.name << " <iso-value>=" << Value.value << "]" << std::endl;
	std::cout << "\t[--" << InterpolationType.name << " <interpolation type>=" << InterpolationType.value << "]" << std::endl;
	for( int i=0 ; i<IsoSurface3D< float >::INTERPOLATE_COUNT ; i++ ) std::cout << "\t\t" << i << "] " << IsoSurface3D< float >::InterpolationNames[i] << std::endl;
	std::cout << "\t[--" << Polygons.name << "]" << std::endl;
	std::cout << "\t[--" << FullCaseTable.name << "]" << std::endl;
	std::cout << "\t[--" << NonManifold.name << "]" << std::endl;
}

typedef VertexFactory::PositionFactory< float , 3 > Factory;
typedef typename Factory::VertexType Vertex;

int main( int argc , char* argv[] )
{
	Miscellany::Timer timer;

	Misha::CmdLineParse( argc-1 , &argv[1] , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	XForm< float , 4 > gridToModel;
	RegularGrid< float , 3 > grid;
	grid.read( In.value , gridToModel );


	std::vector< Vertex > vertices;
	std::vector< std::vector< int > > polygons;
	std::vector< TriangleIndex > triangles;
	if( Polygons.set ) IsoSurface3D< float >::Extract( grid , Value.value , vertices , polygons , FullCaseTable.set , InterpolationType.value );
	else               IsoSurface3D< float >::Extract( grid , Value.value , vertices , triangles , FullCaseTable.set , InterpolationType.value , !NonManifold.set );

	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = gridToModel * vertices[i];

	if( Out.set )
		if( Polygons.set ) PLY::WritePolygons( Out.value , Factory() , vertices , polygons , PLY_BINARY_NATIVE );
		else               PLY::WriteTriangles( Out.value , Factory() , vertices , triangles , PLY_BINARY_NATIVE );

	return EXIT_SUCCESS;
}
