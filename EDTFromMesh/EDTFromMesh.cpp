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
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Miscellany.h"
#include "Misha/SquaredEDT.h"
#include "Misha/RegularGrid.h"
#include "Misha/Rasterizer.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Ply.h"

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< float > Scale( "scale" , 2.f );
Misha::CmdLineParameter< int > Depth( "depth" , 8 ) , Radius( "radius" , -1 );
Misha::CmdLineParameter< int > LockDepth( "lockDepth" , 4 ) , Threads( "threads" , omp_get_num_procs() );
Misha::CmdLineReadable Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In , &Out , &Verbose , &Threads , &Scale , &Depth , &LockDepth , &Radius , NULL
};

void ShowUsage( char* ex )
{
	std::cout << "Usage: " << std::string(ex) << std::endl;
	std::cout << "\t --" << In.name  << " <input image list>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output signed EDT>]" << std::endl;
	std::cout << "\t[--" << Depth.name << " <rasterization depth>=" << Depth.value << "]" << std::endl;
	std::cout << "\t[--" << Scale.name << " <model scale>=" << Scale.value << "]" << std::endl;
	std::cout << "\t[--" << Radius.name << " <distance radius ( less than zero ? Saito : Danielsson )=" << Radius.value << "]" << std::endl;
#ifdef _OPENMP
	std::cout << "\t[--" << Threads.name << " <num threads>=" << Threads.value << "]" << std::endl;
	std::cout << "\t[--" << LockDepth.name << " <lock depth>=" << LockDepth.value << "]" << std::endl;
#endif // _OPENMP
	std::cout << "\t[--" << Verbose.name << "]" << std::endl;
}


RegularGrid< float , 3 > GetEDT( const SimplicialComplex< double , 3 , 2 > &sc , unsigned int depth , unsigned int lockDepth , int radius , double scale , XForm< float , 4 > &gridToModel , bool verbose )
{
	RegularGrid< float , 3 > edt;
	if( radius<0 )
	{
		XForm< double , 4 > _gridToModel;
		RegularGrid< unsigned , 3 > _edt2 = Misha::SquaredEDT< double , 3 >::Saito< int , 2 >( sc , depth , lockDepth , _gridToModel , scale , verbose );
		edt.resize( _edt2.res() );
		for( unsigned int i=0 ; i<edt.res(0) ; i++ ) for( unsigned int j=0 ; j<edt.res(1) ; j++ ) for( unsigned int k=0 ; k<edt.res(2) ; k++ )
			edt(i,j,k) = (float)sqrt( _edt2(i,j,k) );
		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ ) gridToModel(i,j) = (float)_gridToModel(i,j);
	}
	else
	{
		XForm< double , 4 > _gridToModel;
		RegularGrid< double , 3 > _edt2 = Misha::SquaredEDT< double , 3 >::Danielsson< int , 2 >( sc , depth , lockDepth , radius , _gridToModel , scale , verbose );
		edt.resize( _edt2.res() );
		for( unsigned int i=0 ; i<edt.res(0) ; i++ ) for( unsigned int j=0 ; j<edt.res(1) ; j++ ) for( unsigned int k=0 ; k<edt.res(2) ; k++ )
			edt(i,j,k) = (float)sqrt( _edt2(i,j,k) ) * (1<<depth);
		for( int i=0 ; i<4 ; i++ ) for( int j=0 ; j<4 ; j++ ) gridToModel(i,j) = (float)_gridToModel(i,j);
	}

	return edt;
}


int main( int argc , char* argv[] )
{
	double t = Miscellany::Time();

	Misha::CmdLineParse( argc-1 , &argv[1] , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	omp_set_num_threads( Threads.value );


	typedef VertexFactory::PositionFactory< double , 3 > Factory;
	typedef typename Factory::VertexType Vertex;

	std::vector< Vertex > vertices;
	std::vector< TriangleIndex > triangles;

	int file_type;
	PLY::ReadTriangles( In.value , Factory() , vertices , triangles , NULL , file_type );
	std::cout << "Vertices / Triangles: " << vertices.size() << " / " << triangles.size() << std::endl;

	IndexedSimplicialComplex< double , 3 , 2 , int > isc( vertices , triangles );

	XForm< float , 4 > gridToModel;
	RegularGrid< float , 3 > EDT = GetEDT( isc , Depth.value , LockDepth.value , Radius.value , Scale.value , gridToModel, Verbose.set );

	if( Out.set ) EDT.write( Out.value , gridToModel );

	return EXIT_SUCCESS;
}
