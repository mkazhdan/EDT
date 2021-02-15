/*
Copyright (c) 2015, Michael Kazhdan
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
#ifndef ISO_SURFACE_3D_INCLUDED
#define ISO_SURFACE_3D_INCLUDED

#include <omp.h>
#include <unordered_map>
#include <algorithm>
#include "Geometry.h"
#include "MarchingCubes.h"
#include "RegularGrid.h"
#include "Array.h"
#include "Polynomial.h"


template< typename Real >
struct IsoSurface3D
{
	enum
	{
		INTERPOLATE_LINEAR ,
		INTERPOLATE_QUADRATIC ,
		INTERPOLATE_CUBIC ,
		INTERPOLATE_CATMULL_ROM ,
		INTERPOLATE_COUNT
	};
	static const std::string InterpolationNames[];

	static void Extract( const unsigned int res[3] , ConstPointer( Real ) values , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType );
	static void Extract( const unsigned int res[3] , ConstPointer( Real ) values , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< TriangleIndex >& triangles , bool fullCaseTable , int interpolationType , bool manifold );

	static void Extract( const RegularGrid< Real , 3 > &voxelGrid , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType );
	static void Extract( const RegularGrid< Real , 3 > &voxelGrid , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< TriangleIndex >& triangles , bool fullCaseTable , int interpolationType , bool manifold );

protected:
	struct _Vertex
	{
		int dir , idx[3];
		Point3D< Real > p;
		_Vertex( Point3D< Real > _p , int _dir , int x , int y , int z ) : p(_p) , dir(_dir) { idx[0] = x , idx[1] = y , idx[2] = z; }
		static bool CoFacial( const _Vertex &t1 , const _Vertex &t2 );
	};
	static void _Extract( const RegularGrid< Real , 3 > &voxelGrid , Real isoValue , std::vector< _Vertex >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType );
	static void _Extract( const unsigned int res[3] , ConstPointer( Real ) values , Real isoValue , std::vector< _Vertex >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType );

	static Real     _LinearInterpolant( Real x1 , Real x2 , Real isoValue );
	static Real  _QuadraticInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue );
	static Real      _CubicInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue );
	static Real _CatmullRomInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue );

	static void _SetFlags( int resX , int resY , ConstPointer( Real ) values , Real isoValue , Pointer( unsigned char ) flags );
	static void _SetZVertices( int resX , int resY , int z , ConstPointer( Real ) values0 , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , ConstPointer( Real ) values3 , ConstPointer( unsigned char ) flags1 , ConstPointer( unsigned char ) flags2 , Real isoValue , int interpolationType , std::unordered_map< long long , int >& isoVertexMap , std::vector< _Vertex >& vertices );
	static void _SetXYVertices( int resX , int resY , int z , ConstPointer( Real ) values , ConstPointer( unsigned char ) flags , Real isoValue , int interpolationType , std::unordered_map< long long , int >& xIsoVertexMap , std::unordered_map< long long , int >& yIsoVertexMap , std::vector< _Vertex >& vertices );
	static void _SetPolygons( int resX , int resY , int z , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , Real isoValue , bool fullCaseTable , const std::unordered_map< long long , int >& xIsoVertexMap1 , const std::unordered_map< long long , int >& xIsoVertexMap2 , const std::unordered_map< long long , int >& yIsoVertexMap1 , const std::unordered_map< long long , int >& yIsoVertexMap2 , const std::unordered_map< long long , int >& zIsoVertexMap , const std::vector< _Vertex >& vertices , std::vector< std::vector< int > >& polygons );
};

#include "IsoSurface3D.inl"
#endif // ISO_SURFACE_3D_INCLUDED