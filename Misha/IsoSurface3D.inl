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

#include "Polynomial.h"

///////////////////////////
// IsoSurface3D::_Vertex //
///////////////////////////
template< typename Real >
bool IsoSurface3D< Real >::_Vertex::CoFacial( const _Vertex &t1 , const _Vertex &t2 )
{
#define _ABS_( a ) ( (a)<0 ? -(a) : (a) )
	int d[] = { _ABS_( t1.idx[0] - t2.idx[0] ) , _ABS_( t1.idx[1] - t2.idx[1] ) , _ABS_( t1.idx[2] - t2.idx[2] ) };
	if( t1.dir==t2.dir ) return d[t1.dir]==0 && ( ( d[(t1.dir+1)%3]==0 && d[(t1.dir+2)%3]<=1 ) || ( d[(t1.dir+2)%3]==0 && d[(t1.dir+1)%3]<=1 ) );
	else                 return d[ 3 - t1.dir - t2.dir ]==0 && d[t1.dir]<=1 && d[t2.dir]<=1;
#undef _ABS_
}

//////////////////
// IsoSurface3D //
//////////////////
template< typename Real >
const std::string IsoSurface3D< Real >::InterpolationNames[] = { "linear" , "quadratic" , "cubic" , "catmull-rom" };

template< typename Real >
void IsoSurface3D< Real >::Extract( const unsigned int res[3] , ConstPointer( Real ) values , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType )
{
	std::vector< _Vertex > _vertices;
	_Extract( res , values , isoValue , _vertices , polygons , fullCaseTable , interpolationType );

	vertices.resize( _vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = _vertices[i].p;
}
template< typename Real >
void IsoSurface3D< Real >::Extract( const unsigned int res[3] , ConstPointer( Real ) values , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< TriangleIndex >& triangles , bool fullCaseTable , int interpolationType , bool manifold  )
{
	std::vector< _Vertex > _vertices;
	std::vector< std::vector< int > > polygons;
	_Extract( res , values , isoValue , _vertices , polygons , fullCaseTable , interpolationType );
	vertices.resize( _vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = _vertices[i].p;

	MinimalAreaTriangulation< Real , 3 > mat;

	for( int i=0 ; i<polygons.size() ; i++ )
	{
		// To ensure that we have no more than two triangles adjacent on an edge,
		// we avoid creating a minimial area triangulation when it could introduce a new
		// edge that is on a face of a cube
		bool isCofacial = false;
		if( manifold )
			for( int j=0 ; j<(int)polygons[i].size() ; j++ ) for( int k=0 ; k<j ; k++ )
				if( (j+1)%polygons[i].size()!=k && (k+1)%polygons[i].size()!=j )
					if( _Vertex::CoFacial( _vertices[ polygons[i][j] ] , _vertices[ polygons[i][k] ] ) ) isCofacial = true;
		if( isCofacial )
		{
			TriangleIndex triangle;
			Point3D< Real > v;
			for( int j=0 ; j<(int)polygons[i].size() ; j++ ) v += vertices[ polygons[i][j] ];
			v /= (Real)polygons[i].size();
			int cIdx = (int)vertices.size();
			vertices.push_back( v );
			for( int j=0 ; j<(int)polygons[i].size() ; j++ )
			{
				triangle[0] = polygons[i][j];
				triangle[1] = polygons[i][(j+1)%polygons[i].size()];
				triangle[2] = cIdx;
				triangles.push_back( triangle );
			}
		}
		else
		{
			std::vector< Point3D< Real > > polygon( polygons[i].size() );
			std::vector< TriangleIndex > pTriangles;
			for( int j=0 ; j<polygons[i].size() ; j++ ) polygon[j] = vertices[ polygons[i][j] ];
			mat.GetTriangulation( polygon , pTriangles );
			for( int j=0 ; j<pTriangles.size() ; j++ )
			{
				TriangleIndex tri;
				for( int k=0 ; k<3 ; k++ ) tri[k] = polygons[i][ pTriangles[j][k] ];
				triangles.push_back( tri );
			}
		}
	}
}

template< typename Real >
void IsoSurface3D< Real >::Extract( const RegularGrid< Real , 3 > &voxelGrid , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType )
{
	return Extract( voxelGrid.res() , voxelGrid() , isoValue , vertices , polygons , fullCaseTable , interpolationType );
}
template< typename Real >
void IsoSurface3D< Real >::Extract( const RegularGrid< Real , 3 > &voxelGrid , Real isoValue , std::vector< Point3D< Real > >& vertices , std::vector< TriangleIndex >& triangles , bool fullCaseTable , int interpolationType , bool manifold  )
{
	return Extract( voxelGrid.res() , voxelGrid() , isoValue , vertices , triangles , fullCaseTable , interpolationType , manifold );
}

template< typename Real >
void IsoSurface3D< Real >::_Extract( const RegularGrid< Real , 3 > &voxelGrid , Real isoValue , std::vector< _Vertex >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType )
{
	_Extract( voxelGrid.res() , voxelGrid() , isoValue , vertices , polygons , fullCaseTable , interpolationType );
}

template< typename Real >
void IsoSurface3D< Real >::_Extract( const unsigned int res[3] , ConstPointer( Real ) values , Real isoValue , std::vector< _Vertex >& vertices , std::vector< std::vector< int > >& polygons , bool fullCaseTable , int interpolationType )
{
	std::unordered_map< long long , int > xIsoVertexMap[2] , yIsoVertexMap[2] , zIsoVertexMap;
	Pointer( unsigned char ) flags[2];
	flags[0] = NewPointer< unsigned char >( res[0]*res[1] );
	flags[1] = NewPointer< unsigned char >( res[0]*res[1] );

	if( fullCaseTable ) MarchingCubes::SetFullCaseTable();
	else                MarchingCubes::SetCaseTable();

	_SetFlags     ( res[0] , res[1] ,     values , isoValue , flags[0] );
	_SetXYVertices( res[0] , res[1] , 0 , values ,            flags[0] , isoValue , interpolationType , xIsoVertexMap[0] , yIsoVertexMap[0] , vertices );
	for( int z=0 ; z<(int)res[2]-1 ; z++ )
	{
		int z0 = z&1 , z1 = (z+1)&1;
		xIsoVertexMap[z1].clear() , yIsoVertexMap[z1].clear() , zIsoVertexMap.clear();
		_SetFlags     ( res[0] , res[1] ,       values + (z+1)*res[0]*res[1] , isoValue , flags[z1] );
		_SetXYVertices( res[0] , res[1] , z+1 , values + (z+1)*res[0]*res[1] ,            flags[z1] , isoValue , interpolationType , xIsoVertexMap[z1] , yIsoVertexMap[z1] , vertices );
		_SetZVertices ( res[0] , res[1] , z , z>0 ? values + (z-1)*res[0]*res[1] : NullPointer< Real >() , values + z*res[0]*res[1] , values + (z+1)*res[0]*res[1] , z+1<(int)res[2]-1 ? values + (z+2)*res[0]*res[1] : NullPointer< Real >() , flags[z0] , flags[z1] , isoValue , interpolationType , zIsoVertexMap , vertices );
		_SetPolygons  ( res[0] , res[1] , z , values + z*res[0]*res[1] , values+(z+1)*res[0]*res[1] , isoValue , fullCaseTable , xIsoVertexMap[z0] , xIsoVertexMap[z1] , yIsoVertexMap[z0] , yIsoVertexMap[z1] , zIsoVertexMap ,  vertices , polygons );
	}
	DeletePointer( flags[0] );
	DeletePointer( flags[1] );
}

template< typename Real >
Real IsoSurface3D< Real >::_LinearInterpolant( Real x1 , Real x2 , Real isoValue ){ return ( isoValue-x1 ) / ( x2-x1 ); }
template< typename Real >
Real IsoSurface3D< Real >::_QuadraticInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue )
{
	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	if( !x1 ) return 0;
	if( !x2 ) return 1;

	// Estimate the derivatives at x1 and x2
	Real dx1 = (x2-x0) / 2.f , dx2 = (x3-x1) / 2.f;
	// Solve for the quadratic polynomial:
	//		P(x) = a x^2 + b x + c 
	// such that:
	//		P(0) = x1 , P(1) = x2 , and minimizing || P'(0) - dx1 ||^2 + || P'(1) - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || 2*a + b - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || 2*x2 - 2*x1 - b - dx2 ||^2
	//	=>  c = x1 , a = x2 - x1 - b , and minimizing || b - dx1 ||^2 + || b - ( 2*x2 - 2*x1 - dx2 ) ||^2
	//	=>  c = x1 , b = ( 2*x2 - 2*x1 - dx2 + dx1 ) / 2 , a = x2 - x1 - b
	//	=>  c = x1 , b = ( x2 - x1 ) - ( dx2 - dx1 ) / 2 , a = ( dx2 - dx1 ) / 2

	double a = (dx2-dx1)/2.f , b = (dx1-dx2)/2.f + x2 - x1 , c = x1;
	if( !a )
	{
		// Solve b * x + c = 0
		return (Real)( -c / b );
	}
	else
	{
		// Solve a x^2 + b x + c = 0
		b /= a , c /= a;
		double disc = b*b - 4.*c;
		if( disc<0 ) ERROR_OUT( "Negative discriminant: " , disc );
		disc = sqrt( disc );
		double r1 = ( - b - disc ) / 2. , r2 = ( - b + disc ) / 2.;
		if( r2<0 || r1>1 ) ERROR_OUT( "Roots out of bounds: " , r1 , " " , r2 );
		if( r2>1 ) return (Real)r1;
		else       return (Real)r2;
	}
}
template< typename Real >
Real IsoSurface3D< Real >::_CubicInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue )
{
	static bool firstTime = true;
	static Polynomial::Polynomial1D< 3 > lagrangePolynomials[4];
	if( firstTime )
	{
#if 1
		Point< double , 1 > positions[] = { Point< double , 1 >(-1.) , Point< double , 1 >(0.) , Point< double , 1 >(1.) , Point< double , 1 >(2.) };
#else
		double positions[] = { -1. , 0. , 1. , 2. };
#endif
		SquareMatrix< double , 4 > Einv = Polynomial::Polynomial1D< 3 >::EvaluationMatrix( positions ).inverse();
		for( int i=0 ; i<4 ; i++ )
		{
			Point< double , 4 > values;
			values[i] = 1;
			Point< double , 4 > coefficients = Einv * values;
#if 1
			for( int j=0 ; j<4 ; j++ ) lagrangePolynomials[i].coefficient(j) = coefficients[j];
#else
			for( int j=0 ; j<4 ; j++ ) lagrangePolynomials[i][j] = coefficients[j];
#endif
		}
	}

	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	if( !x1 ) return 0;
	if( !x2 ) return 1;

	Polynomial::Polynomial1D< 3 > p = lagrangePolynomials[0] * x0 + lagrangePolynomials[1] * x1 + lagrangePolynomials[2] * x2 + lagrangePolynomials[3] * x3;
	double roots[3] , _roots[3];
#if 1
	int rootCount = Polynomial::Roots( p , roots );
#else
	int rootCount = p.roots( roots );
#endif
	int _rootCount = 0;
	for( int i=0 ; i<rootCount ; i++ ) if( roots[i]>=0 && roots[i]<1 ) _roots[ _rootCount++ ] = roots[i];
	if     ( _rootCount==1 ) return (Real)_roots[0];
	else if( _rootCount==3 ) return (Real)_roots[1];
	else ERROR_OUT( "Unexpected number of roots: " , _rootCount );
	return 0;
}

template< typename Real >
Real IsoSurface3D< Real >::_CatmullRomInterpolant( Real x0 , Real x1 , Real x2 , Real x3 , Real isoValue )
{
	static bool firstTime = true;
	static Polynomial::Polynomial1D< 3 > blendingFunctions[4];
	if( firstTime )
	{
#if 1
		blendingFunctions[0].coefficient(0) =  0.0;
		blendingFunctions[0].coefficient(1) = -0.5;
		blendingFunctions[0].coefficient(2) =  1.0;
		blendingFunctions[0].coefficient(3) = -0.5;

		blendingFunctions[1].coefficient(0) =  1.0;
		blendingFunctions[1].coefficient(1) =  0.0;
		blendingFunctions[1].coefficient(2) = -2.5;
		blendingFunctions[1].coefficient(3) =  1.5;

		blendingFunctions[2].coefficient(0) =  0.0;
		blendingFunctions[2].coefficient(1) =  0.5;
		blendingFunctions[2].coefficient(2) =  2.0;
		blendingFunctions[2].coefficient(3) = -1.5;

		blendingFunctions[3].coefficient(0) =  0.0;
		blendingFunctions[3].coefficient(1) =  0.0;
		blendingFunctions[3].coefficient(2) = -0.5;
		blendingFunctions[3].coefficient(3) =  0.5;
#else
		blendingFunctions[0] = Polynomial::Polynomial1D< 3 >( 0.0 , -0.5 ,  1.0 , -0.5 );
		blendingFunctions[1] = Polynomial::Polynomial1D< 3 >( 1.0 ,  0.0 , -2.5 ,  1.5 );
		blendingFunctions[2] = Polynomial::Polynomial1D< 3 >( 0.0 ,  0.5 ,  2.0 , -1.5 );
		blendingFunctions[3] = Polynomial::Polynomial1D< 3 >( 0.0 ,  0.0 , -0.5 ,  0.5 );
#endif
	}

	// Adjust so that we are looking for a zero-crossing
	x0 -= isoValue , x1 -= isoValue , x2 -= isoValue , x3 -= isoValue;
	if( !x1 ) return 0;
	if( !x2 ) return 1;

	Polynomial::Polynomial1D< 3 > p = blendingFunctions[0] * x0 + blendingFunctions[1] * x1 + blendingFunctions[2] * x2 + blendingFunctions[3] * x3;
	double roots[3] , _roots[3];
#if 1
	int rootCount = Polynomial::Roots( p , roots );
#else
	int rootCount = p.roots( roots );
#endif
	int _rootCount = 0;
	for( int i=0 ; i<rootCount ; i++ ) if( roots[i]>=0 && roots[i]<1 ) _roots[ _rootCount++ ] = roots[i];
	if     ( _rootCount==1 ) return (Real)_roots[0];
	else if( _rootCount==3 ) return (Real)_roots[1];
	else
	{
		std::cout << p << std::endl;
		printf( "Values: %g %g %g %g\n" , x0 , x1 , x2 , x3 );
		printf( "Roots:" ) ; for( int i=0 ; i<rootCount ; i++ ) printf( " %g" , roots[i] ) ; printf( "\n" );
		ERROR_OUT( "Unexpected number of roots: " , _rootCount );
	}
	return 0;
}

template< typename Real >
void IsoSurface3D< Real >::_SetFlags( int resX , int resY , ConstPointer( Real ) values , Real isoValue , Pointer( unsigned char ) flags )
{
#pragma omp parallel for
	for( int i=0 ; i<resX*resY ; i++ ) flags[i] = MarchingCubes::ValueLabel( values[i] , isoValue );
}

template< typename Real >
void IsoSurface3D< Real >::_SetZVertices( int resX , int resY , int z , ConstPointer( Real ) values0 , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , ConstPointer( Real ) values3 , ConstPointer( unsigned char ) flags1 , ConstPointer( unsigned char ) flags2 , Real isoValue , int interpolationType , std::unordered_map< long long , int >& isoVertexMap , std::vector< _Vertex >& vertices )
{
#define INDEX( x , y ) ( x + (y)*resX )
#pragma omp parallel for
	for( int i=0 ; i<resX ; i++ ) for( int j=0 ; j<resY ; j++ )
	{
		int idx = INDEX( i , j );
		if( flags1[idx]!=flags2[idx] )
		{
			Real iso;
			switch( interpolationType )
			{
			case INTERPOLATE_LINEAR:
				iso = _LinearInterpolant( values1[idx] , values2[idx] , isoValue );
				break;
			case INTERPOLATE_QUADRATIC:
				iso = _QuadraticInterpolant( values0 ? values0[idx] : values1[idx] , values1[idx] , values2[idx] , values3 ? values3[idx] : values2[idx] , isoValue );
				break;
			case INTERPOLATE_CUBIC:
				iso = _CubicInterpolant( values0 ? values0[idx] : values1[idx] , values1[idx] , values2[idx] , values3 ? values3[idx] : values2[idx] , isoValue );
				break;
			case INTERPOLATE_CATMULL_ROM:
				iso = _CatmullRomInterpolant( values0 ? values0[idx] : values1[idx] , values1[idx] , values2[idx] , values3 ? values3[idx] : values2[idx] , isoValue );
				break;
			default:
				ERROR_OUT( "Unrecognized interpolation type: " , interpolationType );
			}
			Point3D< Real > p = Point3D< Real >( (Real)i , (Real)j , (Real)z + iso );
			long long key = i + j*(resX);
#pragma omp critical
			{
				isoVertexMap[key] = (int)vertices.size();
				vertices.push_back( _Vertex( p , 2 , i , j , z ) );
			}
		}
	}
#undef INDEX
}

template< typename Real >
void IsoSurface3D< Real >::_SetXYVertices( int resX , int resY , int z , ConstPointer( Real ) values , ConstPointer( unsigned char ) flags , Real isoValue , int interpolationType , std::unordered_map< long long , int >& xIsoVertexMap , std::unordered_map< long long , int >& yIsoVertexMap , std::vector< _Vertex >& vertices )
{
#define INDEX( x , y ) ( x + (y)*resX )
#pragma omp parallel for
	for( int i=0 ; i<resX-1 ; i++ ) for( int j=0 ; j<resY ; j++ )
	{
		int idx1 = INDEX( i , j ) , idx2 = INDEX( i+1 , j );
		if( flags[idx1]!=flags[idx2] )
		{
			Real iso;
			switch( interpolationType )
			{
			case INTERPOLATE_LINEAR:
				iso = _LinearInterpolant( values[idx1] , values[idx2] , isoValue );
				break;
			case INTERPOLATE_QUADRATIC:
				iso = _QuadraticInterpolant( i>0 ? values[ INDEX(i-1,j) ] : values[idx1] , values[idx1] , values[idx2] , i+1<resX-1 ? values[ INDEX(i+2,j) ] : values[idx2] , isoValue );
				break;
			case INTERPOLATE_CUBIC:
				iso = _CubicInterpolant( i>0 ? values[ INDEX(i-1,j) ] : values[idx1] , values[idx1] , values[idx2] , i+1<resX-1 ? values[ INDEX(i+2,j) ] : values[idx2] , isoValue );
				break;
			case INTERPOLATE_CATMULL_ROM:
				iso = _CatmullRomInterpolant( i>0 ? values[ INDEX(i-1,j) ] : values[idx1] , values[idx1] , values[idx2] , i+1<resX-1 ? values[ INDEX(i+2,j) ] : values[idx2] , isoValue );
				break;
			default:
				ERROR_OUT( "Unrecognized interpolation type: " , interpolationType );
			}
			Point3D< Real > p = Point3D< Real >( (Real)i + iso , (Real)j , (Real)z );
			long long key = i + j*(resX);
#pragma omp critical
			{
				xIsoVertexMap[key] = (int)vertices.size();
				vertices.push_back( _Vertex( p , 0 , i , j , z ) );
			}
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<resX ; i++ ) for( int j=0 ; j<resY-1 ; j++ )
	{
		int idx1 = INDEX( i , j ) , idx2 = INDEX( i , j+1 );
		if( flags[idx1]!=flags[idx2] )
		{
			Real iso;
			switch( interpolationType )
			{
			case INTERPOLATE_LINEAR:
				iso = _LinearInterpolant( values[idx1] , values[idx2] , isoValue );
				break;
			case INTERPOLATE_QUADRATIC:
				iso = _QuadraticInterpolant( j>0 ? values[ INDEX(i,j-1) ] : values[idx1] , values[idx1] , values[idx2] , j+1<resY-1 ? values[ INDEX(i,j+2) ] : values[idx2] , isoValue );
				break;
			case INTERPOLATE_CUBIC:
				iso = _CubicInterpolant( j>0 ? values[ INDEX(i,j-1) ] : values[idx1] , values[idx1] , values[idx2] , j+1<resY-1 ? values[ INDEX(i,j+2) ] : values[idx2] , isoValue );
				break;
			case INTERPOLATE_CATMULL_ROM:
				iso = _CatmullRomInterpolant( j>0 ? values[ INDEX(i,j-1) ] : values[idx1] , values[idx1] , values[idx2] , j+1<resY-1 ? values[ INDEX(i,j+2) ] : values[idx2] , isoValue );
				break;
			default:
				ERROR_OUT( "Unrecognized interpolation type: " , interpolationType );
			}
			Point3D< Real > p = Point3D< Real >( (Real)i , (Real)j + iso , (Real)z );
			long long key = i + j*(resX);
#pragma omp critical
			{
				yIsoVertexMap[key] = (int)vertices.size();
				vertices.push_back( _Vertex( p , 1 , i , j , z ) );
			}
		}
	}
#undef INDEX
}

template< typename Real >
void IsoSurface3D< Real >::_SetPolygons( int resX , int resY , int z , ConstPointer( Real ) values1 , ConstPointer( Real ) values2 , Real isoValue , bool fullCaseTable , const std::unordered_map< long long , int >& xIsoVertexMap1 , const std::unordered_map< long long , int >& xIsoVertexMap2 , const std::unordered_map< long long , int >& yIsoVertexMap1 , const std::unordered_map< long long , int >& yIsoVertexMap2 , const std::unordered_map< long long , int >& zIsoVertexMap , const std::vector< _Vertex >& vertices , std::vector< std::vector< int > >& polygons )
{
#define INDEX( x , y ) ( x + (y)*resX )
#pragma omp parallel for
	for( int i=0 ; i<resX-1 ; i++ ) for( int j=0 ; j<resY-1 ; j++ )
	{
		Real _values[Cube::CORNERS];
		for( int cx=0 ; cx<2 ; cx++ ) for( int cy=0 ; cy<2 ; cy++ )
		{
			_values[ Cube::CornerIndex(cx,cy,0) ] = values1[ INDEX(i+cx,j+cy) ];
			_values[ Cube::CornerIndex(cx,cy,1) ] = values2[ INDEX(i+cx,j+cy) ];
		}
		int mcIndex = fullCaseTable ? MarchingCubes::GetFullIndex( _values , isoValue ) : MarchingCubes::GetIndex( _values , isoValue );
		const std::vector< std::vector< int > >& isoPolygons = MarchingCubes::caseTable( mcIndex , fullCaseTable );
		for( int p=0 ; p<isoPolygons.size() ; p++ )
		{
			const std::vector< int >& isoPolygon = isoPolygons[p];
			std::vector< int > polygon( isoPolygon.size() );
			for( int v=0 ; v<isoPolygon.size() ; v++ )
			{
				int orientation , i1 , i2;
				Cube::FactorEdgeIndex( isoPolygon[v] , orientation , i1 , i2 );
				long long key;
				std::unordered_map< long long , int >::const_iterator iter;
				bool success;
				switch( orientation )
				{
				case 0:
					key = (i   ) + (j+i1)*resX;
					if( i2==0 ){ iter = xIsoVertexMap1.find( key ) ; success = iter!=xIsoVertexMap1.end(); }
					else       { iter = xIsoVertexMap2.find( key ) ; success = iter!=xIsoVertexMap2.end(); }
					break;
				case 1:
					key = (i+i1) + (j   )*resX;
					if( i2==0 ){ iter = yIsoVertexMap1.find( key ) ; success = iter!=yIsoVertexMap1.end(); }
					else       { iter = yIsoVertexMap2.find( key ) ; success = iter!=yIsoVertexMap2.end(); }
					break;
				case 2:
					key = (i+i1) + (j+i2)*resX;
					iter = zIsoVertexMap.find( key ) ; success = iter!=zIsoVertexMap.end();
					break;
				}

				if( !success )
				{
					fprintf( stderr , "[ERROR] Couldn't find iso-vertex in map:\n" );
					printf( "\t%d: " , orientation );
					switch( orientation )
					{
					case 0: printf( "%d %d %d\n" , i , j+i1 , z+i2 ) ; break;
					case 1: printf( "%d %d %d\n" , i+i1 , j , z+i2 ) ; break;
					case 2: printf( "%d %d %d\n" , i+i1 , j+i2 , z ) ; break;
					}
					exit( 0 );
				}
				polygon[v] = iter->second;
			}
#pragma omp critical
			polygons.push_back( polygon );
		}
	}
#undef INDEX
}
