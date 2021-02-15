/* -*- C++ -*-
Copyright (c) 2019, Michael Kazhdan
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


#ifndef PLY_INCLUDED
#define PLY_INCLUDED

#define NEW_PLY
#define USE_PLY_FACTORY
#include <vector>
#include <string>
#include <functional>
#include "PlyFile.h"
#include "Geometry.h"
#include "Exceptions.h"

namespace PLY
{
	// Converts from C-type to PLY type
	template< class Real > int Type( void );

	// Converts from C-type to PLY name
	template< typename Integer > struct Traits{ static const std::string name; };

	// A structure representing a face
	template< typename Index >
	struct Face
	{
		unsigned int nr_vertices;
		Index *vertices;

		static PlyProperty Properties[];
	};


	int DefaultFileType( void );

	// PLY read functionality
	void ReadHeader( std::string fileName , const PlyProperty *properties , int propertyNum , bool *readFlags );

	void ReadHeader( std::string fileName , const PlyProperty *properties , int propertyNum , bool *readFlags , int &file_type );

#ifdef USE_PLY_FACTORY
	template< typename VertexFactory , typename Index >
	void Read( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< std::pair< Index , Index > > *edges , std::vector< std::vector< Index > > *polygons , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL );
#else // !USE_PLY_FACTORY
	template< typename Vertex , typename Index >
	void Read( std::string fileName , std::vector< Vertex > &vertices , std::vector< std::pair< Index , Index > > *edges , std::vector< std::vector< Index > > *polygons , const PlyProperty *vertexProperties , bool *vertexPropertiesFlag , int vertexPropertyNum , int &file_type , std::vector< std::string > *comments=NULL );
#endif // USE_PLY_FACTORY

#ifdef USE_PLY_FACTORY
	template< typename VertexFactory >
	void ReadVertices( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL );
#else // !USE_PLY_FACTORY
	template< typename Vertex >
	void ReadVertices( std::string fileName , std::vector< Vertex > &vertices , const PlyProperty *vertexProperties , bool *vertexPropertiesFlag , int vertexPropertyNum , int &file_type , std::vector< std::string > *comments=NULL );
#endif // USE_PLY_FACTORY

#ifdef USE_PLY_FACTORY
	template< typename VertexFactory , typename Index >
	void ReadTriangles( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL );
	template< typename VertexFactory , typename Real , unsigned int Dim , typename Index >
	void ReadTriangles( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , bool *vertexPropertiesFlag , int &file_type , std::function< Point< Real , Dim > ( typename VertexFactory::VertexType ) > VertexToPointFunctor , std::vector< std::string > *comments=NULL );

	template< typename VertexFactory , typename Index >
	void ReadPolygons( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< std::vector< Index > > &polygons ,  bool *readFlags , int &file_type , std::vector< std::string > *comments=NULL );

	template< typename VertexFactory , typename Index >
	void ReadTetrahedra( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< 3 , Index > > &tetrahedra , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL );

	template< class VertexFactory , typename Polygon >
	int ReadPolygons( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType >& vertices , std::vector< Polygon >& polygons , bool *vertexPropertiesFlag , PlyProperty* polygonProperties , bool* polygonPropertiesFlag , int polygonPropertyNum , int& file_type, std::vector< std::string > *comments=NULL );
#else // !USE_PLY_FACTORY
	template< typename Vertex , typename Real , unsigned int Dim , typename Index >
	void ReadTriangles( std::string fileName , const VertexFactory &vFactory , std::vector< Vertex > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , const PlyProperty *vertexProperties , bool *vertexPropertiesFlag , int vertexPropertyNum , int &file_type , std::function< Point< Real , Dim > (Vertex) > VertexToPointFunctor , std::vector< std::string > *comments=NULL );

	template< typename Vertex , typename Index >
	void ReadPolygons( std::string fileName , std::vector< Vertex > &vertices , std::vector< std::vector< Index > > &polygons , const PlyProperty *properties , bool *readFlags , int propertyNum , int &file_type , std::vector< std::string > *comments=NULL );

	template< typename Vertex , typename Index >
	void ReadTetrahedra( std::string fileName , std::vector< Vertex > &vertices , std::vector< SimplexIndex< 3 , Index > > &tetrahedra , const PlyProperty *vertexProperties , bool *vertexPropertiesFlag , int vertexPropertyNum , int &file_type , std::vector< std::string > *comments=NULL );

	template< class Vertex , typename Polygon >
	int ReadPolygons( std::string fileName , std::vector< Vertex >& vertices , std::vector< Polygon >& polygons , PlyProperty*  vertexProperties , bool*  vertexPropertiesFlag , int  vertexPropertyNum , PlyProperty* polygonProperties , bool* polygonPropertiesFlag , int polygonPropertyNum , int &file_type , std::vector< std::string > *comments );
#endif // USE_PLY_FACTORY

	// PLY write functionality
#ifdef USE_PLY_FACTORY
	template< typename VertexFactory >
	void WriteVertices( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const PlyProperty *vertexProperties , int vertexPropertyNum , int file_type , const std::vector< std::string > *comments=NULL );
#else // !USE_PLY_FACTORY
	template< typename Vertex >
	void WriteVertices( std::string fileName , const std::vector< Vertex > &vertices , const PlyProperty *vertexProperties , int vertexPropertyNum , int file_type , const std::vector< std::string > *comments=NULL );
#endif // USE_PLY_FACTORY

#ifdef USE_PLY_FACTORY
	template< typename VertexFactory , typename Index >
	void Write( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< std::pair< Index , Index > > *edges , const std::vector< std::vector< Index > > *polygons , int file_type , const std::vector< std::string > *comments=NULL );

	template< typename VertexFactory , typename Index >
	void WriteTriangles( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles , int file_type , const std::vector< std::string > *comments=NULL );

	template< typename VertexFactory , typename Index >
	void WritePolygons( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< std::vector< Index > > &polygons , int file_type , const std::vector< std::string > *comments=NULL );

	template< typename VertexFactory , typename Index >
	void WriteTetrahedra( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< SimplexIndex< 3 , Index > > &tetrahedra , int file_type , const std::vector< std::string > *comments=NULL );

	template< class VertexFactory , typename Polygon >
	void WritePolygons( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< Polygon > &polygons , PlyProperty *polygonProperties , int polygonPropertyNum , int file_type , const std::vector< std::string > *comments=NULL );
#else // !USE_PLY_FACTORY
	template< typename Vertex , typename Index >
	void Write( std::string fileName , const std::vector< Vertex > &vertices , const std::vector< std::pair< Index , Index > > *edges , const std::vector< std::vector< Index > > *polygons , const PlyProperty *vertexProperties , int vertexPropertyNum , int file_type , const std::vector< std::string > *comments=NULL );

	template< typename Vertex , typename Index >
	void WriteTriangles( std::string fileName , const std::vector< Vertex > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles , const PlyProperty *properties , int propertyNum , int file_type , const std::vector< std::string > *comments=NULL );

	template< typename Vertex , typename Index >
	void WritePolygons( std::string fileName , const std::vector< Vertex > &vertices , const std::vector< std::vector< Index > > &polygons , const PlyProperty *properties , int propertyNum , int file_type , const std::vector< std::string > *comments=NULL );

	template< typename Vertex , typename Index >
	void WriteTetrahedra( std::string fileName , const std::vector< Vertex > &vertices , const std::vector< SimplexIndex< 3 , Index > > &tetrahedra , const PlyProperty *properties , int propertyNum , int file_type , const std::vector< std::string > *comments=NULL );

	template< class Vertex , typename Polygon >
	void WritePolygons( std::string fileName , const std::vector< Vertex > &vertices , const std::vector< Polygon > &polygons , const PlyProperty *properties , int propertyNum , PlyProperty *polygonProperties , int polygonPropertyNum , int file_type , const std::vector< std::string > *comments=NULL );
#endif // USE_PLY_FACTORY

#ifdef USE_PLY_FACTORY
#else // !USE_PLY_FACTORY
	// Read/write specializations for PLY::VertexData< Real , VertexPosition< Real , Dim > ... >
	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void Read( std::string fileName , std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , std::vector< std::pair< Index , Index > > *edges , std::vector< std::vector< Index > > *polygons , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return Read< PlyVertex >( fileName , vertices , edges , polygons , PlyVertex::ReadProperties() , vertexPropertiesFlag , PlyVertex::ReadNum , file_type , []( PlyVertex v ){ return Point< double , Dim >( v.data<0>() ); } , comments );
	}
	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void Write( std::string fileName , const std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , const std::vector< std::pair< Index , Index > > *edges , const std::vector< std::vector< Index > > *polygons , int file_type , const std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return Write( fileName , vertices , edges , polygons , PlyVertex::WriteProperties() , PlyVertex::WriteNum , file_type , comments );
	}

	template< typename Real , unsigned int Dim , typename ... AdditionalVertexData >
	void ReadVertices( std::string fileName , std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return ReadVertices< PlyVertex >( fileName , vertices , PlyVertex::ReadProperties() , vertexPropertiesFlag , PlyVertex::ReadNum , file_type , []( PlyVertex v ){ return Point< double , Dim >( v.data<0>() ); } , comments );
	}
	template< typename Real , unsigned int Dim , typename ... AdditionalVertexData >
	void WriteVertices( std::string fileName , const std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , int file_type , const std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return WriteVertices( fileName , vertices , PlyVertex::WriteProperties() , PlyVertex::WriteNum , file_type , comments );
	}

	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void ReadTriangles( std::string fileName , std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return ReadTriangles< PlyVertex , double , Dim >( fileName , vertices , triangles , PlyVertex::ReadProperties() , vertexPropertiesFlag , PlyVertex::ReadNum , file_type , []( PlyVertex v ){ return Point< double , Dim >( v.data<0>() ); } , comments );
	}
	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void WriteTriangles( std::string fileName , const std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles , int file_type , const std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return WriteTriangles( fileName , vertices , triangles , PlyVertex::WriteProperties() , PlyVertex::WriteNum , file_type , comments );
	}

	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void ReadPolygons( std::string fileName , std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , std::vector< std::vector< Index > > &polygons , bool *readFlags , int &file_type , std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return ReadPolygons( fileName , vertices , polygons , PlyVertex::ReadProperties() , readFlags , PlyVertex::ReadNum , file_type , comments );
	}
	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void WritePolygons( std::string fileName , const std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , const std::vector< std::vector< Index > > &polygons , int file_type , const std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return WritePolygons( fileName , vertices , polygons , PlyVertex::WriteProperties() , PlyVertex::WriteNum , file_type , comments );
	}

	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void ReadTetrahedra( std::string fileName , std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , std::vector< SimplexIndex< 3 , Index > > &tetrahedra , bool *vertexPropertiesFlag , int &file_type , std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return ReadTetrahedra< PlyVertex >( fileName , vertices , tetrahedra , PlyVertex::ReadProperties() , vertexPropertiesFlag , PlyVertex::ReadNum , file_type , []( PlyVertex v ){ return Point< double , Dim >( v.data<0>() ); } , comments );
	}
	template< typename Real , unsigned int Dim , typename Index , typename ... AdditionalVertexData >
	void WriteTetrahedra( std::string fileName , const std::vector< VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > > &vertices , const std::vector< SimplexIndex< 3 , Index > > &tetrahedra , int file_type , const std::vector< std::string > *comments=NULL )
	{
		typedef VertexData< Real , VertexPosition< Real , Dim > , AdditionalVertexData ... > PlyVertex;
		return WriteTetrahedra( fileName , vertices , tetrahedra , PlyVertex::WriteProperties() , PlyVertex::WriteNum , file_type , comments );
	}
#endif // USE_PLY_FACTORY

}
#include "Ply.inl"

#endif // PLY_INCLUDED
