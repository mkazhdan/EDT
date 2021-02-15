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
#include "Misha/SquaredEDT.h"
#include "Misha/RegularGrid.h"
#include "Misha/Image.h"

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineReadable Verbose( "verbose" );
Misha::CmdLineParameter< int > Threads( "threads" , omp_get_num_procs() );
Misha::CmdLineParameterArray< unsigned int , 3 > ID( "id" );

Misha::CmdLineReadable* params[] =
{
	&In , &Out , &Verbose , &Threads , &ID , NULL
};

void ShowUsage( char* ex )
{
	std::cout << "Usage: " << std::string(ex) << std::endl;
	std::cout << "\t --" << In.name  << " <input image list>" << std::endl;
	std::cout << "\t --" << ID.name  << " <color id>" << std::endl;
	std::cout << "\t[--" << Out.name << " <output signed EDT>]" << std::endl;
#ifdef _OPENMP
	std::cout << "\t[--" << Threads.name << " <num threads>=" << Threads.value << "]" << std::endl;
#endif // _OPENMP
	std::cout << "\t[--" << Verbose.name << "]" << std::endl;
}

template< typename InputType >
RegularGrid< float , 3 > GetSignedEDT( const RegularGrid< InputType , 3 > &grid , InputType id , bool verbose )
{
	Miscellany::Timer timer;

	RegularGrid< unsigned char , 3 > mask;
	mask.resize( grid.res() );

	bool foundID = false;
#pragma omp parallel for
	for( int i=0 ; i<(int)grid.res(0) ; i++ ) for( int j=0 ; j<(int)grid.res(1) ; j++ ) for( int k=0 ; k<(int)grid.res(2) ; k++ )
	{
		bool boundary = false;
		if( grid(i,j,k)==id )
		{
			for( int ii=-1 ; ii<=1 ; ii++ ) for( int jj=-1 ; jj<=1 ; jj++ ) for( int kk=-1 ; kk<=1 ; kk++ )
				if( !ii && !jj && !kk ) continue;
				else if( i+ii<0 || i+ii>=(int)grid.res(0) || j+jj<0 || j+jj>=(int)grid.res(1) || k+kk<0 || k+kk>=(int)grid.res(2) ) continue;
				else if( grid( i+ii , j+jj , k+kk )!=id ) boundary = true;
			foundID = true;
		}
		if( boundary ) mask(i,j,k) = 1;
		else           mask(i,j,k) = 0;
	}
	if( !foundID ) ERROR_OUT( "Could not find voxel with target id" );
	std::cout << "Got mask: " << timer.elapsed() << std::endl;

	RegularGrid< unsigned int , 3 > EDT = Misha::SquaredEDT< double , 3 >::Saito( mask , verbose );
	RegularGrid< float , 3 > signedEDT;
	signedEDT.resize( grid.res() );

#pragma omp parallel for
	for( int i=0 ; i<(int)grid.res(0) ; i++ ) for( int j=0 ; j<(int)grid.res(1) ; j++ ) for( int k=0 ; k<(int)grid.res(2) ; k++ )
		if( grid(i,j,k)!=id ) signedEDT(i,j,k) =  (float)sqrt( EDT(i,j,k) );
		else                  signedEDT(i,j,k) = -(float)sqrt( EDT(i,j,k) );

	return signedEDT;
}

int main( int argc , char* argv[] )
{
	double t = Miscellany::Time();

	Misha::CmdLineParse( argc-1 , &argv[1] , params );
	if( !In.set || !ID.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	omp_set_num_threads( Threads.value );

	RegularGrid< unsigned int , 3 > grid;

	Miscellany::Timer timer;
	{
		auto RGBToInt = []( const unsigned char *rgb ){ return (unsigned int)rgb[0] <<16 | (unsigned int)rgb[1]<<8 | (unsigned int)rgb[2]<<0; };

		std::vector< std::string > images = Misha::ReadLines( In.value );

		for( unsigned int k=0 ; k<(unsigned int)images.size() ; k++ )
		{
			unsigned int width , height;
			unsigned char *rgb = ImageReader::ReadColor( images[k].c_str() , width , height );
			if( !k )
			{
				grid.resize( width , height , (unsigned int)images.size() );
				for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int j=0 ; j<height ; j++ )
					grid(i,j,k) = RGBToInt( rgb + 3 *( j*width+i ) );
			}
			else
			{
				if( width!=grid.res(0) ||  height!=grid.res(1) ) ERROR_OUT( "Image[" , k , "] resolution doesn't match: ( " , width , " , " , height , " ) != ( " , grid.res(0) , " , " , grid.res(1) , " )" );
				for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int j=0 ; j<height ; j++ )
					grid(i,j,k) = RGBToInt( rgb + 3 *( j*width+i ) );
			}
			delete[] rgb;
		}
	}
	std::cout << "Read grid: " << timer.elapsed() << std::endl;
	
	RegularGrid< float , 3 > signedEDT = GetSignedEDT< unsigned int >( grid , ID.values[0]<<16 | ID.values[1]<<8 | ID.values[2] , Verbose.set );

	if( Out.set ) signedEDT.template write< float >( Out.value , XForm< float , 4 >::Identity() );
	return EXIT_SUCCESS;
}
