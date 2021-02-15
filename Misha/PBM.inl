#include <stdio.h>
#include <stdlib.h>

inline bool PBMReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	PBMInfo _info;
	_info.fp = fopen( fileName , "rb" );
	if( !_info.fp ) return false;

	static const unsigned int COMMENT_SIZE = 4096;
	char comment[COMMENT_SIZE];
	char temp[512];
	bool haveMagicNumber=false , haveWidth=false , haveHeight=false;
	while( !haveMagicNumber || !haveWidth || !haveHeight )
	{
		if( fscanf( _info.fp , " %s " , temp )!=1 ){ fclose(_info.fp) ; return false; }
		else
		{
			if( temp[0]=='#' ) fgets( comment , COMMENT_SIZE , _info.fp );
			else if( !haveMagicNumber )
			{
				if( temp[0]!='P' ){ fclose(_info.fp) ;  return false; }
				else if( temp[1]=='1' ) _info.binary = false;
				else if( temp[1]=='4' ) _info.binary = true;
				else{ fclose(_info.fp) ; return false; }
				haveMagicNumber = true;
			}
			if( !haveWidth )
			{
				int w;
				if( fscanf( _info.fp , " %d " , &w )!=1 ){ fclose(_info.fp) ; return false; }
				width = w;
				haveWidth = true;
			}
			if( !haveHeight )
			{
				int h;
				if( fscanf( _info.fp , " %d " , &h )!=1 ){ fclose(_info.fp) ; return false; }
				height = h;
				haveHeight = true;
			}
		}
	}
	fclose( _info.fp );
	channels = 1;
	return true;
}

inline PBMReader::PBMReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	_currentRow = 0;
	_info.data = NULL;
	_info.fp = fopen( fileName , "rb" );
	if( !_info.fp ) ERROR_OUT( "PBMInitRead: Failed to open: " , fileName );

	static const unsigned int COMMENT_SIZE = 4096;
	char comment[COMMENT_SIZE];
	char temp[512];
	bool haveMagicNumber=false , haveWidth=false , haveHeight=false;
	while( !haveMagicNumber || !haveWidth || !haveHeight )
	{

		if( fscanf( _info.fp , " %s " , temp )!=1 ) ERROR_OUT( "Failed to read next string" );
		else
		{
			if( temp[0]=='#' ) fgets( comment , COMMENT_SIZE , _info.fp );
			else if( !haveMagicNumber )
			{
				if( temp[0]!='P' ) ERROR_OUT( "Failed to read magic number: " , temp );
				else if( temp[1]=='1' ) _info.binary = false;
				else if( temp[1]=='4' ) _info.binary = true;
				else ERROR_OUT( "Failed to read magic number: " , temp );
				haveMagicNumber = true;
			}
			if( !haveWidth )
			{
				int w;
				if( fscanf( _info.fp , " %d " , &w )!=1 ) ERROR_OUT( "Failed to read width" );
				width = w;
				haveWidth = true;
			}
			if( !haveHeight )
			{
				int h;
				if( fscanf( _info.fp , " %d " , &h )!=1 ) ERROR_OUT( "Failed to read height" );
				height = h;
				haveHeight = true;
			}
		}
	}
	channels = 1;

	_info.width = width;
	if( _info.binary )
	{
		_info.lineLength = (width+7)/8;
		_info.data = new unsigned char[ _info.lineLength ];
	}
	else _info.data = NULL;
}

inline PBMReader::~PBMReader( void ){ fclose( _info.fp ) ; delete[] _info.data; }

inline unsigned int PBMReader::nextRow( unsigned char* row )
{
	if( _info.binary )
	{
		fread( _info.data , sizeof(unsigned char) , _info.lineLength , _info.fp );
		for( unsigned int i=0 ; i<_info.width ; i++ )
		{
			unsigned int _i = i/8;
			if( _info.data[i/8] & ( 1<<( 7 - (i%8) ) ) ) row[i] = 255;
			else row[i] = 0;
		}
	}
	else for( unsigned int i=0 ; i<_info.width ; i++ ) if( fscanf( _info.fp , " %c" , row+i )!=1 ) ERROR_OUT( "Failed to read " , i , "-th character from line" );
	return ++_currentRow;
}

inline PBMWriter::PBMWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int )
{
	if( channels!=1 ) ERROR_OUT( "Only single-channel output supported: " , channels );
	_currentRow = 0;

	_info.fp = fopen( fileName , "wb" );
	if( !_info.fp ) ERROR_OUT( "Failed to open: " , fileName );

	_info.width = width;
	_info.lineLength = (width+7)/8;
	_info.data = new unsigned char[ _info.lineLength ];
	_info.binary = true;

	fprintf( _info.fp , "P4\n" );
	fprintf( _info.fp , "%d %d\n" , width , height );
}
inline PBMWriter::~PBMWriter( void ){ fclose(_info.fp) ; delete[] _info.data; }

inline unsigned int PBMWriter::nextRow( const unsigned char* row )
{
	if( _info.binary )
	{
		memset( _info.data , 0 , sizeof(unsigned char) * _info.lineLength );
		for( unsigned int i=0 ; i<_info.width ; i++ ) if( row[i] ) _info.data[i/8] |= 1<<( 7 - (i%8) );
		fwrite( _info.data , sizeof(unsigned char) , _info.lineLength , _info.fp );
	}
	else
	{
		for( unsigned int i=0 ; i<_info.width ; i++ ) fprintf( _info.fp , " %d" , row[i] ? 1 : 0 );
		fprintf( _info.fp , "\n" );
	}
	return ++_currentRow;
}
inline unsigned int PBMWriter::nextRows( const unsigned char* rows , unsigned int rowNum )
{
	for( unsigned int i=0 ; i<rowNum ; i++ ) nextRow( rows + _info.width*i );
	return _currentRow;
}
