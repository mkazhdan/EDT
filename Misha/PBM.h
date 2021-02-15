#ifndef PBM_INCLUDED
#define PBM_INCLUDED


struct PBMInfo
{
	unsigned char *data;
	FILE *fp;
	bool binary;
	unsigned int width , lineLength;
};

struct PBMReader : public ImageReader
{
	PBMReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	~PBMReader( void );
	unsigned int nextRow( unsigned char* row );
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
protected:
	unsigned int _currentRow;
	PBMInfo _info;
};

struct PBMWriter : public ImageWriter
{
	PBMWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
	~PBMWriter( void );
	unsigned int nextRow( const unsigned char* row );
	unsigned int nextRows( const unsigned char* row , unsigned int rowNum );
protected:
	PBMInfo _info;
	unsigned int _currentRow;
};

#include "PBM.inl"
#endif // PBM_INCLUDED
