/*
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

///////////////////
// Polynomial 0D //
///////////////////
template< unsigned int Degree >
Polynomial< 0 , Degree >::Polynomial( void ){ memset( _coefficients , 0 , sizeof( _coefficients ) ); }

template< unsigned int Degree >
Polynomial< 0 , Degree >::Polynomial( double c ) : Polynomial() { _coefficients[0] = c; }

template< unsigned int Degree >
unsigned int Polynomial< 0 , Degree >::_setCoefficients( const double *coefficients , unsigned int maxDegree )
{
	_coefficients[0] = coefficients[0];
	return 1;
}
template< unsigned int Degree >
unsigned int Polynomial< 0 , Degree >::_getCoefficients( double *coefficients , unsigned int maxDegree ) const
{
	coefficients[0] = _coefficients[0];
	return 1;
}
template< unsigned int Degree >
Point< double , Polynomial< 0 , Degree >::NumCoefficients > Polynomial< 0 , Degree >::coefficients( void ) const
{
	return Point< double , 1 >( _coefficients[0] );
}
template< unsigned int Degree >
Polynomial< 0 , Degree >::Polynomial( Point< double , NumCoefficients > coefficients ){ _setCoefficients( coefficients , 0 ); }

template<unsigned int Degree >
void Polynomial< 0 , Degree >::SetDegrees( unsigned int coefficientIndex , unsigned int degrees[/*0*/] ){}

template< unsigned int Degree >
template< unsigned int _Degree >
Polynomial< 0 , Degree >::Polynomial( const Polynomial< 0 , _Degree > &p )
{
	_coefficients[0] = p._coefficients[0];
}

template< unsigned int Degree >
template< unsigned int _Degree >
Polynomial< 0 , Degree > &Polynomial< 0 , Degree >::operator= ( const Polynomial< 0 , _Degree > &p )
{
	_coefficients[0] = p._coefficients[0];
	return *this;
}

template< unsigned int Degree >
const double &Polynomial< 0 , Degree >::_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const
{
	return _coefficients[0];
}

template< unsigned int Degree >
double &Polynomial< 0 , Degree >::_coefficient( const unsigned int indices[] , unsigned int maxDegree )
{
	return _coefficients[0];
}

template< unsigned int Degree >
double Polynomial< 0 , Degree >::_evaluate( const double coordinates[] , unsigned int maxDegree ) const { return _coefficients[0]; }

template< unsigned int Degree >
template< unsigned int _Dim >
Polynomial< _Dim , Degree > Polynomial< 0 , Degree >::_pullBack( const Matrix< double , _Dim+1 , 0 > &A , unsigned int maxDegree ) const
{
	return Polynomial< _Dim , Degree >( _coefficients[0] );
}

template< unsigned int Degree >
bool Polynomial< 0 , Degree >::_isZero( unsigned int maxDegree ) const { return _coefficients[0]==0; }

template< unsigned int Degree >
bool Polynomial< 0 , Degree >::_isConstant( unsigned int maxDegree ) const { return true; }

template< unsigned int Degree >
const double &Polynomial< 0 , Degree >::coefficient( void ) const { return _coefficients[0]; }

template< unsigned int Degree >
double &Polynomial< 0 , Degree >::coefficient( void ) { return _coefficients[0]; }

template< unsigned int Degree >
double &Polynomial< 0 , Degree >::coefficient( const unsigned int *d ) { return _coefficients[0]; }

template< unsigned int Degree >
const double &Polynomial< 0 , Degree >::coefficient( const unsigned int *d ) const { return _coefficients[0]; }

template< unsigned int Degree >
double &Polynomial< 0 , Degree >::coefficient( unsigned int *d ) { return _coefficients[0]; }

template< unsigned int Degree >
const double &Polynomial< 0 , Degree >::coefficient( unsigned int *d ) const { return _coefficients[0]; }

template< unsigned int Degree >
double Polynomial< 0 , Degree >::operator()( void ) const { return _coefficients[0]; }

template< unsigned int Degree >
double Polynomial< 0 , Degree >::operator()( Point< double , 0 > p ) const { return _coefficients[0]; }

template< unsigned int Degree >
Polynomial< 0 , (Degree>1) ? Degree-1 : 0 > Polynomial< 0 , Degree >::d( unsigned int ) const
{
	return Polynomial< 0 , (Degree>1) ? Degree-1 : 0 >(0);
}

template< unsigned int Degree >
template< unsigned int _Dim >
Polynomial< _Dim-1 , Degree > Polynomial< 0 , Degree >::operator()( const Matrix< double , _Dim , 0 > &A ) const
{
	static_assert( _Dim!=0 , "Dimension cannot be negative" );
	return Polynomial< _Dim-1 , Degree >( _coefficients[0] );
}

template< unsigned int Degree >
double Polynomial< 0 , Degree >::integrateUnitCube( void ) const { return _coefficients[0]; }

template< unsigned int Degree >
double Polynomial< 0 , Degree >::integrateUnitRightSimplex( void ) const { return _coefficients[0]; }

template< unsigned int Degree >
void Polynomial< 0 , Degree >::Scale( double s ){ _coefficients[0] *= s; }

template< unsigned int Degree >
void Polynomial< 0 , Degree >::Add( const Polynomial< 0 , Degree > &p ){ _coefficients[0] += p._coefficients[0]; }

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial< 0 , Degree1 + Degree2 > operator * ( const Polynomial< 0 , Degree1 > &p1 , const Polynomial< 0 , Degree2 > &p2 )
{
	return Polynomial< 0 , Degree1 + Degree2 >( p1._coefficients[0] * p2._coefficients[0] );
}

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial< 0 , Max< Degree1 , Degree2 >::Value > operator + ( const Polynomial< 0 , Degree1 > &p1 , const Polynomial< 0 , Degree2 > &p2 )
{
	return Polynomial< 0 , Max< Degree1 , Degree2 >::Value >( p1._coefficients[0] + p2._coefficients[0] );
}

template< unsigned int Degree1 , unsigned int Degree2 >
Polynomial< 0 , Max< Degree1 , Degree2 >::Value > operator - ( const Polynomial< 0 , Degree1 > &p1 , const Polynomial< 0 , Degree2 > &p2 ){ return p1 + (-p2); }

template< unsigned int Degree >
bool Polynomial< 0 , Degree >::_print( std::ostream &ostream , const std::string varNames[] , bool first ) const
{
	if( _coefficients[0] )
	{
		if( first ) ostream << _coefficients[0];
		else
		{
			if( _coefficients[0]<0 ) ostream << " - " << -_coefficients[0];
			else                     ostream << " + " <<  _coefficients[0];
		}
		first = false;
	}
	return first;
}

template< unsigned int Degree >
bool Polynomial< 0 , Degree >::_print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const
{
	if( _coefficients[0] )
	{
		if( first )
		{
			if     ( _coefficients[0]== 1 ) ostream << suffix;
			else if( _coefficients[0]==-1 ) ostream << "-" << suffix; 
			else                            ostream << _coefficients[0] << "*" << suffix;
		}
		else
		{
			if( _coefficients[0]<0 )
			{
				if( _coefficients[0]==-1 ) ostream << " - " << suffix;
				else                       ostream << " - " << -_coefficients[0] << "*" << suffix;
			}
			else
			{
				if( _coefficients[0]==1 ) ostream << " + " << suffix;
				else                      ostream << " + " <<  _coefficients[0] << "*" << suffix;
			}
		}
		first = false;
	}
	return first;
}

template< unsigned int Degree >
SquareMatrix< double , Polynomial< 0 , Degree >::NumCoefficients > Polynomial< 0 , Degree >::EvaluationMatrix( const Point< double , 0 > positions[NumCoefficients] )
{
	SquareMatrix< double , NumCoefficients > E;
	E(0,0) = 1;
	return E;
}

///////////////////////////////
// Polynomial Dim-dimensions //
///////////////////////////////
template< unsigned int Dim , unsigned int Degree > Polynomial< Dim , Degree >::Polynomial( void ){}

template< unsigned int Dim , unsigned int Degree > Polynomial< Dim , Degree >::Polynomial( double c )
{
	_polynomials[0] = Polynomial< Dim-1 , Degree >( c );
}

template< unsigned int Dim , unsigned int Degree >
unsigned int Polynomial< Dim , Degree >::_getCoefficients( double *coefficients , unsigned int maxDegree ) const
{
	unsigned int offset = 0;
	for( unsigned int d=0 ; d<=maxDegree ; d++ ) offset += _polynomials[d]._getCoefficients( coefficients + offset , maxDegree - d );
	return offset;
}
template< unsigned int Dim , unsigned int Degree >
unsigned int Polynomial< Dim , Degree >::_setCoefficients( const double *coefficients , unsigned int maxDegree )
{
	unsigned int offset = 0;
	for( unsigned int d=0 ; d<=maxDegree ; d++ ) offset += _polynomials[d]._setCoefficients( coefficients + offset , maxDegree - d );
	return offset;
}
template< unsigned int Dim , unsigned int Degree >
Polynomial< Dim , Degree >::Polynomial( Point< double , NumCoefficients > coefficients ){ _setCoefficients( &coefficients[0] , Degree ); }

template< unsigned int Dim , unsigned int Degree >
Point< double , Polynomial< Dim , Degree >::NumCoefficients > Polynomial< Dim , Degree >::coefficients( void ) const
{
	Point< double , NumCoefficients > c;
	_getCoefficients( &c[0] , Degree );
	return c;
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int _D >
typename std::enable_if< (_D<Degree ) >::type Polynomial< Dim , Degree >::_SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] )
{
	if( coefficientIndex < Polynomial< Dim-1 , Degree-_D >::NumCoefficients )
	{
		degrees[0] = _D;
		Polynomial< Dim-1 , Degree-_D >::SetDegrees( coefficientIndex , degrees+1 );
	}
	else _SetDegrees< _D+1 >( coefficientIndex - Polynomial< Dim-1 , Degree-_D >::NumCoefficients , degrees );
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int _D >
typename std::enable_if< (_D==Degree ) >::type Polynomial< Dim , Degree >::_SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] )
{
	if( coefficientIndex < Polynomial< Dim-1 , 0 >::NumCoefficients )
	{
		degrees[0] = _D;
		Polynomial< Dim-1 , Degree-_D >::SetDegrees( coefficientIndex , degrees+1 );
	}
	else ERROR_OUT( "coefficient index too big" );
}

template< unsigned int Dim , unsigned int Degree >
void Polynomial< Dim , Degree >::SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] )
{
	return _SetDegrees< 0 >( coefficientIndex , degrees );
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int _Degree >
Polynomial< Dim , Degree >::Polynomial( const Polynomial< Dim , _Degree > &p )
{
	for( int d=0 ; d<=Degree && d<=_Degree ; d++ ) _polynomials[d] = p._polynomials[d];
	for( int d=_Degree+1 ; d<=Degree ; d++ ) _polynomials[d] = Polynomial< Dim-1 , Degree >();
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int _Degree >
Polynomial< Dim , Degree > &Polynomial< Dim , Degree >::operator= ( const Polynomial< Dim , _Degree > &p )
{
	for( int d=0 ; d<=Degree && d<=_Degree ; d++ ) _polynomials[d] = p._polynomials[d];
	for( int d=_Degree+1 ; d<=Degree ; d++ ) _polynomials[d] = Polynomial< Dim-1 , Degree >();
	return *this;
}

template< unsigned int Dim , unsigned int Degree >
const double &Polynomial< Dim , Degree >::_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const
{
	if( indices[0]>maxDegree ) ERROR_OUT( "degree out of bounds: %d > %d\n" , indices[0] , maxDegree );
	return _polynomials[ indices[0] ]._coefficient( indices+1 , maxDegree-indices[0] );
}

template< unsigned int Dim , unsigned int Degree >
double& Polynomial< Dim , Degree >::_coefficient( const unsigned int indices[] , unsigned int maxDegree )
{
	if( indices[0]>maxDegree ) ERROR_OUT( "degree out of bounds: %d > %d\n" , indices[0] , maxDegree );
	return _polynomials[ indices[0] ]._coefficient( indices+1 , maxDegree-indices[0] );
}

template< unsigned int Dim , unsigned int Degree >
double Polynomial< Dim , Degree >::_evaluate( const double coordinates[] , unsigned int maxDegree ) const
{
	double sum = 0 , tmp = 1;
	for( unsigned int d=0 ; d<=maxDegree ; d++ )
	{
		sum += _polynomials[d]._evaluate( coordinates+1 , maxDegree-d ) * tmp;
		tmp *= coordinates[0];
	}
	return sum;
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int _Dim >
Polynomial< _Dim , Degree > Polynomial< Dim , Degree >::_pullBack( const Matrix< double , _Dim+1 , Dim > &A , unsigned int maxDegree ) const
{
	unsigned int indices[ _Dim>0 ? _Dim : 1 ]; 
	Polynomial< _Dim , 1 > _p;
	Matrix< double , _Dim+1 , Dim-1 > _A;

	for( int i=0 ; i<_Dim+1 ; i++ ) for( int j=1 ; j<Dim ; j++ ) _A(i,j-1) = A(i,j) ;

	for( int i=0 ; i<_Dim ; i++ ) indices[i] = 0;
	_p.coefficient( indices ) = A( _Dim , 0 );
	for( int i=0 ; i<_Dim ; i++ )
	{
		indices[i] = 1;
		_p.coefficient( indices ) = A( i , 0 );
		indices[i] = 0;
	}

	Polynomial< _Dim , Degree > p( 0. ) , __p( 1. );
	for( unsigned int d=0 ; d<=maxDegree ; d++ )
	{
		p += _polynomials[d]._pullBack< _Dim >( _A , maxDegree-d ) * __p;
		__p = __p * _p;
	}
	return p;
}

template< unsigned int Dim , unsigned int Degree >
bool Polynomial< Dim , Degree >::_isZero( unsigned int maxDegree ) const
{
	for( unsigned int d=0 ; d<=maxDegree ; d++ ) if( !_polynomials[d]._isZero( maxDegree-d ) ) return false;
	return true;
}

template< unsigned int Dim , unsigned int Degree >
bool Polynomial< Dim , Degree >::_isConstant( unsigned int maxDegree ) const
{
	if( !_polynomials[0]._isConstant( Degree ) ) return false;
	for( unsigned int d=1 ; d<=maxDegree ; d++ ) if( !_polynomials[d]._isZero( maxDegree-d ) ) return false;
	return true;
}

template< unsigned int Dim , unsigned int Degree >
template< typename ... UnsignedInts >
const double &Polynomial< Dim , Degree >::coefficient( UnsignedInts ... indices ) const
{
	static_assert( sizeof...(indices)==Dim  , "[ERROR] Polynomial< Dim , Degree >::coefficient: Invalid number of indices" );
	unsigned int _indices[] = { indices ... };
	return _coefficient( _indices , Degree );
}

template< unsigned int Dim , unsigned int Degree >
template< typename ... UnsignedInts >
double &Polynomial< Dim , Degree >::coefficient( UnsignedInts ... indices )
{
	static_assert( sizeof...(indices)==Dim , "[ERROR] Polynomial< Dim , Degree >::coefficient: Invalid number of indices" );
	unsigned int _indices[] = { (unsigned int)indices ... };
	return _coefficient( _indices , Degree );
}

template< unsigned int Dim , unsigned int Degree >
const double &Polynomial< Dim , Degree >::coefficient( const unsigned int indices[Dim] ) const
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree >
double &Polynomial< Dim , Degree >::coefficient( const unsigned int indices[Dim] )
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree >
const double &Polynomial< Dim , Degree >::coefficient( unsigned int indices[Dim] ) const
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree >
double &Polynomial< Dim , Degree >::coefficient( unsigned int indices[Dim] )
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree >
template< typename ... Doubles >
double Polynomial< Dim , Degree >::operator()( Doubles ... coordinates ) const
{
	static_assert( sizeof...(coordinates)==Dim , "[ERROR] Polynomial< Dim , Degree >::operator(): Invalid number of coordinates" );
	double _coordinates[] = { coordinates... };
	return _evaluate( _coordinates , Degree );
}

template< unsigned int Dim , unsigned int Degree >
double Polynomial< Dim , Degree >::operator()( Point< double , Dim > p ) const { return _evaluate( &p[0] , Degree ); }

/** This method returns the partial derivative with respect to the prescribed dimension.*/
template< unsigned int Dim , unsigned int Degree >
Polynomial< Dim , (Degree>1) ? Degree-1 : 0 > Polynomial< Dim , Degree >::d( int dim ) const
{
	Polynomial< Dim , (Degree>1) ? Degree-1 : 0 > derivative;
	if( dim==0 ) for( int d=0 ; d<Degree ; d++ ) derivative._polynomials[d] = _polynomials[d+1] * (d+1);
	else         for( int d=0 ; d<Degree ; d++ ) derivative._polynomials[d] = _polynomials[d].d( dim-1 );
	return derivative;
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int _Dim >
Polynomial< _Dim-1 , Degree > Polynomial< Dim , Degree >::operator()( const Matrix< double , _Dim , Dim > &A ) const
{
	static_assert( _Dim!=0 , "Dimension cannot be negative" );
	return _pullBack< _Dim-1 >( A , Degree );
}

template< unsigned int Dim , unsigned int Degree >
double Polynomial< Dim , Degree >::integrateUnitCube( void ) const
{
	// I_d = \int_0^1 ... \int_0^1 x_n^d * P_d(x_1,...,x_{n-1}) dx_n ... dx_1
	//     = 1/(d+1) * \int_0^1 ... \int_0^1 P_d(x_1,...,x_{n-1}) dx_{n-1} ... dx_1
	double integral = 0;
	for( int d=0 ; d<=Degree ; d++ ) integral += 1./(d+1) * _polynomials[d].integrateUnitCube();
	return integral;
}

template< unsigned int Dim , unsigned int Degree >
double Polynomial< Dim , Degree >::integrateUnitRightSimplex( void ) const
{
	// I_d = \int_0^1 \int_0^{1-x_1} ... \int_0^{1-x_1-x_2...-x_{n-1}} x_n^d * P_d(x_1,...,x_{n-1}) dx_n ... dx_1
	//     = 1/(d+1) * \int_0^1 ... \int_0^1 P_d(x_1,...,x_{n-1}) * (1 - x_1 - x_2 - ... - x_{n-1} )^{d+1} dx_{n-1} ... dx_1
	double integral = 0;
	Polynomial< Dim-1 , Degree+1 > p;
	Polynomial< Dim-1 , 1 > _p;
	{
		unsigned int indices[ Dim>1 ? Dim-1 : 1 ];
		for( int d=0 ; d<Dim-1 ; d++ ) indices[d] = 0;
		_p.coefficient( indices ) = 1;
		for( int d=0 ; d<Dim-1 ; d++ )
		{
			indices[d] = 1;
			_p.coefficient( indices ) = -1;
			indices[d] = 0;
		}
	}
	 p = _p;
	for( int d=0 ; d<=Degree ; d++ )
	{
		integral += ( _polynomials[d] * p ).integrateUnitRightSimplex() / ( d+1 );
		p = p * _p;
	}
	return integral;
}

template< unsigned int Dim , unsigned int Degree >
SquareMatrix< double , Polynomial< Dim , Degree >::NumCoefficients > Polynomial< Dim , Degree >::EvaluationMatrix( const Point< double , Dim > positions[NumCoefficients] )
{
	SquareMatrix< double , NumCoefficients > E;
	unsigned int degrees[ Dim ];
	for( unsigned int i=0 ; i<NumCoefficients ; i++ ) for( unsigned int j=0 ; j<NumCoefficients ; j++ )
	{
		SetDegrees( i , degrees );
		double value = 1;
		for( int d=0 ; d<Dim ; d++ ) value *= pow( positions[j][d] , degrees[d] );
		E(i,j) = value;
	}
	return E;
}

template< unsigned int Dim , unsigned int Degree >
bool Polynomial< Dim , Degree >::_print( std::ostream &ostream , const std::string varNames[] , bool first ) const
{
	first = _polynomials[0]._print( ostream , varNames , first );
	for( int d=1 ; d<=Degree ; d++ )
		if( d==1 ) first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] , first );
		else       first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] + "^" + std::to_string( d ) , first );
	return first;
}
template< unsigned int Dim , unsigned int Degree >
bool Polynomial< Dim , Degree >::_print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const
{
	first = _polynomials[0]._print( ostream , varNames , suffix , first );
	for( int d=1 ; d<=Degree ; d++ )
		if( d==1 ) first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] + "*" + suffix , first );
		else       first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] + "^" + std::to_string( d ) + "*" + suffix , first );
	return first;
}

template< unsigned int Dim , unsigned int Degree >
std::ostream &operator << ( std::ostream &stream , const Polynomial< Dim , Degree > &poly )
{
	std::string varNames[Dim];
	if     ( Dim==1 ) varNames[0] = "x";
	else if( Dim==2 ) varNames[0] = "y" , varNames[1] = "x";
	else if( Dim==3 ) varNames[0] = "z" , varNames[1] = "y" , varNames[2] = "x";
	else for( int i=0 ; i<Dim ; i++ ) varNames[i] = std::string( "x" ) + std::string( "_" ) + std::to_string( Dim-i-1 );
	if( poly._print( stream , varNames , true ) ) stream << "0";
	return stream;
}

template< unsigned int Degree >
std::ostream &operator << ( std::ostream &stream , const Polynomial< 0 , Degree > &poly )
{
	std::string varNames[1] = { "" };
	if( poly._print( stream , varNames , true ) ) stream << "0";
	return stream;
}

template< unsigned int Dim , unsigned int Degree >
void Polynomial< Dim , Degree >::Scale( double s )
{
	for( int d=0 ; d<=Degree ; d++ ) _polynomials[d] *= s;
}

template< unsigned int Dim , unsigned int Degree >
void Polynomial< Dim , Degree >::Add( const Polynomial< Dim , Degree > &p )
{
	for( int d=0 ; d<=Degree ; d++ ) _polynomials[d] += p._polynomials[d];
}

template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 >
Polynomial< Dim , Degree1 + Degree2 > operator * ( const Polynomial< Dim , Degree1 > &p1 , const Polynomial< Dim , Degree2 > &p2 )
{
	Polynomial< Dim , Degree1 + Degree2 > p;
	for( int d1=0 ; d1<=Degree1 ; d1++ ) for( int d2=0 ; d2<=Degree2 ; d2++ ) p._polynomials[ d1+d2 ] += p1._polynomials[d1] * p2._polynomials[d2];
	return p;
}

template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 >
Polynomial< Dim , Max< Degree1 , Degree2 >::Value > operator + ( const Polynomial< Dim , Degree1 > &p1 , const Polynomial< Dim , Degree2 > &p2 )
{
	Polynomial< Dim , Max< Degree1 , Degree2 >::Value > p;
	for( int d=0 ; d<=Degree1 ; d++ ) p._polynomials[d] += p1._polynomials[d];
	for( int d=0 ; d<=Degree2 ; d++ ) p._polynomials[d] += p2._polynomials[d];
	return p;
}

template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 >
Polynomial< Dim , Max< Degree1 , Degree2 >::Value > operator - ( const Polynomial< Dim , Degree1 > &p1 , const Polynomial< Dim , Degree2 > &p2 ){ return p1 + (-p2); }

template< unsigned int Degree >
unsigned int Roots( const Polynomial< 1 , Degree > &p , double *r )
{
	ERROR_OUT( "Root functionality not supported for polynomial of degree = %d" , Degree );
	return 0;
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 1 > &p , double *r )
{
	if( p.coefficient(1u)==0 ) return 0;
	else
	{
		r[0] = -p.coefficient(0u) / p.coefficient(1u);
		return 1;
	}
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 2 > &p , double *r )
{
	if( !p.coefficient(2u) ) return Roots( Polynomial< 1 , 1 >( p ) , r );
	double disc = p.coefficient(1u)*p.coefficient(1u) - 4. * p.coefficient(0u) * p.coefficient(2u);
	if( disc<0 ) return 0;
	else if( disc==0 )
	{
		r[0] = - p.coefficient(1u) / ( 2 * p.coefficient(2u) );
		return 1;
	}
	else
	{
		disc = sqrt(disc);
		r[0] = ( -p.coefficient(1u) - disc ) / (2 * p.coefficient(2u) );
		r[1] = ( -p.coefficient(1u) + disc ) / (2 * p.coefficient(2u) );
		return 2;
	}
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 3 > &p , double *r )
{
	if( !p.coefficient(3u) ) return Roots( Polynomial< 1 , 2 >( p ) , r );
	return Poly34::SolveP3( r , p.coefficient(2u)/p.coefficient(3u) , p.coefficient(1u)/p.coefficient(3u) , p.coefficient(0u)/p.coefficient(3u) );
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 4 > &p , double *r )
{
	if( !p.coefficient(4u) ) return Roots( Polynomial< 1 , 3 >( p ) , r );
	return Poly34::SolveP4( r , p.coefficient(3u)/p.coefficient(4u) , p.coefficient(2u)/p.coefficient(4u) , p.coefficient(1u)/p.coefficient(4u) , p.coefficient(0u)/p.coefficient(4u) );
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 5 > &p , double *r )
{
	if( !p.coefficient(5u) ) return Roots( Polynomial< 1 , 4 >( p ) , r );
	return Poly34::SolveP5( r , p.coefficient(4u)/p.coefficient(5u) , p.coefficient(3u)/p.coefficient(5u) , p.coefficient(2u)/p.coefficient(5u) , p.coefficient(1u)/p.coefficient(5u) , p.coefficient(0u)/p.coefficient(5u) );
}
