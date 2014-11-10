 //! Typedefs for the application 
/*!
 * @author diehlpk
 * @date 2012
 */

#ifndef TEMPLATE
#define TEMPLATE

#include <cstdlib>
#include <cassert>

template<typename T>
struct Tuple2
{
public:
	Tuple2()
		: x( 0 ), y( 0 )
	{
	}

	Tuple2( T _x, T _y )
		: x( _x ), y( _y )
	{
	}

	Tuple2<T>& operator =( const Tuple2<T>& rhs )
	{
		x = rhs.x;
		y = rhs.y;
		return *this;
	}

	bool operator ==( const Tuple2<T>& rhs ) const
	{
		return ( x == rhs.x && y == rhs.y );
	}

	bool operator !=( const Tuple2<T>& rhs ) const
	{
		return !( *this == rhs );
	}

	T& operator []( size_t index )
	{
		assert( index < 2 );
		return coord[ index ];
	}

	const T& operator []( size_t index ) const
	{
		assert( index < 2 );
		return coord[ index ];
	}

	Tuple2<T> operator +( const Tuple2<T>& rhs ) const
	{
		return Tuple2<T>( x + rhs.x, y + rhs.y );
	}

	Tuple2<T>& operator +=( const Tuple2<T>& rhs )
	{
		x += rhs.x;
		y += rhs.y;
		return *this;
	}

	Tuple2<T> operator -( const Tuple2<T>& rhs ) const
	{
		return Tuple2<T>( x - rhs.x, y - rhs.y );
	}

	Tuple2<T>& operator -=( const Tuple2<T>& rhs )
	{
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}

public:
	static const Tuple2<T> ZERO;
	static const Tuple2<T> ONE;

	union
	{
		T coord[ 2 ];
		struct
		{
			T x;
			T y;
		};
	};
};

template<typename T>
const Tuple2<T> Tuple2<T>::ZERO( 0, 0 );
template<typename T>
const Tuple2<T> Tuple2<T>::ONE( 1, 1 );

//! Template for a basic array data type.
/*!
 * \param T The datatype or class of the elements in the array.
 * \param size The size of the array.
 */
template < class T, size_t size > class Array
{
private:

  //! The elements of the array.
  T coord[size];
public:

  //! Constructor.
  Array (void)
  {
    for (size_t d = 0; d < size; ++d)
      coord[d] = 0;
  }

  //! Destructor.
  ~Array (void)
  {
  }

  //! Operator returns the element with the index index.
  /*!
   * \param index The index of the element.
   */
  T & operator[](size_t index)
  {
    assert (size > index);
    return coord[index];
  }

  //! Operator returns the element with the index index.
  /*!
   * \param index The index of the element.
   */
  const T & operator[] (size_t index) const
  {
    assert (size > index);
    return coord[index];
  }
};


#endif
