//*************************************************************************
// File:     utility.h
// Author:   Yan LIU
// Date:     01-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#ifndef _UTILITY_H_01_10_2013__
#define _UTILITY_H_01_10_2013__


#include "config.h"

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

#include <set>
#include <vector>
#include <list>
#include <stack>
#include <map>
#include <utility>
#include <set>
#include <bitset>
#include <deque>
#include <string>
#include <queue>
#include <algorithm>
#include <functional>



NS_BEGIN

/******************************************************************************
*                                 Debug                                       *
******************************************************************************/
struct Debug
{
	static int CellID;
};


/******************************************************************************
*                                 Random                                      *
******************************************************************************/

struct   Random 
{
    // Rand()    Generate a random number between 0 and 'choices' - 1.
    static unsigned int Rand(unsigned int choices);
    
    static void         SetSeed(unsigned int seed);
    static unsigned int GetSeed();
    static unsigned int randomseed;

	template<typename Iter>
	static void Shuffle(Iter b,Iter e);
};


/******************************************************************************
*                                 MiniTools                                   * 
******************************************************************************/

template <typename T>
inline const T& Minimum(const T &a, const T &b)
{
    return a<b ? a : b ;
}

template <typename T>
inline const T& Maximum(const T &a, const T &b)
{
    return a>b ? a : b ;
}

template <typename T>
inline  const T& Minimum3(const T &a, const T &b,const T& c)
{
    return Minimum(Minimum(a,b),c);
}

template <typename T>
inline  const T& Minimum4(const T &a, const T &b,const T& c,const T& d)
{
    return Minimum(Minimum(Minimum(a,b),c),d);
}

template <typename T>
inline  const T& Maximum3(const T &a, const T &b,const T& c)
{
    return Maximum(Maximum(a,b),c);
}

template <typename T>
inline  const T& Maximum4(const T &a, const T &b,const T& c,const T& d)
{
    return Maximum(Maximum(Maximum(a,b),c),d);
}

template <typename T>
inline void Swap(T& a, T& b) 
{
    T tmp = a;
    a = b;
    b = tmp;
}

template<typename Iter>
void Random::Shuffle(Iter b,Iter e)
{
	if(b==e) return;
    size_t c = 2;
    Iter i,j;
    for (i=b+1;i!=e;i++,c++)
    {
        j=b+Random::Rand(c);
        Swap(*i,*j);
    }
}

template<typename T>
void Sort2(T& a,T& b)
{
    if(b<a) Swap(a,b);
}

template <typename T>
void Sort3 (T& a, T& b, T& c)
/* Input/Output: a, b, c. */
/* Sorts (&a,&b,&c) with <= 3 comparisons and <= 3 swaps. */
{
    if (b < a)Swap (a, b);
    if (c < b)
    {
        Swap (b, c);
        if (b < a)Swap (a, b);
    }
}


template <typename T>
void Sort4 (T& a, T& b, T& c,T& d)
/* Input/Output: a, b, c, d. */
/* Sorts (&a,&b,&c,&d) with <= 6 = 3 + 3 comparisons
and <= 6 = 3 + 3 swaps (insertion sort). */
{
    /* step1: isort3 */
    if (b < a)
        Swap (a, b);
    if (c < b)
    {
        Swap (b, c);
        if (b<a)
            Swap (a, b);
    }
    /* step2: isort4p */
    if (d < c)
    {
        if (d < a)
        {
            Swap (c, d);
            Swap (b, c);
            Swap (a, b);
        }
        else if (d < b)
        {
            Swap (c, d);
            Swap (b, c);
        }
        else
            Swap (c, d);
    }
}

/******************************************************************************
*                                 Index                                       *
******************************************************************************/
template<typename T>
class  Index2
{
	template<typename Type>
	friend bool operator <  (const Index2<Type>& a, const Index2<Type>& b);
	template<typename Type>
	friend bool operator == (const Index2<Type>& a, const Index2<Type>& b);
	template<typename Type>
	friend bool operator != (const Index2<Type>& a, const Index2<Type>& b);
public:
	Index2(const T& k0 = T(), const T& k1 = T())
	{
		_key[0] = k0;
		_key[1] = k1;
	}

	const T& get(int i) const
	{
		ASSERT(i == 0 || i == 1);
		return _key[i];
	}
	void set(int i, const T& k)
	{
		ASSERT(i == 0 || i == 1);
		_key[i] = k;
	}
	void reset(const T& k0, const T& k1)
	{
		_key[0] = k0;
		_key[1] = k1;
	}
private:
	mutable T   _key[2];
};

template<typename T>
bool operator < (const Index2<T>& a, const Index2<T>& b)
{
	if (a._key[0] < b._key[0]) return true;
	if (b._key[0] < a._key[0]) return false;
	if (a._key[1] < b._key[1]) return true;
	if (b._key[1] < a._key[1]) return false;
	return false;
}

template<typename T>
bool operator == (const Index2<T>& a, const Index2<T>& b)
{
	return !(a<b) && !(b<a);
}

template<typename T>
bool operator != (const Index2<T>& a, const Index2<T>& b)
{
	return !(a == b);
}


template<typename T>
class  Index3
{
	template<typename Type>
	friend bool operator <  (const Index3<Type>& a, const Index3<Type>& b);
	template<typename Type>
	friend bool operator == (const Index3<Type>& a, const Index3<Type>& b);
	template<typename Type>
	friend bool operator != (const Index3<Type>& a, const Index3<Type>& b);
public:
	Index3(const T& k0 = T(), const T& k1 = T(), const T& k2 = T())
	{
		_key[0] = k0;
		_key[1] = k1;
		_key[2] = k2;
	}

	const T& get(int i) const
	{
		ASSERT(0 <= i && i <= 2);
		return _key[i];
	}
	void set(int i, const T& k)
	{
		ASSERT(0 <= i && i <= 2);
		_key[i] = k;
	}
	void reset(const T& k0, const T& k1, const T& k2)
	{
		_key[0] = k0;
		_key[1] = k1;
		_key[2] = k2;
	}

private:
	mutable T   _key[3];
};

template<typename T>
bool operator < (const Index3<T>& a, const Index3<T>& b)
{
	if (a._key[0] < b._key[0]) return true;
	if (b._key[0] < a._key[0]) return false;
	if (a._key[1] < b._key[1]) return true;
	if (b._key[1] < a._key[1]) return false;
	if (a._key[2] < b._key[2]) return true;
	if (b._key[2] < a._key[2]) return false;
	return false;
}

template<typename T>
bool operator == (const Index3<T>& a, const Index3<T>& b)
{
	return !(a<b) && !(b<a);
}

template<typename T>
bool operator != (const Index3<T>& a, const Index3<T>& b)
{
	return !(a == b);
}


/******************************************************************************
*                                 KeyIndex                                       *
******************************************************************************/
template<typename T>
class  KeyIndex2
{
    template<typename Type>
    friend bool operator <  (const KeyIndex2<Type>& a,const KeyIndex2<Type>& b);
    template<typename Type>
    friend bool operator == (const KeyIndex2<Type>& a,const KeyIndex2<Type>& b);
    template<typename Type>
    friend bool operator != (const KeyIndex2<Type>& a,const KeyIndex2<Type>& b);
public:
    KeyIndex2(const T& k0=T(),const T& k1=T())
    {
        _key[0] = k0;
        _key[1] = k1;
        sort();
    }

    const T& get(int i) const
    {
        ASSERT(i==0 || i==1);
        return _key[i];
    }
    void set(int i,const T& k)
    {
        ASSERT(i==0 || i==1);
        _key[i] = k;
        _sorted = false;
    }
    void reset(const T& k0,const T& k1)
    {
        _key[0] = k0;
        _key[1] = k1;
        _sorted = false;
    }
private:
    void sort() const
    {
        Sort2(_key[0],_key[1]);
        _sorted = true;
    }
private:
    mutable bool _sorted;
    mutable T   _key[2];
};

template<typename T>
bool operator < (const KeyIndex2<T>& a,const KeyIndex2<T>& b)
{
    if(a._sorted == false) a.sort();
    if(b._sorted == false) b.sort();
    if(a._key[0] < b._key[0]) return true;
    if(b._key[0] < a._key[0]) return false;
    if(a._key[1] < b._key[1]) return true;
    if(b._key[1] < a._key[1]) return false;
    return false;
}

template<typename T>
bool operator == (const KeyIndex2<T>& a,const KeyIndex2<T>& b)
{
    return !(a<b) && !(b<a);
}

template<typename T>
bool operator != (const KeyIndex2<T>& a,const KeyIndex2<T>& b)
{
    return !(a==b);
}


template<typename T>
class  KeyIndex3
{
    template<typename Type>
    friend bool operator <  (const KeyIndex3<Type>& a,const KeyIndex3<Type>& b);
    template<typename Type>
    friend bool operator == (const KeyIndex3<Type>& a,const KeyIndex3<Type>& b);
    template<typename Type>
    friend bool operator != (const KeyIndex3<Type>& a,const KeyIndex3<Type>& b);
public:
    KeyIndex3(const T& k0=T(),const T& k1=T(),const T& k2=T())
    {
        _key[0] = k0;
        _key[1] = k1;
        _key[2] = k2;
        sort();
    }

    const T& get(int i) const
    {
        ASSERT(0<=i && i<=2);
        return _key[i];
    }
    void set(int i,const T& k)
    {
        ASSERT(0<=i && i<=2);
        _key[i] = k;
        _sorted = false;
    }
    void reset(const T& k0,const T& k1,const T& k2)
    {
        _key[0] = k0;
        _key[1] = k1;
        _key[2] = k2;
        _sorted = false;
    }
private:
    void sort() const
    {
        Sort3(_key[0],_key[1],_key[2]);
        _sorted = true;
    }
private:
    mutable bool _sorted;
    mutable T   _key[3];
};

template<typename T>
bool operator < (const KeyIndex3<T>& a,const KeyIndex3<T>& b)
{
    if(a._sorted == false) a.sort();
    if(b._sorted == false) b.sort();
    if(a._key[0] < b._key[0]) return true;
    if(b._key[0] < a._key[0]) return false;
    if(a._key[1] < b._key[1]) return true;
    if(b._key[1] < a._key[1]) return false;
    if(a._key[2] < b._key[2]) return true;
    if(b._key[2] < a._key[2]) return false;
    return false;
}

template<typename T>
bool operator == (const KeyIndex3<T>& a,const KeyIndex3<T>& b)
{
    return !(a<b) && !(b<a);
}

template<typename T>
bool operator != (const KeyIndex3<T>& a,const KeyIndex3<T>& b)
{
    return !(a==b);
}

/******************************************************************************
*                                 KeyIndexM                                       *
******************************************************************************/

template<typename T>
struct KeyIndexM
{
	KeyIndexM() {}
	KeyIndexM(const std::vector<int>& key)
		:_key(key)
	{
		std::sort(_key.begin(), _key.end());
	}
	bool operator < (const KeyIndexM& b) const;
	std::vector<T> _key;
};

template<typename T>
bool KeyIndexM<T>::operator < (const KeyIndexM<T>& b) const
{
	int size = _key.size() < b._key.size() ? _key.size() : b._key.size();
	for (int i = 0;i < size;i++)
	{
		if (_key[i] < b._key[i])
		{
			return true;
		}
		else if (_key[i] > b._key[i])
		{
			return false;
		}
	}
	if (_key.size() < b._key.size())
	{
		return true;
	}
	else
	{
		return false;
	}
}

NS_END
#endif
