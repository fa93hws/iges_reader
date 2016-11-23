//*************************************************************************
// File:     smartpool.h
// Author:   Yan LIU
// Date:     01-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#ifndef _SMARTPOOL_H_01_10_2013__
#define _SMARTPOOL_H_01_10_2013__

#include "config.h"
#include "utility.h"
#include <vector>

NS_BEGIN

template<typename ObjType, int PAGE_SIZE = 4096>
class StaticPool
{
public:
	StaticPool(int size = 0)
	{
		_pool.reserve(256);
		resize(size);
	}
	StaticPool(const StaticPool& pool)
	{
		_pool.reserve(256);
		operator = (pool);
	}
	StaticPool& operator = (const StaticPool& src)
	{
		if (&src == this) return *this;
		destroy();

		_size     = src._size;
		_capacity = src._capacity;
		int size = src._pool.size();
		_pool.resize(size, NULL);
		for (int i = 0;i < size;i++)
		{
			_pool[i] = allocate();
			for (int j = 0;j < PAGE_SIZE;j++)
			{
				_pool[i][j] = src._pool[i][j];
			}
		}
		return *this;
	}

	~StaticPool()
	{
		destroy();
	}
public:
	void destroy()
	{
		typename std::vector<ObjType*>::iterator iter = _pool.begin();
		for (;iter != _pool.end();iter++)
		{
			delete[](*iter);
		}
		_pool.clear();
		_size = 0;
		_capacity = 0;
	}

	void clear()
	{
		_size = 0;
	}

	void resize(int size)
	{
		destroy();
		_size = size;
		if (size <= 0) return;
		_pool.reserve(size / PAGE_SIZE + 1);
		while (size > 0)
		{
			size -= PAGE_SIZE;
			_pool.push_back(allocate());
		}
		_capacity = _pool.size() * PAGE_SIZE;
	}
	ObjType& operator[] (int index)
	{
		ASSERT(0<=index && index < _size);
		return getObj(index);
	}

	const ObjType& operator[] (int index) const
	{
		ASSERT(0 <= index && index < _size);
		return getObj(index);
	}

	int size() const
	{
		return _size;
	}

	ObjType& push_back(const ObjType& obj)
	{
		if (_size < _capacity)
		{
			ObjType& object = getObj(_size);
			_size++;
			object = obj;
			return object;
		}
		else
		{
			_capacity += PAGE_SIZE;
			_pool.push_back(allocate());
			return push_back(obj);
		}
	}
private:
	const ObjType& getObj(int index) const
	{
		ASSERT(0 <= index && index < _capacity);
		int block = index / PAGE_SIZE;
		return _pool[block][index%PAGE_SIZE];
	}
	ObjType& getObj(int index)
	{
		return const_cast<ObjType&>((static_cast<const StaticPool<ObjType, PAGE_SIZE>*>(this))->getObj(index));
	}
	ObjType* allocate()
	{
		ObjType* page = new ObjType[PAGE_SIZE];
		if (page == NULL)
		{
			std::cout << "memory limited\t" << std::endl;
			ASSERT(false);
		}
		return page;
	}
private:
	std::vector<ObjType*> _pool;
	int _size;
	int _capacity;
};


////functions required to be achieved include append, operator []
////BTW, an extra mark array is recommended control garbage
//template<typename Type,unsigned int PAGE_SIZE=4096>
//class SmartPool
//{
//public:
//	SmartPool(unsigned int size=4096);
//	~SmartPool();
//
//private:
//	SmartPool(const SmartPool& src);
//	SmartPool& operator=(const SmartPool& src);
//
//public:
//	Type*                 produce(const Type&);
//	void                  recycle(Type*);
//	unsigned int          size() const;
//	Type*                 random();
//	void                  reserve(unsigned int);
//
//	void                  begin();
//	bool 				  end();
//	Type*				  elem();
//	void                  move();
//
//	void                  clear();
//
//private:
//	std::vector<Type*>    _depot;
//	std::vector<Type*>    _bin;
//
//	Type*                 _shelf;
//	unsigned int          _pos;
//
//	unsigned int          _d_iter;
//	unsigned int          _s_iter;
//};
//
//template<typename Type,unsigned int PAGE_SIZE>
//SmartPool<Type, PAGE_SIZE>::SmartPool(unsigned int size)
//{
//	_depot.reserve(size);
//	_bin.reserve(size);
//
//	_shelf  = new Type[PAGE_SIZE];
//	ASSERT(_shelf != NULL);
//	_depot.push_back(_shelf);
//	_pos    = 0;
//
//	_d_iter = 0;
//	_s_iter = 0;
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//SmartPool<Type, PAGE_SIZE>::~SmartPool()
//{
//	unsigned int size = _depot.size();
//	for(unsigned int i=0;i<size;i++)
//	{
//		delete [] _depot[i];
//	}
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//Type*  SmartPool<Type, PAGE_SIZE>::produce(const Type& src)
//{
//	Type* elem;
//	if(!_bin.empty())
//	{
//		elem = _bin.back();
//		_bin.pop_back();
//	}else
//	{
//		if(_pos == PAGE_SIZE)
//		{
//			_shelf  = new Type[PAGE_SIZE];
//			ASSERT(_shelf != NULL);
//			_depot.push_back(_shelf);
//			_pos = 0;
//		}
//		elem = &(_shelf[_pos]);
//		_pos++;
//	}
//	*elem = src;
//	return elem;
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//void SmartPool<Type, PAGE_SIZE>::recycle(Type* elem)
//{
//	_bin.push_back(elem);
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//unsigned int SmartPool<Type, PAGE_SIZE>::size() const
//{
//	unsigned int cnt;
//	cnt  = (_depot.size()-1)*PAGE_SIZE + _pos;
//	cnt -= _bin.size();
//	return cnt;
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//void SmartPool<Type, PAGE_SIZE>::reserve(unsigned int size)
//{
//	_depot.reserve(size);
//	_bin.reserve(size);
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//Type* SmartPool<Type, PAGE_SIZE>::random()
//{
//	unsigned int s = (_depot.size()-1)*PAGE_SIZE + _pos;
//	if(s==0) return NULL;
//	unsigned int i = Random::Rand(s);
//	unsigned int p = i/PAGE_SIZE;
//	unsigned int j = i-p*PAGE_SIZE;
//	return &(_depot[p][j]);
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//void SmartPool<Type, PAGE_SIZE>::begin()
//{
//	_d_iter = 0;
//	_s_iter = 0;
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//bool SmartPool<Type, PAGE_SIZE>::end()
//{
//	if(_s_iter==_pos && _d_iter+1==_depot.size())
//	{
//		return true;
//	}else
//	{
//		return false;
//	}
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//Type* SmartPool<Type, PAGE_SIZE>::elem()
//{
//	ASSERT(_d_iter>=0 && _s_iter>=0 && !end());
//	return &(_depot[_d_iter][_s_iter]);
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//void SmartPool<Type, PAGE_SIZE>::move()
//{
//	ASSERT(_s_iter < PAGE_SIZE);
//	_s_iter++;
//	if(_s_iter==PAGE_SIZE)
//	{
//		if(!end())
//		{
//			_d_iter++;
//			_s_iter=0;
//		}
//	}
//}
//
//template<typename Type,unsigned int PAGE_SIZE>
//void SmartPool<Type, PAGE_SIZE>::clear()
//{
//	unsigned int size = _depot.size();
//	for(unsigned int i=1;i<size;i++)
//	{
//		delete [] _depot[i];
//	}
//	_depot.erase(_depot.begin()+1,_depot.end());
//	_bin.clear();
//	_shelf  = _depot[0];
//	_pos    = 0;
//	_d_iter = 0;
//	_s_iter = 0;
//}

NS_END
#endif
