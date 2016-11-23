//*************************************************************************
// File:     config.h
// Author:   Yan LIU
// Date:     01-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************

#ifndef _CONFIG_H_01_10_2013__
#define _CONFIG_H_01_10_2013__

#include <iostream>
#include <stdlib.h>
//macroes

//#define LINUX 
//#define _SELFCHECK
//#define _VS

//namespace 
#define NS UNSW_OCTREE
#define NS_BEGIN namespace UNSW_OCTREE {
#define NS_END }
#define NSN_BEG namespace NURBS {
#define NSI_BEG namespace IGES {
//isnan
#include <math.h>

#ifdef _VS
#define ISNAN					_isnan
#else
#define ISNAN					isnan
#endif

//assert
//#include <assert.h>

//#define ASSERT(arg)                  assert(arg)
#ifdef ASSERT
#undef ASSERT
#define ASSERT(arg)                  if((arg)==false){std::cout<<__FILE__<<": "<<__LINE__<<" Error!"<<std::endl; int c; std::cin>>c; exit(-1);}
#else
#define ASSERT(arg)                  if((arg)==false){std::cout<<__FILE__<<": "<<__LINE__<<" Error!"<<std::endl; int c; std::cin>>c; exit(-1);}
#endif

#ifdef _ASSERTE
#undef _ASSERTE
#define _ASSERTE(expr)               if((expr)==false){std::cout<<__FILE__<<": "<<__LINE__<<" Error!"<<std::endl; int c; std::cin>>c; exit(-1);}
#endif

#endif


NS_BEGIN
NS_END




