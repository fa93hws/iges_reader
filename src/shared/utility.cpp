//*************************************************************************
// File:     utility.cpp
// Author:   Yan LIU
// Date:     01-10-2013
// Email:    liuyan0411@gmail.com
//*************************************************************************


#include "utility.h"
#include <cstdlib>
#include <cstring>


NS_BEGIN

/******************************************************************************
*                                 Debug                                       *
******************************************************************************/
int Debug::CellID = -1;

/******************************************************************************
*                                 Random                                      *
******************************************************************************/

unsigned int Random::randomseed = 0;

void Random::SetSeed(unsigned int seed)
{
    randomseed = seed;
}

unsigned int Random::Rand(unsigned int choices)
{
    unsigned long newrandom;

    if (choices >= 714025l) {
        newrandom = (randomseed * 1366l + 150889l) % 714025l;
        randomseed = (newrandom * 1366l + 150889l) % 714025l;
        newrandom = newrandom * (choices / 714025l) + randomseed;
        if (newrandom >= choices) {
            return newrandom - choices;
        } else {
            return newrandom;
        }
    } else {
        randomseed = (randomseed * 1366l + 150889l) % 714025l;
        return randomseed % choices;
    }
    // Old function.
    // randomseed = (randomseed * 1366l + 150889l) % 714025l;
    // return randomseed / (714025l / choices + 1);
}

unsigned int Random::GetSeed()
{
    return randomseed;
}



NS_END
