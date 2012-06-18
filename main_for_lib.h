/* Copyright 2005-2007 The MathWorks, Inc. */

#ifndef _MAIN_H_
#define _MAIN_H_
#ifndef mclmcrrt_h
/* Defines the proxy layer. */
#include "mclmcrrt.h"
#endif
typedef struct
{
    int ac;
    const char** av;
    int err;
} inputs;

#endif
