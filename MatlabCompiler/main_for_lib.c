/* Copyright 2005-2007 The MathWorks, Inc. */

/*
 * This is the main wrapper for all C library functions.
 * This wrapper has conditional code as is necessary for
 * macintosh platform.
 */

#include "main_for_lib.h" /* for the definition of the structure inputs */

int run_main(int ac, const char  *av[]);

int main(int ac, const char* av[])
{
    mclmcrInitialize();
    return mclRunMain((mclMainFcnType)run_main,ac,av);
}
