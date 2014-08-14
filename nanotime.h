/* Remember to comment this with auto name tag and auto version tag,
   and maybe compilation instructions and auto license tag */


#pragma once
#ifndef NANOTIME_H
#define NANOTIME_H

/* Macro values taken from http://predef.sourceforge.net/ */
#if defined(WIN32) || defined(_WIN32) \
    || defined(__WIN32__) || defined(__WINDOWS__) \
    || defined(__TOS_WIN__)
    #define WINDOWS
    #include <windows.h>
#elif defined(__MACH__) || defined(__APPLE__)
    #define MACOSX
    #include <mach/mach_time.h>
#elif defined(linux) || defined(__linux)
    /* May need to compile with -lrt */
    #define LINUX
    #include <time.h>
/* TODO Add support for other real-time POSIX systems */
#elif defined(sun) || defined(__sun) || defined(_AIX)
    #define RTPOSIX
    #include <sys/time.h>
#else /* Unsupported OS */
    #if defined(_MSC_VER)
        /* MS Vis C++ */
        #pragma message "Unsupported OS - get_nanotime() will return 0"
    #else
        /* #warning May only be supported by GCC, not sure yet... */
        #warning "Unsupported OS - get_nanotime() will return 0"
    #endif
#endif

/* System Header Inclusion */
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Represents time in nanoseconds
 */
typedef uint64_t nanotime_t;


/**
 * Wrapper to get time (in ns) since an arbitrary time
 */
extern inline nanotime_t get_nanotime(void)
{
    #if defined(WINDOWS)
        LARGE_INTEGER time_var, frequency;
        QueryPerformanceCounter(&time_var);
        QueryPerformanceFrequency(&frequency);

        /* Convert to nanoseconds */
        return 1.0e9 * time_var.QuadPart / frequency.QuadPart;
        
    #elif defined(MACOSX)
        uint64_t time_var;
        mach_timebase_info_data_t info;
        
        time_var = mach_absolute_time();
        mach_timebase_info(&info);
        
        /* Convert to nanoseconds */
        return time_var * (info.numer / info.denom);
        
    #elif defined(LINUX)
        timespec time_var;
        clock_gettime(CLOCK_REALTIME, &time_var);
        
        /* Convert to nanoseconds */
        return time_var.tv_sec * 1.0e9 + time_var.tv_nsec;
        
    #elif defined(RTPOSIX)
        return gethrtime();
        
    #endif
    
    /* Unsupported OS */

    /* If compilation is halted (using #error), should maybe assert(0),
       but for now it isn't, because the point of this wrapper is to
       avoid worrying about OS dependant timing and focus on the code
       functionality itself */
    
    /* Zero chosen as error value due to it being impossible in
       practice to have zero time difference after a function call */
    return 0;
}


#ifdef __cplusplus
}
#endif

#endif /* NANOTIME_H */

