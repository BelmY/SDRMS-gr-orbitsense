/*
 * log.h
 *
 *  Created on: Jun 7, 2019
 *      Author: nkaramolegos
 */

#ifndef INCLUDE_ORBITSENSE_LOG_H_
#define INCLUDE_ORBITSENSE_LOG_H_


#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/syscall.h>

#define ORBITSENSE_MESSAGES 0
#define ORBITSENSE_DEBUG_MESSAGES 1

#if ORBITSENSE_MESSAGES
#define ORBITSENSE_LOG_CLASS_INFO(CLASS_DEBUG_ENABLE, M, ...)                               \
        do {                                                                            \
                if(CLASS_DEBUG_ENABLE) {                                                \
                                fprintf(stderr, "[INFO]: " M " \n", ##__VA_ARGS__);     \
                }                                                                       \
        } while(0)

#else
#define ORBITSENSE_LOG_CLASS_INFO(CLASS, M, ...)
#endif

#if ORBITSENSE_MESSAGES
#define ORBITSENSE_LOG_INFO(M, ...)                                                         \
                fprintf(stderr, "[INFO]: " M " \n", ##__VA_ARGS__)

#else
#define ORBITSENSE_LOG_INFO(M, ...)
#endif

#define ORBITSENSE_ERROR(M, ...)                                                    \
        fprintf(stderr, "[ERROR] %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define ORBITSENSE_WARN(M, ...)                                                             \
        fprintf(stderr, "[WARNING] %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#if ORBITSENSE_DEBUG_MESSAGES
#define ORBITSENSE_DEBUG(M, ...)                                                    \
        fprintf(stderr, "[DEBUG]: " M "\n", ##__VA_ARGS__)
#else
#define ORBITSENSE_DEBUG(M, ...)
#endif



#endif /* INCLUDE_ORBITSENSE_LOG_H_ */
