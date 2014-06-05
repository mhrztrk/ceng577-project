/*
 * debug.h
 *
 *  Created on: Jun 2, 2014
 *      Author: mhrztrk
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#define DEBUG

#ifdef DEBUG
# define dbg_printf(...) printf(__VA_ARGS__)
#else
# define dbg_printf(...) do {} while (0)
#endif

#endif /* DEBUG_H_ */
