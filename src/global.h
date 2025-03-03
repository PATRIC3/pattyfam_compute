#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <boost/timer/timer.hpp>

#ifdef DEFINE_GLOBALS
#define G_EXTERN
#else
#define G_EXTERN extern
#endif

namespace boost { namespace program_options { class variables_map; } }

G_EXTERN boost::timer::cpu_timer g_timer;
G_EXTERN boost::program_options::variables_map *g_parameters;

#ifdef BLCR_SUPPORT
#include "libcr.h"
G_EXTERN cr_client_id_t g_cr_client_id;
#endif

#endif /* _GLOBAL_H */
