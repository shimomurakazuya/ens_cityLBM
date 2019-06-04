/**
 * @file option_parser-def.h
 * @brief (File brief)
 *
 * (File explanation)
 *
 * @author Takashi Shimokawabe
 * @date 2010/08/26 Created
 * @version 1.0.0
 *
 * $Id: option_parser-def.h,v 78b6ad3f61e2 2011/01/06 11:59:53 shimokawabe $
 */

#ifndef OPTION_PARSER_HELP
#define OPTION_PARSER_HELP(name , opt, long_opt) /* */
#endif
#ifndef OPTION_PARSER_VERSION
#define OPTION_PARSER_VERSION(name , opt, long_opt, message) /* */
#endif
#ifndef OPTION_PARSER_MESSAGE
#define OPTION_PARSER_MESSAGE(name , opt, long_opt, message, help) /* */
#endif

#ifndef OPTION_PARSER_OPT
#define OPTION_PARSER_OPT(name, opt, long_opt, help) /* */
#endif
#ifndef OPTION_PARSER_OPTARG
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) /* */
#endif
#ifndef OPTION_PARSER_OPTARGS
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) /* */
#endif
#ifndef OPTION_PARSER_OPTCHOICE
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) /* */
#endif


#define G(...) __VA_ARGS__


#include "option_parser-user-def.h"


#undef G

#undef OPTION_PARSER_HELP
#undef OPTION_PARSER_VERSION
#undef OPTION_PARSER_MESSAGE
#undef OPTION_PARSER_OPT
#undef OPTION_PARSER_OPTARG
#undef OPTION_PARSER_OPTARGS
#undef OPTION_PARSER_OPTCHOICE
