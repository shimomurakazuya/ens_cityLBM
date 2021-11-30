/**
 * @file option_parser.cc
 * @brief Option Parser
 *
 * Option Parser
 *
 * @author Takashi Shimokawabe
 * @date 2011/02/08  1.2.2  Support double-precision and other types
 * @date 2011/01/22  1.2.0  Add enable_print() and disable_print()
 * @date 2011/01/22         Fix bugs for parser of integer options
 * @date 2011/01/08         Add function check_narguments()
 * @date 2011/01/06         Add functions for argument()
 * @date 2010/08/26         Created
 * @version 1.2.2
 *
 * $Id: option_parser.cc,v 9410065036a6 2011/02/08 06:42:46 shimokawabe $
 */

#include "option_parser.h"
#include <cstring>
#include <sstream>




OptionParser::OptionParser(const char *prog, const char *args)
    : prog_(prog), args_(args),
      narguments_(0), arguments_(NULL), noprint_(false)
{
    set_default_values();
}


OptionParser::~OptionParser()
{
    delete_arguments();
}

int OptionParser::parse_args(int argc, char **argv, bool *used_argv, bool save_arguments)
{
    for (int i=0; i<argc; i++) { used_argv[i] = false; }

#define OPTION_PARSER_HELP(name , opt, long_opt) \
    if (parse_args_for_opt(argc, argv, opt, long_opt, used_argv)) { \
        print_help();                                               \
        return 2;                                                   \
    }

#define OPTION_PARSER_VERSION(name , opt, long_opt, message)           \
    if (parse_args_for_opt(argc, argv, opt, long_opt, used_argv)) {    \
        print_message(message);                                        \
        return 3;                                                      \
    }

#define OPTION_PARSER_MESSAGE(name , opt, long_opt, message, help) \
    if (parse_args_for_opt(argc, argv, opt, long_opt, used_argv)) {    \
        print_message(message);                                        \
        return 4;                                                      \
    }

#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
    if (parse_args_for_opt(argc, argv, opt, long_opt, used_argv)) { \
        name ## _ = true;                                           \
        name ## _is_given_ = true;                                  \
    }

    const char *opts[2] = {NULL, NULL};
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    opts[0] = opt; opts[1] = long_opt;                                  \
    for (int i=0; i<2; i++) {                                           \
        std::string input_val;                                          \
        int ret = parse_args_for_optarg(argc, argv, opts[i], &name ## _, &input_val, used_argv); \
        if      (ret ==  1) { name ## _is_given_ = true; }              \
        else if (ret == -1) {                                           \
            print_error_for_invalid_arg(opts[i], # type, input_val.c_str()); \
            print_usage(); print_help_for_ ## name(); return 0; }   \
        else if (ret == -2) {                                           \
            print_error_for_requiring_args(opts[i], # type, 1);         \
            print_usage(); print_help_for_ ## name(); return 0; }   \
    }

#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    opts[0] = opt; opts[1] = long_opt;                                  \
    for (int i=0; i<2; i++) {                                           \
        std::string input_vals[n_args];                                 \
        int invalid_val_pos;                                            \
        int ret = parse_args_for_optargs(argc, argv, opts[i], n_args,   \
                                         name ## _, input_vals, used_argv,\
                                         &invalid_val_pos);             \
        if      (ret ==  1) { name ## _is_given_ = true; }              \
        else if (ret == -1) {                                           \
            print_error_for_invalid_arg(opts[i], # type, input_vals[invalid_val_pos].c_str()); \
            print_usage(); print_help_for_ ## name(); return 0; }   \
        else if (ret == -2) {                                           \
            print_error_for_requiring_args(opts[i], # type, n_args);    \
            print_usage(); print_help_for_ ## name(); return 0; }   \
    }

#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice0, default_val, metavar, help) \
    opts[0] = opt; opts[1] = long_opt;                                  \
    for (int i=0; i<2; i++) {                                           \
        std::string input_val;                                          \
        const type choice [] = choice0;                                 \
        int ret = parse_args_for_optchoice(argc, argv, opts[i],         \
                                           n_choice, choice,            \
                                           &name ## _, &input_val, used_argv); \
        if      (ret ==  1) { name ## _is_given_ = true; }              \
        else if (ret == -1) {                                           \
            print_error_for_invalid_choice(opts[i], # type, input_val.c_str(), n_choice, choice); \
            print_usage(); print_help_for_ ## name(); return 0; }   \
        else if (ret == -2) {                                           \
            print_error_for_requiring_choice(opts[i], # type, n_choice, choice); \
            print_usage(); print_help_for_ ## name(); return 0; }   \
    }

#include "option_parser-def.h"



    for (int i=1; i<argc; i++) { if (! used_argv[i]) narguments_++; }

    if (save_arguments) {
        arguments_ = new char *[narguments_];
        if (arguments_ == NULL) return -1;
        int index = 0;
        for (int i=1; i<argc; i++) {
            if (used_argv[i]) continue;

            arguments_[index] = new char[strlen(argv[i]) + 1];
            if (arguments_[index] == NULL) return -1;
            strcpy(arguments_[index], argv[i]);
            index++;
        }
    }

    return 1;
}

int OptionParser::parse_args(int argc, char **argv)
{
    bool *used_argv = new bool[argc];
    if (used_argv == NULL) return -1;

    const bool save_arguments = true;
    int ret = parse_args(argc, argv, used_argv, save_arguments);

    delete [] used_argv; used_argv = NULL;

    return ret;
}


int OptionParser::narguments() const
{
    return narguments_;
}

bool OptionParser::check_narguments(int narguments) const
{
//    std::cout << "narguments, narguments = " << narguments << ", " << narguments_ << std::endl;
    if (narguments == narguments_) return true;

    print_error_for_number_of_positional_arguments(narguments);
    print_usage();

    return false;
}

void OptionParser::delete_arguments()
{
    if (arguments_ == NULL) return;

    for (int i=0; i<narguments_; i++) {
        delete [] arguments_[i]; arguments_[i] = NULL;
    }
    delete [] arguments_; arguments_ = NULL;
}

const char *OptionParser::argument(int index) const
{
    return (arguments_ && index >= 0 && index < narguments_) ? arguments_[index] : NULL;
}

bool OptionParser::parse_args_for_opt(int argc, char **argv,
                                      const char *opt, const char *long_opt,
                                      bool *used_argv)
{
    if (strcmp("", opt) == 0 && strcmp("", long_opt) == 0) return 0;

    for (int i=1; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0 || strcmp(argv[i], long_opt) == 0) {
            used_argv[i] = true;
            return true;
        }
    }
    return false;
}


template <typename T>
int OptionParser::parse_args_for_optarg(int argc, char **argv, const char *opt,
                                        T *val, std::string *input_val, bool *used_argv)
{
    const int n_args = 1;
    T vals[n_args];
    std::string input_vals[n_args];
    int invalid_val_pos;
    int ret = parse_args_for_optargs(argc, argv, opt, n_args,
                                     vals, input_vals, used_argv,
                                     &invalid_val_pos);
    *input_val = input_vals[0];
    if (ret == 1) *val = vals[0];
    return ret;
}

template <typename T>
int OptionParser::parse_args_for_optargs(int argc, char **argv, const char *opt,
                                         int n_args,
                                         T *vals, std::string *input_vals, bool *used_argv,
                                         int *invalid_val_pos)
{
    const int n_choice = 0;
    const T *choice = NULL;
    return parse_args_for_optargs_choice(argc, argv, opt,
                                         n_args,
                                         n_choice, choice,
                                         vals, input_vals, used_argv,
                                         invalid_val_pos);
}

template <typename T>
int OptionParser::parse_args_for_optchoice(int argc, char **argv, const char *opt,
                                           int n_choice, const T *choice,
                                           T *val, std::string *input_val, bool *used_argv)
{
    const int n_args = 1;
    T vals[n_args];
    std::string input_vals[n_args];
    int invalid_val_pos;
    int ret = parse_args_for_optargs_choice(argc, argv, opt,
                                            n_args,
                                            n_choice, choice,
                                            vals, input_vals, used_argv,
                                            &invalid_val_pos);
    *input_val = input_vals[0];
    if (ret == 1) *val = vals[0];
    return ret;
}

template <typename T>
int OptionParser::parse_args_for_optargs_choice(int argc, char **argv, const char *opt,
                                                int n_args,
                                                int n_choice, const T *choice,
                                                T *vals, std::string *input_vals, bool *used_argv,
                                                int *invalid_val_pos)
{
    if (strcmp("", opt) == 0) return 0;

    for (int i=1; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            if (i >= argc - n_args)      return -2; // require args
            if (used_argv[i])            return 0;  // it means this is not option
            for (int j=0; j<n_args; j++) if (used_argv[i+1+j]) return -2; // args are already used

            for (int j=0; j<n_args; j++) input_vals[j] = argv[i+1+j];
            for (int j=0; j<n_args; j++) {
                std::istringstream istr(argv[i+1+j]);
                T val;
                istr >> val;
                if (istr.fail() || istr.eof() != true ||
                    (n_choice > 0 && !check_choice(val, n_choice, choice))) { // for choice
                    *invalid_val_pos = j; return -1; // invalid arg
                }
                vals[j] = val;
            }

            used_argv[i] = true;
            for (int j=0; j<n_args; j++) used_argv[i+1+j] = true;

            return 1;
        }
    }
    return 0; // nothing to do
}

template <typename T>
bool OptionParser::check_choice(const T &val, int n_choice, const T *choice)
{
    for (int i=0; i<n_choice; i++) {
        if (val == choice[i]) return true;
    }
    return false;
}

void OptionParser::print_usage() const
{
    if (noprint_) return ;
    fprintf(stdout, "Usage: %s %s\n", prog_, args_);
}

void OptionParser::print_help() const
{
    if (noprint_) return ;
    print_usage();
    fprintf(stdout, "Options:\n");
#define OPTION_PARSER_HELP(name , opt, long_opt)        \
    print_help_for_ ## name();
#define OPTION_PARSER_VERSION(name , opt, long_opt, message)    \
    print_help_for_ ## name();
#define OPTION_PARSER_MESSAGE(name , opt, long_opt, message, help)      \
    print_help_for_ ## name();
#define OPTION_PARSER_OPT(name, opt, long_opt, help)    \
    print_help_for_ ## name();
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    print_help_for_ ## name();
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    print_help_for_ ## name();
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    print_help_for_ ## name();
#include "option_parser-def.h"

}

void OptionParser::print_message(const char *message) const
{
    if (noprint_) return ;
    fprintf(stdout, "%s: %s\n", prog_, message);
}

void OptionParser::enable_print()
{
    noprint_ = false;
}

void OptionParser::disable_print()
{
    noprint_ = true;
}

// print_help_for_name
void OptionParser::print_help_for_opt(const char *opt, const char *long_opt,
                                      const char *help) const
{
    if (noprint_) return ;
    int l_opt  = strlen(opt);
    int l_lopt = strlen(long_opt);
    int l_dim  = 0;

    fprintf(stdout, "  ");
    if      (l_opt > 0 && l_lopt > 0) { l_dim = 2; fprintf(stdout, "%s, %-*s", opt, 24-l_dim-l_opt, long_opt); }
    else if (l_opt > 0)               {            fprintf(stdout, "%-24s", opt); }
    else if (l_lopt > 0)              {            fprintf(stdout, "%-24s", long_opt); }

    if (l_opt + l_dim + l_lopt >= 24) { fprintf(stdout, "\n  %24s", ""); }
    fprintf(stdout, "%s\n", help);
}

#define OPTION_PARSER_HELP(name , opt, long_opt) \
void OptionParser::print_help_for_ ## name() const   \
{                                                    \
    return print_help_for_opt(opt, long_opt, "Show this help message and exit");  \
}
#define OPTION_PARSER_VERSION(name , opt, long_opt, message) \
void OptionParser::print_help_for_ ## name() const   \
{                                                    \
    return print_help_for_opt(opt, long_opt, "Display the version number and exit"); \
}
#define OPTION_PARSER_MESSAGE(name , opt, long_opt, message, help) \
void OptionParser::print_help_for_ ## name() const   \
{                                                    \
    return print_help_for_opt(opt, long_opt, help);  \
}
#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
void OptionParser::print_help_for_ ## name() const   \
{                                                    \
    return print_help_for_opt(opt, long_opt, help);  \
}
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
void OptionParser::print_help_for_ ## name() const      \
{                                                       \
    const char *optarg  = strlen(opt)      > 0          \
        ? opt " " metavar                               \
        : opt;                                          \
    const char *loptarg = strlen(long_opt) > 0          \
        ? long_opt " " metavar                          \
        : long_opt;                                     \
    return print_help_for_opt(optarg, loptarg, help);   \
}
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
void OptionParser::print_help_for_ ## name() const                      \
{                                                                       \
    const char *mvars[] = metavars;                                     \
    std::string mvstr(" ");                                             \
    mvstr += mvars[0];                                                  \
    for (int i=1; i<n_args; i++) { mvstr += " "; mvstr += mvars[i]; }   \
    std::string optargs(opt), loptargs(long_opt);                       \
    if (strlen(opt)      > 0) { optargs  += mvstr; }                    \
    if (strlen(long_opt) > 0) { loptargs += mvstr; }                    \
    return print_help_for_opt(optargs.c_str(), loptargs.c_str(), help); \
}
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help)
#include "option_parser-def.h"



void OptionParser::print_error_for_number_of_positional_arguments(int narguments) const
{
    if (noprint_) return ;
    if (narguments == 0 || narguments == 1)
        fprintf(stdout, "%s: requires %d positional argument.\n", prog_, narguments);
    else
        fprintf(stdout, "%s: requires %d positional arguments.\n", prog_, narguments);
}

void OptionParser::print_error_for_requiring_args(const char *opt, const char *type,
                                                  int n_args) const
{
    if (noprint_) return ;

    fprintf(stdout, "%s: option %s requires ", prog_, opt);
    if (n_args > 1) fprintf(stdout, "%d ", n_args);

    if      (strcmp(type, "int")         == 0) fprintf(stdout, "integer");
    else if (strcmp(type, "float")       == 0) fprintf(stdout, "floating");
    else if (strcmp(type, "double")      == 0) fprintf(stdout, "floating");
    else if (strcmp(type, "std::string") == 0) fprintf(stdout, "string");
    else                                       fprintf(stdout, "%s", type);

    if (n_args > 1) fprintf(stdout, " value.\n");
    else            fprintf(stdout, " values.\n");
}

void OptionParser::print_error_for_invalid_arg(const char *opt, const char *type,
                                               const char *input_val) const
{
    if (noprint_) return ;

    fprintf(stdout, "%s: option %s: invalid ", prog_, opt);

    if      (strcmp(type, "int")    == 0) fprintf(stdout, "integer");
    else if (strcmp(type, "float")  == 0) fprintf(stdout, "floating");
    else if (strcmp(type, "double") == 0) fprintf(stdout, "floating");
    else                                  fprintf(stdout, "%s", type);

    fprintf(stdout, " value '%s'.\n", input_val);
}

template <typename T>
void OptionParser::print_error_for_requiring_choice(const char *opt, const char *type,
                                                    int n_choice, const T *choice) const
{
    const int n_args = 1;
    print_error_for_requiring_args(opt, type, n_args);
    print_error_for_candidates_of_choice(n_choice, choice);
}

template <typename T>
void OptionParser::print_error_for_invalid_choice(const char *opt, const char *type,
                                                  const char *input_val,
                                                  int n_choice, const T *choice) const
{
    if (noprint_) return ;

    fprintf(stdout, "%s: option %s: invalid choice '%s'.\n", prog_, opt, input_val);
    print_error_for_candidates_of_choice(n_choice, choice);
}

template <typename T>
void OptionParser::print_error_for_candidates_of_choice(int n_choice, const T *choice) const
{
    if (noprint_) return ;

    std::ostringstream ostr;
    ostr << " Choose from '" << choice[0] << "'";
    for (int i=1; i<n_choice; i++)
        ostr << ", '" << choice[i] << "'";

    fprintf(stdout, "%s.\n", ostr.str().c_str());
}


void OptionParser::set_default_values()
{

#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
    name ## _ = false;                                                  \
    name ## _is_given_ = false;
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    name ## _ = default_val;                                            \
    name ## _is_given_ = false;
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    type name ## dvals[] =  default_vals;                                          \
    for (int i=0; i<n_args; i++) { name ## _[i] = name ## dvals[i]; }   \
    name ## _is_given_ = false;
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    name ## _ = default_val;                                            \
    name ## _is_given_ = false;
#include "option_parser-def.h"

}

void OptionParser::print_all_values() const
{
    if (noprint_) return ;

#define OPTION_PARSER_OPT(name, opt, long_opt, help)                    \
    {                                                                   \
        std::ostringstream ostr;                                        \
        ostr << name ## _;                                              \
        fprintf(stdout, "%-20s: %-3s %-20s = '%s'\n", # name, opt, long_opt, ostr.str().c_str()); \
    }
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    OPTION_PARSER_OPT(name, opt, long_opt, help)
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    {                                                                   \
        std::ostringstream ostr;                                        \
        ostr << "'" << name ## _[0] << "'";                             \
        for (int i=1; i<n_args; i++) {                                  \
            ostr << ", '" << name ## _[i] << "'";                       \
        }                                                               \
        fprintf(stdout, "%-20s: %-3s %-20s = [%s]\n", # name, opt, long_opt, ostr.str().c_str()); \
    }
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    OPTION_PARSER_OPT(name, opt, long_opt, help)
#include "option_parser-def.h"

}



