/**
 * @file option_parser.h
 * @brief Option Parser
 *
 * Option Parser
 *
 * @author Takashi Shimokawabe
 * @date 2011/02/08  1.2.2  Support double-precision and other types
 * @date 2011/01/22  1.2.0  Add enable_print() and disable_print()
 * @date 2011/01/08         Add function check_narguments()
 * @date 2011/01/06         Add functions for argument()
 * @date 2010/08/26         Created
 * @version 1.2.0
 *
 * $Id: option_parser.h,v 9410065036a6 2011/02/08 06:42:46 shimokawabe $
 */

#pragma once
#ifndef OPTION_PARSER_H
#define OPTION_PARSER_H

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdint>

#define OPTION_PARSER_LIB_VERSION "1.2.2"

/**
 * @brief Option Parser class
 *
 * Option Parser class
 */
class OptionParser {
public:
    /**
     * @brief constructor
     */
    OptionParser(const char *prog, const char *args);
    ~OptionParser();


    /**
     * @brief Parse arguments
     * @param[in] argc           The Number of arguments passed to the program.
     * @param[in] argv           Pointer to an array of character strings
                                 that contain the arguments.
     * @param[in] used_argv      Pointer to an array of used_array flags.
     * @param[in] save_arguments Positional arguments are saved into this class and
                                 these are accessible by function arguments()
                                 when this flag is true.
     * @param[out] used_argv     Pointer to an array of used_array flags.
                                 An element contains true after parsing options,
                                 when the corresponding argument were used for options.
     * @retval 1 Succeeded to parse
     * @retval 2 Print help message
     * @retval 3 Print version number
     * @retval 4 Print message text
     * @retval 0 Failed to parse
     */
    int parse_args(int argc, char **argv, bool *used_argv, bool save_arguments = false);

    /**
     * @brief Parse arguments
     *
     * Parse arguments.
     * Positional arguments are accessible by function arguments() after calling this function.
     *
     * @param[in] argc           The Number of arguments passed to the program.
     * @param[in] argv           Pointer to an array of character strings
     *                           that contain the arguments.
     * @retval 1 Succeeded to parse
     * @retval 2 Print help message
     * @retval 3 Print version number
     * @retval 4 Print message text
     * @retval 0 Failed to parse
     */
    int parse_args(int argc, char **argv);

    void print_usage() const;
    void print_help() const;
    void print_message(const char *message) const;

    void enable_print();
    void disable_print();

    void print_all_values() const;

    /**
     * @brief Get the number of positional arguments.
     * @return the number of positional arguments
     */
    int narguments() const;

    /**
     * @brief Check the number of positional arguments.
     * @param[in] narguments the number of positional arguments
     *                       the program expects to receive
     * @retval true  the number of positional arguments is narguments
     * @retval false the number of positional arguments is not narguments and
     *               print error messages.
     */
    bool check_narguments(int narguments) const;

    /**
     * @brief Get the index-th positional argument
     *
     * Return positional arguments when parse_args() was called with save_arguments = true.
     * The index should be the number from 0 to narguments() - 1.
     * When parse_args() was called with save_arguments = false, NULL was returned.
     *
     * @return positional argument
     */
    const char *argument(int index) const;

    void delete_arguments();

#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
    bool name() const { return name ## _; }
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    type name() const { return name ## _; }
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    type name(int i) const { return name ## _[i]; }
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    type name() const { return name ## _; }
#include "option_parser-def.h"

#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
    bool name ## _is_given() const { return name ## _is_given_; }
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    bool name ## _is_given() const { return name ## _is_given_; }
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    bool name ## _is_given() const { return name ## _is_given_; }
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    bool name ## _is_given() const { return name ## _is_given_; }
#include "option_parser-def.h"

private:
    bool parse_args_for_opt(int argc, char **argv,
                            const char *opt, const char *long_opt,
                            bool *used_argv);
    template <typename T>
    int parse_args_for_optarg(int argc, char **argv, const char *opt,
                              T *val, std::string *input_val, bool *used_argv);
    template <typename T>
    int parse_args_for_optargs(int argc, char **argv, const char *opt,
                               int n_args,
                               T *vals, std::string *input_vals, bool *used_argv,
                               int *invalid_val_pos);
    template <typename T>
    int parse_args_for_optchoice(int argc, char **argv, const char *opt,
                                 int n_choice, const T *choice,
                                 T *val, std::string *input_val, bool *used_argv);

    template <typename T>
    int parse_args_for_optargs_choice(int argc, char **argv, const char *opt,
                                      int n_args,
                                      int n_choice, const T *choice,
                                      T *vals, std::string *input_vals, bool *used_argv,
                                      int *invalid_val_pos);
    template <typename T>
    bool check_choice(const T &val, int n_choice, const T *choice);


    void print_help_for_opt(const char *opt, const char *long_opt, const char *help) const;

#define OPTION_PARSER_HELP(name , opt, long_opt)        \
    void print_help_for_ ## name() const ;
#define OPTION_PARSER_VERSION(name , opt, long_opt, message) \
    void print_help_for_ ## name() const ;
#define OPTION_PARSER_MESSAGE(name , opt, long_opt, message, help) \
    void print_help_for_ ## name() const ;
#define OPTION_PARSER_OPT(name, opt, long_opt, help)    \
    void print_help_for_ ## name() const ;
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    void print_help_for_ ## name() const ;
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    void print_help_for_ ## name() const ;
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    void print_help_for_ ## name() const ;
#include "option_parser-def.h"

    void print_error_for_number_of_positional_arguments(int narguments) const;
    void print_error_for_requiring_args(const char *opt, const char *type, int n_args) const;
    void print_error_for_invalid_arg(const char *opt, const char *type, const char *input_val) const;
    template <typename T>
    void print_error_for_requiring_choice(const char *opt, const char *type,
                                          int n_choice, const T *choice) const;
    template <typename T>
    void print_error_for_invalid_choice(const char *opt, const char *type,
                                        const char *input_val,
                                        int n_choice, const T *choice) const;
    template <typename T>
    void print_error_for_candidates_of_choice(int n_choice, const T *choice) const;

    void set_default_values();

private:
    const char *prog_;
    const char *args_;

    int narguments_;
    char **arguments_;

    bool noprint_;

#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
    bool name ## _;
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    type name ## _;
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    type name ## _[n_args];
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    type name ## _;
#include "option_parser-def.h"

#define OPTION_PARSER_OPT(name, opt, long_opt, help) \
    bool name ## _is_given_;
#define OPTION_PARSER_OPTARG(name, opt, long_opt, type, default_val, metavar, help) \
    bool name ## _is_given_;
#define OPTION_PARSER_OPTARGS(name, opt, long_opt, type, n_args, default_vals, metavars, help) \
    bool name ## _is_given_;
#define OPTION_PARSER_OPTCHOICE(name, opt, long_opt, type, n_choice, choice, default_val, metavar, help) \
    bool name ## _is_given_;
#include "option_parser-def.h"


    // Disallow the constructor with no argument
    OptionParser();

    // Disallow the copy constructor and assignment operator
    OptionParser(const OptionParser&);
    void operator=(const OptionParser&);
};


#endif /* OPTION_PARSER_H */
