/**
 * @file option_parser-user-def.h
 * @brief (File brief)
 *
 * (File explanation)
 *
 * @author Takashi Shimokawabe
 * @date 2010/08/26 Created
 * @version 0.1.0
 *
 * $Id: option_parser-user-def.h,v 8c9c367349ca 2011/01/06 10:25:12 shimokawabe $
 */


/*
  Available type:
      int, float, double, std::string.
      You may also use other types. However, you can not use "char *".

  Available options:

  For Help.
      OPTION_PARSER_HELP   ( help    , "-h" , "--help")

  For Version information.
      OPTION_PARSER_VERSION( version , "-v" , "--version", "Version 10.0.2")

  Only print message.
      OPTION_PARSER_MESSAGE( author  , "-a" , "--author" , "Author: Takashi Shimokawabe", "Show the author")

  For a flag (required no argument).
      OPTION_PARSER_OPT    ( flag   , "-f" , "--flag"    ,                                "flag option" )

  For a option required an argument.
      OPTION_PARSER_OPTARG ( intopt , ""   , "--intopt"  , int        , 2000     , "N"  , "int option" )

  For a option required more than one argument.
  Use the keyword "G" for more than one argument.
      OPTION_PARSER_OPTARGS( intopts, ""   , "--intopts", int        , 2, G({10, 20}), G({"N1", "N2"}), "int options" )

  For a option required more than one specified keyword.
      OPTION_PARSER_OPTCHOICE( choice, "-c" , "--choice" , int        , 3, G({10, 11, 12}), 10, "C", "Choose from 10, 11, 12")

*/


// default
OPTION_PARSER_HELP   ( help    , "-h" , "--help")
OPTION_PARSER_VERSION( version , "-v" , "--version", "Version 10.0.2")
OPTION_PARSER_MESSAGE( author  , "-a" , "--author" , "Author: Takashi Shimokawabe", "Show the author")

//OPTION_PARSER_OPT    ( flag    , "-f" , "--flag"   ,                                "flag option" )
//OPTION_PARSER_OPTARG ( intopt  , "-i" , "--intopt" , int        , 3000     , "N"  , "Requires one argument" )
//OPTION_PARSER_OPTARG ( dopt    , "-d" , "--dopt"   , double     , 100.0    , "N"  , "Requires one argument" )
//OPTION_PARSER_OPTARG ( uopt    , "-u" , "--uopt"   , unsigned int, 20     , "N"  , "Requires one argument" )
//OPTION_PARSER_OPTARG ( str     , "-s" , "--str"    , std::string, "string" , "STR", "Requires one argument" )
//OPTION_PARSER_OPTARGS( intopts , "-i2", "--i2"     , int        , 2, G({10, 20}), G({"N1", "N2"}), "Requires two arguments" )
//
//OPTION_PARSER_OPTCHOICE( choice, "-c" , "--choice" , int        , 3, G({10, 11, 12}), 10, "C", "Choose from 10, 11, 12")


// user define //
// number of GPUs per node //
OPTION_PARSER_OPTARG (
	gpu_per_node,
	 "-gpu_per_node",
	"--gpu_per_node",
	int,
	3,
	"number of GPUs per node",
	"Requires one argument" )


OPTION_PARSER_OPTARGS (
	number_of_grid_point,
   	 "-number_of_grid_point",
   	"--number_of_grid_point",
	int, 2,
	G({16, 0}),
	G({"nx", "dim_dir" }),
	"Requires two arguments" )


OPTION_PARSER_OPTARGS (
	domain_min,
   	 "-domain_min",
   	"--domain_min",
	double, 3,
	G({0.0, 0.0, 0.0}),
	G({"xmin", "ymin", "zmin" }),
	"Requires three arguments" )


OPTION_PARSER_OPTARGS (
	domain_length,
   	 "-domain_length",
   	"--domain_length",
	double, 3,
	G({1.0, 1.0, 1.0}),
	G({"xlength", "ylength", "zlength" }),
	"Requires three arguments" )


OPTION_PARSER_OPTARG (
        time_end,
        "-time_end",
        "--time_end",
        double,
        10.0,
        "time_end",
        "Requires one argument" )


OPTION_PARSER_OPTARGS (
	velocity_lbm,
   	 "-velocity_lbm",
   	"--velocity_lbm",
	double, 2,
	G({10.0, 0.1}),
	G({"vel_ref", "cfl at vel_ref" }),
	"Requires two arguments" )


OPTION_PARSER_OPTARGS (
	cfr_steps,
   	 "-cfr_steps",
   	"--cfr_steps",
	int, 2,
	G({100, 200}),
	G({"cout_step", "fout_step" }),
	"Requires two arguments" )


OPTION_PARSER_OPTARGS (
	cfr_flags,
   	 "-cfr_flags",
   	"--cfr_flags",
	bool, 2,
	G({1, 1}),
	G({"cout_flag", "fout_flag" }),
	"Requires two arguments" )


OPTION_PARSER_OPTARGS (
	restart_flags_and_step,
   	 "-restart_flags_and_step",
   	"--restart_flags_and_step",
	int, 2,
	G({0, 1}),
	G({"restart_flag", "restart_step" }),
	"Requires two arguments" )

