#pragma once
#ifndef PARSER_H_
#define PARSER_H_


#include "option_parser.h"
#include <mpi.h>


class  Parser {
private:
    OptionParser* optionParser_;


private:
    Parser () {}
public:
    Parser (int  argc, char* argv[])
    {
        const char  program_name[] = "main";
        const char  program_args[] = "[options...]";
        optionParser_ = new OptionParser(program_name, program_args);

//        std::cout << "argc " << argc << std::endl;
//        for (int i=0; i<argc; i++) {
//            std::cout << "argv[] = " << argv[i] << std::endl;
//        }


        const int   ret = optionParser_->parse_args(argc, argv);
        if (ret != 1)                               { exit(2); }
        if (! optionParser_->check_narguments(0))   { exit(0); }

    }

    ~Parser () { delete optionParser_; }


public:
    const OptionParser*  optionParser ()    const   { return  optionParser_; }
};


#endif
