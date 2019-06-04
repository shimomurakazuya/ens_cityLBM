#include "CreateTreeMPI.h"
#include "Parser.h"


namespace CreateTreeMPI {


void create_tree_data(int argc, char* argv[])
{
    Parser  parser(argc, argv);

    Field  field(parser.optionParser());
    field.preset_field();

    field.write_field(0);
}


void read_tree_data(int argc, char* argv[])
{
    Parser  parser(argc, argv);

    Field  field(parser.optionParser());
    field.read_field(0);

    field.write_field(1);
}


};
