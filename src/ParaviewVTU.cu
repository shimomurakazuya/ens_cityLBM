#include "ParaviewVTU.h"


// OutputVTU //
void ParaviewVTU::OutputVTUFiles(int step)
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    const int64_t np = number_of_points_;
    const int64_t nc = number_of_cells_;

    // vtu //
    std::ofstream fout;
    fout.open(io_folder() + "/" + fname_vtu(step, rank_));

    HeaderVTUFile(fout);

    HeaderData(fout, "UnstructuredGrid");
    HeaderPieceFile(fout);

    HeaderData(fout, "Points"); // xyz //
    SetDataArray(fout, "geometry", VTKDataType::Float32, 3, (BYTE*)xyz_, np*3);
    FooterData(fout, "Points");

    HeaderData(fout, "Cells"); // connections, offsets, types //
    SetDataArray(fout, "connectivity", VTKDataType::Int64, 1, (BYTE*)connectivity_, nc*8);
    SetDataArray(fout, "offsets",      VTKDataType::Int64, 1, (BYTE*)offsets_,      nc);
    SetDataArray(fout, "types",        VTKDataType::UInt8, 1, (BYTE*)types_,        nc);

    FooterData(fout, "Cells");


    HeaderData(fout, "PointData");
    for (auto& elem : pointData_) { SetDataArray(fout, elem.plot_name, elem.dataType, elem.num_components, elem.val, np * elem.num_components); }
    FooterData(fout, "PointData");

    HeaderData(fout, "CellData");
    for (auto& elem : cellData_) { SetDataArray(fout, elem.plot_name, elem.dataType, elem.num_components, elem.val, nc * elem.num_components); }
    FooterData(fout, "CellData");

    FooterPieceFile(fout);
    FooterData(fout, "UnstructuredGrid");

    // appended data//
    HeaderAppendedData(fout);
    WriteAppendedData(fout);
    FooterAppendedData(fout);

    FooterVTUFile(fout);

    fout.close();
}


void ParaviewVTU::OutputPVTUFiles(int step)
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    // vtu //
    std::ofstream fout;
    fout.open(io_folder() + "/" + fname_pvtu(step));

    HeaderPVTUFile(fout);

    HeaderData(fout, "PUnstructuredGrid GhostLevel=\"0\"");

    HeaderData(fout, "PPoints"); // xyz //
    SetPDataArray(fout, "geometry", VTKDataType::Float32, 3);
    FooterData(fout, "PPoints");

    HeaderData(fout, "PCells"); // connections, offsets, types //
    SetPDataArray(fout, "connectivity", VTKDataType::Int64, 1);
    SetPDataArray(fout, "offsets",      VTKDataType::Int64, 1);
    SetPDataArray(fout, "types",        VTKDataType::UInt8, 1);
    FooterData(fout, "PCells");


    HeaderData(fout, "PPointData");
    for (auto& elem : pointData_) { SetPDataArray(fout, elem.plot_name, elem.dataType, elem.num_components); }
    FooterData(fout, "PPointData");

    HeaderData(fout, "PCellData");
    for (auto& elem : cellData_) { SetPDataArray(fout, elem.plot_name, elem.dataType, elem.num_components); }
    FooterData(fout, "PCellData");

    int  nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    SetPieceSource(fout, step, nprocs);

    FooterData(fout, "PUnstructuredGrid");

    FooterVTUFile(fout);

    fout.close();
}


// vtu //
void ParaviewVTU::SetDataArray(std::ofstream& fout, std::string name, VTKDataType type, int number_of_components, BYTE* val, uint64_t total_elem)
{
    const uint64_t total_byte = total_elem * get_sizeof(type);

    fout << " ";
    fout << "<DataArray" << " "
         << "NumberOfComponents=\"" << number_of_components << "\"" << " "
         << "type=\"" << cout_type(type) << "\"" << " "
         << "Name=\"" << name << "\"" << " "
         << "format=\"appended\"" << " "
         << "offset=\"" << offset_ << "\">"
         << "</DataArray>"
         << std::endl;

    offset_ += get_sizeof(header_type_) + total_byte;


    add_appended_data<uint64_t>( data_, total_byte );
    add_appended_data<BYTE>    ( data_, val, total_byte );
}


void ParaviewVTU::HeaderVTUFile(std::ofstream& fout)
{
    fout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"" << cout_type(header_type_) << "\"> \n";
}


void ParaviewVTU::FooterVTUFile(std::ofstream& fout)
{
    fout << "</VTKFile>\n";
}


void ParaviewVTU::HeaderPieceFile(std::ofstream& fout)
{
    fout << "<Piece NumberOfPoints=\"" << number_of_points_ << "\" NumberOfCells=\"" << number_of_cells_ << "\"> \n";
}


void ParaviewVTU::FooterPieceFile(std::ofstream& fout)
{
    fout << "</Piece>";
}


void ParaviewVTU::HeaderAppendedData(std::ofstream& fout)
{
    fout << "<AppendedData encoding=\"raw\">\n";
    fout << "_";
}


void ParaviewVTU::FooterAppendedData(std::ofstream& fout)
{
    fout << "</AppendedData>\n";
}


void ParaviewVTU::HeaderData(std::ofstream& fout, std::string str)
{
    fout << "<" << str << ">" << std::endl;
}


void ParaviewVTU::FooterData(std::ofstream& fout, std::string str)
{
    fout << "</" << str << ">" << std::endl;
}


void ParaviewVTU::WriteAppendedData(std::ofstream& fout)
{
    fout.write( ( char * )data_.data(), data_.size() * sizeof(BYTE)/sizeof(char) );
}


// pvtu //
void ParaviewVTU::SetPDataArray(std::ofstream& fout, std::string name, VTKDataType type, int number_of_components)
{
    fout << " ";
    fout << "<PDataArray" << " "
         << "NumberOfComponents=\"" << number_of_components << "\"" << " "
         << "type=\"" << cout_type(type) << "\"" << " "
         << "Name=\"" << name << "\"/>"
         << std::endl;
}


void ParaviewVTU::HeaderPVTUFile(std::ofstream& fout)
{
    fout << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"> \n";
}


void ParaviewVTU::SetPieceSource(std::ofstream& fout, int step, int nprocs)
{
    for (int i=0; i<nprocs; i++) {
        fout << " <Piece Source=\"" << fname_vtu(step, i) << "\"/>" << std::endl;
    }
}

