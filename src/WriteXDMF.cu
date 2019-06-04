#if 0

#include "WriteXDMF.h"


namespace  WriteXDMF {


void
write_xdmf_header (
    FILE*   xmf
    )
{
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
}


void
write_xdmf_footer (
    FILE*   xmf
    )
{
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
}


void
write_xdmf_datas (
          FILE*     xmf,
    const std::vector<WriteFormat>& geometries,
    const std::vector<WriteFormat>& attributes,
    const std::vector<WriteFormat>& lbm_attributes,
    const std::vector<WriteFormat>& grid_attributes
    )
{
    constexpr int dim_geometry  = 3;
    const int  number_of_leaves   = geometries.size()      / dim_geometry;
    const int  dim_attribute      = attributes.size()      / number_of_leaves;
    const int  dim_lbm_attribute  = lbm_attributes.size()  / number_of_leaves;
    const int  dim_grid_attribute = grid_attributes.size() / number_of_leaves;

    for (int ii=0; ii<number_of_leaves; ii++) {
        const int  index_geometry       = ii*dim_geometry;
        const int  index_attribute      = ii*dim_attribute;
        const int  index_lbm_attribute  = ii*dim_lbm_attribute;
        const int  index_grid_attribute = ii*dim_grid_attribute;

        // output xmf //
        fprintf(xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n", geometries[index_geometry].nameList.grid_name.c_str());

        write_xdmf_geometry_x_y_z(xmf, geometries,      index_geometry,       dim_geometry);
        write_xdmf_attribute     (xmf, attributes,      index_attribute,      dim_attribute);
        write_xdmf_attribute     (xmf, lbm_attributes,  index_lbm_attribute,  dim_lbm_attribute);
        write_xdmf_grid_attribute(xmf, grid_attributes, index_grid_attribute, dim_grid_attribute);

        fprintf(xmf, "   </Grid>\n");
    }
}


void
write_xdmf_datas_for_visualization (
          FILE*                     xmf,
    const std::vector<WriteFormat>& geometries,
    const std::vector<WriteFormat>& attributes,
    const std::vector<WriteFormat>& grid_attributes,
    const std::vector<bool>&        flag_visualization
    )
{
    constexpr int dim_geometry  = 3;
    const int  number_of_leaves   = geometries.size()      / dim_geometry;
    const int  dim_attribute      = attributes.size()      / number_of_leaves;
    const int  dim_grid_attribute = grid_attributes.size() / number_of_leaves;

    for (int ii=0; ii<number_of_leaves; ii++) {
        if ( !flag_visualization[ii] ) { continue; }

        const int  index_geometry       = ii*dim_geometry;
        const int  index_attribute      = ii*dim_attribute;
        const int  index_grid_attribute = ii*dim_grid_attribute;

        // output xmf //
        fprintf(xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n", geometries[index_geometry].nameList.grid_name.c_str());

        write_xdmf_geometry_x_y_z(xmf, geometries,      index_geometry,       dim_geometry);
        write_xdmf_attribute     (xmf, attributes,      index_attribute,      dim_attribute);
        write_xdmf_grid_attribute(xmf, grid_attributes, index_grid_attribute, dim_grid_attribute);

        fprintf(xmf, "   </Grid>\n");
    }
}


void
write_xdmf_datas_for_lv (
          FILE*                     xmf,
    const std::vector<WriteFormat>& geometries,
    const std::vector<WriteFormat>& attributes,
    const std::vector<WriteFormat>& grid_attributes,
    const std::vector<int>&         lv_leaf,
    const int                       lv
    )
{
    constexpr int dim_geometry  = 3;
    const int  number_of_leaves   = geometries.size()      / dim_geometry;
    const int  dim_attribute      = attributes.size()      / number_of_leaves;
    const int  dim_grid_attribute = grid_attributes.size() / number_of_leaves;

    for (int ii=0; ii<number_of_leaves; ii++) {
        if ( lv_leaf[ii] != lv ) { continue; }

        const int  index_geometry       = ii*dim_geometry;
        const int  index_attribute      = ii*dim_attribute;
        const int  index_grid_attribute = ii*dim_grid_attribute;

        // output xmf //
        fprintf(xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n", geometries[index_geometry].nameList.grid_name.c_str());

        write_xdmf_geometry_x_y_z(xmf, geometries,      index_geometry,       dim_geometry);
        write_xdmf_attribute     (xmf, attributes,      index_attribute,      dim_attribute);
        write_xdmf_grid_attribute(xmf, grid_attributes, index_grid_attribute, dim_grid_attribute);

        fprintf(xmf, "   </Grid>\n");
    }
}


void
write_xdmf_geometry_x_y_z (
          FILE*                     xmf,
    const std::vector<WriteFormat>& vec,
    const int                       index,
    const int                       dim
    )
{
    auto  is_header = [dim](const int i){ return  (i%dim == 0    ); };
    auto  is_footer = [dim](const int i){ return  (i%dim == dim-1); };

    for (int ii=0; ii<dim; ii++) {
        auto  elem = vec[index + ii];

        if (is_header(ii)) {
            fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", (int)elem.dims[0], (int)elem.dims[1], (int)elem.dims[2]);
            fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
        }

        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"%s\" Precision=\"%d\" Format=\"HDF\">\n", (int)elem.dims[0], (int)elem.dims[1], (int)elem.dims[2], elem.numberType().c_str(), elem.precision());
        fprintf(xmf, "        %s.h5:%s\n", elem.nameList.file_name.c_str(), elem.nameList.dataset_name().c_str());
        fprintf(xmf, "       </DataItem>\n");

        if (is_footer(ii)) {
            fprintf(xmf, "     </Geometry>\n");
        }
    }
}


void
write_xdmf_attribute (
          FILE*                     xmf,
    const std::vector<WriteFormat>& vec,
    const int                       index,
    const int                       dim
    )
{
    for (int ii=0; ii<dim; ii++) {
        auto  elem = vec[index + ii];

        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"%s\">\n", elem.nameList.attribute_name.c_str(), elem.nameList.pos_name.c_str());
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"%s\" Precision=\"%d\" Format=\"HDF\">\n", (int)elem.dims[0], (int)elem.dims[1], (int)elem.dims[2], elem.numberType().c_str(), elem.precision());
        fprintf(xmf, "        %s.h5:/%s\n", elem.nameList.file_name.c_str(), elem.nameList.dataset_name().c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
    }
}


void
write_xdmf_grid_attribute (
          FILE*                     xmf,
    const std::vector<WriteFormat>& vec,
    const int                       index,
    const int                       dim
    )
{
    for (int ii=0; ii<dim; ii++) {
        auto  elem = vec[index + ii];

        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"%s\">\n", elem.nameList.attribute_name.c_str(), elem.nameList.pos_name.c_str());
        fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"%s\" Precision=\"%d\" Format=\"HDF\">\n", (int)elem.dims[0], elem.numberType().c_str(), elem.precision());
        fprintf(xmf, "        %s.h5:/%s\n", elem.nameList.file_name.c_str(), elem.nameList.dataset_name().c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
    }
}


};

#endif
