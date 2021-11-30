#include "Field.h"
#include "FuncLoop.h"
#include "FuncMath.h"
#include "FuncMapData.h"
#include "IOData.h"
#include <cmath>
#include "runtime_error.hpp"
#include "range.hpp"

// public //
void  Field::
init_ptrs()
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        itp_meshval_ [i] = 0;
        itpn_meshval_[i] = 1;

//        std::cout << i << " : itp, itpn = " << itp_meshval_[i] << ", " << itpn_meshval_[i] << std::endl;
    }

    tree_ptr_[0] = &tree0_;
    tree_ptr_[1] = &tree1_;
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        meshValues_ptr_[0][i] = &meshValues0_[i];
        meshValues_ptr_[1][i] = &meshValues1_[i];
    }

}


//void  Field::
//preset_field()
//{
//    if (rank_ == 0) { std::cout << "start :  preset_tree\n"; }
//
//    // grid information //
//    init_grid(grids_);
//
//    // connenctions //
//    preset_tree(tree0_);
//
//    // taskID //
//    init_taskID(taskID_, tree0_);
//
//    // elapsedTimeInfo //
//    init_elapsedTimeInfo(parameters_);
//
//    // values //
//    init_meshValues(grids_, tree0_, meshValues0_);
//    init_ValueStat(tree0_);
//    init_valueBuff(grids_, tree0_, valueBuff_);
//
//    // copy //
//    const bool  is_reset = true;
//    copy_tree(tree1_, tree0_, is_reset);
//
//
////    init_meshValues(grids_, tree0_, meshValues1_);
//    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);
//
//    // postprocessmonitor_ //
//    postprocessmonitor_.setupdata(rank_, tree0_, meshValues0_);
//
//    if (rank_ == 0) { std::cout << "end : preset_tree\n"; }
//}


//void  Field::
//init_field()
//{
//    // grid information //
//    init_grid(grids_);
//
//    // connenctions //
//    init_tree(tree0_);
//
//    // taskID //
//    init_taskID(taskID_, tree0_);
//
//    // elapsedTimeInfo //
//    init_elapsedTimeInfo(parameters_);
//
//    // values //
//    init_meshValues(grids_, tree0_, meshValues0_);
//    init_ValueStat(tree0_);
//    init_valueBuff(grids_, tree0_, valueBuff_);
//
//    // copy //
//    const bool  is_reset = true;
//    copy_tree(tree1_, tree0_, is_reset);
//
//    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);
//}


void Field::
init_field_wo_map() {
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    auto&& map = MapData();
    init_grid(grids_);
    init_tree_with_map(tree0_, map);
    init_others();
}

void Field::
init_field_with_map(const MapData& map) {
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    init_grid(grids_);
    init_tree_with_map(tree0_, map);
    init_others();
}

void Field::
init_others() {
    // taskID //
    init_taskID(taskID_, tree0_);

    // elapsedTimeInfo //
    init_elapsedTimeInfo(parameters_);

    // values //
    init_meshValues(grids_, tree0_, meshValues0_);
    init_ValueStat(tree0_);
    init_valueBuff(grids_, tree0_, valueBuff_);
    init_valuePBVR(valuePBVR_);
    init_ValueTimeAverage(tree0_, meshValues0_);

    // copy //
    const bool  is_reset = true;
    copy_tree(tree1_, tree0_, is_reset);


//    init_meshValues(grids_, tree0_, meshValues1_);
    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);

    // postprocessmonitor_ //
    postprocessmonitor_.setupdata(tree0_, meshValues0_);

    if (comm_.is_rank0()) { std::cout << "end : preset_tree\n"; }
}


void  Field::
read_field_full(const int step)
{
    const int rstep = read_latest_step(step);
    runtime_assert(step == rstep, "InvalidConfiguration: restart timestep mismatched");
    runtime_assert(rstep % 2 == 0, "InvalidConfiguration: invalid restart timestep; even number required"); // save coherency of double buffering

    if (comm_.is_rank0()) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "rstep = " << rstep << std::endl;
    }

    const auto rank = comm_.world().rank();

    // parameter //
    read_parameters(rstep, parameters_);

    // grid information //
    read_grid(rstep, grids_, rank);

    // connenctions //
    read_tree(rstep, tree0_, rank);

    // valueBuff //
    init_valueBuff(grids_, tree0_, valueBuff_);

    // taskID //
    read_taskID(rstep, taskID_, rank);

    // elapsedTimeInfo //
    read_elapsedTimeInfo(parameters_, rank);

    // postprocessmonitor_ //
    postprocessmonitor_.setupdata(tree0_, meshValues0_);

    // values //
    init_ValueTimeAverage(tree0_, meshValues0_);

    read_meshValues(rstep, grids_, tree0_, meshValues0_, valueTimeAverage1min_, rank);
    init_ValueStat(tree0_);
    init_valuePBVR(valuePBVR_);

    // copy //
    const bool  is_reset = true;
    copy_tree(tree1_, tree0_, is_reset);

    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);
}

void  Field::
read_field_params(const int step)
{
    const int rstep = read_latest_step(step);
    runtime_assert(step == rstep, "InvalidConfiguration: restart timestep mismatched");
    runtime_assert(rstep % 2 == 0, "InvalidConfiguration: invalid restart timestep; even number required"); // save coherency of double buffering

    if (comm_.is_rank0()) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "rstep = " << rstep << std::endl;
    }

    // parameter //
    read_parameters(rstep, parameters_);
}

void  Field::
read_field_meshval(const int step)
{
    const int rstep = read_latest_step(step);
    runtime_assert(step == rstep, "InvalidConfiguration: restart timestep mismatched");
    runtime_assert(rstep % 2 == 0, "InvalidConfiguration: invalid restart timestep; even number required"); // save coherency of double buffering

    if (comm_.is_rank0()) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "rstep = " << rstep << std::endl;
    }

    // parameter //
    read_parameters(rstep, parameters_);

    // values //
    //read_meshValues(rstep, grids_, tree0_, meshValues0_, valueTimeAverage1min_, rank_);
    //copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);

    //for(const auto& i: util::irange(DefAMR::LV_MAX)) {
    //    const auto& nn_max = DefAMR::NN_LEAF*tree0_.number_of_nodes_lv(i);
    //    meshValues0_[i].init(nn_max);
    //    meshValues1_[i].init(nn_max);
    //}
    IOData  ioData(comm_, step);
    ioData.readBinaries(grids_, tree0_, meshValues0_, nullptr, valueTimeAverage1min_);
    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, true);
}


void  Field::
write_field(const int step)
{

#ifdef NO_FIELD_WRITEDAT
#ifdef NO_IOTHREAD
    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT: " << step << std::endl; }
#endif
#else


    // oklahoma //
    //if ( (step >= 5 && step <= 119) || (step >= 180 && step <= 239) || (step >= 300 && step <= 359) ) { return; }
    #ifdef REGRESSION_TEST
    //if( step == ) { return; }
    #endif

    #if !defined(REGRESSION_TEST)
    const int skip_iofield = optionParser().iofield_freq();
    if(step % skip_iofield != 0) {
        if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip io field @ t=" << step << std::endl; }
        return;
    }
    #endif

    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": do t=" << step << std::endl; }

    const int  wstep = step;
    const auto rank = comm_.world().rank();

    write_latest_step(wstep, rank);

    // parameter //
    write_parameters(wstep, parameters_io_, rank);

    // grid information //
    write_grid(wstep, grids_io_, rank);

    // connenctions //
    write_tree(wstep, tree_io_, rank);

    // taskID //
    write_taskID(wstep, taskID_io_, rank);

    // values //
    write_meshValues(wstep, grids_io_, tree0_, parameters_io_, meshValues_io_, rank);
#endif // NO_FIELD_WRITEDAT
}


void  Field::
write_monitor(const int step)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  wstep = step;

    // postprocessmonitor_ //
    #ifdef NO_FIELD_WRITEDAT
    const MeshValue* meshValues_ptr[DefAMR::LV_MAX];
    for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
        meshValues_ptr[lv] = meshValues_ptr_[itpn_meshval_[lv]][lv]; // &(meshValues_new(lv))
    }
    postprocessmonitor_.OutputMonitorData(wstep, tree0_, parameters_, meshValues_ptr
        #ifdef USE_VALUE_STAT
        , valueStat_
        #endif // USE_VALUE_STAT
        );
    #else
    postprocessmonitor_.OutputMonitorData(wstep, tree0_, parameters_io_, meshValues_io_
        #ifdef USE_VALUE_STAT
        , valueStat_io_
        #endif // USE_VALUE_STAT
        );
    #endif
//    postprocessmonitor_.OutputMonitorData(wstep, tree0_, parameters_, meshValues0_, rank_);
}


void  Field::
copy_io_field()
{
    const auto rank = comm_.world().rank();
#ifdef NO_FIELD_WRITEDAT
    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT" << std::endl; }
#else

    static bool io_is_allocated = false;
    const bool  is_reset = (!io_is_allocated) ? true : false;

    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    // parameters //
    parameters_io_.copy(parameters_);

    // grid information //
    for (int i=0; i<DefAMR::LV_MAX; i++) { grids_io_[i].copy(grids_[i]); }

    // connenctions //
    copy_tree(tree_io_, tree0_, is_reset);
//    copy_tree(tree_io_, tree0_, true);

    // taskID //
    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": copyTaskID" << std::endl; }
    taskID_io_.copyTaskID(taskID_, is_reset);

    // values //
    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": copy_meshValues" << std::endl; }
    update_meshValues_io(grids_io_, tree_io_, is_reset);

    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": copy_valueStat" << std::endl; }
    copy_valueStat_io();
    for (auto& v: valueStat_) { v.zeroset(); } // reset when output

    copy_valueTimeAverage_io();

    if (!io_is_allocated) { io_is_allocated = true; }
#endif
}


void  Field::
swap_meshValue(const int lv)
{
    const int tmp = itp_meshval_[lv];
    itp_meshval_ [lv] = itpn_meshval_[lv];
    itpn_meshval_[lv] = tmp;
}


void  Field::
parameter_time(const int t)
{
    int  step = t;
    real dt   = parameters_.dt0();
    real time = t*dt;

    parameters_.time_step_update(step, time);
}


// copy //
void  Field::
copy_tree(Tree&  tree_new, const Tree&  tree, const bool  is_reset)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//    std::cout << "----------" << std::endl;
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    tree_new.copy_stNodeValArray (tree.number_of_nodes(), tree.stnodeValArray(), is_reset);
    tree_new.copy_node_structures(tree, is_reset);
    tree_new.copyMPIPutGetInfo(tree.mpiPutGetInfo(), is_reset);

//    std::cout << "----------" << std::endl;
}


void  Field::
copy_meshValues(
          MeshValue*    meshValue_new,
    const MeshValue*    meshValue,
    const Grid*         grids,
    const Tree&         tree,
    const bool          is_reset
    )
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//    std::cout << "----------" << std::endl;
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    for (int i=0; i<DefAMR::LV_MAX; i++) {
        meshValue_new[i].copy_MeshValue(
            meshValue[i],
            grids[i],
            tree,
            is_reset
            );
    }

//    std::cout << "----------" << std::endl;
}

void Field::
copy_valueStat_io() 
{
#ifdef NO_FIELD_WRITEDAT
    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT" << std::endl; }
#else
    for(int lv=0; lv<DefAMR::LV_MAX; lv++) {
         valueStat_io_[lv].clone_from(valueStat_[lv]);
    }
#endif
}

void Field::
copy_valueTimeAverage_io() 
{
#ifdef NO_FIELD_WRITEDAT
    if (rank== 0) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT" << std::endl; }
#else
    for(int lv=0; lv<DefAMR::LV_MAX; lv++) {
         valueTimeAverage1min_io_[lv].copy(valueTimeAverage1min_[lv].nn_max(), valueTimeAverage1min_[lv]);
    }
#endif
}

//void  Field::
//update_meshValues(
//    const Grid*         grids,
//    const Tree&         tree
//    )
//{
////    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//
//    const bool is_reset = false;
//
//    for (int i=0; i<DefAMR::LV_MAX; i++) {
//        meshValue(i).copy_MeshValue(
//            meshValue_new(i),
//            grids[i],
//            tree,
//            is_reset
//            );
//    }
//
////    std::cout << "----------" << std::endl;
//} = delete;

void  Field::
update_meshValues_io(
    const Grid*         grids,
    const Tree&         tree,
    const bool          is_reset
    )
{
#ifndef NO_FIELD_WRITEDAT
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    for (int i=0; i<DefAMR::LV_MAX; i++) {
        meshValues_io_[i].copy_MeshValue(meshValue_new(i), grids[i], tree, is_reset);
        #ifdef MESHVALUES_IO_MASK_HALO
        meshValues_io_[i].mask_meshValue_wo_halo(tree, i);
        meshValues_io_[i].mask_meshValue_wo_obj (tree, i);
        #endif
    }
#endif
}


//// private //
//void  Field::
//preset_tree(Tree&   tree)
//{
//    tree.preset_tree_data(grids());
//}


// init //
void  Field::
init_grid(Grid*  grids)
{
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    // optionParser //
    const int  nx_tmp  = optionParser_.number_of_grid_point(0);
    const int  dim_dir = optionParser_.number_of_grid_point(1);
    if ( !(dim_dir == 0 || dim_dir == 1 || dim_dir == 2) ) {
        std::cout << "error : " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    // arrange domain length //
    const double _domain_length_tmp[] = { static_cast<double>( optionParser_.domain_length(0) ),
                                          static_cast<double>( optionParser_.domain_length(1) ),
                                          static_cast<double>( optionParser_.domain_length(2) ) };

    const double dx_leaf_lv0 = _domain_length_tmp[dim_dir] / nx_tmp;
#if defined(VVTARGET_TestOklahoma)
#ifdef EXPECT_DX_MESH
    const double dx_mesh = dx_leaf_lv0 * std::pow(0.5, DefAMR::LV_MAX-1) / DefAMR::NX_LEAF;
    std::cout << "dx_mesh, EXPECT_DX_MESH = " << dx_mesh << ", " << EXPECT_DX_MESH << std::endl;
    runtime_assert(fabs(dx_mesh - EXPECT_DX_MESH) < 1e-5 * EXPECT_DX_MESH, "InvalidConfiguration: unexpected mesh resolution");
#endif
#endif

    const Vector3d<int>  nx0( static_cast<int>( ceil( _domain_length_tmp[0]/dx_leaf_lv0 ) ),
                              static_cast<int>( ceil( _domain_length_tmp[1]/dx_leaf_lv0 ) ),
                              static_cast<int>( ceil( _domain_length_tmp[2]/dx_leaf_lv0 ) )  );

    // domain //
    const Vector3d<double> offset( optionParser_.domain_min(0),
                                   optionParser_.domain_min(1),
                                   optionParser_.domain_min(2) );

    const double  dx0 = dx_leaf_lv0;

    // update //
    for (int l=0; l<DefAMR::LV_MAX; l++) {
        const int  ratio = pow(2, l);
        const Vector3d<int>  nx = nx0*ratio;
        const double         dx = dx0/ratio;

        grids[l].init(nx, offset, dx);
    }
}


//void  Field::
//init_tree(Tree&   tree)
//{
//    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//
//    tree.init_tree_data(grids());
//}

void Field::
init_tree_with_map(Tree& tree, const MapData& map) {
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    tree.init_tree_data_with_map(grids_, map);
}


void  Field::
init_elapsedTimeInfo(const Parameters& parameter)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    allTimeInfo_ .InitElapsedTimeInfo("all",  parameters_);
    funcTimeInfo_.InitElapsedTimeInfo("func", parameters_);
    mpiTimeInfo_ .InitElapsedTimeInfo("mpi",  parameters_);
}


void  Field::
init_parameters(Parameters& parameters)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    parameters.init(optionParser_);
    if (comm_.is_rank0()) { parameters.coutParameters(); }
}


void  Field::
init_taskID(TaskID& taskID, const Tree& tree)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    taskID.set_taskID(tree);
}


void  Field::
init_taskID_opt(TaskID& taskID, const Tree& tree, const MeshValue* meshValues)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    taskID.set_taskID_opt(tree, meshValues);
}


void  Field::
init_meshValues(const Grid* grids, const Tree& tree, MeshValue*  meshValues)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    const int n_scalars = optionParser_.number_of_scalars();

    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each( [&grids, &tree, &meshValues, n_scalars](const int lv){
        const int  nn_max = DefAMR::NN_LEAF*tree.number_of_nodes_lv(lv);
        meshValues[lv].init( nn_max, n_scalars );
        meshValues[lv].set_coordinate( lv, grids[lv], tree );
        } );

}

void Field::
init_ValueStat(const Tree& tree) {
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    for(int lv=0; lv<DefAMR::LV_MAX; lv++) {
        const std::intptr_t nn_max = DefAMR::NN_LEAF * tree.number_of_nodes_lv(lv);
        valueStat_   [lv].init(nn_max);
        #ifndef NO_FIELD_WRITEDAT
        valueStat_io_[lv].init(nn_max);
        #endif
    }
}

void Field::
init_valueBuff(const Grid* grids, const Tree& tree, ValueBuff& valueBuff) {
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    auto nn_max = DefAMR::NN_LEAF*tree.number_of_nodes_lv(0);
    for(int lv=1; lv<DefAMR::LV_MAX; lv++) {
        nn_max = std::max(nn_max, DefAMR::NN_LEAF*tree.number_of_nodes_lv(lv));
    }
    valueBuff.init(nn_max);
}

void Field::
init_valuePBVR(ValuePBVR& valuePBVR) {
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    valuePBVR.init(this);
}

void Field::
init_ValueTimeAverage(const Tree& tree, const MeshValue* meshValues) {
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    for(int lv=0; lv<DefAMR::LV_MAX; lv++) {
        const std::intptr_t nn_max = DefAMR::NN_LEAF * tree.number_of_nodes_lv(lv);
        // 1min
        valueTimeAverage1min_[lv].init(1.0, nn_max);
        valueTimeAverage1min_[lv].copy_from_ValueNS(nn_max, meshValues[lv].valueNS());
        valueTimeAverage1min_io_[lv].init(1.0, nn_max);

//        // 15min
//        valueTimeAverage15min_[lv].init(15.0, nn_max);
//        valueTimeAverage15min_[lv].copy_from_ValueNS(nn_max, meshValues[lv].valueNS());
    }
}


// read //
int  Field::
read_latest_step(const int step)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    std::string  filename = Foldernames::io_folder + "/" + Filenames::info_name0 + ".dat";

    int  tmp_step;

    std::ifstream  fout;
    fout.open(filename);
    fout >> tmp_step;
    fout.close();

    return  std::min(step, tmp_step);
}


void  Field::
read_parameters(const int step, Parameters&  parameters)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    // read //
    const std::string  filename = parameters.filename(step);

    parameters.readParameters(optionParser(), filename);
}


void  Field::
read_grid(const int step, Grid*  grid, const int rank)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each(
            [rank, step, &grid](const int i)
            {
                const std::string  filename = grid[i].filename(i, rank, step);

                grid[i].readGrid(filename);
            }
            );
}


void  Field::
read_tree(const int step, Tree&  tree, const int rank)
{
    const std::string  filename = tree.filename(rank, step);

    std::ifstream  fin;
    fin.open(filename);

    tree.read_stNodeValArray(fin, filename);
    tree.read_node_structures(fin, filename);
    tree.read_MPIPutGetInfo(fin, filename);

    fin.close();

    // init //
    tree.createMPIPutInformations();
}


void  Field::
read_taskID(const int step, TaskID& taskID, const int rank)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const std::string  filename = taskID.filename(rank, step);

    taskID.readTaskID(filename);
}


void  Field::
read_elapsedTimeInfo(const Parameters& parameter, const int rank)
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    allTimeInfo_ .ReadElapsedTimeInfo("all",  parameters_);
    funcTimeInfo_.ReadElapsedTimeInfo("func", parameters_);
    mpiTimeInfo_ .ReadElapsedTimeInfo("mpi",  parameters_);
}


void  Field::
read_meshValues(
    const int           step,
    const Grid*         grids,
    const Tree&         tree,
          MeshValue*    meshValues,
          ValueTimeAverage* valueTimeAverage,
    const int           rank__unused__
    )
const
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    auto  nn_max = [&tree](const int lv){ return  DefAMR::NN_LEAF*tree.number_of_nodes_lv(lv); };
    const int n_scalars = optionParser_.number_of_scalars();

    // allocate //
    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each( [nn_max, n_scalars, &meshValues, &tree](const int lv){
        meshValues[lv].init(nn_max(lv), n_scalars);
    } );

    // read //
    IOData  ioData(comm_, step);

//    ioData.readHDF5Values(grids, tree, meshValues);
    ioData.readBinaries(grids, tree, meshValues, nullptr, valueTimeAverage);
}


// write //
void  Field::
write_latest_step(const int step, const int rank)
const
{
#ifdef NO_FIELD_WRITEDAT
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT" << step << std::endl; }
#else
    if ( rank != 0 ) { return; }

    std::string  filename = Foldernames::io_folder + "/" + Filenames::info_name0 + ".dat";

    std::ofstream  fout;
    fout.open(filename);
    fout << step;
    fout.close();
#endif
}


void  Field::
write_parameters (const int step, const Parameters&  parameters, const int rank)
const
{
#if defined(NO_FIELD_WRITEDAT) || defined(NO_IODATA_RESTART)
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT/NO_IODATA_RESTART" << step << std::endl; }
#else
    if (rank != 0) { return; }

    const std::string  filename = parameters.filename(step);

    parameters.writeParameters(filename, step);
#endif
}


void  Field::
write_grid(const int step, const Grid*  grid, const int rank)
const
{
#if defined(NO_FIELD_WRITEDAT) || defined(NO_IODATA_RESTART)
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT/NO_IODATA_RESTART" << step << std::endl; }
#else
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each(
            [rank, step, &grid](const int i)
            {
                const std::string  filename = grid[i].filename(i, rank, step);

                grid[i].writeGrid(filename);
            }
            );
#endif
}


void  Field::
write_tree(const int step, const Tree&  tree, const int rank)
const
{
#if defined(NO_FIELD_WRITEDAT) || defined(NO_IODATA_RESTART)
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT/NO_IODATA_RESTART" << step << std::endl; }
#else
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const std::string  filename = tree.filename(rank, step);

    std::ofstream  fout;
    fout.open(filename);

    tree.write_stNodeValArray(fout, filename);
    tree.write_node_structures(fout, filename, tree.make_vectorIONode());
    tree.write_MPIPutGetInfo(fout, filename);

    fout.close();
#endif
}


void  Field::
write_taskID(const int step, const TaskID& taskID, const int rank)
const
{
#if defined(NO_FIELD_WRITEDAT) || defined(NO_IODATA_RESTART)
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT/NO_IODATA_RESTART" << step << std::endl; }
#else
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const std::string  filename = taskID.filename(rank, step);

    taskID.writeTaskID(filename);
#endif
}


void  Field::
write_meshValues(
    const int           step,
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
    const MeshValue*    meshValues,
    const int           rank__unused__
    )
const
{
#ifdef NO_FIELD_WRITEDAT
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT" << step << std::endl; }
#else
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    IOData  ioData(comm_, step);

    // vtu //
    #ifdef PARAVIEW_FULL
    ioData.writeVTKFile   (grids, tree, parameters, meshValues, valueStat_io_, VTKOutputScale::Full);
    #endif
    #ifdef PARAVIEW_DOWNSIZE_2
    ioData.writeVTKFile   (grids, tree, parameters, meshValues, valueStat_io_, VTKOutputScale::Downsize2);
    #endif
    #ifdef PARAVIEW_DOWNSIZE_4
    ioData.writeVTKFile   (grids, tree, parameters, meshValues, valueStat_io_, VTKOutputScale::Downsize4);
    #endif

    #ifdef PARAVIEW_SLICE
    ioData.writeVTKFile   (grids, tree, parameters, meshValues, valueStat_io_, VTKOutputScale::Slice_Full, 0,
        std::vector<real>{2., 50., 100., 1000.}, //z (alititude meter)
        std::vector<real>{0.}, //y (center)
        std::vector<real>{0.}  //x (center)
    );
    #endif

//    // hdf5 //
//    ioData.writeHDF5Values(grids, tree, meshValues);
    #ifdef REGRESSION_TEST // write anyway
    ioData.writeBinaries(grids, tree, meshValues, nullptr, valueTimeAverage1min_io_);
    #else
    if(step % 2 == 0) { // coherently safe on double buffering
        ioData.writeBinaries(grids, tree, meshValues, nullptr, valueTimeAverage1min_io_);
    }
    #endif

#endif // NO_FIELD_WRITEDAT
}

void Field::
write_ValueStatCsv(int step) {
#ifdef NO_FIELD_WRITEDAT
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skip @ NO_FIELD_WRITEDAT" << step << std::endl; }
#else
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    // csv //
    IOData(comm_, step).writeCsv_valueStat_integral_dtdxdy(tree_io_, parameters_io_, meshValues_io_, valueStat_io_);
#endif
}

