#include "Field.h"
#include "FuncLoop.h"
#include "FuncMath.h"
#include "IOData.h"
#include <cmath>


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


void  Field::
preset_field()
{
    if (rank_ == 0) { std::cout << "start :  preset_tree\n"; }

    // grid information //
    init_grid(grids_);

    // connenctions //
    preset_tree(tree0_);

    // parameters //
    init_parameters(parameters_);

    // taskID //
    init_taskID(taskID_, tree0_);

    // elapsedTimeInfo //
    init_elapsedTimeInfo(parameters_);

    // values //
    init_meshValues(grids_, tree0_, meshValues0_);

    // copy //
    const bool  is_reset = true;
    copy_tree(tree1_, tree0_, is_reset);


//    init_meshValues(grids_, tree0_, meshValues1_);
    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);

    // postprocessmonitor_ //
    postprocessmonitor_.setupdata(rank_, tree0_, meshValues0_);

    if (rank_ == 0) { std::cout << "end : preset_tree\n"; }
}


void  Field::
init_field()
{
    // grid information //
    init_grid(grids_);

    // connenctions //
    init_tree(tree0_);

    // parameters //
    init_parameters(parameters_);

    // taskID //
    init_taskID(taskID_, tree0_);

    // elapsedTimeInfo //
    init_elapsedTimeInfo(parameters_);

    // values //
    init_meshValues(grids_, tree0_, meshValues0_);

    // copy //
    const bool  is_reset = true;
    copy_tree(tree1_, tree0_, is_reset);

    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);
}


void  Field::
read_field(const int step)
{
    const int rstep = read_latest_step(step);

    if (rank_ == 0) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "rstep = " << rstep << std::endl;
    }

    // parameter //
    read_parameters(rstep, parameters_);

    // grid information //
    read_grid(rstep, grids_, rank_);

    // connenctions //
    read_tree(rstep, tree0_, rank_);

    // taskID //
    read_taskID(rstep, taskID_, rank_);

    // elapsedTimeInfo //
    read_elapsedTimeInfo(parameters_, rank_);


    // values //
    read_meshValues(rstep, grids_, tree0_, meshValues0_, rank_);

    // copy //
    const bool  is_reset = true;
    copy_tree(tree1_, tree0_, is_reset);

    copy_meshValues(meshValues1_, meshValues0_, grids_, tree0_, is_reset);
}


void  Field::
write_field(const int step)
const
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  wstep = step;

    write_latest_step(wstep, rank_);

    // parameter //
    write_parameters(wstep, parameters_io_, rank_);

    // grid information //
    write_grid(wstep, grids_io_, rank_);

    // connenctions //
    write_tree(wstep, tree_io_, rank_);

    // taskID //
    write_taskID(wstep, taskID_io_, rank_);

    // values //
    write_meshValues(wstep, grids_io_, tree0_, parameters_io_, meshValues_io_, rank_);
}


void  Field::
write_monitor(const int step)
const
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  wstep = step;

    // postprocessmonitor_ //
    postprocessmonitor_.OutputMonitorData(wstep, tree0_, parameters_io_, meshValues_io_, rank_);
//    postprocessmonitor_.OutputMonitorData(wstep, tree0_, parameters_, meshValues0_, rank_);
}


void  Field::
copy_io_field()
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    // parameters //
    parameters_io_.copy(parameters_);

    // grid information //
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        grids_io_[i].copy(grids_[i]);
    }

    // connenctions //
    copy_tree(tree_io_, tree0_, true);

    // taskID //
    taskID_io_.copyTaskID(taskID_);

    // values //
    copy_meshValues(meshValues_io_, meshValues0_, grids_io_, tree_io_, true);
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
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//    std::cout << "----------" << std::endl;
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    tree_new.copy_stNodeValArray (tree.number_of_nodes(), tree.stnodeValArray());
    tree_new.copy_node_structures(tree, is_reset);
    tree_new.copyMPIPutGetInfo(tree.mpiPutGetInfo());

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
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
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


void  Field::
update_meshValues(
    const Grid*         grids,
    const Tree&         tree
    )
{
//    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const bool is_reset = false;

    for (int i=0; i<DefAMR::LV_MAX; i++) {
        meshValue(i).copy_MeshValue(
            meshValue_new(i),
            grids[i],
            tree,
            is_reset
            );
    }

//    std::cout << "----------" << std::endl;
}


// private //
void  Field::
preset_tree(Tree&   tree)
{
    tree.preset_tree_data(grids());
}


// init //
void  Field::
init_grid(Grid*  grids)
{
    // optionParser //
    const int  nx_tmp  = optionParser_->number_of_grid_point(0);
    const int  dim_dir = optionParser_->number_of_grid_point(1);
    if ( !(dim_dir == 0 || dim_dir == 1 || dim_dir == 2) ) {
        std::cout << "error : " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    // arrange domain length //
    const real _domain_length_tmp[] = { static_cast<real>( optionParser_->domain_length(0) ),
                                        static_cast<real>( optionParser_->domain_length(1) ),
                                        static_cast<real>( optionParser_->domain_length(2) ) };

    const real dx_leaf_lv0 = _domain_length_tmp[dim_dir] / nx_tmp;

    const Vector3d<int>  nx0( static_cast<int>( ceil( _domain_length_tmp[0]/dx_leaf_lv0 ) ),
                              static_cast<int>( ceil( _domain_length_tmp[1]/dx_leaf_lv0 ) ),
                              static_cast<int>( ceil( _domain_length_tmp[2]/dx_leaf_lv0 ) )  );

    // domain //
    const Vector3d<real> offset( optionParser_->domain_min(0),
                                 optionParser_->domain_min(1),
                                 optionParser_->domain_min(2) );

    const real  dx0 = dx_leaf_lv0;

    // update //
    auto  init_grid_i = [nx0, offset, dx0, grids](const int i) {
                                const int  ratio = pow(2, i);
                                const Vector3d<int>  nx = nx0*ratio;
                                const real           dx = dx0/ratio;

                                grids[i].init(nx, offset, dx);
                            };

    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each( init_grid_i );
}


void  Field::
init_tree(Tree&   tree)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    tree.init_tree_data(grids());
}


void  Field::
init_elapsedTimeInfo(const Parameters& parameter)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    allTimeInfo_ .InitElapsedTimeInfo("all",  parameters_);
    funcTimeInfo_.InitElapsedTimeInfo("func", parameters_);
    mpiTimeInfo_ .InitElapsedTimeInfo("mpi",  parameters_);
}


void  Field::
init_parameters(Parameters& parameters)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    parameters.init(optionParser_);
}


void  Field::
init_taskID(TaskID& taskID, const Tree& tree)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    taskID.set_taskID(tree);
}


void  Field::
init_taskID_opt(TaskID& taskID, const Tree& tree, const MeshValue* meshValues)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    taskID.set_taskID_opt(tree, meshValues);
}


void  Field::
init_meshValues(const Grid* grids, const Tree& tree, MeshValue*  meshValues)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each( [&grids, &tree, &meshValues](const int lv){
        const int  nn_max = DefAMR::NN_LEAF*tree.number_of_nodes_lv(lv);
        meshValues[lv].init( nn_max );
        meshValues[lv].set_coordinate( lv, grids[lv], tree );
        } );

}


// read //
int  Field::
read_latest_step(const int step)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

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
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    // read //
    const std::string  filename = parameters.filename(step);

    parameters.readParameters(optionParser(), filename);
}


void  Field::
read_grid(const int step, Grid*  grid, const int rank)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

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
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const std::string  filename = taskID.filename(rank, step);

    taskID.readTaskID(filename);
}


void  Field::
read_elapsedTimeInfo(const Parameters& parameter, const int rank)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

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
    const int           rank
    )
const
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    auto  nn_max = [&tree](const int i){ return  DefAMR::NN_LEAF*tree.number_of_nodes_lv(i); };

    // allocate //
    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each( [nn_max, &meshValues](const int i){ meshValues[i].init(nn_max(i)); } );

    // read //
    IOData  ioData(rank, step);

//    ioData.readHDF5Values(grids, tree, meshValues);
}


// write //
void  Field::
write_latest_step(const int step, const int rank)
const
{
    if ( rank != 0 ) { return; }

    std::string  filename = Foldernames::io_folder + "/" + Filenames::info_name0 + ".dat";

    std::ofstream  fout;
    fout.open(filename);
    fout << step;
    fout.close();
}


void  Field::
write_parameters (const int step, const Parameters&  parameters, const int rank)
const
{
    if (rank != 0) { return; }

    const std::string  filename = parameters.filename(step);

    parameters.writeParameters(filename, step);
}


void  Field::
write_grid(const int step, const Grid*  grid, const int rank)
const
{
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    FuncLoop::Loop1d  loop1d(0, DefAMR::LV_MAX);
    loop1d.for_each(
            [rank, step, &grid](const int i)
            {
                const std::string  filename = grid[i].filename(i, rank, step);

                grid[i].writeGrid(filename);
            }
            );
}


void  Field::
write_tree(const int step, const Tree&  tree, const int rank)
const
{
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const std::string  filename = tree.filename(rank, step);

    std::ofstream  fout;
    fout.open(filename);

    tree.write_stNodeValArray(fout, filename);
    tree.write_node_structures(fout, filename, tree.make_vectorIONode());
    tree.write_MPIPutGetInfo(fout, filename);

    fout.close();
}


void  Field::
write_taskID(const int step, const TaskID& taskID, const int rank)
const
{
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const std::string  filename = taskID.filename(rank, step);

    taskID.writeTaskID(filename);
}


void  Field::
write_meshValues(
    const int           step,
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
    const MeshValue*    meshValues,
    const int           rank
    )
const
{
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    IOData  ioData(rank, step);

    // vtu //
//    ioData.writeVTKFile   (grids, tree, parameters, meshValues, VTKOutputScale::Full);
    ioData.writeVTKFile   (grids, tree, parameters, meshValues, VTKOutputScale::Downsize2);
    ioData.writeVTKFile   (grids, tree, parameters, meshValues, VTKOutputScale::Downsize4);

//    // hdf5 //
//    ioData.writeHDF5Values(grids, tree, meshValues);

}
