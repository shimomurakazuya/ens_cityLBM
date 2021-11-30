#ifndef CITYLBM_MPICOMMENSEMBLE_H_
#define CITYLBM_MPICOMMENSEMBLE_H_

#include <mpi.h>
#include "mpi_wrapper.hpp"

// definition:
// - the whole matrix of the ensemble simluation is defined as:
//   X = (x^(1), x^(2), ..., x^(Ne))
//   - each column vector (x^(k)) is the state vector of k-th ensemble member
//   - each row vector is the ensemble samples of i-th state variable

class MPICommEnsemble {
    // premitives
    MPI_Comm comm_col_vector_, comm_row_vector_;
    const int ofs_ensemble_idx_;

public:
    // getter
    util::mpi world     () const { return util::mpi(MPI_COMM_WORLD  ); }
    util::mpi col_vector() const { return util::mpi(comm_col_vector_); }
    util::mpi row_vector() const { return util::mpi(comm_row_vector_); }

    bool is_rank0() const { return 0 == this->world().rank(); }
    int col_id_wo_offset() const { return row_vector().rank(); }
    int col_id() const { return col_id_wo_offset() + ofs_ensemble_idx_; }
    int row_id() const { return col_vector().rank(); }

    // getter alias
    auto n_cols() const { return row_vector().size(); }
    auto n_rows() const { return col_vector().size(); }
    auto ensemble_size() const { return n_cols(); }
    auto ensemble_id() const { return col_id(); }

    // ctor
    MPICommEnsemble() = delete;
    MPICommEnsemble(const MPICommEnsemble&) = default;
    MPICommEnsemble(MPICommEnsemble&&) = default;
    MPICommEnsemble(const std::size_t ensemble_size, const int ofs_ensemble_idx): ofs_ensemble_idx_(ofs_ensemble_idx) {
        const auto rank = world().rank();
        const auto size = world().size();

        // mpisize --> n_rows (grid parallelization)
        runtime_assert(size % ensemble_size == 0, "InvalidConfiguration: mpi size vs. ensemble size mismatched");
        const auto n_rows = size / ensemble_size;

        // comm split: each ensemble member (column vector) //
        const auto k_col = rank / n_rows;
        MPI_Comm_split(MPI_COMM_WORLD, k_col, 0, &comm_col_vector_);

        // comm split: common grid lacation over ensmembers (row vector) //
        const auto i_row = rank % n_rows;
        MPI_Comm_split(MPI_COMM_WORLD, i_row, 0, &comm_row_vector_);
    }

    // debug
    void cout_info() const {
        world().for_each_rank([=]() {
            std::cout << "this rank is " 
                << world     ().rank() << "/" << world     ().size() << " (world) >>> " 
                << col_vector().rank() << "/" << col_vector().size() << " (@ens member idx k = " << col_id() << ") >>> " 
                << row_vector().rank() << "/" << row_vector().size() << " (@grid parallelization idx i = " << row_id() << ")"
                << std::endl;
        });
    }

};

#endif

