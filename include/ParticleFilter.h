#ifndef PARTICLEFILTER_H__
#define PARTICLEFILTER_H__

#include <iostream>
#include <valarray>
#include <numeric>
#include <fstream>
#include <string>
#include "defineFilenames.h"
#include "mt19937r.hpp"
#include "MPICommEnsemble.h"
#include "mpi_wrapper.hpp"


class ParticleFilter {
private:
    using UniformDist = std::uniform_real_distribution<double>;
    using NormalDist  = std::normal_distribution<double>;
    using Engine = util::mt19937r;
    const MPICommEnsemble comm_;

    int n_particles_;
    std::valarray<double> err_; // y - fx
    std::valarray<double> x_; // require restart !!!!!
    std::valarray<double> xre_;
    std::valarray<double> w_; // weight

    double sigmax_;
    double sigmay_;

    double x_min_;
    double x_max_;

    Engine engine_{0};

public:
    ParticleFilter(const MPICommEnsemble comm): comm_(comm){};
    ~ParticleFilter(){}

public:
    const std::valarray<double> x()   const { return x_; }
    const std::valarray<double> xre() const { return xre_; }

public:
    void  init(int n_particles, double xre, double sigmax, double sigmay, double x_min, double x_max)
    {
        n_particles_ = n_particles;
        sigmax_ = sigmax;
        sigmay_ = sigmay;

        x_min_ = x_min;
        x_max_ = x_max;

        err_ = std::valarray<double>(0.0, n_particles_);
        x_   = std::valarray<double>(0.0, n_particles_);
        xre_ = std::valarray<double>(xre, n_particles_); // initial value
        w_   = std::valarray<double>(0.0, n_particles_);

        auto&& v = uni_valarray(n_particles_, sigmax_);
        x_   = xre_ + v;
    }

    void  update(const double* err, int n)
    {
        std::valarray<double> _err = std::valarray<double>(err, n);
        update(_err);
    }

    void  update(const std::valarray<double> err)
    {
        if(err.size() != err_.size()) { std::cout << __PRETTY_FUNCTION__ << ", input error : err.size of input is not equal to err_.size of class memeber!" << std::endl; }
        if(comm_.is_rank0()) {std::cout << __PRETTY_FUNCTION__ << ": update "<< std::endl;}

        // update observation model by x_
        // set y<-obs, and fx<-x 
        // err = y - fx
        err_ = err;

        w_ = norm_likelihood(err_, sigmay_*sigmay_) + 1.0e-7;
        auto&& wsum = w_.sum();
        w_ = w_/wsum;

        #if 1 // particle filter

        auto&& k_list = resampling(w_);
        xre_ = copy_by_list(x_, k_list);
        std::sort(std::begin(xre_), std::end(xre_));

        #else // merging particle filter (mpf)

        std::valarray<double> xre_av = std::valarray<double>(0.0, n_particles_);
        const std::valarray<double> weights = { 3.0/4.0, (sqrt(13.0) + 1.0)/8.0, -(sqrt(13.0) - 1.0)/8.0 };
//        const std::valarray<double> weights = { 19.0/20.0, (sqrt(77.0) + 1.0)/40.0, -(sqrt(77.0) - 1.0)/40.0 };
        for (auto&& weight : weights) {
            std::valarray<double> xre_tmp = std::valarray<double>(0.0, n_particles_);

            auto&& k_list = resampling(w_);
            xre_tmp = copy_by_list(x_, k_list);
            std::sort(std::begin(xre_tmp), std::end(xre_tmp));

            xre_av = xre_av + weight*xre_tmp;
        }
        xre_ = xre_av;
        std::sort(std::begin(xre_), std::end(xre_));

        #endif

        auto&& v = rand_valarray(n_particles_, NormalDist(0.0, sigmax_), engine_);
        const auto vlim = 2.0*sigmax_;
        limit_min_max(v, -vlim, vlim);

        x_ = xre_ + v;
        std::sort(std::begin(x_), std::end(x_));

        for(int i=0; i<n_particles_; i++) {
            if(x_[i] <= x_min_) { x_[i] = x_min_; }
            if(x_[i] >= x_max_) { x_[i] = x_max_; }
        }
    }

    void output_pf(int t) const
    {
        output(t, "x",   x_);
        output(t, "xre", xre_);
    }

    void save_ParticleFiter(const int t) const
    {
        save_x(t);
        save_engine(t);
    }

    void load_ParticleFiter(const int t) 
    {
        load_x(t);
        load_engine(t);
    }

private:
    std::string filename_x(const int t) const
    {
        const std::string prefix = Foldernames::io_folder + "/restart/" + std::to_string(t);
        const std::string suffix = "particle_x" + std::to_string(comm_.world().rank()) + ".dat";
        return prefix + "/" + suffix;
    }

    std::string filename_engine(const int t) const
    {
        const std::string prefix = Foldernames::io_folder + "/restart/" + std::to_string(t);
        const std::string suffix = "rand_rank" + std::to_string(comm_.world().rank()) + ".dat";
        return prefix + "/" + suffix;
    }

    void save_x(const int t) const
    {
        if(comm_.is_rank0()) {std::cout << __PRETTY_FUNCTION__ << ":  save time =" << t << std::endl;}    

        FILE* fp = fopen(filename_x(t).c_str(), "wb");
        if(comm_.is_rank0()) {std::cout << __PRETTY_FUNCTION__ << ":  filename =" << filename_x(t).c_str() << std::endl;}    
        runtime_assert(fp != NULL, "IOError");
        size_t size = x_.size();
        fwrite(&size, sizeof(size_t), 1, fp);
        for(const double& v: x_) {
            fwrite(&v, sizeof(double), 1, fp);
        }
        fclose(fp);
    }

    void load_x(const int t)
    {
        if(comm_.is_rank0()) {std::cout << __PRETTY_FUNCTION__ << ":  load time =" << t << std::endl;}    

        FILE* fp = fopen(filename_x(t).c_str(), "rb");
        runtime_assert(fp != NULL, "IOError");
        size_t size;
        fread(&size, sizeof(size_t), 1, fp);
        x_.resize(size);
        for(double& v: x_) {
            fread(&v, sizeof(double), 1, fp);
        }
        fclose(fp);
    }

    void save_engine(const int t) const
    {
        if(comm_.is_rank0()) {std::cout << __PRETTY_FUNCTION__ << ":  save time =" << t << std::endl;}    

        engine_.save(filename_engine(t));
    }

    void load_engine(const int t) 
    {
        if(comm_.is_rank0()) {std::cout << __PRETTY_FUNCTION__ << ":  restart time =" << t << std::endl;}    

        engine_.load(filename_engine(t));
    }

    template<class DIST, class ENGINE>
    std::valarray<typename DIST::result_type>  rand_valarray(size_t size, DIST&& dist, ENGINE&& engine)
    {
        std::valarray<typename DIST::result_type> ret(size);
        for(auto& r: ret) {
            r = dist(engine);
        }
        return ret;
    }

    std::valarray<double> uni_valarray(size_t size, double sigma_x, double offset = 0.0)
    {
        std::valarray<double> uni(size);

        const double uni_dx = 2.0*sigma_x/n_particles_;
        for(int i=0; i<n_particles_; i++) {
            uni[i] = -sigma_x + (0.5 + i)*uni_dx + offset;
        }

        return uni;
    }

    template<typename T>
    void  limit_min_max(std::valarray<T>& val, T val_min, T val_max) const
    {
        for (int i=0; i<val.size(); i++) {
            const T val_tmp = val[i];

            val[i] =  (val_tmp < val_min) ? val_min
                    : (val_tmp > val_max) ? val_max
                    :                       val_tmp;
        }
    }

    template<typename T>
    std::valarray<T>  cumsum(const std::valarray<T>& val) const
    {
        std::valarray<T> cum(val.size());
        std::partial_sum(std::begin(val), std::end(val), std::begin(cum));

        return cum;
    }

    std::valarray<double>  norm_likelihood(const std::valarray<double>& err, double s2) const // s2 : sigma^2
    {
        return std::exp( -(err*err)/(2.0*s2) ) / sqrt(2.0*M_PI*s2);
    }

    int  F_inv(const std::valarray<double>& w_cumsum, const std::valarray<int>& idx, double u)
    {
//        if(std::valarray<bool> b = (w_cumsum < u);  b.sum() == false) { // <-C++17;  python, if np.any(w_cumsum < u) == False: 
        std::valarray<bool> b = (w_cumsum < u);
        if(b.sum() == false) { // if np.any(w_cumsum < u) == False:
            return 0;
        }
        else {
            std::valarray<int> idx_k = idx[w_cumsum < u];
            int k = idx_k.max();
            return k+1;
        }
    }

    std::valarray<int>  resampling(const std::valarray<double>& weights)
    {
        auto&& w_cumsum = cumsum<double>(weights); // cumulative sum (accumulation sum)
        std::valarray<int> k_list(0, n_particles_); // sampled k lists
        std::valarray<int> idx   (0, n_particles_);
        for(int i=0; i<n_particles_; i++) { idx[i] = i; }

        #if 0 // uniform [0 ~ 1]
            std::valarray<double> uni(n_particles_);
            for(int i=0; i<n_particles_; i++) {
                const double uni_dx = 1.0/n_particles_;
                uni[i] = (0.5 + i)*uni_dx;
            }
        #else // rand [0 ~ 1]
            auto&& uni = rand_valarray(n_particles_, UniformDist(0., 1.), engine_);
        #endif
        std::sort(std::begin(uni), std::end(uni));
        for(int i=0; i<n_particles_; i++) {
            k_list[i] = F_inv(w_cumsum, idx, uni[i]);
        }

        // check out of range
        for(int i=0; i<n_particles_; i++) {
            if(k_list[i] < 0 || k_list[i] >= n_particles_) {
                std::cout << __PRETTY_FUNCTION__ << " : error k_list[" << i << "]" << std::endl;
            }
        }

        return k_list;
    }

    std::valarray<double>  copy_by_list(std::valarray<double>& x, std::valarray<int>& list) const
    {
        std::valarray<double> xn(x.size());
        for(int i=0; i<x.size(); i++) {
            xn[i] = x[list[i]];
        }
        return xn;
    }

    void  output(int t, const std::string ss, const std::valarray<double>& val) const
    {
        std::ofstream fout;
        std::string filename = "output/" + ss + "-" + std::to_string(t) + ".txt";
        fout.open(filename, std::ios::out);
        for(auto& v: val) {
            fout << v << std::endl;
        }
        fout.close();
    }

}; // class ParticleFilter


#endif // PARTICLEFILTER_H
