#pragma once
#ifndef NODECALFLAGS_H_
#define NODECALFLAGS_H_


class  NodeCalFlags {
private:
    bool  Cal_,
          L2C_, C2L_,
          L2F_, F2L_,
          putMPI_, getMPI_;

public:
    NodeCalFlags ()
    {
        reset();
    }

    ~NodeCalFlags () {}

    void  reset()
    {
        Cal_    = false;
        L2C_    = false;
        C2L_    = false;
        L2F_    = false;
        F2L_    = false;
        putMPI_ = false;
        getMPI_ = false;
    }

public:
    bool  Cal() const { return  Cal_; }

    bool  L2C() const { return  L2C_; }
    bool  C2L() const { return  C2L_; }

    bool  L2F() const { return  L2F_; }
    bool  F2L() const { return  F2L_; }

    bool  putMPI()    const { return  putMPI_; }
    bool  getMPI()    const { return  getMPI_; }

public:
    void  set_flag_Cal    (const bool  flag)    { Cal_ = flag; }
    void  set_flag_L2C    (const bool  flag)    { L2C_ = flag; }
    void  set_flag_C2L    (const bool  flag)    { C2L_ = flag; }
    void  set_flag_L2F    (const bool  flag)    { L2F_ = flag; }
    void  set_flag_F2L    (const bool  flag)    { F2L_ = flag; }
    void  set_flag_putMPI (const bool  flag)    { putMPI_    = flag; }
    void  set_flag_getMPI (const bool  flag)    { getMPI_    = flag; }

private:

};


#endif
