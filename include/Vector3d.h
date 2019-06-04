#pragma once
#ifndef VECTOR3D_H__
#define VECTOR3D_H__


#include <iostream>
#include <cmath>
#include <string>


template <typename TYPE>
class Vector3d {
private:
    union {
        TYPE  xi_[3];
        struct {
            TYPE  x_, y_, z_;
        };
    };

public:
    TYPE x() const { return x_; }
    TYPE y() const { return y_; }
    TYPE z() const { return z_; }

    void set_xyz(const TYPE x, const TYPE y, const TYPE z) { set_x(x); set_y(y); set_z(z); }

    void set_x(const TYPE x) { x_ = x; }
    void set_y(const TYPE y) { y_ = y; }
    void set_z(const TYPE z) { z_ = z; }

    const TYPE xi(const int i) const {
        if (i == 0 || i == 1 || i == 2) { return xi_[i]; }
        else                            { std::cout << __PRETTY_FUNCTION__ << " : error xi"; exit(-1); }
    }


public:
    // Numerical constructor
    Vector3d(const TYPE &x, const TYPE &y, const TYPE &z): x_(x), y_(y), z_(z) {}; // For 3-rotations
    Vector3d(): x_(TYPE()), y_(TYPE()), z_(TYPE()) {};


    // Copy constructor and assignment
    Vector3d(const Vector3d &q): x_(q.x_), y_(q.y_), z_(q.z_) {};
    Vector3d& operator=(const Vector3d &q) { x_=q.x_; y_=q.y_; z_=q.z_; return *this; }

    // Unary operators
    Vector3d operator-() const { return Vector3d(-x_, -y_, -z_); }

    // In-place operators
    Vector3d& operator+=(const TYPE &r) { x_ += r; y_ += r; z_ += r; return *this; }
    Vector3d& operator-=(const TYPE &r) { x_ -= r; y_ -= r; z_ -= r; return *this; }

    Vector3d& operator+=(const Vector3d &q) { x_ += q.x_; y_ += q.y_; z_ += q.z_; return *this; }
    Vector3d& operator-=(const Vector3d &q) { x_ -= q.x_; y_ -= q.y_; z_ -= q.z_; return *this; }

    Vector3d& operator*=(const TYPE &r) { x_ *= r; y_ *= r; z_ *= r; return *this; }
    Vector3d& operator/=(const TYPE &r) { x_ /= r; y_ /= r; z_ /= r; return *this; }

    Vector3d& operator*=(const Vector3d &q) { x_ *= q.x_; y_ *= q.y_; z_ *= q.z_; return *this; }
    Vector3d& operator/=(const Vector3d &q) { x_ /= q.x_; y_ /= q.y_; z_ /= q.z_; return *this; }

    // Binary operators based on in-place operators
    Vector3d operator+(const TYPE &r)     const { return Vector3d(*this) += r; }
    Vector3d operator-(const TYPE &r)     const { return Vector3d(*this) -= r; }

    Vector3d operator+(const Vector3d &q) const { return Vector3d(*this) += q; }
    Vector3d operator-(const Vector3d &q) const { return Vector3d(*this) -= q; }

    Vector3d operator*(const TYPE &r)     const { return Vector3d(*this) *= r; }
    Vector3d operator/(const TYPE &r)     const { return Vector3d(*this) /= r; }

    Vector3d operator*(const Vector3d &q) const { return Vector3d(*this) *= q; }
    Vector3d operator/(const Vector3d &q) const { return Vector3d(*this) /= q; }

    // Comparison operators, as much as they make sense
    bool operator==(const Vector3d &q) const { return (x_ == q.x_) && (y_ == q.y_) && (z_ == q.z_); }
    bool operator!=(const Vector3d &q) const { return !operator==(q); }

    // cast ??? //
    template<typename T>
    const Vector3d<T>  castVector3d() const { return Vector3d<T>{ (T)x_, (T)y_, (T)z_ }; }

    // function //
    TYPE sum    ()      const { return x_ + y_ + z_; }
    TYPE product()      const { return x_ * y_ * z_; }

    TYPE normSquared()  const { return x_*x_ + y_*y_ + z_*z_; }

public:
    void  cout(const std::string  strings) const { std::cout << strings << " : " << x() << ", " << y() << ", " << z() << std::endl; }

};


#endif
