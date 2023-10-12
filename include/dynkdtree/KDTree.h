#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cwchar>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <variant>
#include <vector>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

template <int Dim, class Scalar>
using cref = const Eigen::Ref<const Eigen::Matrix<Scalar, Dim, 1>> &;
template <int Dim, class Scalar>
using ref = Eigen::Ref<Eigen::Matrix<Scalar, Dim, 1>>;

namespace dynotree {
template <typename T, typename Scalar>
void choose_split_dimension_default(const T &lb, const T &ub, int &ii,
                                    Scalar &width) {
  for (std::size_t i = 0; i < lb.size(); i++) {
    Scalar dWidth = ub[i] - lb[i];
    if (dWidth > width) {
      ii = i;
      width = dWidth;
    }
  }
}

// template <int Dim, typename Scalar>
// void sample_uniform_default(cref<Dim,Scalar> lb, cref<Dim,Scalar> ub,
// ref<Dim,Scalar> x) { x = Eigen::Matrix<Scalar, T::RowsAtCompileTime,
// 1>::Random(lb.size()); x.setRandom(); x = lb + (ub -
// lb).cwiseProduct(x.array() + 1.) / 2.;
// }

template <typename Scalar, int Dimensions = -1> struct RnL1 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;

  Eigen::Matrix<Scalar, Dimensions, 1> lb;
  Eigen::Matrix<Scalar, Dimensions, 1> ub;

  void set_bounds(cref_t lb_, cref_t ub_) {
    lb = lb_;
    ub = ub_;
  }

  inline void sample_uniform(ref_t x) const {
    x.setRandom();
    x.array() += 1;
    x /= .2;
    x = lb + (ub - lb).cwiseProduct(x);
  }

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    out = from + t * (to - from);
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  inline Scalar distance_to_rectangle(cref_t &location1, cref_t &lb,
                                      cref_t &ub) const {

    Scalar d = 0;
    Scalar dist = 0;

    if constexpr (Dimensions == Eigen::Dynamic) {

      assert(location1.size());
      assert(ub.size());
      assert(lb.size());
      assert(location1.size() == ub.size());
      assert(location1.size() == lb.size());

      for (size_t i = 0; i < location1.size(); i++) {
        Scalar xx = std::max(lb(i), std::min(ub(i), location1(i)));
        Scalar dif = xx - location1(i);
        dist += std::abs(dif);
      }
    } else {
      for (size_t i = 0; i < Dimensions; i++) {
        Scalar xx = std::max(lb(i), std::min(ub(i), location1(i)));
        Scalar dif = xx - location1(i);
        dist += std::abs(dif);
      }
    }
    return dist;
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    assert(location1.size());
    assert(location2.size());
    return (location1 - location2).cwiseAbs().sum();
  }
};

template <typename Scalar> struct SO2 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 1, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 1, 1>>;

  void sample_uniform(ref_t x) const {

    x(0) = (double(rand()) / (RAND_MAX + 1.)) * 2. * M_PI - M_PI;
  }

  void set_bounds(cref_t lb_, cref_t ub_) {
    std::stringstream ss;
    ss << "so2 has no bounds " << __FILE__ << ":" << __LINE__;
    throw std::runtime_error(ss.str());
  }

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {

    assert(t >= 0);
    assert(t <= 1);

    Eigen::Matrix<Scalar, 1, 1> d = to - from;

    if (d(0) > M_PI) {
      d(0) -= 2 * M_PI;
    } else if (d(0) < -M_PI) {
      d(0) += 2 * M_PI;
    }

    out = from + t * d;

    if (out(0) > M_PI) {
      out(0) -= 2 * M_PI;
    } else if (out(0) < -M_PI) {
      out(0) += 2 * M_PI;
    }
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  inline Scalar distance_to_rectangle(cref_t location1, cref_t lb,
                                      cref_t ub) const {

    assert(location1(0) >= -M_PI);
    assert(location1(0) <= M_PI);

    assert(lb(0) >= -M_PI);
    assert(lb(0) <= M_PI);

    assert(ub(0) >= -M_PI);
    assert(ub(0) <= M_PI);

    if (location1(0) >= lb(0) && location1(0) <= ub(0)) {
      return 0;
    } else if (location1(0) > ub(0)) {
      Scalar d1 = location1(0) - ub(0);
      Scalar d2 = lb(0) - (location1(0) - 2 * M_PI);
      assert(d2 >= 0);
      assert(d1 >= 0);
      return std::min(d1, d2);
    } else if (location1(0) < lb(0)) {
      Scalar d1 = lb(0) - location1(0);
      Scalar d2 = (location1(0) + 2 * M_PI) - ub(0);
      assert(d2 >= 0);
      assert(d1 >= 0);
      return std::min(d1, d2);
    } else {
      assert(false);
      return 0;
    }
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    assert(location1(0) >= -M_PI);
    assert(location2(0) >= -M_PI);

    assert(location1(0) <= M_PI);
    assert(location2(0) <= M_PI);

    Scalar dif = location1(0) - location2(0);
    if (dif > M_PI) {
      dif -= 2 * M_PI;
    } else if (dif < -M_PI) {
      dif += 2 * M_PI;
    }
    Scalar out = std::abs(dif);
    return out;
  }
};

template <typename Scalar> struct SO2Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 1, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 1, 1>>;

  Eigen::Matrix<Scalar, 1, 1> lb;
  Eigen::Matrix<Scalar, 1, 1> ub;

  SO2<Scalar> so2;

  inline void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    so2.interpolate(from, to, t, out);
  }

  inline void sample_uniform(ref_t x) const { so2.sample_uniform(x); }

  inline void set_bounds(cref_t lb_, cref_t ub_) {
    std::stringstream ss;
    ss << "so2 has no bounds " << __FILE__ << ":" << __LINE__;
    throw std::runtime_error(ss.str());
  }

  inline void choose_split_dimension(cref_t lb, cref_t ub, int &ii,
                                     Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  inline Scalar distance_to_rectangle(cref_t location1, cref_t lb,
                                      cref_t ub) const {

    Scalar d = so2.distance_to_rectangle(location1, lb, ub);
    return d * d;
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    Scalar d = so2.distance(location1, location2);
    return d * d;
  }
};

template <typename Scalar, int Dimensions = -1> struct RnSquared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;

  Eigen::Matrix<Scalar, Dimensions, 1> lb;
  Eigen::Matrix<Scalar, Dimensions, 1> ub;

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    out = from + t * (to - from);
  }

  void set_bounds(cref_t lb_, cref_t ub_) {
    lb = lb_;
    ub = ub_;
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii,
                              Scalar &width) const {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  void sample_uniform(ref_t x) const {

    x.setRandom();
    x.array() += 1;
    x /= .2;
    x = lb + (ub - lb).cwiseProduct(x);

    //
    // x.setRandom();
    // x.array() += 1;
    // x /= .2;
    // x = lb + (ub - lb).cwiseProduct(x);
    // x.setRandom();
    // x = lb + (ub - lb).cwiseProduct(x.array() + 1.) / 2.;
  }

  inline Scalar distance_to_rectangle(cref_t location1, cref_t lb,
                                      cref_t ub) const {

    Scalar d = 0;
    Scalar dist = 0;

    if constexpr (Dimensions == Eigen::Dynamic) {

      assert(location1.size());
      assert(ub.size());
      assert(lb.size());
      assert(location1.size() == ub.size());
      assert(location1.size() == lb.size());

      for (size_t i = 0; i < location1.size(); i++) {
        Scalar xx = std::max(lb(i), std::min(ub(i), location1(i)));
        Scalar dif = xx - location1(i);
        dist += dif * dif;
      }
    } else {
      for (size_t i = 0; i < Dimensions; i++) {
        Scalar xx = std::max(lb(i), std::min(ub(i), location1(i)));
        Scalar dif = xx - location1(i);
        dist += dif * dif;
      }
    }

    // Issue: memory allocation
    // Eigen::Matrix<Scalar, Dimensions, 1> bb = location1;
    // dist = (location1 - bb.cwiseMax(lb).cwiseMin(ub)).squaredNorm();

    return dist;
  }

  inline Scalar distance(cref_t &location1, cref_t &location2) const {

    assert(location1.size());
    assert(location2.size());

    // if constexpr (Dimensions == Eigen::Dynamic) {
    //   std::cout << "dynamic" << std::endl;
    // } else {
    //   std::cout << "" << std::endl;
    // }

    return (location1 - location2).squaredNorm();
  }
};

template <typename Scalar, int Dimensions = -1> struct Rn {

  using RnSquared = RnSquared<Scalar, Dimensions>;
  RnSquared rn_squared;

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;

  Eigen::Matrix<Scalar, Dimensions, 1> lb;
  Eigen::Matrix<Scalar, Dimensions, 1> ub;

  void set_bounds(cref_t lb_, cref_t ub_) {

    std::cout << "setting bounds " << std::endl;
    std::cout << lb_.transpose() << std::endl;
    std::cout << ub_.transpose() << std::endl;
    lb = lb_;
    ub = ub_;
  }

  inline void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    out = from + t * (to - from);
  }

  inline void sample_uniform(ref_t x) const {

    std::cout << "bounds are" << std::endl;
    std::cout << lb.transpose() << std::endl;
    std::cout << ub.transpose() << std::endl;

    x.setRandom();
    x.array() += 1;
    x /= 2.;
    x = lb + (ub - lb).cwiseProduct(x);

    // x.setRandom();
    // x = lb + (ub - lb).cwiseProduct(x.array() + 1.) / 2.;
  }

  inline void choose_split_dimension(cref_t lb, cref_t ub, int &ii,
                                     Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  inline Scalar distance_to_rectangle(cref_t &location1, cref_t &lb,
                                      cref_t &ub) const {

    Scalar d = rn_squared.distance_to_rectangle(location1, lb, ub);
    return std::sqrt(d);
  };

  inline Scalar distance(cref_t &location1, cref_t &location2) const {
    Scalar d = rn_squared.distance(location1, location2);
    return std::sqrt(d);
  }
};

template <typename Scalar> struct R2SO2 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 3, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 3, 1>>;

  using cref2_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 2, 1>> &;
  using ref2_t = Eigen::Ref<Eigen::Matrix<Scalar, 2, 1>>;

  // Eigen::Matrix<Scalar, 3, 1> lb;
  // Eigen::Matrix<Scalar, 3, 1> ub;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar angular_weight = 1.0;

  Rn<Scalar, 2> l2;
  SO2<Scalar> so2;

  void set_bounds(cref2_t lb_, cref2_t ub_) { l2.set_bounds(lb_, ub_); }

  inline void sample_uniform(ref_t x) const {
    l2.sample_uniform(x.template head<2>());
    so2.sample_uniform(x.template tail<1>());
  }

  inline void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    l2.interpolate(from.template head<2>(), to.template head<2>(), t,
                   out.template head<2>());
    so2.interpolate(from.template tail<1>(), to.template tail<1>(), t,
                    out.template tail<1>());
  }

  inline Scalar distance_to_rectangle(cref_t location1, cref_t lb,
                                      cref_t ub) const {

    Scalar d1 =
        l2.distance_to_rectangle(location1.template head<2>(),
                                 lb.template head<2>(), ub.template head<2>());
    Scalar d2 =
        so2.distance_to_rectangle(location1.template tail<1>(),
                                  lb.template tail<1>(), ub.template tail<1>());
    return d1 + angular_weight * d2;
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    Scalar d1 =
        l2.distance(location1.template head<2>(), location2.template head<2>());
    Scalar d2 = so2.distance(location1.template tail<1>(),
                             location2.template tail<1>());
    return d1 + angular_weight * d2;
  };
};

template <typename Scalar> struct R2SO2Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 3, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 3, 1>>;

  using cref2_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 2, 1>> &;
  using ref2_t = Eigen::Ref<Eigen::Matrix<Scalar, 2, 1>>;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar angular_weight = 1.0;

  RnSquared<Scalar, 2> rn_squared;
  SO2Squared<Scalar> so2squared;

  void set_bounds(cref2_t lb_, cref2_t ub_) { rn_squared.set_bounds(lb_, ub_); }

  void sample_uniform(ref_t x) const {
    rn_squared.sample_uniform(x.template head<2>());
    so2squared.sample_uniform(x.template tail<1>());
  }

  inline Scalar distance_to_rectangle(cref_t location1, cref_t lb,
                                      cref_t ub) const {

    Scalar d1 = rn_squared.distance_to_rectangle(location1.template head<2>(),
                                                 lb.template head<2>(),
                                                 ub.template head<2>());
    Scalar d2 = so2squared.distance_to_rectangle(location1.template tail<1>(),
                                                 lb.template tail<1>(),
                                                 ub.template tail<1>());
    return d1 + angular_weight * d2;
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    Scalar d1 = rn_squared.distance(location1.template head<2>(),
                                    location2.template head<2>());
    Scalar d2 = so2squared.distance(location1.template tail<1>(),
                                    location2.template tail<1>());
    return d1 + angular_weight * d2;
  };
};

template <typename Scalar> struct SO3Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 4, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 4, 1>>;

  RnSquared<Scalar, 4> rn_squared;

  void sample_uniform(ref_t x) const {
    x = Eigen::Quaterniond().UnitRandom().coeffs();
  }

  void set_bounds(cref_t lb_, cref_t ub_) {
    std::stringstream ss;
    ss << "so3 has no bounds " << __FILE__ << ":" << __LINE__;
    throw std::runtime_error(ss.str());
  }

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    throw std::runtime_error("not implemented interpolate in so3");
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  inline Scalar distance_to_rectangle(cref_t &location1, cref_t &lb,
                                      cref_t &ub) const {

    assert(std::abs(location1.norm() - 1) < 1e-6);

    Scalar d1 = rn_squared.distance_to_rectangle(location1, lb, ub);
    Scalar d2 = rn_squared.distance_to_rectangle(-1. * location1, lb, ub);
    return std::min(d1, d2);
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    assert(location1.size() == 4);
    assert(location2.size() == 4);

    assert(std::abs(location1.norm() - 1) < 1e-6);
    assert(std::abs(location2.norm() - 1) < 1e-6);

    Scalar d1 = rn_squared.distance(location1, location2);
    Scalar d2 = rn_squared.distance(-location1, location2);
    return std::min(d1, d2);
  };
};

template <typename Scalar> struct SO3 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 4, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 4, 1>>;

  SO3Squared<Scalar> so3squared;

  void sample_uniform(ref_t x) const { so3squared.sample_uniform(x); }

  void set_bounds(cref_t lb_, cref_t ub_) {
    std::stringstream ss;
    ss << "so3 has no bounds " << __FILE__ << ":" << __LINE__;
    throw std::runtime_error(ss.str());
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    throw std::runtime_error("not implemented interpolate in so3");
  }

  inline Scalar distance_to_rectangle(cref_t &location1, cref_t &lb,
                                      cref_t &ub) const {

    return std::sqrt(so3squared.distance_to_rectangle(location1, lb, ub));
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    return std::sqrt(so3squared.distance(location1, location2));
  };
};

// Rigid Body: Pose and Velocities
template <typename Scalar> struct R9SO3Squared {};

// SE3
template <typename Scalar> struct R3SO3Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 7, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 7, 1>>;
  using cref3_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 3, 1>> &;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  RnSquared<Scalar, 3> l2;
  SO3Squared<Scalar> so3;

  void set_bounds(cref3_t lb_, cref3_t ub_) { l2.set_bounds(lb_, ub_); }

  inline void sample_uniform(cref3_t lb, cref3_t ub, ref_t x) const {
    l2.sample_uniform(x.template head<3>());
    so3.sample_uniform(x.template tail<4>());
  }

  inline Scalar distance_to_rectangle(cref_t &location1, cref_t &lb,
                                      cref_t &ub) const {

    Scalar d1 =
        l2.distance_to_rectangle(location1.template head<3>(),
                                 lb.template head<3>(), ub.template head<3>());

    Scalar d2 =
        so3.distance_to_rectangle(location1.template tail<4>(),
                                  lb.template tail<4>(), ub.template tail<4>());

    return d1 + d2;
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    Scalar d1 =
        l2.distance(location1.template head<3>(), location2.template head<3>());
    Scalar d2 = so3.distance(location1.template tail<4>(),
                             location2.template tail<4>());
    return d1 + d2;
  };
};

template <typename Scalar> struct R3SO3 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 7, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 7, 1>>;

  using cref3_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 3, 1>> &;
  using ref3_t = Eigen::Ref<Eigen::Matrix<Scalar, 3, 1>>;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {

    std::stringstream error_msg;
    error_msg << "not implemented interpolate in " << __PRETTY_FUNCTION__;
    throw std::runtime_error(error_msg.str());
  }

  Rn<Scalar, 3> l2;
  SO3<Scalar> so3;

  void set_bounds(cref3_t lb_, cref3_t ub_) { l2.set_bounds(lb_, ub_); }

  void sample_uniform(ref_t x) const {
    l2.sample_uniform(x.template head<3>());
    so3.sample_uniform(x.template tail<4>());
  }

  inline Scalar distance_to_rectangle(cref_t &location1, cref_t &lb,
                                      cref_t &ub) const {

    Scalar d1 =
        l2.distance_to_rectangle(location1.template head<3>(),
                                 lb.template head<3>(), ub.template head<3>());

    Scalar d2 =
        so3.distance_to_rectangle(location1.template tail<4>(),
                                  lb.template tail<4>(), ub.template tail<4>());

    return d1 + d2;
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    Scalar d1 =
        l2.distance(location1.template head<3>(), location2.template head<3>());
    Scalar d2 = so3.distance(location1.template tail<4>(),
                             location2.template tail<4>());
    return d1 + d2;
  };
};

enum class DistanceType {
  RnL1,
  Rn,
  RnSquared,
  SO2,
  SO2Squared,
  SO3,
  SO3Squared
};

inline bool starts_with(const std::string &str, const std::string &prefix) {
  return str.size() >= prefix.size() &&
         str.compare(0, prefix.size(), prefix) == 0;
}

// get the substring after : as an integer number
inline int get_number(const std::string &str) {
  std::string delimiter = ":";
  size_t pos = 0;
  std::string token;
  // while ((
  pos = str.find(delimiter);
  if (pos == std::string::npos) {
    throw std::runtime_error("delimiter not found");
  }

  token = str.substr(pos + delimiter.length(), str.size());
  int out = std::stoi(token);
  std::cout << "out " << out << std::endl;
  return out;
}

template <typename Scalar> struct Combined {
  // TODO: test this!! How i am going to give this as input? -- it is not a
  // static function anymore...
  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, -1, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, -1, 1>>;

  using Space =
      std::variant<RnL1<Scalar>, Rn<Scalar>, RnSquared<Scalar>, SO2<Scalar>,
                   SO2Squared<Scalar>, SO3<Scalar>, SO3Squared<Scalar>>;
  std::vector<Space> spaces;
  std::vector<int> dims; // TODO: remove this and get auto from spaces

  Combined(const std::vector<Space> &spaces, const std::vector<int> &dims)
      : spaces(spaces), dims(dims) {
    assert(spaces.size() == dims.size());
  }

  Combined(const std::vector<std::string> &spaces_str) {
    // throw std::runtime_error("not implemented " + __PRETTY_FUNCTION__);

    for (size_t i = 0; i < spaces_str.size(); i++) {

      if (spaces_str.at(i) == "SO2") {
        spaces.push_back(SO2<Scalar>());
        dims.push_back(1);
      } else if (spaces_str.at(i) == "SO2Squared") {
        spaces.push_back(SO2Squared<Scalar>());
        dims.push_back(1);
      } else if (spaces_str.at(i) == "SO3") {
        spaces.push_back(SO3<Scalar>());
        dims.push_back(4);
      } else if (spaces_str.at(i) == "SO3Squared") {
        spaces.push_back(SO3Squared<Scalar>());
      } else if (starts_with(spaces_str.at(i), "RnL1")) {
        spaces.push_back(RnL1<Scalar>());
        int dim = get_number(spaces_str.at(i));
        dims.push_back(dim);
      } else if (starts_with(spaces_str.at(i), "Rn") &&
                 !starts_with(spaces_str.at(i), "RnSquared")) {
        spaces.push_back(Rn<Scalar>());
        int dim = get_number(spaces_str.at(i));
        dims.push_back(dim);
      } else if (starts_with(spaces_str.at(i), "RnSquared")) {
        spaces.push_back(RnSquared<Scalar>());
        int dim = get_number(spaces_str.at(i));
        dims.push_back(dim);
      } else {
        std::stringstream error_msg;
        error_msg << "Unknown space " << spaces_str.at(i) << " in "
                  << __PRETTY_FUNCTION__ << std::endl
                  << __FILE__ << ":" << __LINE__;
        throw std::runtime_error(error_msg.str());
      }
    }
    assert(spaces.size() == dims.size());
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  void set_bounds(const std::vector<Eigen::VectorXd> &lbs,
                  const std::vector<Eigen::VectorXd> &ubs) {

    assert(lbs.size() == ubs.size());

    for (size_t i = 0; i < lbs.size(); i++) {
      std::visit(
          [&](auto &obj) {
            if (lbs[i].size())
              obj.set_bounds(lbs[i], ubs[i]);
          },
          spaces[i]);
    }
  }

  void sample_uniform(ref_t x) const {

    assert(spaces.size() == dims.size());
    assert(spaces.size());
    int counter = 0;
    for (size_t i = 0; i < spaces.size(); i++) {
      std::visit(
          [&](const auto &obj) {
            obj.sample_uniform(x.segment(counter, dims[i]));
          },
          spaces[i]);
      counter += dims[i];
    }
  }

  // void sample_uniform(std::vector<cref_t> lb, std::vector<cref_t> ub, ref_t
  // x) const {
  //   assert(spaces.size() == dims.size());
  //   assert(spaces.size());
  //   int counter = 0;
  //   for (size_t i = 0; i < spaces.size(); i++) {
  //     std::visit(
  //         [&](const auto &obj) {
  //           obj.sample_uniform(x.segment(counter, dims[i]));
  //         },
  //         spaces[i]);
  //     counter += dims[i];
  //   }
  //
  // }

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {

    assert(spaces.size() == dims.size());
    assert(spaces.size());
    Scalar d = 0;
    int counter = 0;
    for (size_t i = 0; i < spaces.size(); i++) {
      std::visit(
          [&](const auto &obj) {
            obj.interpolate(from.segment(counter, dims[i]),
                            to.segment(counter, dims[i]), t,
                            out.segment(counter, dims[i]));
          },
          spaces[i]);
      counter += dims[i];
    }
  }

  inline Scalar distance(cref_t location1, cref_t location2) const {

    assert(spaces.size() == dims.size());
    assert(spaces.size());
    Scalar d = 0;
    int counter = 0;
    for (size_t i = 0; i < spaces.size(); i++) {
      auto caller = [&](const auto &obj) {
        return obj.distance(location1.segment(counter, dims[i]),
                            location2.segment(counter, dims[i]));
      };

      d += std::visit(caller, spaces[i]);
      counter += dims[i];
    }
    return d;
  }

  inline Scalar distance_to_rectangle(cref_t location1, cref_t lb,
                                      cref_t ub) const {

    assert(spaces.size() == dims.size());
    assert(spaces.size());

    Scalar d = 0;
    int counter = 0;
    for (size_t i = 0; i < spaces.size(); i++) {

      auto caller = [&](const auto &obj) {
        return obj.distance_to_rectangle(location1.segment(counter, dims[i]),
                                         lb.segment(counter, dims[i]),
                                         ub.segment(counter, dims[i]));
      };

      d += std::visit(caller, spaces[i]);
      counter += dims[i];
    }

    return d;
  }
  //   switch (spaces[i]) {
  //   case DistanceType::RnL1: {
  //     d += RnL1<Scalar, -1>::distance_to_rectangle(
  //         location1.template segment(counter, dims[i]),
  //         lb.template segment(counter, dims[i]),
  //         ub.template segment(counter, dims[i]));
  //   } break;
  //   case DistanceType::Rn: {
  //     d += Rn<Scalar, -1>::distance_to_rectangle(
  //         location1.template segment(counter, dims[i]),
  //         lb.template segment(counter, dims[i]),
  //         ub.template segment(counter, dims[i]));
  //   } break;
  //   case DistanceType::RnSquared: {
  //     d += RnSquared<Scalar, -1>::distance_to_rectangle(
  //         location1.template segment(counter, dims[i]),
  //         lb.template segment(counter, dims[i]),
  //         ub.template segment(counter, dims[i]));
  //   } break;
  //   case DistanceType::SO2: {
  //     assert(dims[i] == 1);
  //     d += SO2<Scalar>::distance_to_rectangle(
  //         location1.template segment<1>(counter),
  //         lb.template segment<1>(counter), ub.template
  //         segment<1>(counter));
  //   } break;
  //   case DistanceType::SO2Squared: {
  //     assert(dims[i] == 1);
  //     d += SO2Squared<Scalar>::distance_to_rectangle(
  //         location1.template segment<1>(counter),
  //         lb.template segment<1>(counter), ub.template
  //         segment<1>(counter));
  //   } break;
  //   case DistanceType::SO3: {
  //     assert(dims[i] == 4);
  //     d += SO3<Scalar>::distance_to_rectangle(
  //         location1.template segment<4>(counter),
  //         lb.template segment<4>(counter), ub.template
  //         segment<4>(counter));
  //   } break;
  //   case DistanceType::SO3Squared: {
  //     assert(dims[i] == 4);
  //     d += SO3Squared<Scalar>::distance_to_rectangle(
  //         location1.template segment<4>(counter),
  //         lb.template segment<4>(counter), ub.template
  //         segment<4>(counter));
  //   } break;
  //   default:
  //     assert(false);
  //   }
  //   counter += dims[i];
  // }
  // return d;

  // protected:
  //   std::vector<DistanceType> spaces;
  //   std::vector<int> dims;
};

template <class Id, int Dimensions, std::size_t BucketSize = 32,
          typename Scalar = double,
          typename StateSpace = RnSquared<Scalar, Dimensions>>
class KDTree {
private:
  struct Node;
  std::vector<Node> m_nodes;
  std::set<std::size_t> waitingForSplit;
  StateSpace state_space;

public:
  using scalar_t = Scalar;
  using id_t = Id;
  using point_t = Eigen::Matrix<Scalar, Dimensions, 1>;
  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;
  int m_dimensions = Dimensions;
  static const std::size_t bucketSize = BucketSize;
  // TODO: I also want Dimensions at runtime!!
  using tree_t = KDTree<Id, Dimensions, BucketSize, Scalar, StateSpace>;

  // void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
  //   state_space.interpolate(from, to, t, out);
  // }
  //
  StateSpace &getStateSpace() { return state_space; }

  KDTree(int runtime_dimension = -1,
         const StateSpace &state_space = StateSpace())
      : state_space(state_space) {

    if constexpr (Dimensions == Eigen::Dynamic) {
      assert(runtime_dimension > 0);
      m_dimensions = runtime_dimension;
      m_nodes.emplace_back(BucketSize, m_dimensions);
    } else {
      m_nodes.emplace_back(BucketSize, -1);
    }
  }

  size_t size() const { return m_nodes[0].m_entries; }

  void addPoint(const point_t &location, const Id &payload,
                bool autosplit = true) {
    std::size_t addNode = 0;

    assert(m_dimensions > 0);
    while (m_nodes[addNode].m_splitDimension != m_dimensions) {
      m_nodes[addNode].expandBounds(location);
      if (location[m_nodes[addNode].m_splitDimension] <
          m_nodes[addNode].m_splitValue) {
        addNode = m_nodes[addNode].m_children.first;
      } else {
        addNode = m_nodes[addNode].m_children.second;
      }
    }
    m_nodes[addNode].add(LocationId{location, payload});

    if (m_nodes[addNode].shouldSplit() &&
        m_nodes[addNode].m_entries % BucketSize == 0) {
      if (autosplit) {
        split(addNode);
      } else {
        waitingForSplit.insert(addNode);
      }
    }
  }

  void splitOutstanding() {
    std::vector<std::size_t> searchStack(waitingForSplit.begin(),
                                         waitingForSplit.end());
    waitingForSplit.clear();
    while (searchStack.size() > 0) {
      std::size_t addNode = searchStack.back();
      searchStack.pop_back();
      if (m_nodes[addNode].m_splitDimension == m_dimensions &&
          m_nodes[addNode].shouldSplit() && split(addNode)) {
        searchStack.push_back(m_nodes[addNode].m_children.first);
        searchStack.push_back(m_nodes[addNode].m_children.second);
      }
    }
  }

  struct DistanceId {
    Scalar distance;
    Id payload;
    inline bool operator<(const DistanceId &dp) const {
      return distance < dp.distance;
    }
  };

  std::vector<DistanceId> searchKnn(const point_t &location,
                                    std::size_t maxPoints) const {
    return searcher().search(location, std::numeric_limits<Scalar>::max(),
                             maxPoints, state_space);
  }

  std::vector<DistanceId> searchBall(const point_t &location,
                                     Scalar maxRadius) const {
    return searcher().search(location, maxRadius,
                             std::numeric_limits<std::size_t>::max(),
                             state_space);
  }

  std::vector<DistanceId>
  searchCapacityLimitedBall(const point_t &location, Scalar maxRadius,
                            std::size_t maxPoints) const {
    return searcher().search(location, maxRadius, maxPoints, state_space);
  }

  DistanceId search(const point_t &location) const {
    DistanceId result;
    result.distance = std::numeric_limits<Scalar>::infinity();

    if (m_nodes[0].m_entries > 0) {
      std::vector<std::size_t> searchStack;
      searchStack.reserve(
          1 +
          std::size_t(1.5 * std::log2(1 + m_nodes[0].m_entries / BucketSize)));
      searchStack.push_back(0);

      while (searchStack.size() > 0) {
        std::size_t nodeIndex = searchStack.back();
        searchStack.pop_back();
        const Node &node = m_nodes[nodeIndex];
        if (result.distance > node.pointRectDist(location, state_space)) {
          if (node.m_splitDimension == m_dimensions) {
            for (const auto &lp : node.m_locationId) {
              Scalar nodeDist = state_space.distance(location, lp.location);
              if (nodeDist < result.distance) {
                result = DistanceId{nodeDist, lp.payload};
              }
            }
          } else {
            node.queueChildren(location, searchStack);
          }
        }
      }
    }
    return result;
  }

  class Searcher {
  public:
    Searcher(const tree_t &tree) : m_tree(tree) {}
    Searcher(const Searcher &searcher) : m_tree(searcher.m_tree) {}

    // NB! this method is not const. Do not call this on same instance from
    // different threads simultaneously.
    const std::vector<DistanceId> &search(const point_t &location,
                                          Scalar maxRadius,
                                          std::size_t maxPoints,
                                          const StateSpace &state_space) {
      // clear results from last time
      m_results.clear();

      // reserve capacities
      m_searchStack.reserve(
          1 + std::size_t(1.5 * std::log2(1 + m_tree.m_nodes[0].m_entries /
                                                  BucketSize)));
      if (m_prioqueueCapacity < maxPoints &&
          maxPoints < m_tree.m_nodes[0].m_entries) {
        std::vector<DistanceId> container;
        container.reserve(maxPoints);
        m_prioqueue = std::priority_queue<DistanceId, std::vector<DistanceId>>(
            std::less<DistanceId>(), std::move(container));
        m_prioqueueCapacity = maxPoints;
      }

      m_tree.searchCapacityLimitedBall(location, maxRadius, maxPoints,
                                       m_searchStack, m_prioqueue, m_results,
                                       state_space);

      m_prioqueueCapacity = std::max(m_prioqueueCapacity, m_results.size());
      return m_results;
    }

  private:
    const tree_t &m_tree;

    std::vector<std::size_t> m_searchStack;
    std::priority_queue<DistanceId, std::vector<DistanceId>> m_prioqueue;
    std::size_t m_prioqueueCapacity = 0;
    std::vector<DistanceId> m_results;
  };

  // NB! returned class has no const methods. Get one instance per thread!
  Searcher searcher() const { return Searcher(*this); }

private:
  struct LocationId {
    point_t location;
    Id payload;
  };
  std::vector<LocationId> m_bucketRecycle;

  void searchCapacityLimitedBall(
      const point_t &location, Scalar maxRadius, std::size_t maxPoints,
      std::vector<std::size_t> &searchStack,
      std::priority_queue<DistanceId, std::vector<DistanceId>> &prioqueue,
      std::vector<DistanceId> &results, const StateSpace &state_space) const {
    std::size_t numSearchPoints = std::min(maxPoints, m_nodes[0].m_entries);

    if (numSearchPoints > 0) {
      searchStack.push_back(0);
      while (searchStack.size() > 0) {
        std::size_t nodeIndex = searchStack.back();
        searchStack.pop_back();
        const Node &node = m_nodes[nodeIndex];
        Scalar minDist = node.pointRectDist(location, state_space);
        if (maxRadius > minDist && (prioqueue.size() < numSearchPoints ||
                                    prioqueue.top().distance > minDist)) {
          if (node.m_splitDimension == m_dimensions) {
            node.searchCapacityLimitedBall(location, maxRadius, numSearchPoints,
                                           prioqueue, state_space);
          } else {
            node.queueChildren(location, searchStack);
          }
        }
      }

      results.reserve(prioqueue.size());
      while (prioqueue.size() > 0) {
        results.push_back(prioqueue.top());
        prioqueue.pop();
      }
      std::reverse(results.begin(), results.end());
    }
  }

  bool split(std::size_t index) {
    if (m_nodes.capacity() < m_nodes.size() + 2) {
      m_nodes.reserve((m_nodes.capacity() + 1) * 2);
    }
    Node &splitNode = m_nodes[index];
    splitNode.m_splitDimension = m_dimensions;
    Scalar width(0);
    state_space.choose_split_dimension(splitNode.m_lb, splitNode.m_ub,
                                       splitNode.m_splitDimension, width);

    if (splitNode.m_splitDimension == m_dimensions) {
      return false;
    }

    std::vector<Scalar> splitDimVals;
    splitDimVals.reserve(splitNode.m_entries);
    for (const auto &lp : splitNode.m_locationId) {
      splitDimVals.push_back(lp.location[splitNode.m_splitDimension]);
    }
    std::nth_element(splitDimVals.begin(),
                     splitDimVals.begin() + splitDimVals.size() / 2 + 1,
                     splitDimVals.end());
    std::nth_element(splitDimVals.begin(),
                     splitDimVals.begin() + splitDimVals.size() / 2,
                     splitDimVals.begin() + splitDimVals.size() / 2 + 1);
    splitNode.m_splitValue = (splitDimVals[splitDimVals.size() / 2] +
                              splitDimVals[splitDimVals.size() / 2 + 1]) /
                             Scalar(2);

    splitNode.m_children = std::make_pair(m_nodes.size(), m_nodes.size() + 1);
    std::size_t entries = splitNode.m_entries;
    m_nodes.emplace_back(m_bucketRecycle, entries, m_dimensions);
    Node &leftNode = m_nodes.back();
    m_nodes.emplace_back(entries, m_dimensions);
    Node &rightNode = m_nodes.back();

    for (const auto &lp : splitNode.m_locationId) {
      if (lp.location[splitNode.m_splitDimension] < splitNode.m_splitValue) {
        leftNode.add(lp);
      } else {
        rightNode.add(lp);
      }
    }

    if (leftNode.m_entries ==
        0) // points with equality to splitValue go in rightNode
    {
      splitNode.m_splitValue = 0;
      splitNode.m_splitDimension = m_dimensions;
      splitNode.m_children = std::pair<std::size_t, std::size_t>(0, 0);
      std::swap(rightNode.m_locationId, m_bucketRecycle);
      m_nodes.pop_back();
      m_nodes.pop_back();
      return false;
    } else {
      splitNode.m_locationId.clear();
      // if it was a standard sized bucket, recycle the memory to reduce
      // allocator pressure otherwise clear the memory used by the bucket
      // since it is a branch not a leaf anymore
      if (splitNode.m_locationId.capacity() == BucketSize) {
        std::swap(splitNode.m_locationId, m_bucketRecycle);
      } else {
        std::vector<LocationId> empty;
        std::swap(splitNode.m_locationId, empty);
      }
      return true;
    }
  }

  struct Node {
    Node(std::size_t capacity, int runtime_dimension = -1) {
      init(capacity, runtime_dimension);
    }

    Node(std::vector<LocationId> &recycle, std::size_t capacity,
         int runtime_dimension) {
      std::swap(m_locationId, recycle);
      init(capacity, runtime_dimension);
    }

    void init(std::size_t capacity, int runtime_dimension) {

      if constexpr (Dimensions == Eigen::Dynamic) {
        assert(runtime_dimension > 0);
        m_lb.resize(runtime_dimension);
        m_ub.resize(runtime_dimension);
        m_splitDimension = runtime_dimension;
      }

      m_lb.setConstant(std::numeric_limits<Scalar>::max());
      m_ub.setConstant(std::numeric_limits<Scalar>::lowest());
      m_locationId.reserve(std::max(BucketSize, capacity));
    }

    void expandBounds(const point_t &location) {
      m_lb = m_lb.cwiseMin(location);
      m_ub = m_ub.cwiseMax(location);
      m_entries++;
    }

    void add(const LocationId &lp) {
      expandBounds(lp.location);
      m_locationId.push_back(lp);
    }

    bool shouldSplit() const { return m_entries >= BucketSize; }

    void searchCapacityLimitedBall(const point_t &location, Scalar maxRadius,
                                   std::size_t K,
                                   std::priority_queue<DistanceId> &results,
                                   const StateSpace &state_space) const {

      std::size_t i = 0;

      // this fills up the queue if it isn't full yet
      for (; results.size() < K && i < m_entries; i++) {
        const auto &lp = m_locationId[i];
        Scalar distance = state_space.distance(location, lp.location);
        if (distance < maxRadius) {
          results.emplace(DistanceId{distance, lp.payload});
        }
      }

      // this adds new things to the queue once it is full
      for (; i < m_entries; i++) {
        const auto &lp = m_locationId[i];
        Scalar distance = state_space.distance(location, lp.location);
        if (distance < maxRadius && distance < results.top().distance) {
          results.pop();
          results.emplace(DistanceId{distance, lp.payload});
        }
      }
    }

    void queueChildren(const point_t &location,
                       std::vector<std::size_t> &searchStack) const {
      if (location[m_splitDimension] < m_splitValue) {
        searchStack.push_back(m_children.second);
        searchStack.push_back(m_children.first); // left is popped first
      } else {
        searchStack.push_back(m_children.first);
        searchStack.push_back(m_children.second); // right is popped first
      }
    }

    Scalar pointRectDist(const point_t &location,
                         const StateSpace &distance) const {
      return distance.distance_to_rectangle(location, m_lb, m_ub);
    }

    std::size_t m_entries = 0; /// size of the tree, or subtree

    int m_splitDimension = Dimensions; /// split dimension of this node
    Scalar m_splitValue = 0;           /// split value of this node

    // struct Range {
    //   Scalar min, max;
    // };

    // std::array<Range, Dimensions> m_bounds; /// bounding box of this node
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    Eigen::Matrix<Scalar, Dimensions, 1> m_lb;
    Eigen::Matrix<Scalar, Dimensions, 1> m_ub;

    std::pair<std::size_t, std::size_t>
        m_children; /// subtrees of this node (if not a leaf)
    std::vector<LocationId> m_locationId; /// data held in this node (if a leaf)
  };
};

// template class myclass<int>;
//
//
// template
//
// class

// template <class Id, std::size_t Dimensions, std::size_t BucketSize =
// 32,
//           class Distance = RnSquared, typename Scalar = double>

// template class KDTree<size_t, 4, 32, double, SO3Squared<double>>;
// template class KDTree<size_t, 2, 32>;
// template class KDTree<size_t, -1, 32>;
// TODO:...
// Create a Python Interface :)

} // namespace dynotree
//
//
//
// <class Id, std::size_t Dimensions, std::size_t BucketSize = 32,
//           class Distance = RnSquared, typename Scalar = double>
//
