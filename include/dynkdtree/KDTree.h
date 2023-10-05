#pragma once

#include <cwchar>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>

template <class T> using cref = const Eigen::Ref<const T> &;

/**
 * KDTree.h by Julian Kent
 * A C++11 KD-Tree with the following features:
 *     single file
 *     header only
 *     high performance K Nearest Neighbor and ball searches
 *     dynamic insertions
 *     simple API
 *     depends only on the STL
 *     templatable on your custom data type to store in the leaves. No need to
 * keep a separate data structure! templatable on double, float etc templatable
 * on L1, L2Squared or custom distance functor templated on number of dimensions
 * for efficient inlining
 *
 * -------------------------------------------------------------------
 *
 * This Source Code Form is subject to the terms of the Mozilla Public License,
 * v. 2.0. If a copy of the MPL was not distributed with this  file, you can
 * obtain one at http://mozilla.org/MPL/2.0/.
 *
 * A high level explanation of MPLv2: You may use this in any software provided
 * you give attribution. You *must* make available any changes you make to the
 * source code of this file to anybody you distribute your software to.
 *
 * Upstreaming features and bugfixes are highly appreciated via:
 *
 * https://github.com/jkflying/bucket-pr-kdtree/tree/master/C%2B%2B
 *
 * For additional licensing rights, feature requests or questions, please
 * contact Julian Kent <jkflying@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * Example usage:
 *
 * // setup
 * using tree_t = jk::tree::KDTree<std::string, 2>;
 * using point_t = std::array<double, 2>;
 * tree_t tree;
 * tree.addPoint(point_t{{1, 2}}, "George");
 * tree.addPoint(point_t{{1, 3}}, "Harold");
 * tree.addPoint(point_t{{7, 7}}, "Melvin");
 *
 * // KNN search
 * point_t lazyMonsterLocation{{6, 6}}; // this monster will always try to eat
 * the closest people const std::size_t monsterHeads = 2; // this monster can
 * eat two people at once auto lazyMonsterVictims =
 * tree.searchKnn(lazyMonsterLocation, monsterHeads); for (const auto& victim :
 * lazyMonsterVictims)
 * {
 *     std::cout << victim.payload << " closest to lazy monster, with distance "
 * << sqrt(victim.distance) << "!"
 *               << std::endl;
 * }
 *
 * // ball search
 * point_t stationaryMonsterLocation{{8, 8}}; // this monster doesn't move, so
 * can only eat people that are close const double neckLength = 6.0; // it can
 * only reach within this range auto potentialVictims =
 * tree.searchBall(stationaryMonsterLocation, neckLength * neckLength); //
 * metric is L2Squared std::cout << "Stationary monster can reach any of " <<
 * potentialVictims.size() << " people!" << std::endl;
 *
 * // hybrid KNN/ball search
 * auto actualVictims
 *     = tree.searchCapacityLimitedBall(stationaryMonsterLocation, neckLength *
 * neckLength, monsterHeads); std::cout << "The stationary monster will try to
 * eat "; for (const auto& victim : actualVictims)
 * {
 *     std::cout << victim.payload << " and ";
 * }
 * std::cout << "nobody else." << std::endl;
 *
 * Output:
 *
 * Melvin closest to lazy monster, with distance 1.41421!
 * Harold closest to lazy monster, with distance 5.83095!
 * Stationary monster can reach any of 1 people!
 * The stationary monster will try to eat Melvin and nobody else.
 *
 * -------------------------------------------------------------------
 *
 * Tuning tips:
 *
 * If you need to add a lot of points before doing any queries, set the optional
 * `autosplit` parameter to false, then call splitOutstanding(). This will
 * reduce temporaries and result in a better balanced tree.
 *
 * Set the bucket size to be at least twice the K in a typical KNN query. If you
 * have more dimensions, it is better to have a larger bucket size. 32 is a good
 * starting point. If possible use powers of 2 for the bucket size.
 *
 * If you experience linear search performance, check that you don't have a
 * bunch of duplicate point locations. This will result in the tree being unable
 * to split the bucket the points are in, degrading search performance.
 *
 * The tree adapts to the parallel-to-axis dimensionality of the problem. Thus,
 * if there is one dimension with a much larger scale than the others, most of
 * the splitting will happen on this dimension. This is achieved by trying to
 * keep the bounding boxes of the data in the buckets equal lengths in all axes.
 *
 * Random data performs worse than 'real world' data with structure. This is
 * because real world data has tighter bounding boxes, meaning more branches of
 * the tree can be eliminated sooner.
 *
 * On pure random data, more than 7 dimensions won't be much faster than linear.
 * However, most data isn't actually random. The tree will adapt to any locally
 * reduced dimensionality, which is found in most real world data.
 *
 * Hybrid ball/KNN searches are faster than either type on its own, because
 * subtrees can be more aggresively eliminated.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <vector>

namespace jk {
namespace tree {

template <typename T, typename Scalar>
void choose_split_dimension_default(T lb, T ub, int &ii, Scalar &width) {
  for (std::size_t i = 0; i < lb.size(); i++) {
    Scalar dWidth = ub[i] - lb[i];
    if (dWidth > width) {
      ii = i;
      width = dWidth;
    }
  }
}

template <typename Scalar, int Dimensions> struct L1 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    double d = 0;
    double dist = 0;

    if constexpr (Dimensions == Eigen::Dynamic) {

      assert(location1.size());
      assert(ub.size());
      assert(lb.size());
      assert(location1.size() == ub.size());
      assert(location1.size() == lb.size());

      for (size_t i = 0; i < location1.size(); i++) {
        double xx = std::max(lb(i), std::min(ub(i), location1(i)));
        double dif = xx - location1(i);
        dist += std::abs(dif);
      }
    } else {
      for (size_t i = 0; i < Dimensions; i++) {
        double xx = std::max(lb(i), std::min(ub(i), location1(i)));
        double dif = xx - location1(i);
        dist += std::abs(dif);
      }
    }
    return dist;
  }

  Scalar distance(cref_t &location1, cref_t &location2) {

    assert(location1.size());
    assert(location2.size());
    return (location1 - location2).cwiseAbs().sum();
  }
};

template <typename Scalar> struct SO2 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 1, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 1, 1>>;

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

  Scalar distanceToRect(cref_t location1, cref_t lb, cref_t ub) const {

    assert(location1(0) >= -M_PI);
    assert(location1(0) <= M_PI);

    assert(lb(0) >= -M_PI);
    assert(lb(0) <= M_PI);

    assert(ub(0) >= -M_PI);
    assert(ub(0) <= M_PI);

    if (location1(0) >= lb(0) && location1(0) <= ub(0)) {
      return 0;
    } else if (location1(0) > ub(0)) {
      double d1 = location1(0) - ub(0);
      double d2 = lb(0) - (location1(0) - 2 * M_PI);
      assert(d2 >= 0);
      assert(d1 >= 0);
      return std::min(d1, d2);
    } else if (location1(0) < lb(0)) {
      double d1 = lb(0) - location1(0);
      double d2 = (location1(0) + 2 * M_PI) - ub(0);
      assert(d2 >= 0);
      assert(d1 >= 0);
      return std::min(d1, d2);
    } else {
      assert(false);
    }
  }

  Scalar distance(cref_t location1, cref_t location2) const {

    assert(location1(0) >= -M_PI);
    assert(location2(0) >= -M_PI);

    assert(location1(0) <= M_PI);
    assert(location2(0) <= M_PI);

    double dif = location1(0) - location2(0);
    if (dif > M_PI) {
      dif -= 2 * M_PI;
    } else if (dif < -M_PI) {
      dif += 2 * M_PI;
    }
    double out = std::abs(dif);
    return out;
  }
};

template <typename Scalar> struct SO2Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 1, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 1, 1>>;

  SO2<Scalar> so2;

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    so2.interpolate(from, to, t, out);
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar distanceToRect(cref_t location1, cref_t lb, cref_t ub) const {

    double d = so2.distanceToRect(location1, lb, ub);
    return d * d;
  }

  Scalar distance(cref_t location1, cref_t location2) const {

    double d = so2.distance(location1, location2);
    return d * d;
  }
};

template <typename Scalar, int Dimensions> struct L2Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    out = from + t * (to - from);
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    double d = 0;
    double dist = 0;

    if constexpr (Dimensions == Eigen::Dynamic) {

      assert(location1.size());
      assert(ub.size());
      assert(lb.size());
      assert(location1.size() == ub.size());
      assert(location1.size() == lb.size());

      for (size_t i = 0; i < location1.size(); i++) {
        double xx = std::max(lb(i), std::min(ub(i), location1(i)));
        double dif = xx - location1(i);
        dist += dif * dif;
      }
    } else {
      for (size_t i = 0; i < Dimensions; i++) {
        double xx = std::max(lb(i), std::min(ub(i), location1(i)));
        double dif = xx - location1(i);
        dist += dif * dif;
      }
    }
    return dist;
  }

  Scalar distance(cref_t &location1, cref_t &location2) const {

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

template <typename Scalar, int Dimensions> struct L2 {

  using L2squared = L2Squared<Scalar, Dimensions>;
  L2squared l2squared;

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    out = from + t * (to - from);
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    double d = l2squared.distanceToRect(location1, lb, ub);
    return std::sqrt(d);
  };

  Scalar distance(cref_t &location1, cref_t &location2) const {
    double d = l2squared.distance(location1, location2);
    return std::sqrt(d);
  }
};

template <typename Scalar> struct R2SO2 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 3, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 3, 1>> ;



  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  double angular_weight = 1.0;

  L2<Scalar, 2> l2;
  SO2<Scalar> so2;


  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    assert(t >= 0);
    assert(t <= 1);
    l2.interpolate(from.template head<2>(), to.template head<2>(), t,
                   out.template head<2>());
    so2.interpolate(from.template tail<1>(), to.template tail<1>(), t,
                    out.template tail<1>());
    std::cout << out.transpose() << std::endl;

  }




  Scalar distanceToRect(cref_t location1, cref_t lb, cref_t ub) const {

    double d1 = l2.distanceToRect(location1.template head<2>(),
                                  lb.template head<2>(), ub.template head<2>());
    double d2 =
        so2.distanceToRect(location1.template tail<1>(), lb.template tail<1>(),
                           ub.template tail<1>());
    return d1 + angular_weight * d2;
  }

  Scalar distance(cref_t location1, cref_t location2) const {

    double d1 =
        l2.distance(location1.template head<2>(), location2.template head<2>());
    double d2 = so2.distance(location1.template tail<1>(),
                             location2.template tail<1>());
    return d1 + angular_weight * d2;
  };
};

template <typename Scalar> struct R2SO2Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 3, 1>> &;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  double angular_weight = 1.0;

  L2Squared<Scalar, 2> l2squared;
  SO2Squared<Scalar> so2squared;
  Scalar distanceToRect(cref_t location1, cref_t lb, cref_t ub) const {

    double d1 =
        l2squared.distanceToRect(location1.template head<2>(),
                                 lb.template head<2>(), ub.template head<2>());
    double d2 =
        so2squared.distanceToRect(location1.template tail<1>(),
                                  lb.template tail<1>(), ub.template tail<1>());
    return d1 + angular_weight * d2;
  }

  Scalar distance(cref_t location1, cref_t location2) const {

    double d1 = l2squared.distance(location1.template head<2>(),
                                   location2.template head<2>());
    double d2 = so2squared.distance(location1.template tail<1>(),
                                    location2.template tail<1>());
    return d1 + angular_weight * d2;
  };
};

template <typename Scalar> struct SO3Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 4, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, 4, 1>>;

  L2Squared<Scalar, 4> l2squared;

  void interpolate(cref_t from, cref_t to, Scalar t, ref_t out) const {
    throw std::runtime_error("not implemented interpolate in so3");
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    assert(std::abs(location1.norm() - 1) < 1e-6);

    double d1 = l2squared.distanceToRect(location1, lb, ub);
    double d2 = l2squared.distanceToRect(-1. * location1, lb, ub);
    return std::min(d1, d2);
  }

  double distance(cref_t location1, cref_t location2) const {

    assert(location1.size() == 4);
    assert(location2.size() == 4);

    assert(std::abs(location1.norm() - 1) < 1e-6);
    assert(std::abs(location2.norm() - 1) < 1e-6);

    double d1 = l2squared.distance(location1, location2);
    double d2 = l2squared.distance(-location1, location2);
    return std::min(d1, d2);
  };
};

template <typename Scalar> struct SO3 {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 4, 1>> &;

  SO3Squared<Scalar> so3squared;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    return std::sqrt(so3squared.distanceToRect(location1, lb, ub));
  }

  double distance(cref_t location1, cref_t location2) const {

    return std::sqrt(so3squared.distance(location1, location2));
  };
};

// Rigid Body: Pose and Velocities
template <typename Scalar> struct R9SO3Squared {};

// SE3
template <typename Scalar> struct R3SO3Squared {

  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, 7, 1>> &;

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  L2Squared<Scalar, 3> l2;
  SO3Squared<Scalar> so3;

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    double d1 = l2.distanceToRect(location1.template head<3>(),
                                  lb.template head<3>(), ub.template head<3>());

    double d2 =
        so3.distanceToRect(location1.template tail<4>(), lb.template tail<4>(),
                           ub.template tail<4>());

    return d1 + d2;
  }

  double distance(cref_t location1, cref_t location2) const {

    double d1 =
        l2.distance(location1.template head<3>(), location2.template head<3>());
    double d2 = so3.distance(location1.template tail<4>(),
                             location2.template tail<4>());
    return d1 + d2;
  };
};

enum class DistanceType { L1, L2, L2Squared, SO2, SO2Squared, SO3, SO3Squared };

template <typename Scalar> struct Combined {
  // TODO: test this!! How i am going to give this as input? -- it is not a
  // static function anymore...
  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, -1, 1>> &;

  std::vector<void *> m_distances;
  Combined(std::vector<DistanceType> types, std::vector<int> dims)
      : types(types), dims(dims) {

    assert(types.size() == dims.size());
  }

  void choose_split_dimension(cref_t lb, cref_t ub, int &ii, Scalar &width) {
    choose_split_dimension_default(lb, ub, ii, width);
  }

  double distance(cref_t location1, cref_t location2) const {

    double d = 0;
    int counter = 0;
    for (size_t i = 0; i < dims.size(); i++) {
      switch (types[i]) {
      case DistanceType::L1: {
        d += L1<Scalar, -1>::distance(
            location1.template segment(counter, dims[i]),
            location2.template segment(counter, dims[i]));
      } break;
      case DistanceType::L2: {
        d += L2<Scalar, -1>::distance(
            location1.template segment(counter, dims[i]),
            location2.template segment(counter, dims[i]));
      } break;
      case DistanceType::L2Squared: {
        d += L2Squared<Scalar, -1>::distance(
            location1.template segment(counter, dims[i]),
            location2.template segment(counter, dims[i]));
      } break;
      case DistanceType::SO2: {
        assert(dims[i] == 1);
        d += SO2<Scalar>::distance(location1.template segment<1>(counter),
                                   location2.template segment<1>(counter));
      } break;
      case DistanceType::SO2Squared: {
        assert(dims[i] == 1);
        d += SO2Squared<Scalar>::distance(
            location1.template segment<1>(counter),
            location2.template segment<1>(counter));
      } break;
      case DistanceType::SO3: {
        assert(dims[i] == 4);
        d += SO3<Scalar>::distance(location1.template segment<4>(counter),
                                   location2.template segment<4>(counter));
      } break;
      case DistanceType::SO3Squared: {
        assert(dims[i] == 4);
        d += SO3Squared<Scalar>::distance(
            location1.template segment<4>(counter),
            location2.template segment<4>(counter));
      } break;
      default:
        assert(false);
      }
      counter += dims[i];
    }
    return d;
  }

  Scalar distanceToRect(cref_t &location1, cref_t &lb, cref_t &ub) const {

    double d = 0;
    int counter = 0;
    for (size_t i = 0; i < dims.size(); i++) {
      switch (types[i]) {
      case DistanceType::L1: {
        d += L1<Scalar, -1>::distanceToRect(
            location1.template segment(counter, dims[i]),
            lb.template segment(counter, dims[i]),
            ub.template segment(counter, dims[i]));
      } break;
      case DistanceType::L2: {
        d += L2<Scalar, -1>::distanceToRect(
            location1.template segment(counter, dims[i]),
            lb.template segment(counter, dims[i]),
            ub.template segment(counter, dims[i]));
      } break;
      case DistanceType::L2Squared: {
        d += L2Squared<Scalar, -1>::distanceToRect(
            location1.template segment(counter, dims[i]),
            lb.template segment(counter, dims[i]),
            ub.template segment(counter, dims[i]));
      } break;
      case DistanceType::SO2: {
        assert(dims[i] == 1);
        d += SO2<Scalar>::distanceToRect(location1.template segment<1>(counter),
                                         lb.template segment<1>(counter),
                                         ub.template segment<1>(counter));
      } break;
      case DistanceType::SO2Squared: {
        assert(dims[i] == 1);
        d += SO2Squared<Scalar>::distanceToRect(
            location1.template segment<1>(counter),
            lb.template segment<1>(counter), ub.template segment<1>(counter));
      } break;
      case DistanceType::SO3: {
        assert(dims[i] == 4);
        d += SO3<Scalar>::distanceToRect(location1.template segment<4>(counter),
                                         lb.template segment<4>(counter),
                                         ub.template segment<4>(counter));
      } break;
      case DistanceType::SO3Squared: {
        assert(dims[i] == 4);
        d += SO3Squared<Scalar>::distanceToRect(
            location1.template segment<4>(counter),
            lb.template segment<4>(counter), ub.template segment<4>(counter));
      } break;
      default:
        assert(false);
      }
      counter += dims[i];
    }
    return d;
  }

protected:
  std::vector<DistanceType> types;
  std::vector<int> dims;
};

template <class Payload, int Dimensions, std::size_t BucketSize = 32,
          typename Scalar = double,
          typename Distance = L2Squared<Scalar, Dimensions>>
class KDTree {
private:
  struct Node;
  std::vector<Node> m_nodes;
  std::set<std::size_t> waitingForSplit;
  Distance distance_fun;

public:
  using distance_t = Distance;
  using scalar_t = Scalar;
  using payload_t = Payload;
  using point_t = Eigen::Matrix<Scalar, Dimensions, 1>;
  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;
  int m_dimensions = Dimensions;
  static const std::size_t bucketSize = BucketSize;
  // TODO: I also want Dimensions at runtime!!
  using tree_t = KDTree<Payload, Dimensions, BucketSize, Scalar, Distance>;

  void interpolate(cref_t from, cref_t to, Scalar t,
                   ref_t out) const {
    distance_fun.interpolate(from, to, t, out);
    std::cout << "out is 2 " << out.transpose() << std::endl;
  }

  Distance &getDistanceFun() { return distance_fun; }

  KDTree(int runtime_dimension = -1, Distance distance_fun = Distance())
      : distance_fun(distance_fun) {

    if constexpr (Dimensions == Eigen::Dynamic) {
      assert(runtime_dimension > 0);
      m_dimensions = runtime_dimension;
      m_nodes.emplace_back(BucketSize, m_dimensions);
    } else {
      m_nodes.emplace_back(BucketSize, -1);
    }
  }

  size_t size() const { return m_nodes[0].m_entries; }

  void addPoint(const point_t &location, const Payload &payload,
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
    m_nodes[addNode].add(LocationPayload{location, payload});

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

  struct DistancePayload {
    Scalar distance;
    Payload payload;
    bool operator<(const DistancePayload &dp) const {
      return distance < dp.distance;
    }
  };

  std::vector<DistancePayload> searchKnn(const point_t &location,
                                         std::size_t maxPoints) const {
    return searcher().search(location, std::numeric_limits<Scalar>::max(),
                             maxPoints, distance_fun);
  }

  std::vector<DistancePayload> searchBall(const point_t &location,
                                          Scalar maxRadius) const {
    return searcher().search(location, maxRadius,
                             std::numeric_limits<std::size_t>::max(),
                             distance_fun);
  }

  std::vector<DistancePayload>
  searchCapacityLimitedBall(const point_t &location, Scalar maxRadius,
                            std::size_t maxPoints) const {
    return searcher().search(location, maxRadius, maxPoints, distance_fun);
  }

  DistancePayload search(const point_t &location) const {
    DistancePayload result;
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
        if (result.distance > node.pointRectDist(location, distance_fun)) {
          if (node.m_splitDimension == m_dimensions) {
            for (const auto &lp : node.m_locationPayloads) {
              Scalar nodeDist = distance_fun.distance(location, lp.location);
              if (nodeDist < result.distance) {
                result = DistancePayload{nodeDist, lp.payload};
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
    const std::vector<DistancePayload> &search(const point_t &location,
                                               Scalar maxRadius,
                                               std::size_t maxPoints,
                                               Distance distance_fun) {
      // clear results from last time
      m_results.clear();

      // reserve capacities
      m_searchStack.reserve(
          1 + std::size_t(1.5 * std::log2(1 + m_tree.m_nodes[0].m_entries /
                                                  BucketSize)));
      if (m_prioqueueCapacity < maxPoints &&
          maxPoints < m_tree.m_nodes[0].m_entries) {
        std::vector<DistancePayload> container;
        container.reserve(maxPoints);
        m_prioqueue =
            std::priority_queue<DistancePayload, std::vector<DistancePayload>>(
                std::less<DistancePayload>(), std::move(container));
        m_prioqueueCapacity = maxPoints;
      }

      m_tree.searchCapacityLimitedBall(location, maxRadius, maxPoints,
                                       m_searchStack, m_prioqueue, m_results,
                                       distance_fun);

      m_prioqueueCapacity = std::max(m_prioqueueCapacity, m_results.size());
      return m_results;
    }

  private:
    const tree_t &m_tree;

    std::vector<std::size_t> m_searchStack;
    std::priority_queue<DistancePayload, std::vector<DistancePayload>>
        m_prioqueue;
    std::size_t m_prioqueueCapacity = 0;
    std::vector<DistancePayload> m_results;
  };

  // NB! returned class has no const methods. Get one instance per thread!
  Searcher searcher() const { return Searcher(*this); }

private:
  struct LocationPayload {
    point_t location;
    Payload payload;
  };
  std::vector<LocationPayload> m_bucketRecycle;

  void searchCapacityLimitedBall(
      const point_t &location, Scalar maxRadius, std::size_t maxPoints,
      std::vector<std::size_t> &searchStack,
      std::priority_queue<DistancePayload, std::vector<DistancePayload>>
          &prioqueue,
      std::vector<DistancePayload> &results, Distance distance_fun) const {
    std::size_t numSearchPoints = std::min(maxPoints, m_nodes[0].m_entries);

    if (numSearchPoints > 0) {
      searchStack.push_back(0);
      while (searchStack.size() > 0) {
        std::size_t nodeIndex = searchStack.back();
        searchStack.pop_back();
        const Node &node = m_nodes[nodeIndex];
        Scalar minDist = node.pointRectDist(location, distance_fun);
        if (maxRadius > minDist && (prioqueue.size() < numSearchPoints ||
                                    prioqueue.top().distance > minDist)) {
          if (node.m_splitDimension == m_dimensions) {
            node.searchCapacityLimitedBall(location, maxRadius, numSearchPoints,
                                           prioqueue, distance_fun);
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
    distance_fun.choose_split_dimension(splitNode.m_lb, splitNode.m_ub,
                                        splitNode.m_splitDimension, width);

    if (splitNode.m_splitDimension == m_dimensions) {
      return false;
    }

    std::vector<Scalar> splitDimVals;
    splitDimVals.reserve(splitNode.m_entries);
    for (const auto &lp : splitNode.m_locationPayloads) {
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

    for (const auto &lp : splitNode.m_locationPayloads) {
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
      std::swap(rightNode.m_locationPayloads, m_bucketRecycle);
      m_nodes.pop_back();
      m_nodes.pop_back();
      return false;
    } else {
      splitNode.m_locationPayloads.clear();
      // if it was a standard sized bucket, recycle the memory to reduce
      // allocator pressure otherwise clear the memory used by the bucket
      // since it is a branch not a leaf anymore
      if (splitNode.m_locationPayloads.capacity() == BucketSize) {
        std::swap(splitNode.m_locationPayloads, m_bucketRecycle);
      } else {
        std::vector<LocationPayload> empty;
        std::swap(splitNode.m_locationPayloads, empty);
      }
      return true;
    }
  }

  struct Node {
    Node(std::size_t capacity, int runtime_dimension = -1) {
      init(capacity, runtime_dimension);
    }

    Node(std::vector<LocationPayload> &recycle, std::size_t capacity,
         int runtime_dimension) {
      std::swap(m_locationPayloads, recycle);
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
      m_locationPayloads.reserve(std::max(BucketSize, capacity));
    }

    void expandBounds(const point_t &location) {
      m_lb = m_lb.cwiseMin(location);
      m_ub = m_ub.cwiseMax(location);
      m_entries++;
    }

    void add(const LocationPayload &lp) {
      expandBounds(lp.location);
      m_locationPayloads.push_back(lp);
    }

    bool shouldSplit() const { return m_entries >= BucketSize; }

    void
    searchCapacityLimitedBall(const point_t &location, Scalar maxRadius,
                              std::size_t K,
                              std::priority_queue<DistancePayload> &results,
                              Distance distance_fun) const {

      std::size_t i = 0;

      // this fills up the queue if it isn't full yet
      for (; results.size() < K && i < m_entries; i++) {
        const auto &lp = m_locationPayloads[i];
        Scalar distance = distance_fun.distance(location, lp.location);
        if (distance < maxRadius) {
          results.emplace(DistancePayload{distance, lp.payload});
        }
      }

      // this adds new things to the queue once it is full
      for (; i < m_entries; i++) {
        const auto &lp = m_locationPayloads[i];
        Scalar distance = distance_fun.distance(location, lp.location);
        if (distance < maxRadius && distance < results.top().distance) {
          results.pop();
          results.emplace(DistancePayload{distance, lp.payload});
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

    Scalar pointRectDist(const point_t &location, Distance distance) const {
      return distance.distanceToRect(location, m_lb, m_ub);
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
    std::vector<LocationPayload>
        m_locationPayloads; /// data held in this node (if a leaf)
  };
};

// template class myclass<int>;
//
//
// template
//
// class

// template <class Payload, std::size_t Dimensions, std::size_t BucketSize =
// 32,
//           class Distance = L2Squared, typename Scalar = double>

template class KDTree<size_t, 4, 32, double, SO3Squared<double>>;
template class KDTree<size_t, 2, 32>;
template class KDTree<size_t, -1, 32>;
// TODO:...
// Create a Python Interface :)

} // namespace tree
} // namespace jk

//
//
//
// <class Payload, std::size_t Dimensions, std::size_t BucketSize = 32,
//           class Distance = L2Squared, typename Scalar = double>
//
