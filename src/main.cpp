#include "nigh/nigh_forward.hpp"
#define BOOST_TEST_MODULE test_0
#define BOOST_TEST_DYN_LINK
#include "ompl/base/ScopedState.h"
#include "ompl/base/spaces/SO3StateSpace.h"
#include "ompl/datastructures/NearestNeighborsGNAT.h"
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <iostream>

#include "dynkdtree/KDTree.h"

// Include the header file for the space that we use.  Euclidean space
// has an L^2 metric, so we include a space file for general L^p
// metric spaces, which includes the L2Space alias.
#include <nigh/lp_space.hpp>
#include <nigh/so3_space.hpp>

// Nigh support several nearest neighbor strategies, we will use the
// one that benchmarks the fastest and supports concurrent queries and
// inserts.
#include <nigh/kdtree_batch.hpp>

// Our actual points use vectors from the Eigen matrix library, since
// include the necessary header for it.
#include <Eigen/Dense>

// int main() {

// TODO: all the tests from the original code.
// Try SO2.
// Try SO3 -- just using the rectangle seems enough?

BOOST_AUTO_TEST_CASE(t_hello) {

  {
    using tree_t = jk::tree::KDTree<std::string, 2>;
    using point_t = Eigen::Vector2d;
    tree_t tree;
    tree.addPoint(point_t(Eigen::Vector2d(1, 2)), "George");
    tree.addPoint(point_t(Eigen::Vector2d(1, 3)), "Harold");
    tree.addPoint(point_t(Eigen::Vector2d(7, 7)), "Melvin");

    // KNN search
    point_t lazyMonsterLocation(Eigen::Vector2d(
        6, 6)); // this monster will always try to eat the closest people
    const std::size_t monsterHeads =
        2; // this monster can eat two people at once
    auto lazyMonsterVictims = tree.searchKnn(lazyMonsterLocation, monsterHeads);
    for (const auto &victim : lazyMonsterVictims) {
      std::cout << victim.payload << " closest to lazy monster, with distance "
                << sqrt(victim.distance) << "!" << std::endl;
    }
  }
}

// SO2

BOOST_AUTO_TEST_CASE(t_hello2) {
  using tree_t =
      jk::tree::KDTree<std::string, 1, 32, double, jk::tree::SO2<double>>;

  using point_t = Eigen::VectorXd;
  using V1d = Eigen::Matrix<double, 1, 1>;
  tree_t tree;
  tree.addPoint(point_t(V1d(0)), "George");
  tree.addPoint(point_t(V1d(3.)), "Harold");
  tree.addPoint(point_t(V1d(-3.)), "Melvin");

  // KNN search
  point_t lazyMonsterLocation(V1d(3.1));
  // 6, 6)); // this monster will always try to eat the closest people
  const std::size_t monsterHeads = 2; // this monster can eat two people at once
  auto lazyMonsterVictims = tree.searchKnn(lazyMonsterLocation, monsterHeads);
  for (const auto &victim : lazyMonsterVictims) {
    std::cout << victim.payload << " closest to lazy monster, with distance "
              << victim.distance << "!" << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(t_hello3)

{
  std::srand(0);
  using tree_t = jk::tree::KDTree<size_t, 1, 4, double, jk::tree::SO2<double>>;

  using point_t = Eigen::VectorXd;
  using V1d = Eigen::Matrix<double, 1, 1>;
  tree_t tree;

  std::vector<point_t> points;
  size_t N = 100;
  for (size_t ii = 0; ii < N; ++ii) {
    double i = M_PI * (-1 + 2 * static_cast<double>(rand()) /
                                static_cast<double>(RAND_MAX));
    points.push_back(point_t(V1d(i)));
    tree.addPoint(point_t(V1d(i)), ii);
  }
  double radius = 1;

  for (size_t jj = 0; jj < 100; ++jj) {
    double j = M_PI * (-1 + 2 * static_cast<double>(rand()) /
                                static_cast<double>(RAND_MAX));
    point_t lazyMonsterLocation((V1d(j)));
    auto lazyMonsterVictims = tree.searchBall(lazyMonsterLocation, radius);
    double min_distance = std::numeric_limits<double>::max();
    int counter = 0;
    for (size_t ii = 0; ii < points.size(); ++ii) {
      double d =
          tree.getDistanceFun().distance(points[ii], lazyMonsterLocation);
      if (d < radius) {
        counter++;
      }
    }

    std::cout << "counter/size: " << counter << " " << lazyMonsterVictims.size()
              << std::endl;
    if (counter != lazyMonsterVictims.size()) {
      std::cout << "Error: " << counter << " " << lazyMonsterVictims.size()
                << std::endl;
      throw std::runtime_error("Error");
    }
  }
}

BOOST_AUTO_TEST_CASE(t_hello4) {
  std::srand(0);
  using tree_t =
      jk::tree::KDTree<size_t, 3, 4, double, jk::tree::R2SO2Squared<double>>;

  using point_t = Eigen::Vector3d;
  tree_t tree;

  std::vector<point_t> points;
  size_t N = 100;
  for (size_t ii = 0; ii < N; ++ii) {
    point_t v = M_PI * point_t::Random();
    points.push_back(point_t(v));
    tree.addPoint(point_t(v), ii);
  }

  double radius = 1.;

  for (size_t jj = 0; jj < 100; ++jj) {

    point_t vv = M_PI * point_t::Random();
    auto lazyMonsterVictims = tree.searchBall(vv, radius);
    double min_distance = std::numeric_limits<double>::max();
    int counter = 0;
    for (size_t ii = 0; ii < points.size(); ++ii) {
      double d = tree.getDistanceFun().distance(points[ii], vv);
      // double d = 0;
      if (d < radius) {
        counter++;
      }
    }

    std::cout << "counter/size: " << counter << " " << lazyMonsterVictims.size()
              << std::endl;
    if (counter != lazyMonsterVictims.size()) {
      std::cout << "Error: " << counter << " " << lazyMonsterVictims.size()
                << std::endl;
      throw std::runtime_error("Error");
    }
  }
}

BOOST_AUTO_TEST_CASE(bench_run_vs_compile) {
  std::srand(0);
  std::cout << "benchmark in c++" << std::endl;

  using TreeRX = jk::tree::KDTree<int, -1>;
  using TreeR4 = jk::tree::KDTree<int, 4>;

  TreeRX treex(4);
  TreeR4 tree4(-1);

  Eigen::MatrixXd X = Eigen::MatrixXd::Random(10000, 4);

  for (size_t i = 0; i < X.rows(); ++i) {
    treex.addPoint(X.row(i), i);
    tree4.addPoint(X.row(i).head<4>(), i);
  }
  int num_neighs = 10;

  Eigen::Vector4d x = Eigen::Vector4d::Random();

  {
    auto t0 = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < 10; i++)
      tree4.searchBall(x, .5);
    // treex.searchKnn(x, num_neighs);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "dt X:" << dt << std::endl;
  }

  {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 10; i++)
      tree4.searchBall(x, .5);
    // tree4.searchKnn(x, num_neighs);

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "dt 4:" << dt << std::endl;
  }
}

// incremental benchmark

BOOST_AUTO_TEST_CASE(bench_run_vs_compile2) {
  std::srand(0);
  std::cout << "benchmark in c++" << std::endl;

  using TreeRX = jk::tree::KDTree<int, -1>;
  using TreeR4 = jk::tree::KDTree<int, 4>;

  TreeRX treex(4);
  TreeR4 tree4(-1);

  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::Vector4d x4 = x;

  auto X = Eigen::MatrixXd::Random(4, 10000);
  int num_neighs = 10;

  {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < X.cols(); ++i) {
      treex.searchKnn(X.col(i), num_neighs);
      treex.addPoint(X.col(i), i);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "dt X:" << dt << std::endl;
  }

  {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < X.cols(); ++i) {
      tree4.searchKnn(X.col(i), num_neighs);
      tree4.addPoint(X.col(i), i);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "dt X:" << dt << std::endl;
  }
}

struct MyNode {
  int id;
  Eigen::Vector4d point_;

  MyNode(const int &name, const Eigen::Vector4d &pt) : id(name), point_(pt) {}
};

struct MyNodeKey {
  // The functor should return a constant reference to the data
  // member or a copy of it.
  const Eigen::Vector4d &operator()(const MyNode &node) const {
    return node.point_;
  }
};

struct MyNodeQuat {
  int id;
  Eigen::Quaterniond point_;
  MyNodeQuat(const int &name, const Eigen::Quaterniond &pt)
      : id(name), point_(pt) {}
};

struct MyNodeKeyQuat {
  // The functor should return a constant reference to the data
  // member or a copy of it.
  const Eigen::Quaterniond &operator()(const MyNodeQuat &node) const {
    return node.point_;
  }
};

using namespace unc::robotics;

BOOST_AUTO_TEST_CASE(t_against_nigh) {

  nigh::Nigh<MyNode, nigh::L2Space<double, 4>, MyNodeKey, nigh::NoThreadSafety,
             nigh::KDTreeBatch<32>>
      nn;

  std::srand(0);
  using TreeR4 = jk::tree::KDTree<int, 4>;

  TreeR4 tree4(-1);

  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::Vector4d x4 = x;

  Eigen::MatrixXd X = Eigen::MatrixXd::Random(4, 100000);
  int num_neighs = 10;

  for (size_t i = 0; i < X.cols(); ++i) {
    nn.insert(MyNode(i, X.col(i)));
    tree4.addPoint(X.col(i), i);
  }

  int best = -1;
  {

    std::cout << "linear nn" << std::endl;

    double dist = std::numeric_limits<double>::max();
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < X.cols(); ++i) {
      double d = (X.col(i) - x).norm();
      if (d < dist) {
        dist = d;
        best = i;
      }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    std::cout << " linear time " << dt << std::endl;
    std::cout << "best is " << best << " with dist " << dist << "dist2 "
              << dist * dist << std::endl;
  }

  int num_experiments = 1;
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    int num_neighs = 10;
    std::vector<std::pair<MyNode, double>> nbh;
    for (size_t i = 0; i < num_experiments; i++) {
      nn.nearest(nbh, x4, num_neighs, 1e8);
    }
    // tree4.searchBall(x, .5);
    // treex.searchKnn(x, num_neighs);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    for (size_t j = 0; j < nbh.size(); ++j) {
      std::cout << nbh[j].first.id << " " << nbh[j].second << std::endl;
    }

    std::cout << "dt nigh:" << dt << std::endl;
  }
  std::vector<TreeR4::DistancePayload> nnt;
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_experiments; i++) {
      nnt = tree4.searchKnn(x, num_neighs);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    for (size_t j = 0; j < nnt.size(); ++j) {
      std::cout << nnt[j].payload << " " << nnt[j].distance << std::endl;
    }
    std::cout << "dt tree:" << dt << std::endl;
  }

  BOOST_TEST(best, nnt[0].payload);

  {

    using TreeR4 = jk::tree::KDTree<int, 4>;
    TreeR4 tree4(-1);
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < X.cols(); ++i) {
      auto out = tree4.searchKnn(X.col(i), num_neighs);
      tree4.addPoint(X.col(i), i);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "dt nigh inc:" << dt << std::endl;
  }
  {

    using TreeR4 = jk::tree::KDTree<int, 4>;
    TreeR4 tree4(-1);
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<MyNode, double>> nbh;
    for (size_t i = 0; i < X.cols(); ++i) {
      nn.insert(MyNode(i, X.col(i)));
      nn.nearest(nbh, x4, num_neighs, 1e8);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "dt tree inc:" << dt << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(t_against_nigh_so3) {
  // compare the two SO3 Implementations

  nigh::Nigh<MyNodeQuat, nigh::SO3Space<double>, MyNodeKeyQuat,
             nigh::NoThreadSafety, nigh::KDTreeBatch<32>>
      nn;

  std::srand(0);
  using TreeR4 = jk::tree::KDTree<int, 4, 32, double, jk::tree::SO3<double>>;

  TreeR4 tree4(-1);

  Eigen::VectorXd x = Eigen::VectorXd::Random(4);
  Eigen::Vector4d x4 = x;
  x4.normalize();

  Eigen::MatrixXd X = Eigen::MatrixXd::Random(4, 100000);

  bool real_part_positive = true;
  int index_real_part = 3;

  if (real_part_positive) {
    if (x4(index_real_part) < 0) {
      x4 *= -1;
    }
  }

  for (size_t i = 0; i < X.cols(); ++i) {
    X.col(i).normalize();
    if (real_part_positive) {
      if (X.col(i)(index_real_part) < 0) {
        X.col(i) *= -1;
      }
    }
  }

  int num_neighs = 10;

  for (size_t i = 0; i < X.cols(); ++i) {
    Eigen::Vector4d q = X.col(i);
    nn.insert(MyNodeQuat(i, Eigen::Quaterniond(q)));
    tree4.addPoint(q, i);
  }

  int best = -1;
  {

    std::cout << "linear nn" << std::endl;

    double dist = std::numeric_limits<double>::max();
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < X.cols(); ++i) {
      double d = std::min((X.col(i) - x).norm(), (X.col(i) + x).norm());
      if (d < dist) {
        dist = d;
        best = i;
      }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    std::cout << " linear time " << dt << std::endl;
    std::cout << "best is " << best << " with dist " << dist << "dist2 "
              << dist * dist << std::endl;
  }

  std::cout << "quaternion" << std::endl;
  int num_experiments = 1;
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<MyNodeQuat, double>> nbh;
    for (size_t i = 0; i < num_experiments; i++) {
      nn.nearest(nbh, x4, num_neighs, 1e8);
    }
    // tree4.searchBall(x, .5);
    // treex.searchKnn(x, num_neighs);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    for (size_t j = 0; j < nbh.size(); ++j) {
      std::cout << nbh[j].first.id << " "
                << nbh[j].first.point_.coeffs().transpose() << " "
                << nbh[j].second << std::endl;
    }
    std::cout << "dt nigh:" << dt << std::endl;
  }
  std::vector<TreeR4::DistancePayload> nnt;
  {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_experiments; i++) {
      nnt = tree4.searchKnn(x4, num_neighs);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    for (size_t j = 0; j < nnt.size(); ++j) {
      std::cout << nnt[j].payload << " " << nnt[j].distance << std::endl;
    }

    std::cout << "dt tree:" << dt << std::endl;
  }

  std::cout << "best: " << best << " " << nnt[0].payload << std::endl;
  BOOST_TEST(best == nnt[0].payload);

  // ompl
  {
    ompl::base::StateSpacePtr space(new ompl::base::SO3StateSpace());

    std::vector<ompl::base::State *> states;

    for (size_t i = 0; i < X.cols(); ++i) {
      ompl::base::State *state = space->allocState();

      auto ptr = state->as<ompl::base::SO3StateSpace::StateType>();
      ptr->x = X.col(i)(0);
      ptr->y = X.col(i)(1);
      ptr->z = X.col(i)(2);
      ptr->w = X.col(i)(3);
      states.push_back(state);
    }

    // state->setY(0.2);
    // state->setYaw(0.0);

    ompl::NearestNeighbors<ompl::base::State *> *tt =
        new ompl::NearestNeighborsGNAT<ompl::base::State *>();

    tt->setDistanceFunction(
        [&](auto &a, auto &b) { return space->distance(a, b); });

    for (auto &s : states)
      tt->add(s);

    ompl::base::State *query = space->allocState();

    auto ptr = query->as<ompl::base::SO3StateSpace::StateType>();
    ptr->x = x4(0);
    ptr->y = x4(1);
    ptr->z = x4(2);
    ptr->w = x4(3);

    std::vector<ompl::base::State *> nbh;

    // free memory...

    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_experiments; i++) {
      tt->nearestK(query, num_neighs, nbh);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dt =
        1.e-9 *
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    std::cout << "nns " << std::endl;
    for (size_t j = 0; j < nbh.size(); ++j) {
      space->printState(nbh[j]);
    }
    //
    std::cout << "ompl tree:" << dt << std::endl;

    for (size_t i = 0; i < states.size(); ++i) {
      space->freeState(states[i]);
    }

    space->freeState(query);

    // t

    // ompl::base::SO3StateSpace space;
    //
    // using state_t = ompl::base::ScopedState<ompl::base::SO3StateSpace> ;
  }

  // TODO: try R3 x SO(3) and R9 x SO(3)
}
