

#include <iostream>

#include "KDTree.h"

int main() {

  // TODO: all the tests from the original code.
  // Try SO2.
  // Try SO3 -- just using the rectangle seems enough?

  std::cout << "Hello World!\n";

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

  // SO2

  {
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
    const std::size_t monsterHeads =
        2; // this monster can eat two people at once
    auto lazyMonsterVictims = tree.searchKnn(lazyMonsterLocation, monsterHeads);
    for (const auto &victim : lazyMonsterVictims) {
      std::cout << victim.payload << " closest to lazy monster, with distance "
                << victim.distance << "!" << std::endl;
    }
  }

  {
    std::srand(0);
    using tree_t =
        jk::tree::KDTree<size_t, 1, 4, double, jk::tree::SO2<double>>;

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
            jk::tree::SO2<double>::distance(points[ii], lazyMonsterLocation);
        if (d < radius) {
          counter++;
        }
      }

      std::cout << "counter/size: " << counter << " "
                << lazyMonsterVictims.size() << std::endl;
      if (counter != lazyMonsterVictims.size()) {
        std::cout << "Error: " << counter << " " << lazyMonsterVictims.size()
                  << std::endl;
        throw std::runtime_error("Error");
      }
    }
  }

  // SO3

#if 0
  {

    std::srand(0);
    using tree_t = jk::tree::KDTree<size_t, 4, 4, double, jk::tree::SO3<double>>;

    using point_t = Eigen::Vector4d;
    tree_t tree;

    std::vector<point_t> points;
    size_t N = 100;
    for (size_t ii = 0; ii < N; ++ii) {
      Eigen::Vector4d v = Eigen::Vector4d::Random().normalized();
      points.push_back(point_t(v));
      tree.addPoint(point_t(v), ii);
    }

    double radius = .2;

    for (size_t jj = 0; jj < 100; ++jj) {

      Eigen::Vector4d vv = Eigen::Vector4d::Random().normalized();
      point_t lazyMonsterLocation(vv);
      auto lazyMonsterVictims = tree.searchBall(lazyMonsterLocation, radius);
      double min_distance = std::numeric_limits<double>::max();
      int counter = 0;
      for (size_t ii = 0; ii < points.size(); ++ii) {
        double d = jk::tree::SO3<double>::distance(points[ii], lazyMonsterLocation);
        // double d = jk::tree::SquaredL2::distance(points[ii],
        // lazyMonsterLocation);
        if (d < radius) {
          counter++;
        }
      }

      std::cout << "counter/size: " << counter << " "
                << lazyMonsterVictims.size() << std::endl;
      if (counter != lazyMonsterVictims.size()) {
        std::cout << "Error: " << counter << " " << lazyMonsterVictims.size()
                  << std::endl;
        throw std::runtime_error("Error");
      }
    }
  }
#endif

  // SE2
  {
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
        double d =
            jk::tree::R2SO2Squared<double>::distance(points[ii], vv);
        // double d = 0;
        if (d < radius) {
          counter++;
        }
      }

      std::cout << "counter/size: " << counter << " "
                << lazyMonsterVictims.size() << std::endl;
      if (counter != lazyMonsterVictims.size()) {
        std::cout << "Error: " << counter << " " << lazyMonsterVictims.size()
                  << std::endl;
        throw std::runtime_error("Error");
      }
    }
  }
  return 0;
}
