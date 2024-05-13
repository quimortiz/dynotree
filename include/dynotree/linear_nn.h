
#include "KDTree.h"

namespace dynotree {
template <class Id, int Dimensions, typename Scalar = double,
          typename StateSpace = Rn<Scalar, Dimensions>>
class LinearKNN {

public:
  using scalar_t = Scalar;
  using id_t = Id;
  using point_t = Eigen::Matrix<Scalar, Dimensions, 1>;
  using cref_t = const Eigen::Ref<const Eigen::Matrix<Scalar, Dimensions, 1>> &;
  using ref_t = Eigen::Ref<Eigen::Matrix<Scalar, Dimensions, 1>>;
  int m_dimensions = Dimensions;

  StateSpace &getStateSpace() { return state_space; }

  LinearKNN(int runtime_dimension = -1,
            const StateSpace &state_space = StateSpace())
      : state_space(state_space) {

    if constexpr (Dimensions == Eigen::Dynamic) {
      assert(runtime_dimension > 0);
      m_dimensions = runtime_dimension;
    }
  }

  size_t size() const { return m_points.size(); }

  void addPoint(const point_t &x, const Id &id, bool dummy) {
    (void)dummy; // unused, only to match interface of kdtree

    std::size_t addNode = 0;
    m_points.emplace_back(PointId{x, id});
  }

  struct DistanceId {
    Scalar distance;
    Id id;
    inline bool operator<(const DistanceId &dp) const {
      return distance < dp.distance;
    }
  };

  std::vector<DistanceId> searchKnn(const point_t &x,
                                    std::size_t maxPoints) const {

    std::priority_queue<DistanceId> max_heap;

    for (const auto &point : m_points) {
      Scalar distance = state_space.distance(x, point.point);
      if (max_heap.size() < maxPoints) {
        max_heap.push({distance, point.id});
      } else if (distance < max_heap.top().distance) {
        max_heap.pop();
        max_heap.push({distance, point.id});
      }
    }

    std::vector<DistanceId> neighbors;
    while (!max_heap.empty()) {
      neighbors.push_back(max_heap.top());
      max_heap.pop();
    }

    std::reverse(neighbors.begin(), neighbors.end());

    return neighbors;
  }

  std::vector<DistanceId> searchBall(const point_t &x, Scalar maxRadius) const {

    std::vector<DistanceId> out;

    for (auto &point : m_points) {
      Scalar distance = state_space.distance(x, point.point);
      if (distance < maxRadius) {
        out.emplace_back(DistanceId{distance, point.id});
      }
    }
    return out;
  }

  DistanceId searchNN(const point_t &x) const {

    double min_dist = std::numeric_limits<double>::max();
    DistanceId out;
    out.distance = std::numeric_limits<double>::max();
    for (auto &point : m_points) {
      Scalar distance = state_space.distance(x, point.point);
      if (distance < out.distance) {
        out.id = point.id;
        out.distance = distance;
      }
    }
    return out;
  }

private:
  struct PointId {
    point_t point;
    Id id;
  };

  std::vector<PointId> m_points;
  StateSpace state_space;
};
} // namespace dynotree
