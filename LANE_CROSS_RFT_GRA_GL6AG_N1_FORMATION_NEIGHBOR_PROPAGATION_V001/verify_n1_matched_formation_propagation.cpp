#include <gmpxx.h>

#include <array>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using Integer = mpz_class;
using Vector = std::vector<Integer>;

namespace {

constexpr int kLinks = 16;
constexpr int kDimension = 1 << kLinks;
constexpr int kMaxOrder = 16;
constexpr int kPairs[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
constexpr int kE[6][2] = {
    {1, 1}, {-1, 0}, {0, -1}, {0, -1}, {-1, 0}, {1, 1}};

long choose(int n, int k) {
  long value = 1;
  for (int j = 1; j <= k; ++j) value = value * (n - k + j) / j;
  return value;
}

struct Replay {
  std::vector<int> diagonal;
  std::array<bool, kLinks> formed{};
  std::vector<Vector> blank_powers;

  Replay(const std::array<int, 4>& source_pattern, bool shared_bridges)
      : diagonal(kDimension) {
    for (int a = 0; a < 4; ++a) formed[a] = source_pattern[a] != 0;
    for (int cell = 1; cell < 4; ++cell) {
      for (int a = 0; a < 4; ++a) formed[4 * cell + a] = true;
    }

    std::vector<std::pair<int, int>> interactions;
    for (int cell = 0; cell < 4; ++cell) {
      for (int pair = 0; pair < 6; ++pair) {
        interactions.push_back(
            {4 * cell + kPairs[pair][0], 4 * cell + kPairs[pair][1]});
      }
    }
    if (shared_bridges) {
      for (int cell = 0; cell < 4; ++cell) {
        for (int other = cell + 1; other < 4; ++other) {
          interactions.push_back({4 * cell + other, 4 * other + cell});
        }
      }
    }
    const std::size_t expected = shared_bridges ? 30 : 24;
    if (interactions.size() != expected) {
      throw std::runtime_error("interaction census");
    }

    for (int word = 0; word < kDimension; ++word) {
      for (const auto& edge : interactions) {
        diagonal[word] +=
            2 * ((word >> edge.first) & 1) * ((word >> edge.second) & 1);
      }
    }

    blank_powers.emplace_back(kDimension);
    blank_powers[0][0] = 1;
    for (int order = 0; order < kMaxOrder; ++order) {
      blank_powers.push_back(apply_h(blank_powers.back()));
    }
  }

  Vector apply_h(const Vector& in) const {
    Vector out(kDimension);
    for (int word = 0; word < kDimension; ++word) {
      out[word] += diagonal[word] * in[word];
      for (int q = 0; q < kLinks; ++q) {
        if (formed[q]) out[word] -= in[word ^ (1 << q)];
      }
    }
    return out;
  }

  std::array<Integer, 2> e_raw(int cell, int order) const {
    std::array<Integer, 2> answer{};
    for (int split = 0; split <= order; ++split) {
      std::array<Integer, 2> matrix_element{};
      const Vector& left = blank_powers[order - split];
      const Vector& right = blank_powers[split];
      for (int word = 0; word < kDimension; ++word) {
        std::array<int, 4> z{};
        for (int port = 0; port < 4; ++port) {
          z[port] = 1 - 2 * ((word >> (4 * cell + port)) & 1);
        }
        for (int pair = 0; pair < 6; ++pair) {
          const Integer amplitude = left[word] * right[word] *
              z[kPairs[pair][0]] * z[kPairs[pair][1]];
          for (int coordinate = 0; coordinate < 2; ++coordinate) {
            matrix_element[coordinate] +=
                kE[pair][coordinate] * amplitude;
          }
        }
      }
      Integer weight = choose(order, split);
      if (split & 1) weight = -weight;
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        answer[coordinate] += weight * matrix_element[coordinate];
      }
    }
    return answer;
  }
};

void require(bool condition, const char* message, int& checks) {
  if (!condition) throw std::runtime_error(message);
  ++checks;
}

}  // namespace

int main() {
  int checks = 0;
  using OrderTable =
      std::array<std::array<std::array<Integer, kMaxOrder + 1>, 2>, 4>;
  std::array<OrderTable, 16> tables{};

  auto replay_mask = [&](int mask, bool shared_bridges) {
    std::array<int, 4> pattern{};
    for (int port = 0; port < 4; ++port) pattern[port] = (mask >> port) & 1;
    Replay replay(pattern, shared_bridges);
    OrderTable table{};
    for (int order = 0; order <= kMaxOrder; ++order) {
      for (int cell = 0; cell < 4; ++cell) {
        const auto value = replay.e_raw(cell, order);
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          table[cell][coordinate][order] = value[coordinate];
        }
      }
    }
    return table;
  };

  tables[0] = replay_mask(0, true);
  for (int port = 0; port < 4; ++port) {
    tables[1 << port] = replay_mask(1 << port, true);
  }

  for (int a = 0; a < 4; ++a) {
    for (int b = a + 1; b < 4; ++b) {
      const int mask = (1 << a) | (1 << b);
      tables[mask] = replay_mask(mask, true);
      int pair_index = 0;
      while (kPairs[pair_index][0] != a || kPairs[pair_index][1] != b) {
        ++pair_index;
      }
      for (int order = 0; order < 4; ++order) {
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          require(tables[mask][0][coordinate][order] ==
                      tables[0][0][coordinate][order],
                  "local earlier coefficient", checks);
        }
      }
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        const Integer local = tables[mask][0][coordinate][4] -
                              tables[0][0][coordinate][4];
        require(local == 96 * kE[pair_index][coordinate],
                "local order-four pair coefficient", checks);
      }

      for (int cell = 1; cell < 4; ++cell) {
        int receiver_pair = 0;
        while (kPairs[receiver_pair][0] != 0 ||
               kPairs[receiver_pair][1] != cell) {
          ++receiver_pair;
        }
        for (int order = 0; order < 12; ++order) {
          for (int coordinate = 0; coordinate < 2; ++coordinate) {
            require(tables[mask][cell][coordinate][order] ==
                        tables[0][cell][coordinate][order],
                    "receiver earlier coefficient", checks);
          }
        }
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          const Integer remote = tables[mask][cell][coordinate][12] -
                                 tables[0][cell][coordinate][12];
          const Integer expected = -Integer(63371264) *
              ((mask >> cell) & 1) * kE[receiver_pair][coordinate];
          require(remote == expected,
                  "receiver order-twelve matched coefficient", checks);
          require(tables[mask][cell][coordinate][13] ==
                      tables[0][cell][coordinate][13],
                  "receiver odd order thirteen", checks);
        }

        for (int order = 0; order < 16; ++order) {
          for (int coordinate = 0; coordinate < 2; ++coordinate) {
            const Integer mobius = tables[mask][cell][coordinate][order] -
                tables[1 << a][cell][coordinate][order] -
                tables[1 << b][cell][coordinate][order] +
                tables[0][cell][coordinate][order];
            require(mobius == 0,
                    "remote pair Mobius earlier than sixteen", checks);
          }
        }
      }
    }
  }

  for (int mask = 1; mask < 16; ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) > 2) {
      tables[mask] = replay_mask(mask, true);
    }
  }

  for (int mask = 0; mask < 16; ++mask) {
    for (int cell = 1; cell < 4; ++cell) {
      int receiver_pair = 0;
      while (kPairs[receiver_pair][0] != 0 ||
             kPairs[receiver_pair][1] != cell) {
        ++receiver_pair;
      }
      for (int order = 0; order < 12; ++order) {
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          require(tables[mask][cell][coordinate][order] ==
                      tables[0][cell][coordinate][order],
                  "all-pattern receiver earlier coefficient", checks);
        }
      }
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        const Integer remote = tables[mask][cell][coordinate][12] -
                               tables[0][cell][coordinate][12];
        const Integer expected = -Integer(63371264) *
            ((mask >> cell) & 1) * kE[receiver_pair][coordinate];
        require(remote == expected,
                "all-pattern receiver order-twelve coefficient", checks);
        require(tables[mask][cell][coordinate][13] ==
                    tables[0][cell][coordinate][13],
                "all-pattern receiver odd order thirteen", checks);
      }
    }
  }

  for (int port = 0; port < 4; ++port) {
    const int mask = 1 << port;
    for (int cell = 1; cell < 4; ++cell) {
      int receiver_pair = 0;
      while (kPairs[receiver_pair][0] != 0 ||
             kPairs[receiver_pair][1] != cell) {
        ++receiver_pair;
      }
      for (int order = 0; order < 12; ++order) {
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          require(tables[mask][cell][coordinate][order] ==
                      tables[0][cell][coordinate][order],
                  "single receiver earlier coefficient", checks);
        }
      }
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        const Integer remote = tables[mask][cell][coordinate][12] -
                               tables[0][cell][coordinate][12];
        const Integer expected = -Integer(63371264) *
            (port == cell) * kE[receiver_pair][coordinate];
        require(remote == expected,
                "single receiver order-twelve coefficient", checks);
      }
    }
  }

  bool found_pair_mobius_at_sixteen = false;
  for (int a = 0; a < 4; ++a) {
    for (int b = a + 1; b < 4; ++b) {
      const int mask = (1 << a) | (1 << b);
      for (int cell = 1; cell < 4; ++cell) {
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          const Integer mobius = tables[mask][cell][coordinate][16] -
              tables[1 << a][cell][coordinate][16] -
              tables[1 << b][cell][coordinate][16] +
              tables[0][cell][coordinate][16];
          int receiver_pair = 0;
          while (kPairs[receiver_pair][0] != 0 ||
                 kPairs[receiver_pair][1] != cell) {
            ++receiver_pair;
          }
          const Integer expected = Integer(123422773248) *
              ((mask >> cell) & 1) * kE[receiver_pair][coordinate];
          require(mobius == expected,
                  "remote pair Mobius order-sixteen coefficient", checks);
          if (mobius != 0) found_pair_mobius_at_sixteen = true;
        }
      }
    }
  }
  require(found_pair_mobius_at_sixteen,
          "a genuine remote pair Mobius coefficient opens at sixteen",
          checks);

  std::array<OrderTable, 16> bridge_off{};
  bridge_off[0] = replay_mask(0, false);
  for (int port = 0; port < 4; ++port) {
    bridge_off[1 << port] = replay_mask(1 << port, false);
  }
  for (int a = 0; a < 4; ++a) {
    for (int b = a + 1; b < 4; ++b) {
      const int mask = (1 << a) | (1 << b);
      bridge_off[mask] = replay_mask(mask, false);
    }
  }
  for (int mask = 0; mask < 16; ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) > 2) continue;
    for (int cell = 1; cell < 4; ++cell) {
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        for (int order = 0; order <= kMaxOrder; ++order) {
          require(bridge_off[mask][cell][coordinate][order] ==
                      bridge_off[0][cell][coordinate][order],
                  "bridge-off matched receiver coefficient", checks);
        }
      }
    }
  }

  std::cout << "PASS GL6AG exact matched-formation checks " << checks << "/"
            << checks << "\n";
  return 0;
}
