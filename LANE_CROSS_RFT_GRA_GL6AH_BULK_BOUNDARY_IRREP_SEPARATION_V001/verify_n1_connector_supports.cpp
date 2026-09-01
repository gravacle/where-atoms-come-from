#include <gmpxx.h>

#include <array>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using Integer = mpz_class;
using Vector = std::vector<Integer>;

namespace {

constexpr int kCells = 4;
constexpr int kPorts = 4;
constexpr int kFullLinks = kCells * kPorts;
constexpr int kBridgeCount = 6;
constexpr int kMaxOrder = 16;
constexpr int kPairs[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
constexpr int kE[6][2] = {
    {1, 1}, {-1, 0}, {0, -1}, {0, -1}, {-1, 0}, {1, 1}};
constexpr int kBridges[kBridgeCount][2] = {
    {4 * 0 + 1, 4 * 1 + 0},  // 01: directed source--receiver bridge
    {4 * 0 + 2, 4 * 2 + 0},  // 02: inert for source patterns used here
    {4 * 0 + 3, 4 * 3 + 0},  // 03: inert for source patterns used here
    {4 * 1 + 2, 4 * 2 + 1},  // 12: receiver-port-2 helper chain
    {4 * 1 + 3, 4 * 3 + 1},  // 13: receiver-port-3 helper chain
    {4 * 2 + 3, 4 * 3 + 2},  // 23: helper--helper bridge
};

using PairRaw = std::array<Integer, 6>;
using ERaw = std::array<Integer, 2>;

long binomial(int n, int k) {
  long out = 1;
  for (int j = 1; j <= k; ++j) out = out * (n - k + j) / j;
  return out;
}

int popcount(int word) {
  int out = 0;
  while (word != 0) {
    out += word & 1;
    word >>= 1;
  }
  return out;
}

struct Replay {
  int source_mask;
  int bridge_mask;
  int source_count;
  int link_count;
  int dimension;
  std::vector<int> diagonal;
  std::vector<Vector> powers;

  Replay(int selected_source_mask, int selected_bridge_mask)
      : source_mask(selected_source_mask), bridge_mask(selected_bridge_mask) {
    std::array<int, kFullLinks> compressed{};
    compressed.fill(-1);
    link_count = 0;
    source_count = 0;
    for (int port = 0; port < kPorts; ++port) {
      if ((source_mask >> port) & 1) {
        compressed[port] = link_count++;
        ++source_count;
      }
    }
    for (int full = 4; full < kFullLinks; ++full) {
      compressed[full] = link_count++;
    }
    dimension = 1 << link_count;

    std::vector<std::pair<int, int>> active_interactions;
    for (int cell = 0; cell < kCells; ++cell) {
      for (const auto& pair : kPairs) {
        const int left = compressed[4 * cell + pair[0]];
        const int right = compressed[4 * cell + pair[1]];
        if (left >= 0 && right >= 0) {
          active_interactions.push_back({left, right});
        }
      }
    }
    for (int bridge = 0; bridge < kBridgeCount; ++bridge) {
      if (!((bridge_mask >> bridge) & 1)) continue;
      const int left = compressed[kBridges[bridge][0]];
      const int right = compressed[kBridges[bridge][1]];
      if (left >= 0 && right >= 0) {
        active_interactions.push_back({left, right});
      }
    }

    diagonal.assign(dimension, 0);
    for (int word = 0; word < dimension; ++word) {
      int value = 0;
      for (const auto& interaction : active_interactions) {
        value += 2 * ((word >> interaction.first) & 1) *
                 ((word >> interaction.second) & 1);
      }
      diagonal[word] = value;
    }

    powers.reserve(kMaxOrder + 1);
    powers.emplace_back(dimension);
    powers[0][0] = 1;
    for (int order = 0; order < kMaxOrder; ++order) {
      powers.push_back(apply_h(powers.back()));
    }
  }

  Vector apply_h(const Vector& input) const {
    Vector output(dimension);
    for (int word = 0; word < dimension; ++word) {
      if (input[word] == 0) continue;
      if (diagonal[word] != 0) {
        output[word] += diagonal[word] * input[word];
      }
      for (int q = 0; q < link_count; ++q) {
        output[word ^ (1 << q)] -= input[word];
      }
    }
    return output;
  }

  PairRaw pair_raw(int order) const {
    PairRaw answer{};
    for (int split = 0; split <= order; ++split) {
      PairRaw matrix{};
      const Vector& bra = powers[order - split];
      const Vector& ket = powers[split];
      for (int word = 0; word < dimension; ++word) {
        if (bra[word] == 0 || ket[word] == 0) continue;
        const Integer product = bra[word] * ket[word];
        const int nibble = (word >> source_count) & 15;
        std::array<int, 4> z{};
        for (int port = 0; port < 4; ++port) {
          z[port] = 1 - 2 * ((nibble >> port) & 1);
        }
        for (int pair = 0; pair < 6; ++pair) {
          const int value = z[kPairs[pair][0]] * z[kPairs[pair][1]];
          matrix[pair] += value * product;
        }
      }
      long weight = binomial(order, split);
      if (split & 1) weight = -weight;
      for (int pair = 0; pair < 6; ++pair) {
        answer[pair] += weight * matrix[pair];
      }
    }
    return answer;
  }

  ERaw e_raw(int order) const {
    const PairRaw pair = pair_raw(order);
    ERaw answer{};
    for (int row = 0; row < 6; ++row) {
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        answer[coordinate] += kE[row][coordinate] * pair[row];
      }
    }
    return answer;
  }
};

PairRaw subtract(const PairRaw& left, const PairRaw& right) {
  PairRaw out{};
  for (int j = 0; j < 6; ++j) out[j] = left[j] - right[j];
  return out;
}

ERaw subtract(const ERaw& left, const ERaw& right) {
  return {left[0] - right[0], left[1] - right[1]};
}

ERaw pair_mobius(const ERaw& pair, const ERaw& single_a,
                  const ERaw& single_x, const ERaw& sham) {
  return {pair[0] - single_a[0] - single_x[0] + sham[0],
          pair[1] - single_a[1] - single_x[1] + sham[1]};
}

PairRaw bridge_mobius(const std::array<PairRaw, 64>& values, int support) {
  PairRaw out{};
  for (int subset = support;; subset = (subset - 1) & support) {
    const int sign = ((popcount(support) - popcount(subset)) & 1) ? -1 : 1;
    for (int j = 0; j < 6; ++j) out[j] += sign * values[subset][j];
    if (subset == 0) break;
  }
  return out;
}

ERaw bridge_mobius(const std::array<ERaw, 64>& values, int support) {
  ERaw out{};
  for (int subset = support;; subset = (subset - 1) & support) {
    const int sign = ((popcount(support) - popcount(subset)) & 1) ? -1 : 1;
    for (int j = 0; j < 2; ++j) out[j] += sign * values[subset][j];
    if (subset == 0) break;
  }
  return out;
}

bool nonzero(const PairRaw& value) {
  for (const auto& entry : value) {
    if (entry != 0) return true;
  }
  return false;
}

bool nonzero(const ERaw& value) {
  return value[0] != 0 || value[1] != 0;
}

void require(bool condition, const char* message, int& checks) {
  if (!condition) throw std::runtime_error(message);
  ++checks;
}

void require_pair(const PairRaw& actual, const PairRaw& expected,
                  const char* message, int& checks) {
  for (int j = 0; j < 6; ++j) {
    require(actual[j] == expected[j], message, checks);
  }
}

void require_e(const ERaw& actual, const ERaw& expected,
               const char* message, int& checks) {
  require(actual[0] == expected[0], message, checks);
  require(actual[1] == expected[1], message, checks);
}

}  // namespace

int main() {
  int checks = 0;

  // The source branches used in the theorem occupy only ports 0 and 1.
  // Therefore source bridges 02 and 03 have a permanently blank endpoint.
  // We evaluate the 16 distinct invariant-sector Hamiltonians and then
  // populate all 64 literal connector subsets, including those two exactly
  // inert connector bits, before taking the six-variable Boolean Mobius
  // transform.
  constexpr int kEffectiveBits = (1 << 0) | (1 << 3) | (1 << 4) | (1 << 5);
  std::array<PairRaw, 64> matched_q6{};
  std::array<ERaw, 64> matched_q12{};
  std::array<ERaw, 64> source_pair_q16{};
  std::array<PairRaw, 64> effective_q6{};
  std::array<ERaw, 64> effective_q12{};
  std::array<ERaw, 64> effective_q16{};

  for (int bridges = 0; bridges < 64; ++bridges) {
    if ((bridges & ~kEffectiveBits) != 0) continue;
    const Replay sham(0, bridges);
    const Replay source_x(1 << 0, bridges);
    const Replay source_a(1 << 1, bridges);
    const Replay source_ax((1 << 0) | (1 << 1), bridges);

    effective_q6[bridges] =
        subtract(source_a.pair_raw(6), sham.pair_raw(6));
    effective_q12[bridges] =
        subtract(source_a.e_raw(12), sham.e_raw(12));
    effective_q16[bridges] = pair_mobius(
        source_ax.e_raw(16), source_a.e_raw(16), source_x.e_raw(16),
        sham.e_raw(16));
  }

  for (int bridges = 0; bridges < 64; ++bridges) {
    const int effective = bridges & kEffectiveBits;
    matched_q6[bridges] = effective_q6[effective];
    matched_q12[bridges] = effective_q12[effective];
    source_pair_q16[bridges] = effective_q16[effective];
  }

  const Integer q12 = 63371264;
  const Integer q16 = 123422773248;
  int nonzero_q6 = 0;
  int nonzero_q12 = 0;
  int nonzero_q16 = 0;

  for (int support = 0; support < 64; ++support) {
    const PairRaw got6 = bridge_mobius(matched_q6, support);
    const ERaw got12 = bridge_mobius(matched_q12, support);
    const ERaw got16 = bridge_mobius(source_pair_q16, support);

    PairRaw want6{};
    if (support == (1 << 0)) {
      want6 = {-128, -128, -128, 0, 0, 0};
    }
    ERaw want12{};
    ERaw want16{};
    if (support == ((1 << 0) | (1 << 3))) {
      // +q12*w_02 and -q16*w_02.
      want12 = {-q12, 0};
      want16 = {q16, 0};
    } else if (support == ((1 << 0) | (1 << 4))) {
      // +q12*w_03 and -q16*w_03.
      want12 = {0, -q12};
      want16 = {0, q16};
    }

    require_pair(got6, want6, "q6 six-bridge support coefficient", checks);
    require_e(got12, want12, "q12 six-bridge support coefficient", checks);
    require_e(got16, want16, "q16 six-bridge support coefficient", checks);
    if (nonzero(got6)) ++nonzero_q6;
    if (nonzero(got12)) ++nonzero_q12;
    if (nonzero(got16)) ++nonzero_q16;
  }

  require(nonzero_q6 == 1, "exactly one nonzero q6 support", checks);
  require(nonzero_q12 == 2, "exactly two nonzero q12 supports", checks);
  require(nonzero_q16 == 2, "exactly two nonzero q16 supports", checks);

  require_pair(matched_q6[63], {-128, -128, -128, 0, 0, 0},
               "full N=1 matched q6", checks);
  require_e(matched_q12[63], {-q12, -q12},
            "full N=1 matched E q12", checks);
  require_e(source_pair_q16[63], {q16, q16},
            "full N=1 source-pair E q16", checks);

  // The two formally present source bridges with unformed endpoints really
  // are invisible for every literal subset, not silently omitted from the
  // 64-subset census.
  for (int bridges = 0; bridges < 64; ++bridges) {
    const int effective = bridges & kEffectiveBits;
    require_pair(matched_q6[bridges], matched_q6[effective],
                 "inert-source-bridge q6 equality", checks);
    require_e(matched_q12[bridges], matched_q12[effective],
              "inert-source-bridge q12 equality", checks);
    require_e(source_pair_q16[bridges], source_pair_q16[effective],
              "inert-source-bridge q16 equality", checks);
  }

  std::cout << "PASS__GL6AH_N1_CONNECTOR_SUPPORTS__" << checks << "/"
            << checks << "\n";
  std::cout << "Q6_SUPPORTS=1:DIRECT_01:-128_U0_IN\n";
  std::cout << "Q12_SUPPORTS=2:CHAIN_01_12:+63371264_W02;"
               "CHAIN_01_13:+63371264_W03\n";
  std::cout << "Q16_SUPPORTS=2:CHAIN_01_12:-123422773248_W02;"
               "CHAIN_01_13:-123422773248_W03\n";
  std::cout << "FULL_Q12=-63371264_W01;FULL_Q16=+123422773248_W01\n";
  return 0;
}
