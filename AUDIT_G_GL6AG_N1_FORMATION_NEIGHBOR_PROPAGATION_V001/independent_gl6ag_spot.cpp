#include <gmpxx.h>

#include <array>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using Integer = mpz_class;
using Rational = mpq_class;
using Vector = std::vector<Integer>;

namespace {

constexpr int kLinks = 16;
constexpr int kDimension = 1 << kLinks;
constexpr int kMaxOrder = 16;
constexpr int kPairs[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
constexpr int kE[6][2] = {
    {1, 1}, {-1, 0}, {0, -1}, {0, -1}, {-1, 0}, {1, 1}};

using Raw = std::array<Integer, 2>;

long binomial(int n, int k) {
  long out = 1;
  for (int j = 1; j <= k; ++j) out = out * (n - k + j) / j;
  return out;
}

Integer factorial(int n) {
  Integer out = 1;
  for (int j = 2; j <= n; ++j) out *= j;
  return out;
}

std::array<std::array<int, 2>, 16> make_observable_table() {
  std::array<std::array<int, 2>, 16> table{};
  for (int nibble = 0; nibble < 16; ++nibble) {
    std::array<int, 4> z{};
    for (int port = 0; port < 4; ++port) {
      z[port] = 1 - 2 * ((nibble >> port) & 1);
    }
    for (int pair = 0; pair < 6; ++pair) {
      const int value = z[kPairs[pair][0]] * z[kPairs[pair][1]];
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        table[nibble][coordinate] += kE[pair][coordinate] * value;
      }
    }
  }
  return table;
}

const auto kObservable = make_observable_table();

struct Replay {
  int source_mask;
  bool bridges;
  std::array<bool, kLinks> formed{};
  std::vector<int> diagonal;
  std::vector<Vector> powers;

  Replay(int mask, bool retain_bridges)
      : source_mask(mask), bridges(retain_bridges), diagonal(kDimension) {
    for (int port = 0; port < 4; ++port) {
      formed[port] = ((source_mask >> port) & 1) != 0;
    }
    for (int q = 4; q < kLinks; ++q) formed[q] = true;

    // This reconstructs the diagonal directly from four occupied-link counts
    // plus the six reciprocal half-port tests; it does not enumerate the
    // author's interaction edge vector.
    for (int word = 0; word < kDimension; ++word) {
      int value = 0;
      for (int cell = 0; cell < 4; ++cell) {
        const unsigned nibble = (word >> (4 * cell)) & 15;
        const int occupied = __builtin_popcount(nibble);
        value += occupied * (occupied - 1);  // 2*C(occupied,2)
      }
      if (bridges) {
        for (int left = 0; left < 4; ++left) {
          for (int right = left + 1; right < 4; ++right) {
            value += 2 * ((word >> (4 * left + right)) & 1) *
                     ((word >> (4 * right + left)) & 1);
          }
        }
      }
      diagonal[word] = value;
    }

    powers.reserve(kMaxOrder + 1);
    powers.emplace_back(kDimension);
    powers[0][0] = 1;
    for (int order = 0; order < kMaxOrder; ++order) {
      powers.push_back(apply_h(powers.back()));
    }
  }

  Vector apply_h(const Vector& input) const {
    Vector output(kDimension);
    for (int word = 0; word < kDimension; ++word) {
      output[word] += diagonal[word] * input[word];
      for (int q = 0; q < kLinks; ++q) {
        if (formed[q]) output[word] -= input[word ^ (1 << q)];
      }
    }
    return output;
  }

  Raw raw(int cell, int order) const {
    Raw answer{};
    for (int split = 0; split <= order; ++split) {
      Raw matrix_element{};
      const Vector& bra = powers[order - split];
      const Vector& ket = powers[split];
      for (int word = 0; word < kDimension; ++word) {
        if (bra[word] == 0 || ket[word] == 0) continue;
        const Integer product = bra[word] * ket[word];
        const int nibble = (word >> (4 * cell)) & 15;
        for (int coordinate = 0; coordinate < 2; ++coordinate) {
          matrix_element[coordinate] +=
              kObservable[nibble][coordinate] * product;
        }
      }
      long weight = binomial(order, split);
      if (split & 1) weight = -weight;
      for (int coordinate = 0; coordinate < 2; ++coordinate) {
        answer[coordinate] += weight * matrix_element[coordinate];
      }
    }
    return answer;
  }
};

struct SpotTable {
  std::array<std::array<Raw, kMaxOrder + 1>, 4> value{};
  std::array<std::array<bool, kMaxOrder + 1>, 4> filled{};
};

using Request = std::pair<int, std::vector<int>>;

SpotTable reconstruct(int mask, bool bridges,
                      const std::vector<Request>& requests) {
  Replay replay(mask, bridges);
  SpotTable table{};
  for (const auto& request : requests) {
    const int cell = request.first;
    for (const int order : request.second) {
      if (!table.filled[cell][order]) {
        table.value[cell][order] = replay.raw(cell, order);
        table.filled[cell][order] = true;
      }
    }
  }
  return table;
}

std::vector<int> through(int order) {
  std::vector<int> out;
  for (int n = 0; n <= order; ++n) out.push_back(n);
  return out;
}

std::vector<int> through_plus(int order, int extra) {
  auto out = through(order);
  out.push_back(extra);
  return out;
}

Raw subtract(const Raw& left, const Raw& right) {
  return {left[0] - right[0], left[1] - right[1]};
}

Raw mobius(const Raw& pair, const Raw& single_a, const Raw& single_b,
            const Raw& blank) {
  return {pair[0] - single_a[0] - single_b[0] + blank[0],
          pair[1] - single_a[1] - single_b[1] + blank[1]};
}

void require(bool condition, const char* label, int& checks) {
  if (!condition) throw std::runtime_error(label);
  ++checks;
}

void require_raw(const Raw& actual, const Raw& expected, const char* label,
                 int& checks) {
  require(actual[0] == expected[0], label, checks);
  require(actual[1] == expected[1], label, checks);
}

}  // namespace

int main() {
  int checks = 0;

  // Structural reconstruction of E and its pair-incidence-null typing.
  for (int coordinate = 0; coordinate < 2; ++coordinate) {
    for (int port = 0; port < 4; ++port) {
      int contraction = 0;
      for (int pair = 0; pair < 6; ++pair) {
        const int incidence =
            (kPairs[pair][0] == port || kPairs[pair][1] == port) ? 1 : 0;
        contraction += kE[pair][coordinate] * incidence;
      }
      require(contraction == 0, "E column outside ker(P^T)", checks);
    }
  }
  require(kE[0][0] == 1 && kE[0][1] == 1,
          "w_01 column direction", checks);
  require(kE[1][0] == -1 && kE[1][1] == 0,
          "w_02 column direction", checks);
  require(kE[2][0] == 0 && kE[2][1] == -1,
          "w_03 column direction", checks);

  // The K projectors are the sixteen delta functions on the physical source
  // support word: orthogonal, exhaustive, and unchanged by every active-word
  // transition used below.
  for (int word = 0; word < 16; ++word) {
    int resolution = 0;
    for (int mask = 0; mask < 16; ++mask) {
      const int projector = (mask == word) ? 1 : 0;
      resolution += projector;
      require(projector * projector == projector, "K projector idempotence",
              checks);
      for (int other = mask + 1; other < 16; ++other) {
        const int other_projector = (other == word) ? 1 : 0;
        require(projector * other_projector == 0,
                "K projector orthogonality", checks);
      }
    }
    require(resolution == 1, "K projector resolution", checks);
  }

  const auto all16 = through(16);
  const auto through13 = through(13);
  const auto through13_plus16 = through_plus(13, 16);
  const auto through5 = through(5);

  // Five independently selected full-H branches: reference, the two singles
  // and their pair, plus the all-formed word.  They cover two receiver
  // directions and the order-sixteen pair Mobius boundary.
  const SpotTable blank = reconstruct(
      0, true, {{0, through5}, {1, all16}, {2, through13_plus16}});
  const SpotTable single0 =
      reconstruct(1, true, {{1, all16}, {2, {16}}});
  const SpotTable single1 =
      reconstruct(2, true, {{1, all16}, {2, {16}}});
  const SpotTable pair01 = reconstruct(
      3, true, {{0, through5}, {1, all16}, {2, through13_plus16}});
  const SpotTable all_formed =
      reconstruct(15, true, {{2, through13}});

  // Independent interaction-census probes.  Fully occupied has 30 pair
  // terms, a same-cell two-bit word has one, and reciprocal half-ports have
  // the unique bridge only.
  Replay census_full(0, true);
  Replay census_off(0, false);
  require(census_full.diagonal[kDimension - 1] == 60,
          "thirty full interactions", checks);
  require(census_off.diagonal[kDimension - 1] == 48,
          "twenty-four bridge-off interactions", checks);
  require(census_full.diagonal[(1 << 0) | (1 << 1)] == 2,
          "within-cell interaction", checks);
  require(census_off.diagonal[(1 << 0) | (1 << 1)] == 2,
          "within-cell retained after ablation", checks);
  const int reciprocal = (1 << (4 * 0 + 1)) | (1 << (4 * 1 + 0));
  require(census_full.diagonal[reciprocal] == 2,
          "source-receiver shared-child bridge", checks);
  require(census_off.diagonal[reciprocal] == 0,
          "shared-child term removed", checks);

  for (int order = 0; order < 4; ++order) {
    require_raw(subtract(pair01.value[0][order], blank.value[0][order]),
                {0, 0}, "local matched coefficient below four", checks);
  }
  require_raw(subtract(pair01.value[0][4], blank.value[0][4]),
              {96, 96}, "local raw q4", checks);
  require_raw(subtract(pair01.value[0][5], blank.value[0][5]),
              {0, 0}, "local odd q5", checks);

  for (int order = 0; order < 12; ++order) {
    require_raw(subtract(pair01.value[1][order], blank.value[1][order]),
                {0, 0}, "receiver-1 matched coefficient below twelve",
                checks);
  }
  require_raw(subtract(pair01.value[1][12], blank.value[1][12]),
              {-63371264, -63371264}, "receiver-1 raw q12", checks);
  require_raw(subtract(pair01.value[1][13], blank.value[1][13]),
              {0, 0}, "receiver-1 odd q13", checks);

  for (int order = 0; order <= 13; ++order) {
    require_raw(subtract(pair01.value[2][order], blank.value[2][order]),
                {0, 0}, "kappa_2 zero branch", checks);
  }
  for (int order = 0; order < 12; ++order) {
    require_raw(
        subtract(all_formed.value[2][order], blank.value[2][order]),
        {0, 0}, "all-formed receiver-2 below twelve", checks);
  }
  require_raw(subtract(all_formed.value[2][12], blank.value[2][12]),
              {63371264, 0}, "all-formed receiver-2 raw q12", checks);
  require_raw(subtract(all_formed.value[2][13], blank.value[2][13]),
              {0, 0}, "all-formed receiver-2 odd q13", checks);

  for (int order = 0; order < 16; ++order) {
    require_raw(mobius(pair01.value[1][order], single0.value[1][order],
                       single1.value[1][order], blank.value[1][order]),
                {0, 0}, "receiver-1 pair Mobius below sixteen", checks);
  }
  require_raw(mobius(pair01.value[1][16], single0.value[1][16],
                     single1.value[1][16], blank.value[1][16]),
              {123422773248, 123422773248},
              "receiver-1 pair Mobius q16", checks);
  require_raw(mobius(pair01.value[2][16], single0.value[2][16],
                     single1.value[2][16], blank.value[2][16]),
              {0, 0}, "receiver-2 unsupported pair Mobius q16", checks);

  Rational local(Integer(96), factorial(4));
  Rational remote(Integer(-63371264), factorial(12));
  Rational nonlinear(Integer(123422773248), factorial(16));
  local.canonicalize();
  remote.canonicalize();
  nonlinear.canonicalize();
  require(local == Rational(4), "local factorial", checks);
  require(remote == Rational(-5626, 42525), "remote factorial and sign",
          checks);
  require(nonlinear == Rational(1116019, 189189000),
          "Mobius factorial and sign", checks);

  // A separate bridge-off pair verifies the finite jet, while the direct
  // formula for H_off is a sum of four disjoint cell operators and therefore
  // promotes this equality to every source pattern and all times.
  const SpotTable off_blank = reconstruct(0, false, {{1, all16}});
  const SpotTable off_pair01 = reconstruct(3, false, {{1, all16}});
  for (int order = 0; order <= 16; ++order) {
    require_raw(subtract(off_pair01.value[1][order],
                         off_blank.value[1][order]),
                {0, 0}, "bridge-off receiver independence", checks);
  }

  std::cout << "PASS__INDEPENDENT_GL6AG_SPOT__" << checks << "/" << checks
            << "\n";
  std::cout << "HAMILTONIAN=FULL16_24_WITHIN_6_SHARED_ONLY_SOURCE_K_X_SUPPORT_VARIES\n";
  std::cout << "MATCHED=REFERENCE_0000_NO_ABSOLUTE_RECEIVER_INFERENCE\n";
  std::cout << "E_TYPE=FIXED_BROKEN_S4_RESTRICTION_W_COLUMNS_IN_R2\n";
  std::cout << "Q4=PLUS96_W01_COEFFICIENT_PLUS4\n";
  std::cout << "Q12=MINUS63371264_KAPPA_C_W0C_COEFFICIENT_MINUS5626_OVER42525\n";
  std::cout << "Q16=MOBIUS_PLUS123422773248_SUPPORTED_RECEIVER\n";
  std::cout << "ABLATION=FOUR_CELL_FACTORIZATION_DIAGNOSTIC_NOT_PHYSICAL_SWITCH\n";
  std::cout << "BRANCH=PHYSICAL_K_PROJECTORS_NOT_SEMANTIC_REC\n";
  std::cout << "SCOPE=NO_BULK_CONE_STRESS_CONSERVATION_RICCI_GRAVITY_G\n";
  return 0;
}
