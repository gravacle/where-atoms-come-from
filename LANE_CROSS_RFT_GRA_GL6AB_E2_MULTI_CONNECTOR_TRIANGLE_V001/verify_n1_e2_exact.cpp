#include <gmpxx.h>

#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

using Integer = mpz_class;
using Vector = std::vector<Integer>;

namespace {

long choose(int n, int k) {
  long value = 1;
  for (int j = 1; j <= k; ++j) value = value * (n - k + j) / j;
  return value;
}

const std::array<Integer, 15> kExpected19 = {
    Integer("-67530641899520"), Integer("-135061283799040"),
    Integer("-165289008267264"), Integer("-158213815304192"),
    Integer("-132193339604992"), Integer("-105585215864832"),
    Integer("-87422498291712"), Integer("-77413660606464"),
    Integer("-70744951750656"), Integer("-62882751594496"),
    Integer("-51947031625728"), Integer("-38655920619520"),
    Integer("-25051091828736"), Integer("-13213688233984"),
    Integer("-4594932678656")};

int replay_block(int source_basis, int read_basis) {
  constexpr int kLinks = 16;
  constexpr int kDimension = 1 << kLinks;
  constexpr int kMaxTotal = 20;
  constexpr int kPairs[6][2] = {
      {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  constexpr int kE[2][6] = {
      {1, -1, 0, 0, -1, 1}, {1, 0, -1, -1, 0, 1}};

  std::vector<std::pair<int, int>> interactions;
  for (int cell = 0; cell < 4; ++cell) {
    for (int a = 0; a < 4; ++a) {
      for (int b = a + 1; b < 4; ++b) {
        interactions.push_back({4 * cell + a, 4 * cell + b});
      }
    }
  }
  for (int cell = 0; cell < 4; ++cell) {
    for (int other = cell + 1; other < 4; ++other) {
      interactions.push_back({4 * cell + other, 4 * other + cell});
    }
  }
  if (interactions.size() != 30) throw std::runtime_error("interaction census");

  std::vector<int> diagonal(kDimension), source(kDimension), read(kDimension);
  for (int word = 0; word < kDimension; ++word) {
    for (const auto& edge : interactions) {
      diagonal[word] +=
          2 * ((word >> edge.first) & 1) * ((word >> edge.second) & 1);
    }
    for (int pair = 0; pair < 6; ++pair) {
      const int a = kPairs[pair][0];
      const int b = kPairs[pair][1];
      const int za = 1 - 2 * ((word >> a) & 1);
      const int zb = 1 - 2 * ((word >> b) & 1);
      source[word] += kE[source_basis][pair] * za * zb;
      const int wa = 1 - 2 * ((word >> (4 + a)) & 1);
      const int wb = 1 - 2 * ((word >> (4 + b)) & 1);
      read[word] += kE[read_basis][pair] * wa * wb;
    }
  }

  auto apply_h = [&](const Vector& in) {
    Vector out(kDimension);
    for (int word = 0; word < kDimension; ++word) {
      out[word] += diagonal[word] * in[word];
      for (int q = 0; q < kLinks; ++q) out[word] -= in[word ^ (1 << q)];
    }
    return out;
  };

  std::vector<Vector> blank_powers;
  blank_powers.emplace_back(kDimension);
  blank_powers[0][0] = 1;
  for (int q = 0; q <= kMaxTotal; ++q) {
    blank_powers.push_back(apply_h(blank_powers.back()));
  }

  std::map<std::tuple<int, int, int>, Integer> table;
  for (int c = 0; c <= kMaxTotal; ++c) {
    Vector work(kDimension);
    for (int word = 0; word < kDimension; ++word) {
      work[word] = source[word] * blank_powers[c][word];
    }
    for (int b = 0; b + c <= kMaxTotal; ++b) {
      for (int a = 0; a + b + c <= kMaxTotal; ++a) {
        Integer value = 0;
        for (int word = 0; word < kDimension; ++word) {
          value += blank_powers[a][word] * read[word] * work[word];
        }
        table[{a, b, c}] = value;
      }
      if (b + c < kMaxTotal) work = apply_h(work);
    }
  }

  auto mixed = [&](int p, int r) {
    Integer value = 0;
    for (int l = 0; l <= p; ++l) {
      for (int j = 0; j <= r; ++j) {
        long coefficient = choose(p, l) * choose(r, j);
        if ((l + j) & 1) coefficient = -coefficient;
        value += Integer(coefficient) *
                 (table[{p - l + r - j, j, l}] -
                  table[{l + j, r - j, p - l}]);
      }
    }
    return value;
  };

  int checks = 0;
  for (int total = 0; total <= 18; ++total) {
    for (int r = 1; r <= total; ++r) {
      if (mixed(total - r, r) != 0) throw std::runtime_error("early nonzero");
      ++checks;
    }
  }
  for (int index = 0; index < 15; ++index) {
    const int p = 16 - index;
    const int r = 3 + index;
    if (mixed(p, r) != kExpected19[index]) {
      throw std::runtime_error("order-19 coefficient mismatch");
    }
    ++checks;
  }
  if (mixed(18, 1) != 0 || mixed(17, 2) != 0 ||
      mixed(1, 18) != 0 || mixed(0, 19) != 0) {
    throw std::runtime_error("order-19 endpoint mismatch");
  }
  checks += 4;
  for (int r = 1; r <= 20; ++r) {
    if (mixed(20 - r, r) != 0) throw std::runtime_error("even total mismatch");
    ++checks;
  }
  return checks;
}

}  // namespace

int main() {
  int checks = 0;
  for (int source = 0; source < 2; ++source) {
    for (int read = 0; read < 2; ++read) {
      checks += replay_block(source, read);
    }
  }
  std::cout << "PASS GL6AB exact full-N1 checks " << checks << "/" << checks
            << "\n";
  return 0;
}
