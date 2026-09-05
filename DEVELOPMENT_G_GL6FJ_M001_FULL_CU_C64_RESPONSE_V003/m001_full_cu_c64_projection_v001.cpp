// GL6FJ V001 fixed-m001 polarization-complete scalar response wrapper.
//
// The numerical dynamics and owner-once scalar path ledger are inherited
// byte-for-byte from the independently audited GL6FA implementation.  This
// bounded successor only supplies a canonical 300-ray polarization design
// at the non-self-conjugate L4 character m=(0,0,1).  Direct long-run output
// is subordinate raw material and is never production evidence.
// The inherited checkpoint ledger retains stay_contact, writer_contact,
// stay_spectral, writer_spectral, stay_writer_interference, and disconnected
// exactly once; this wrapper never recomputes or substitutes an owner.

#include "../DEVELOPMENT_G_GL6FA_DIRECTIONALLY_CONTRACTED_CU_STAY_ENGINE_V002/v001_engine_with_renamed_entrypoint.inc"

#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gl6fj {

constexpr const char* TABLE_PATH =
    "DEVELOPMENT_G_GL6ER_COMPILED_FULL_CU_H6_PATH_SCORE_ENGINE_V001/"
    "LOCAL_OWNER_TABLES.bin";
constexpr const char* RAW_TOKEN =
    "GL6FJ_V001_SUBORDINATE_300_RAY_RAW_NOT_PRODUCTION";

constexpr int L = 4;
constexpr int CHANNELS = gl6er::CHANNELS;
constexpr int RAYS = gl6er::UPPER;
constexpr int POPULATION = 1024;
constexpr int BURN = 2560;
constexpr int WINDOW_LENGTH = 256;
constexpr int WINDOWS = 6;
constexpr uint64_t STAY_UNITS = 4096;
constexpr uint64_t WRITER_UNITS = 64;
constexpr uint64_t RANDOM_SEED = 96511001;

using Direction = gl6fa::Direction;

struct Ray {
  int index{};
  int left{};
  int right{};
  Direction direction{};
};

std::array<Ray, RAYS> make_rays() {
  std::array<Ray, RAYS> rays{};
  int cursor = 0;
  for (int left = 0; left < CHANNELS; ++left) {
    for (int right = left; right < CHANNELS; ++right) {
      Ray ray;
      ray.index = cursor;
      ray.left = left;
      ray.right = right;
      if (left == right) {
        ray.direction[left] = 1.0;
      } else {
        const double value = 1.0 / std::sqrt(2.0);
        ray.direction[left] = value;
        ray.direction[right] = value;
      }
      if (cursor != gl6er::upper_index(left, right))
        throw std::runtime_error("GL6FJ ray/upper-index mismatch");
      rays[cursor++] = ray;
    }
  }
  if (cursor != RAYS) throw std::runtime_error("GL6FJ ray census");
  return rays;
}

const std::array<Ray, RAYS>& rays() {
  static const std::array<Ray, RAYS> value = make_rays();
  return value;
}

constexpr gl6er::Character character() { return {{0, 0, 1}}; }

int parse_bounded_integer(const char* text, int lower, int upper,
                          const char* label) {
  const std::string value(text);
  std::size_t used = 0;
  const int answer = std::stoi(value, &used);
  if (used != value.size() || answer < lower || answer > upper)
    throw std::runtime_error(std::string("GL6FJ ") + label + " domain");
  return answer;
}

bool close(double left, double right, double tolerance = 3e-13) {
  return std::abs(left - right) <= tolerance;
}

void catalog_selftest() {
  const auto& catalog = rays();
  std::set<std::pair<int, int>> pairs;
  for (const Ray& ray : catalog) {
    double norm = 0.0;
    int support = 0;
    for (double value : ray.direction) {
      norm += value * value;
      if (value != 0.0) ++support;
    }
    if (!close(norm, 1.0) ||
        support != (ray.left == ray.right ? 1 : 2) ||
        !pairs.emplace(ray.left, ray.right).second)
      throw std::runtime_error("GL6FJ canonical ray failure");
  }
  if (pairs.size() != RAYS)
    throw std::runtime_error("GL6FJ unique ray failure");

  // Deterministic synthetic polarization replay.  It checks the exact
  // reconstruction convention used by the Python postprocessor.
  std::array<double, CHANNELS> first{};
  std::array<double, RAYS> second{};
  std::array<double, RAYS> q{};
  std::array<double, RAYS> linear{};
  for (int index = 0; index < CHANNELS; ++index)
    first[index] = (index % 7) - 3.0 + 0.125 * index;
  for (int left = 0; left < CHANNELS; ++left)
    for (int right = left; right < CHANNELS; ++right)
      second[gl6er::upper_index(left, right)] =
          (left + 1) * 0.25 - (right + 2) * 0.0625;
  for (const Ray& ray : catalog) {
    for (int index = 0; index < CHANNELS; ++index)
      linear[ray.index] += ray.direction[index] * first[index];
    for (int left = 0; left < CHANNELS; ++left) {
      for (int right = left; right < CHANNELS; ++right) {
        double value = second[gl6er::upper_index(left, right)]
                       * ray.direction[left] * ray.direction[right];
        if (left != right) value *= 2.0;
        q[ray.index] += value;
      }
    }
  }
  double residual = 0.0;
  for (const Ray& ray : catalog) {
    double recovered = 0.0;
    if (ray.left == ray.right) {
      recovered = q[ray.index];
      residual = std::max(
          residual, std::abs(linear[ray.index] - first[ray.left]));
    } else {
      const double diagonal_left =
          q[gl6er::upper_index(ray.left, ray.left)];
      const double diagonal_right =
          q[gl6er::upper_index(ray.right, ray.right)];
      recovered = q[ray.index]
                  - 0.5 * (diagonal_left + diagonal_right);
      const double expected_linear =
          (first[ray.left] + first[ray.right]) / std::sqrt(2.0);
      residual = std::max(
          residual, std::abs(linear[ray.index] - expected_linear));
    }
    residual = std::max(residual,
                        std::abs(recovered - second[ray.index]));
  }
  if (residual > 2e-13)
    throw std::runtime_error("GL6FJ polarization reconstruction failure");
  std::cout << std::setprecision(17)
            << "{\"schema\":\"GL6FJ_CANONICAL_300_RAY_SELFTEST_V001\""
            << ",\"channels\":24,\"symmetric_entries\":300"
            << ",\"unit_rays\":24,\"normalized_pair_rays\":276"
            << ",\"maximum_reconstruction_residual\":" << residual
            << ",\"result\":\"PASS\"}\n";
}

void character_selftest() {
  const gl6er::Model model = gl6er::build_model(L);
  const gl6er::Character m001 = character();
  const gl6er::Character minus_m001{{0, 0, -1}};
  bool self_conjugate = true;
  for (int axis = 0; axis < 3; ++axis) {
    self_conjugate = self_conjugate &&
        gl6er::mod(-m001.m[axis], L) == gl6er::mod(m001.m[axis], L);
  }
  if (self_conjugate)
    throw std::runtime_error("GL6FJ character conjugacy classification");
  const uint32_t p0 = model.p_node({0, 0, 0});
  const uint32_t c0 = model.c_node({0, 0, 0});
  const uint32_t pz = model.p_node({0, 0, 1});
  const gl6er::Profile pp0 = gl6er::profile(model, p0, 0, m001);
  const gl6er::Profile pc0 = gl6er::profile(model, c0, 0, m001);
  const gl6er::Profile ppz = gl6er::profile(model, pz, 0, m001);
  const gl6er::Profile mc0 = gl6er::profile(model, c0, 0, minus_m001);
  if (!close(pp0.value[0], 2.0) || !close(pp0.value[1], 0.0) ||
      !close(ppz.value[0], 0.0) || !close(ppz.value[1], -2.0) ||
      !close(pc0.value[0], 2.0 * std::cos(M_PI / 8.0)) ||
      !close(pc0.value[1], 2.0 * std::sin(M_PI / 8.0)) ||
      !close(mc0.value[0], 2.0 * std::cos(M_PI / 8.0)) ||
      !close(mc0.value[1], -2.0 * std::sin(M_PI / 8.0)))
    throw std::runtime_error("GL6FJ m001 P/C centering phases");
  for (uint32_t node = 0; node < model.node_count; ++node) {
    for (int pair = 0; pair < 6; ++pair) {
      const gl6er::Profile plus = gl6er::profile(model, node, pair, m001);
      const gl6er::Profile minus =
          gl6er::profile(model, node, pair, minus_m001);
      if (plus.channel != minus.channel ||
          !close(plus.value[0], minus.value[0], 2e-12) ||
          !close(plus.value[1], -minus.value[1], 2e-12))
        throw std::runtime_error("GL6FJ m001 signed conjugate transport");
    }
  }
  std::cout
      << "{\"schema\":\"GL6FJ_CHARACTER_SELFTEST_V001\""
      << ",\"L\":4,\"m001_self_conjugate\":false"
      << ",\"m001_minus_integer_lift\":[0,0,-1]"
      << ",\"m001_signed_conjugate_profiles_checked\":true"
      << ",\"m001_P_z1_phase\":\"exp(+i*pi/2)\""
      << ",\"m001_C_plus_phase\":\"exp(-i*pi/8)\""
      << ",\"m001_C_minus_phase\":\"exp(+i*pi/8)\""
      << ",\"translation_forced_anomalous_expectation_zero\":true"
      << ",\"finite_run_anomalous_zero_imposed\":false"
      << ",\"real_source_chart_dimension\":24"
      << ",\"real_source_chart_unaliased\":true"
      << ",\"integer_lift_centering_phases_retained\":true"
      << ",\"result\":\"PASS\"}\n";
}

void run_source(const Direction& direction, const Ray& ray, bool smoke) {
  const gl6er::Model model = gl6er::build_model(L);
  const gl6er::LocalTables tables(TABLE_PATH);
  const gl6er::Character m = character();
  gl6fa::SingleDirectionalCache cache(model, tables, m, direction);
  const gl6er::Walker seed =
      gl6er::initial_walker(model, gl6er::SeedKind::GL6CCDense);
  const std::vector<uint8_t> seed_digits = cache.digits(seed);
  const gl6fa::ScalarJet seed_stay = cache.cold(seed, seed_digits);

  const int population_size = smoke ? 32 : POPULATION;
  const int burn = smoke ? 4 : BURN;
  const int window_length = smoke ? 4 : WINDOW_LENGTH;
  const int windows = smoke ? 1 : WINDOWS;
  const uint64_t random_seed = smoke ? 96510001 : RANDOM_SEED;
  std::vector<gl6fa::StationaryParticle> population;
  population.reserve(population_size);
  for (int index = 0; index < population_size; ++index)
    population.push_back(
        {seed, seed_digits, seed_stay, {}, static_cast<uint32_t>(index)});
  std::vector<gl6fa::StationaryParticle> next = population;

  gl6fa::RunDiagnostics burn_diagnostics;
  int global_step = 0;
  for (; global_step < burn; ++global_step)
    gl6fa::stationary_step(
        model, cache, population, next, STAY_UNITS, WRITER_UNITS,
        random_seed, global_step, false, burn_diagnostics);

  const std::array<int, 3> checkpoints{{
      window_length / 4, window_length / 2, window_length}};
  std::cout << std::setprecision(17)
            << "{\"schema\":\""
            << "GL6FJ_SCALAR_PROJECTION_RAW_V001"
            << "\""
            << ",\"mode\":\"" << (smoke ? "smoke" : "subordinate_raw")
            << "\",\"L\":4,\"volume\":64"
            << ",\"component\":\"GL6CC_SELECTED_DENSE_H6_COMPONENT\""
            << ",\"action_order\":\"through_h6\""
            << ",\"momentum_label\":\"m001\""
            << ",\"character\":[" << m.m[0] << ',' << m.m[1] << ','
            << m.m[2] << ']'
            << ",\"channels\":24,\"symmetric_entries\":300";
  std::cout << ",\"design_kind\":\"canonical_full_response\""
            << ",\"ray_index\":" << ray.index
            << ",\"ray_left\":" << ray.left
            << ",\"ray_right\":" << ray.right
            << ",\"ray_kind\":\""
            << (ray.left == ray.right ? "unit" : "normalized_pair_sum")
            << "\"";
  std::cout << ",\"population\":" << population_size
            << ",\"burn\":" << burn
            << ",\"window_length\":" << window_length
            << ",\"windows\":" << windows
            << ",\"stay_units\":" << STAY_UNITS
            << ",\"writer_units\":" << WRITER_UNITS
            << ",\"stay_weight\":64"
            << ",\"random_seed\":" << random_seed
            << ",\"support_cache_dense_entries\":" << cache.dense_entries()
            << ",\"burn_stays\":" << burn_diagnostics.stays
            << ",\"burn_writers\":" << burn_diagnostics.writers
            << ",\"burn_digest\":" << burn_diagnostics.digest
            << ",\"checkpoints\":[";
  bool first = true;
  for (int window = 0; window < windows; ++window) {
    for (int index = 0; index < population_size; ++index) {
      population[index].ledger = {};
      population[index].ancestor = static_cast<uint32_t>(index);
      next[index] = population[index];
    }
    gl6fa::RunDiagnostics diagnostics;
    int checkpoint_index = 0;
    for (int step = 1; step <= window_length; ++step, ++global_step) {
      gl6fa::stationary_step(
          model, cache, population, next, STAY_UNITS, WRITER_UNITS,
          random_seed, global_step, true, diagnostics);
      if (step == checkpoints[checkpoint_index]) {
        if (!first) std::cout << ',';
        first = false;
        gl6fa::print_scalar_checkpoint(population, diagnostics, window, step,
                                       model.volume);
        ++checkpoint_index;
      }
    }
    if (checkpoint_index != 3)
      throw std::runtime_error("GL6FJ checkpoint census");
  }
  double sampled_cold_residual = 0.0;
  for (int index = 0; index < std::min(population_size, 8); ++index) {
    const gl6fa::ScalarJet direct =
        cache.cold(population[index].walker, population[index].digits);
    sampled_cold_residual = std::max(
        sampled_cold_residual,
        std::max(std::abs(direct.first - population[index].stay.first),
                 std::abs(direct.second - population[index].stay.second)));
  }
  if (sampled_cold_residual > 3e-7)
    throw std::runtime_error("GL6FJ final cached-stay residual");
  std::cout << "]"
            << ",\"support_cache_hits\":" << cache.support_hits()
            << ",\"support_cache_misses\":" << cache.support_misses()
            << ",\"writer_cache_hits\":" << cache.writer_hits()
            << ",\"writer_cache_misses\":" << cache.writer_misses()
            << ",\"sampled_cold_residual\":" << sampled_cold_residual
            << ",\"production_execution\":false"
            << ",\"physical_ward_rank_claimed\":false"
            << ",\"translation_forced_anomalous_expectation_zero\":true"
            << ",\"finite_run_anomalous_zero_imposed\":false"
            << ",\"real_source_chart_dimension\":24"
            << ",\"real_source_chart_unaliased\":true"
            << ",\"ceiling\":\"SELECTED_COMPONENT_L4_H6_M001__FULL_REAL_"
               "RESPONSE_RECONSTRUCTIBLE_ONLY_AFTER_ALL_300_SAME_STATE_"
               "RAYS__NO_RANGE_QUOTIENT_OR_WARD_CLAIM\"}\n";
}

void run_projection(int ray_index, bool smoke) {
  const Ray& ray = rays().at(ray_index);
  run_source(ray.direction, ray, smoke);
}

}  // namespace gl6fj

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "catalog-selftest") {
      gl6fj::catalog_selftest();
      return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "character-selftest") {
      gl6fj::character_selftest();
      return 0;
    }
    if (argc == 3 && std::string(argv[1]) == "smoke-projection") {
      const int ray = gl6fj::parse_bounded_integer(
          argv[2], 0, gl6fj::RAYS - 1, "ray");
      gl6fj::run_projection(ray, true);
      return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "pilot-projection-raw") {
      if (std::string(argv[3]) != gl6fj::RAW_TOKEN)
        throw std::runtime_error("GL6FJ subordinate raw token");
      const int ray = gl6fj::parse_bounded_integer(
          argv[2], 0, gl6fj::RAYS - 1, "ray");
      gl6fj::run_projection(ray, false);
      return 0;
    }
    std::cerr
        << "usage: m001_full_cu_c64_projection_v001 catalog-selftest\n"
        << "   or: ... character-selftest\n"
        << "   or: ... smoke-projection RAY_INDEX\n"
        << "   or: ... pilot-projection-raw RAY_INDEX "
           "SUBORDINATE_TOKEN\n";
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "ERROR=" << error.what() << '\n';
    return 1;
  }
}
