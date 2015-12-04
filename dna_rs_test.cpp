#include "catch.hpp"
#include "debug.hpp"
#include "dna_rs.hpp"
#include "dna_bv_rs.hpp"
#include "utility.hpp"
#include <chrono>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/next_prior.hpp>

using namespace cosmo::index;

/*
TEST_CASE("Vectorized popcount", "[bithack][vector]") {
  auto popcount_v = vectorized_popcount<popcounter<3>>();
  vector<uint64_t> v{0,0,0,0,0,0,0};
  REQUIRE(popcount_v(v, 0) == v.size() * 21);
  REQUIRE(popcount_v(v, 1) == 0);
  set_symbol<3>(v[1], 0, 1);
  set_symbol<3>(v[3], 20, 1);
  REQUIRE(popcount_v(v, 1) == 2);
}

TEST_CASE("Block manipulation", "[bithack]") {
  typedef uint64_t block_t;
  block_t x = 0;
  
  SECTION("Set symbols", "[set]") {
    // TODO: add test for throwing
    x = set_symbol<3>(x, 0, 2);
    x = set_symbol<3>(x, 1, 7);
    x = set_symbol<3>(x, 4, 3);
    x = set_symbol<3>(x, 2, 3);
    x = set_symbol<3>(x, 7, 5);
    x = set_symbol<3>(x, 9, 5);
 
    SECTION("Get symbols", "[get]") {
      REQUIRE(get_symbol<3>(x, 0) == 2);
      REQUIRE(get_symbol<3>(x, 1) == 7);
      REQUIRE(get_symbol<3>(x, 2) == 3);
      REQUIRE(get_symbol<3>(x, 3) == 0);
      REQUIRE(get_symbol<3>(x, 4) == 3);
      REQUIRE(get_symbol<3>(x, 5) == 0);
      REQUIRE(get_symbol<3>(x, 7) == 5);
      REQUIRE(get_symbol<3>(x, 9) == 5);
    }

    SECTION("Update symbols", "[update]") {
      x = set_symbol<3>(x, 1, 4);
      SECTION("Get updated symbols", "[get]") {
        REQUIRE(get_symbol<3>(x, 0) == 2);
        REQUIRE(get_symbol<3>(x, 1) == 4);
        REQUIRE(get_symbol<3>(x, 2) == 3);
        REQUIRE(get_symbol<3>(x, 3) == 0);
        REQUIRE(get_symbol<3>(x, 4) == 3);
        REQUIRE(get_symbol<3>(x, 5) == 0);
      }
    }

    SECTION("Popcount symbols", "[popcount]") {
      auto popcount = popcounter<3>();
      REQUIRE(popcount(x,0) == 15);
      REQUIRE(popcount(x,1) == 0);
      REQUIRE(popcount(x,2) == 1);
      REQUIRE(popcount(x,3) == 2);
      REQUIRE(popcount(x,4) == 0);
      REQUIRE(popcount(x,5) == 2);
      REQUIRE(popcount(x,6) == 0);
      REQUIRE(popcount(x,7) == 1);
    }

    SECTION("Rank symbols", "[rank]") {
      cout << bitset<64>(x) << endl;

      SECTION("Rank at 0 is always 0") {
        for (size_t c = 0; c < 9; ++c) {
          REQUIRE(rank<3>(x,0,c) == 0);
        }
      }
      REQUIRE(rank<3>(x,1,2) == 1);
    }
  }
}
*/

namespace cosmo {
namespace internal {
template <typename begin, typename end, typename F>
struct static_for_each {
  static void call() {
    typedef typename boost::mpl::deref<begin>::type current_type;
    typedef typename boost::mpl::next<begin>::type next_type;

    F::template call<current_type>();
    static_for_each<next_type, end, F>::call();
  }
};

template <typename end, typename F>
struct static_for_each<end, end, F> { static void call() {} };
} // namespace internal

template <typename sequence, typename F>
void static_for_each() {
  typedef typename boost::mpl::begin<sequence>::type begin;
  typedef typename boost::mpl::end<sequence>::type   end;

  internal::static_for_each<begin, end, F>::call();
}

template<typename ... types> using type_list = boost::mpl::list<types...>;
} // namespace cosmo


#define TYPED_TEST_CASE(name, tags, types) \
  struct COSMO_##__LINE__ { template <typename T> static void call(); }; \
  TEST_CASE(name, tags) { cosmo::static_for_each<types, COSMO_##__LINE__>(); } \
  template <typename T> void COSMO_##__LINE__::call()


TEST_CASE("Large Query", "[benchmark]") {
  using namespace sdsl;
  size_t n = 10e5;
  size_t m = 10e5;
  auto input   = cosmo::random_string("$acgtACGT", n);
  int_vector<8> temp(input.size());
  for (size_t i = 0; i < input.size(); ++i) temp[i] = input[i];

  auto rnk_idxs    = cosmo::random_uints(0, n,   m);
  auto acc_idxs    = cosmo::random_uints(0, n-1, m);
  auto sel_idxs    = cosmo::random_uints(1, n,   m);
  auto queries = cosmo::random_string("$acgtACGT", m);

  // Construct
  typedef wt_huff<rrr_vector<63>> wt_t;
  typedef dna_bv_rs<sdsl::hyb_vector<>> rs_t;
  //typedef dna_bv_rs<> rs_t;

  auto start = std::chrono::steady_clock::now();
  rs_t rs(input);
  auto end = std::chrono::steady_clock::now();
  std::cout << "RS Construct average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;

  start = std::chrono::steady_clock::now();
  wt_t wt;
  construct_im(wt, temp);
  end = std::chrono::steady_clock::now();
  std::cout << "WT Construct average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;

  // Access
  start = std::chrono::steady_clock::now();
  for (auto i : acc_idxs) {
    auto x = rs[i];
  }
  end = std::chrono::steady_clock::now();
  std::cout << "RS Access average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;

  start = std::chrono::steady_clock::now();
  for (auto i : acc_idxs) {
    auto x = wt[i];
  }
  end = std::chrono::steady_clock::now();
  std::cout << "WT Access average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;


  // Rank
  start = std::chrono::steady_clock::now();
  for (auto i:rnk_idxs) {
    auto x = rs.rank(i, queries[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "RS Rank average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;

  start = std::chrono::steady_clock::now();
  for (auto i:rnk_idxs) {
    auto x = wt.rank(i, queries[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "WT Rank average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;

  /*
  // Select
  start = std::chrono::steady_clock::now();
  for (auto i : sel_idxs) {
    size_t max = rs.rank(rs.size(), queries[i]);
    auto x = rs.select(std::min(i,max), queries[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "RS Rank average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;

  start = std::chrono::steady_clock::now();
  for (auto i : sel_idxs) {
    size_t max = wt.rank(wt.size(), queries[i]);
    auto x = wt.select(std::min(i, max), queries[i]);
  }
  end = std::chrono::steady_clock::now();
  std::cout << "WT Rank average time per element: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/double(m) << " ns" << std::endl;
  */
}

typedef cosmo::type_list<dna_bv_rs<>> rank_types;
TYPED_TEST_CASE("DNA index queries are answered", "[index][dna_rs]", rank_types) {
  using namespace sdsl;
  const std::string input  = "acgt$acgt$ACGTacgt";
  int_vector<8> temp(input.size());
  for (size_t i = 0; i < input.size(); ++i) temp[i] = input[i];
  wt_blcd<> wt;
  construct_im(wt, temp);

  T x(input);

  REQUIRE(x.size() == input.length());

  // Access
  SECTION("original elements are accessible", "[access]") {
    for (size_t i = 0; i < x.size(); ++i) {
      REQUIRE(x[i] == input[i]);
    }
  }

  // Rank
  SECTION("ranks are computed for each symbol over [0, i)", "[rank]") {
    size_t c_i = 0;
    for (char c : std::string("$acgtACGT")) {
      for (size_t i = 0; i <= x.size(); ++i) {
        REQUIRE(x.rank(i,c) == wt.rank(i,c));
      }
      ++c_i;
    }
  }

  // Select
  SECTION("Select is computed for each symbol over [1, m]", "[select]") {
    size_t c_i = 0;
    for (char c : std::string("$acgtACGT")) {
      for (size_t i = 1; i <= x.rank(x.size(), c); ++i) {
        REQUIRE(x.select(i,c) == wt.select(i,c));
      }
      ++c_i;
    }
  }
}
