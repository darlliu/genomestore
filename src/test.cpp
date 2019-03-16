#define CATCH_CONFIG_MAIN
#include "base.hpp"
#include <catch2/catch.hpp>

TEST_CASE("TESTING BASE DB IO FUNCTIONS") {
  auto bdb = basedb("testdb", "test.db");
  SECTION("TESTING DB INITIALIZATION") {
    REQUIRE(bdb.getldb().cursor != nullptr);
  }
  SECTION("TESTING DB PUT STRING KEY AND VAL") {
    bdb.setdb("testkey", "testval");
    auto msg = bdb.getdb("testkey");
    REQUIRE(msg == "testval");
  }
}
