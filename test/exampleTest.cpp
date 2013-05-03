

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE example_test
#include <boost/test/unit_test.hpp>

//________________________________________________________________//


BOOST_AUTO_TEST_CASE( universeInOrder )
{
  BOOST_REQUIRE( 4 == 4 );
}


//________________________________________________________________//

// EOF
