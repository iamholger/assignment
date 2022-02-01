#include <gtest/gtest.h>
#include <stdexcept>

#include "matlib.h"


TEST(MatlibTest, BasicAssertions) {

  Matrix<3, 3, double> A;
  EXPECT_EQ(A.isCorrectlyInitialised(), false);
  
  A = 1,2,3,4,5,6,7,8,9;
  EXPECT_EQ(A.determinant(), 0);
  EXPECT_EQ(A.isCorrectlyInitialised(), true);
  
  EXPECT_EQ(A.isCorrectlyInitialised(), true);
}

TEST(MatlibTest, operations)
{

  Matrix<3, 4, double> A;
  try
  {

    transpose(A);
    FAIL() << "Matrix not properly initialised";
  }
  catch(std::logic_error const & err) {
      EXPECT_EQ(err.what(),std::string("Matrix not properly initialised"));
  }
  catch(...) {
      FAIL() << "Expected Matrix not properly initialised";
  }
}

TEST(MatlibTest, cofactor)
{

  Matrix<3,3,double> B;
  B = 2,-1,0,0,1,2,1,1,0;


  EXPECT_EQ(B.cofactor(0,0), -2);
  EXPECT_EQ(B.cofactor(0,1),  2);
  EXPECT_EQ(B.cofactor(0,2), -1);
  EXPECT_EQ(B.cofactor(1,0),  0);
  EXPECT_EQ(B.cofactor(1,1),  0);
  EXPECT_EQ(B.cofactor(1,2), -3);
  EXPECT_EQ(B.cofactor(2,0), -2);
  EXPECT_EQ(B.cofactor(2,1), -4);
  EXPECT_EQ(B.cofactor(2,2),  2);
  
  try
  {
    B.cofactor(9,9);
    FAIL() << "Trying to access element out of range";
  }
  catch(std::out_of_range const & err) {
      EXPECT_EQ(err.what(),std::string("Trying to access element out of range"));
  }
  catch(...) {
      FAIL() << "Trying to access element out of range";
  }
}
