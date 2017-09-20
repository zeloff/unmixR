library(unmixR)
context("Stirling numbers")

test_that("Test that s(0, 0) == 1", {
					expect_equal(stirling1(0, 0), matrix(1))
})

test_that("Test that negative numbers aren't allowed", {
					expect_error(stirling1(-1, 0))
					expect_error(stirling1(-1, 1))
					expect_error(stirling1(0, -1))
					expect_error(stirling1(1, -1))
	})

test_that("Test that s(9,9) is correct", {
	expect_equal(
		stirling1(9, 9),
		matrix(c( 1,  0,  0,  0,  0,   0,    0,     0,      0,       0,
						 NA,  1, -1 , 2, -6,  24, -120,   720,  -5040,   40320,
						 NA, NA,  1, -3, 11, -50 , 274, -1764,  13068, -109584,
						 NA, NA, NA,  1, -6 , 35, -225,  1624, -13132,  118124,
						 NA, NA, NA, NA,  1, -10,   85,  -735,   6769,  -67284,
						 NA, NA, NA, NA, NA,   1,  -15,   175,  -1960,   22449,
						 NA, NA, NA, NA, NA,  NA,    1,   -21,    322,   -4536,
						 NA, NA, NA, NA, NA,  NA,   NA,     1,    -28,     546,
						 NA, NA, NA, NA, NA,  NA,   NA,    NA,      1,     -36,
						 NA, NA, NA, NA, NA,  NA,   NA,    NA,     NA,       1),
					 nrow = 10, dimnames = list(n = 0:9, k = 0:9)))
	})

