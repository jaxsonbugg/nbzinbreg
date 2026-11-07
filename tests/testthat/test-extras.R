test_that("provisional extras work", {
  skip("planned: extras phase, after the core matches the references")
  # Planned cases:
  # - confint(): the profile interval for theta contains theta-hat and stays
  #   positive;
  # - anova() and overdispersion_test() against pscl::odTest();
  # - vuong_test() against pscl::vuong();
  # - simulate() returns nsim response vectors of the right length;
  # - expected_counts() sums of observed and expected agree with n;
  # - tidy() and glance() return one row per parameter / one row per model.
})
