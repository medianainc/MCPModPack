context("target_dose_estimation")

test_that("Target dose estimation for the linear dose response function", {

	delta = seq(0, 1, 0.1)
	coef = c(0.2, 1.2)
	dose = rep(0, length(delta))
	value = rep(0, length(delta))

	for (i in 1:length(delta)) {
		dose[i] = TestFindTargetDose(delta[i], 1, coef)
     	if (dose[i] >= 0) {
     		value[i] = coef[2] * dose[i]
     	} else {
     		value[i] = -1
     		delta[i] = -1
     	}
	}

    expect_equal(delta, value)

})

test_that("Target dose estimation for the quadratic dose response function", {

	delta = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, -1)
	dose = rep(0, length(delta))
	value = rep(0, length(delta))

	for (i in 1:length(delta)) {
		dose[i] = TestFindTargetDose(delta[i], 2, coef)
     	if (dose[i] >= 0) {
     		value[i] = coef[2] * dose[i] + coef[3] * dose[i]^2
     	} else {
     		value[i] = -1
     		delta[i] = -1
     	}
	}

    expect_equal(delta, value)

})

test_that("Target dose estimation for the exponential dose response function", {

	delta = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, 1)
	dose = rep(0, length(delta))
	value = rep(0, length(delta))

	for (i in 1:length(delta)) {
		dose[i] = TestFindTargetDose(delta[i], 3, coef)
     	if (dose[i] >= 0) {
     		value[i] = coef[2] * (exp(dose[i] / coef[3]) - 1)
     	} else {
     		value[i] = -1
     		delta[i] = -1
     	}
	}

    expect_equal(delta, value)

})


test_that("Target dose estimation for the Emax dose response function", {

	delta = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, 1)
	dose = rep(0, length(delta))
	value = rep(0, length(delta))

	for (i in 1:length(delta)) {
		dose[i] = TestFindTargetDose(delta[i], 4, coef)
     	if (dose[i] >= 0) {
     		value[i] = coef[2] * dose[i] / (coef[3] + dose[i]) 
     	} else {
     		value[i] = -1
     		delta[i] = -1
     	}
	}

    expect_equal(delta, value)

})

test_that("Target dose estimation for the logistic dose response function", {

	delta = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, 0.1, 1)
	dose = rep(0, length(delta))
	value = rep(0, length(delta))

	for (i in 1:length(delta)) {
		dose[i] = TestFindTargetDose(delta[i], 5, coef)
     	if (dose[i] >= 0) {
     		value[i] = coef[2] / (1 + exp((coef[3] - dose[i]) / coef[4])) - coef[2] / (1 + exp((coef[3]) / coef[4]))
     	} else {
     		value[i] = -1
     		delta[i] = -1
     	}
	}

    expect_equal(delta, value)

})

test_that("Target dose estimation for the sigEmax dose response function", {

	delta = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, 0.1, 1)
	dose = rep(0, length(delta))
	value = rep(0, length(delta))

	for (i in 1:length(delta)) {
		dose[i] = TestFindTargetDose(delta[i], 6, coef)
     	if (dose[i] >= 0) {
     		value[i] = coef[2] * dose[i]^coef[4] / (dose[i]^coef[4] + coef[3]^coef[4])
     	} else {
     		value[i] = -1
     		delta[i] = -1
     	}
	}

    expect_equal(delta, value)

})
