context("dose_response_function")

test_that("Evaluate the linear dose response function", {

	x = seq(0, 1, 0.1)
	coef = c(0.2, 1.2)
	input = rep(0, length(x))
	output = rep(0, length(x))

	for (i in 1:length(x)) {
		input[i] = coef[1] + coef[2] * x[i]
		output[i] = TestDoseResponseFunction(x[i], 1, coef)
	}

    expect_equal(input, output)

})

test_that("Evaluate the quadratic dose response function", {

	x = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, -1)
	input = rep(0, length(x))
	output = rep(0, length(x))

	for (i in 1:length(x)) {
		input[i] = coef[1] + coef[2] * x[i] + coef[3] * x[i]^2
		output[i] = TestDoseResponseFunction(x[i], 2, coef)
	}

    expect_equal(input, output)

})

test_that("Evaluate the exponential dose response function", {

	x = seq(0, 1, 0.1)
	coef = c(0.2, 1.2, 1)
	input = rep(0, length(x))
	output = rep(0, length(x))

	for (i in 1:length(x)) {
		input[i] = coef[1] + coef[2] * (exp(x[i] / coef[3]) - 1) 
		output[i] = TestDoseResponseFunction(x[i], 3, coef)
	}

    expect_equal(input, output)

})

test_that("Evaluate the Emax dose response function", {

    x = seq(0, 1, 0.1)
    coef = c(0.2, 1.2, 1)
    input = rep(0, length(x))
    output = rep(0, length(x))

    for (i in 1:length(x)) {
        input[i] = coef[1] + coef[2] * x[i] / (coef[3] + x[i]) 
        output[i] = TestDoseResponseFunction(x[i], 4, coef)
    }

    expect_equal(input, output)

})

test_that("Evaluate the logistic dose response function", {

    x = seq(0, 1, 0.1)
    coef = c(0.2, 1.2, 0.1, 1)
    input = rep(0, length(x))
    output = rep(0, length(x))

    for (i in 1:length(x)) {
        input[i] = coef[1] + coef[2] / (1 + exp((coef[3] - x[i]) / coef[4]))
        output[i] = TestDoseResponseFunction(x[i], 5, coef)
    }

    expect_equal(input, output)

})


test_that("Evaluate the sigEmax dose response function", {

    x = seq(0, 1, 0.1)
    coef = c(0.2, 1.2, 0.1, 1)
    input = rep(0, length(x))
    output = rep(0, length(x))

    for (i in 1:length(x)) {
        input[i] = coef[1] + coef[2] * x[i]^coef[4] / (x[i]^coef[4] + coef[3]^coef[4])
        output[i] = TestDoseResponseFunction(x[i], 6, coef)
    }

    expect_equal(input, output)

})



