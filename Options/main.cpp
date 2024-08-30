#include <iostream>

#include "options.h"

BlackScholesParams bs_params
{
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	OptionType::Call,
	0.0
};

MonteCarloParams mtc_params
{
	1'000'000,
	0.0,
	0.0,
	0.0,
	0.0,
	0.0,
	OptionType::Call,
	0.0
};

int main(int agrc, char** argv)
{
	FinancialCalculator financial_calculator;
	std::cout << financial_calculator.calculateBlackScholes(bs_params) << std::endl;
	std::cout << financial_calculator.calculateMonteCarlo(mtc_params) << std::endl;

	return EXIT_SUCCESS;
}
