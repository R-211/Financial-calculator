#include <iostream>

#include "options.h"

int main(int agrc, char** argv)
{
	BlackScholesParams bs_params
	{
		0.2,
		100.0,
		105.0,
		0.5,
		0.3,
		OptionType::Call,
		5.0
	};

	MonteCarloParams mtc_params
	{
		1'000'000,
		0.2,
		100.0,
		105.0,
		0.5,
		0.3,
		OptionType::Call,
		5.0
	};

	FinancialCalculator financial_calculator;
	std::cout << financial_calculator.calculateBlackScholes(bs_params) << std::endl;
	std::cout << financial_calculator.calculateMonteCarlo(mtc_params) << std::endl;

	return EXIT_SUCCESS;
}