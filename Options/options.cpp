#include "options.h"

[[nodiscard]] double FinancialCalculator::calculateBlackScholes(const BlackScholesParams& params) const
{
	const auto d1 = (log(params.underlying_price / params.strike_price) + (params.interest_rate + (params.volatility * params.volatility) / 2.0) * params.time) / (params.volatility * std::sqrt(params.time));
	const auto d2 = d1 - params.volatility * std::sqrt(params.time);

	switch (params.option_type)
	{
	case OptionType::Call:
		return params.underlying_price * normalCdf(d1) - params.strike_price * std::exp(-params.interest_rate * params.time) * normalCdf(d2);
	case OptionType::Put:
		return params.strike_price * std::exp(-params.interest_rate * params.time) * normalCdf(-d2) - params.underlying_price * normalCdf(-d1);
	default:
		throw std::invalid_argument("[!] Invalid option type: must be Call or Put.");
	}
}

[[nodiscard]] inline double FinancialCalculator::calculateFutures(const FuturesParams& params) const
{
	return params.present_value * std::pow(1 + params.interest_rate, params.time);
}

[[nodiscard]] double FinancialCalculator::calculateGreeks(const GreeksParams& params, const Greeks greek) const
{
	if (params.time <= 0) throw std::runtime_error("[!] Time must be positive");
	if (params.volatility <= 0) throw std::runtime_error("[!] Volatility must be positive");

	const auto d1 = (log(params.underlying_price / params.strike_price) + (params.interest_rate + (params.volatility * params.volatility) / 2.0) * params.time) / (params.volatility * std::sqrt(params.time));
	const auto d2 = d1 - params.volatility * std::sqrt(params.time);

	const auto discount = std::exp(-params.interest_rate * params.time);			// Measures the present value of $1 to be received at a future date
	const auto dividend_discount = std::exp(-params.dividend_yield * params.time);  // It's similar to discount, but it's based on the dividend yield instead of the interest rate

	switch (greek)
	{
		case Greeks::Delta :
		{
			return params.option_type == OptionType::Call ?
				dividend_discount * normalCdf(d1) :
				dividend_discount * (normalCdf(d1) - 1);
		}
		case Greeks::Gamma :
		{
			return (dividend_discount * normalPdf(d1)) / (params.underlying_price * params.volatility * std::sqrt(params.time));
		}
		case Greeks::Theta : 
		{
			const double theta_part1 = -(params.underlying_price * params.volatility * dividend_discount * normalPdf(d1)) / (2 * std::sqrt(params.time));
			const double theta_part2 = params.option_type == OptionType::Call ?
				-params.interest_rate * params.strike_price * discount * normalCdf(d2) + params.dividend_yield * params.underlying_price * dividend_discount * normalCdf(d1) :
				params.interest_rate * params.strike_price * discount * normalCdf(-d2) - params.dividend_yield * params.underlying_price * dividend_discount * normalCdf(-d1);
			
			return theta_part1 + theta_part2;
		}
		case Greeks::Vega :
		{
			return params.underlying_price * dividend_discount * normalPdf(d1) * std::sqrt(params.time);
		}
		case Greeks::Rho :
		{
			return params.option_type == OptionType::Call ?
				params.strike_price * params.time * discount * normalCdf(d2) :
				-params.strike_price * params.time * discount * normalCdf(-d2);
		}
		default:
			throw std::invalid_argument("[!] Invalid Greek specified");
	}
}

/*
	@calculateMonteCarlo: Gets the Monte Carlo pricing

	Calculate the relation between one day and the time passed. For example 1 day represents 0.00274 of a year.

	The underlying asset price is modeled using GBM - dS = μS dt + σS dW
		-. Where:
			-.S is the stock price
			-.μ is the drift (expected return)
			-.σ is the volatility
			-.dW is a Wiener process

	Each day:
		-. Box-Muller Transform will be calculated
			-. Two random uniform variables u1 and u2 will be generated
			-. Transforms them using the following formula Z = sqrt(-2 * ln(u1)) * cos(2π * u2)
			-. This basically converts uniform randomness into normally distributed randomness. In an uniform distribution [0,1] all the 
			values between 0 and 1 have the same chance of appearing, but with the normal distribution, where the values that are closer
			to the mean have more chances of appearing. (the normal distribution is often used as an approximation or simplification when modeling stock)
		
		-. Geometric Brownian Motion will be calculated to simulate assets movement
			-. Drift: The drif basically represents the expected change in the asset's price over time
			-. Diffusion: Attempts to model the fluctuations in the asset's price that are random
			-. Exp is used in order to not get negative values for the underlying price

	After each day:
		-. Calculate the payoff and add it to the total payoff

	After all the iterations:
		-. Calculate the average payoff
		-. Calculate the discount to get the present value of the option. Essentially we're adjusting a future payoff to its present-day equivalent,
		considering what that money could earn if invested at the risk-free rate instead.
*/
[[nodiscard]] double FinancialCalculator::calculateMonteCarlo(const MonteCarloParams& params) const
{
	RandomGenerator<double> random_generator_uniform(0.0, 1.0);

	Price total_payoff{ 0.0 };
	Price underlying_price = params.underlying_price;

	const std::size_t total_days = static_cast<std::size_t>(params.time * 365.0);
	const double time_step = params.time / total_days;

	for (std::size_t i{ 0 }; i < params.number_of_simulations; ++i)
	{
		underlying_price = params.underlying_price; // Reset underlying price for each simulation

		for (std::size_t day{ 0 }; day < total_days; ++day)
		{
			// Box-Muller Transform to generate standard normal random variable
			const double u1 = random_generator_uniform.getRandomValue();
			const double u2 = random_generator_uniform.getRandomValue();
			const double random_normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * MATH_PI * u2);

			// Geometric Brownian Motion
			const auto drift = (params.interest_rate - 0.5 * std::pow(params.volatility, 2)) * time_step;
			const auto diffusion = params.volatility * std::sqrt(time_step) * random_normal;

			underlying_price *= std::exp(drift + diffusion);
		}

		double payoff{ 0.0 };
		if (params.option_type == OptionType::Call)
		{
			payoff = std::max(underlying_price - params.strike_price, 0.0);
		}
		else if (params.option_type == OptionType::Put)
		{
			payoff = std::max(params.strike_price - underlying_price, 0.0);
		}

		total_payoff += payoff;
	}

	const auto average_payoff = total_payoff / params.number_of_simulations;
	const auto discount_factor = std::exp(-params.interest_rate * params.time);

	return average_payoff * discount_factor;
}

[[nodiscard]] inline double normalCdf(double x) noexcept
{
	return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

[[nodiscard]] inline double normalPdf(double x) noexcept
{
	return (1.0 / std::sqrt(2.0 * MATH_PI)) * std::exp(-0.5 * x * x);
}

Option::Option(Price strike, Price premium, OptionType option_type) : strike_(strike), premium_(premium), option_type_(option_type) {};

// @calculatePayoff: Simply calculates the option's payoff at expiration or current value if exercised immediately, minus the premium paid

[[nodiscard]] inline double Option::calculatePayoff(Price spotPrice) const
{
	if (option_type_ == OptionType::Call)
	{
		return std::max(spotPrice - strike_, 0.0) - premium_;
	}
	
	return std::max(strike_ - spotPrice, 0.0) - premium_;
}

// @getPutSpread: A put spread is buying a put option with a higher strike price and selling a put option with a lower strike price
[[nodiscard]] inline StrategyPayoff CalculateStrategy::getPutSpread(const Option& long_put, const Option& short_put, Price spotPrice) const
{
	if (long_put.getStrike() <= short_put.getStrike())
	{
		throw std::runtime_error("[!] Long put strike should be higher than short put strike");
	}
	
	return long_put.calculatePayoff(spotPrice) - short_put.calculatePayoff(spotPrice);
}

// @getCallSpread: A call spread is buying a put option with a lower strike price and selling a put option with a higher strike price
[[nodiscard]] inline StrategyPayoff CalculateStrategy::getCallSpread(const Option& long_call, const Option& short_call, Price spotPrice) const
{
	if (long_call.getStrike() >= short_call.getStrike())
	{
		throw std::runtime_error("[!] Long call strike should be lower than short call strike");
	}
	
	return long_call.calculatePayoff(spotPrice) - short_call.calculatePayoff(spotPrice);
}

/*
	@getButterfly: A butterfly is:
		1-. Buying a call option with a lower strike price (which is the wing1)
		2-. Selling two call options with a middle strike price (which is the body)
		3-. Finally buying a call option with a higher strike price (which is the wing2).
*/
[[nodiscard]] inline StrategyPayoff CalculateStrategy::getButterfly(const Option& wing1, const Option& body, const Option& wing2, Price spotPrice) const
{
	if (wing1.getStrike() >= body.getStrike() || body.getStrike() >= wing2.getStrike())
	{
		throw std::runtime_error("[!] Strikes should be in ascending order");
	}
	
	return wing1.calculatePayoff(spotPrice) - 2.0 * body.calculatePayoff(spotPrice) + wing2.calculatePayoff(spotPrice);
}

// @getStrangle: A strangle is buying a put option with a lower strike price and buying a call option with a higher strike price.
[[nodiscard]] inline StrategyPayoff CalculateStrategy::getStrangle(const Option& put, const Option& call, Price spotPrice) const
{
	if (put.getStrike() >= call.getStrike())
	{
		throw std::runtime_error("[!] Put strike should be lower than Call strike");
	}
	
	return put.calculatePayoff(spotPrice) + call.calculatePayoff(spotPrice);
}

// @getStraddle : A straddle is buying a put option and a call option with the same strike price
[[nodiscard]] inline StrategyPayoff CalculateStrategy::getStraddle(const Option& put, const Option& call, Price spotPrice) const
{
	if (put.getStrike() != call.getStrike()) 
	{
		throw std::runtime_error("For Straddle, Put and Call strikes should be the same");
	}
	return put.calculatePayoff(spotPrice) + call.calculatePayoff(spotPrice);
}