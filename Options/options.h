#pragma once

#include <cmath>
#include <stdexcept>
#include <random>

constexpr const double MATH_PI = 3.14159265358979323846;

enum class CalculationType
{
	Futures,	// How much an investment or asset will be worth at a specific point in the future		
	BlackScholes,	// Theoretical estimate for the price of European-style options.
	MonteCarlo,	// Calculates the price of a European call option uisng Monte Carlo simulation
	Greeks,			
};

// The greeks are measures which are used in order to assess the sensitivity of an option's price to changes
enum class Greeks
{
	Delta,  // Δ : How much the price of an option is expected to move per unit change in the underlying asset's price
	Gamma,  // Γ : Rate of change of delta over time or the rate of change in the underlying asset's price (Higher gamma means the option is more sensitive to price changes in the underlying asset)
	Theta,  // Θ : Rate of decline in the value of an option as it approaches its expiration date
	Vega,   // ν : Shows how much the price of an option is expected to change with a 1 % change in implied volatility (Higher vega means the option is more sensitive to volatility changes)
	Rho     // ρ : Reflects how much the price of an option would change with a 1% change in interest rates
};

/*
*	Delta is related to underlying price (how much an option's price changes if the underlying price changes)
*	Gamma is related to underlying (rate of change in delta as the underlying price changes)
*	Theta is related to time (how much an option's price decreases as time passes)
*	Vega is related to volatility (how much an option's price changes when the implied volatility of the underlying asset changes)
*	Rho is related to interest rate (how much an option's price changes when interest rates change)
* 
*	Theta (Θ) is typically negative for both calls and puts
*	Vega (ν) is typically highest for at-the-money options
*	Rho (ρ) is positive for calls and negative for puts
*/

enum class OptionType
{
	Call,
	Put
};

// Aliases for better readability
using Price = double;
using Time = double;
using InterestRate = double;
using PriceTheo = double;
using RateChange = double;
using StrategyPayoff = double;
using DividendYield = double;
using Volatility = double;

// Struct to hold the parameters needed for the Black-Scholes calculations
struct BlackScholesParams
{
	InterestRate interest_rate{};	//  r          : Represents the interest rate (is the theoretical rate of return on an investment with zero risk)
	Price underlying_price{};       //  S          : Represents the current underlying price (current market price of the asset)
	Price strike_price{};           //  K          : Represents the current strike price (fixed price at which the owner of a call option can buy, or the owner of a put option can sell, the underlying asset).
	Time time{};                    //  T          : Represents the time remaning until the option's expiration (Years).
	Volatility volatility{};        //  σ          : Represents the volatility (variation in the price of the underlying asset over time).
	OptionType option_type{};       //  Call / Put : Represents the type (can either be a Call or a Put).
	Price paid_price{};             //  $          : Represents the paid price
};

// Struct to hold the parameters needed for the futures calculations
struct FuturesParams
{
	Price present_value{};		// Initial investment amount or current value of the asset
	InterestRate interest_rate{};   // Annual interest rate or rate of return on the investment
	Time time{};			// Time horizon over which the investment will grow (Years)		
};

// Struct to hold the parameters needed for the greeks calculations
struct GreeksParams
{
	InterestRate interest_rate{};
	Price underlying_price{};    
	Price strike_price{};  
	Time time{};
	Volatility volatility{};       
	OptionType option_type{};
	Price paid_price{};            
	DividendYield dividend_yield{};	/* Annual dividend payment expressed as a percentage of the stock's price. For example, if a stock pays 1$ in dividends
									per year and the stock price is 100$, then the dividend yield is (1/100) * 100*/
};

// Struct to hold the parameters needed for the Monte Carlo calculations
struct MonteCarloParams
{
	std::size_t number_of_simulations{};
	InterestRate interest_rate{};	
	Price underlying_price{};      
	Price strike_price{};          
	Time time{};    
	Volatility volatility{};
	OptionType option_type{}; 
	Price paid_price{};
};

// Small class used to generate values under certain type (only numeric values allowed) and bound constraints
template<typename RandomType>
class RandomGenerator
{
private:
	std::random_device rd{};
	std::mt19937_64 gen{ rd() };

	static_assert(std::is_arithmetic_v<RandomType>, "[!] The type must be a numeric type");
	typename std::conditional<std::is_floating_point_v<RandomType>,
		std::uniform_real_distribution<RandomType>,
		std::uniform_int_distribution<RandomType>>::type distribution;

public:
	RandomGenerator(const RandomType left_limit, const RandomType right_limit) noexcept
		: distribution(std::min(left_limit, right_limit), std::max(left_limit, right_limit)) {}

	[[nodiscard]] inline RandomType getRandomValue() noexcept { return distribution(gen); }
};

class FinancialCalculator
{
public:
	FinancialCalculator() = default;

	[[nodiscard]] Price calculateBlackScholes(const BlackScholesParams& params) const;
	[[nodiscard]] inline double calculateFutures(const FuturesParams& params) const;
	[[nodiscard]] Price calculateGreeks(const GreeksParams& params, const Greeks greek) const;
	[[nodiscard]] Price calculateMonteCarlo(const MonteCarloParams& params) const;
};

// @normalCdf : Calculates the cumulative distribution function (CDF) of the standard normal distribution (median 0 and variance 1)
[[nodiscard]] inline double normalCdf(double x) noexcept;

// @normalPdf : Calculates the probability density function (PDF) of the standard normal distribution
[[nodiscard]] inline double normalPdf(double x) noexcept;

class Option
{
private:
	Price strike_{};
	Price premium_{};
	OptionType option_type_{};
public:
	Option(Price strike, Price premium, OptionType type);

	[[nodiscard]] inline Price calculatePayoff(Price spotPrice) const;

	[[nodiscard]] inline Price getStrike() const noexcept { return strike_; }
	[[nodiscard]] inline Price getPremium()  const noexcept { return premium_; }
	[[nodiscard]] inline OptionType getType() const noexcept { return option_type_; }
};

class CalculateStrategy
{
public:
	CalculateStrategy() = default;

	[[nodiscard]] inline StrategyPayoff getPutSpread(const Option& long_put, const Option& short_put, Price spotPrice) const;
	[[nodiscard]] inline StrategyPayoff getCallSpread(const Option& long_call, const Option& short_call, Price spotPrice) const;
	[[nodiscard]] inline StrategyPayoff getButterfly(const Option& wing1, const Option& body, const Option& wing2, Price spotPrice) const;
	[[nodiscard]] inline StrategyPayoff getStrangle(const Option& put, const Option& call, Price spotPrice) const;
	[[nodiscard]] inline StrategyPayoff getStraddle(const Option& put, const Option& call, Price spotPrice) const;
};
