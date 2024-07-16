#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>


struct Contract
{
    double premium;
    double strike;
    int dte;
    double delta;
    double gamma;
    double theta;
    double vega;
    double rho;
    double implied_volatility;
    double intrinsic_value;
    double market_price; 
};

// Error function approximation
double erf(double x) {
    const double A1 = 0.254829592;
    const double A2 = -0.284496736;
    const double A3 = 1.421413741;
    const double A4 = -1.453152027;
    const double A5 = 1.061405429;
    const double P = 0.3275911;

    // Save the sign of x
    int sign = (x >= 0) ? 1 : -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + P * x);
    double y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) * t * exp(-x * x);

    return sign * y;
}

// Cumulative standard normal density function
double cumulativeStandardNormal(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

// Black-Scholes option price calculation
double blackScholesPrice(double S0, double K, double r, double sigma, double T, bool isCallOption) {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    if (isCallOption) {
        return S0 * cumulativeStandardNormal(d1) - K * exp(-r * T) * cumulativeStandardNormal(d2);
    } else {
        return K * exp(-r * T) * cumulativeStandardNormal(-d2) - S0 * cumulativeStandardNormal(-d1);
    }
}

// Newton-Raphson method to find implied volatility
double findImpliedVolatility(double marketPrice, double S0, double K, double r, double T, bool isCallOption) {
    double sigma = 0.2; // initial guess
    double tolerance = 1e-5;
    int maxIterations = 100;
    for (int i = 0; i < maxIterations; ++i) {
        double price = blackScholesPrice(S0, K, r, sigma, T, isCallOption);
        double vega = S0 * sqrt(T) * exp(-0.5 * sigma * sigma * T) * cumulativeStandardNormal(0.5 * (log(S0 / K) / (sigma * sqrt(T)) + sigma * sqrt(T)));
        double priceDifference = marketPrice - price;
        if (fabs(priceDifference) < tolerance) {
            break;
        }
        sigma += priceDifference / vega; // Newton-Raphson step
    }
    return sigma;
}

// Black-Scholes Model for option pricing
std::vector<Contract> blackScholesOptionPricing(double S0, double K, double r, double sigma, double T, bool isCallOption) {
    std::vector<Contract> chain;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 0.5); 
    
    for (int i = -5; i < 5; i++) {
        Contract con;
        int days_till_expiry = T * 365.2425;
        con.dte = days_till_expiry;
        con.strike = K + i;
        double d1 = (log(S0 / (K + i)) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);

        if (isCallOption) {
            con.premium = S0 * cumulativeStandardNormal(d1) - (K + i) * exp(-r * T) * cumulativeStandardNormal(d2);
            con.delta = cumulativeStandardNormal(d1);
            con.intrinsic_value = std::max(S0 -(K + i), 0.0);
        } else {
            con.premium = (K + i) * exp(-r * T) * cumulativeStandardNormal(-d2) - S0 * cumulativeStandardNormal(-d1);
            con.delta = cumulativeStandardNormal(d1) - 1;
            con.intrinsic_value = std::max((K + i) - S0, 0.0);
        }

        con.gamma = cumulativeStandardNormal(d1) / (S0 * sigma * sqrt(T));
        con.theta = (-(S0 * cumulativeStandardNormal(d1) * sigma) / (2 * sqrt(T))) - (r * (K + i) * exp(-r * T) * cumulativeStandardNormal(d2));
        con.vega = S0 * cumulativeStandardNormal(d1) * sqrt(T);
        con.rho = (K + i) * T * exp(-r * T) * cumulativeStandardNormal(d2);
        
        // Simulate a market price with some noise
        double market_noise = distribution(generator);
        con.market_price = con.premium + market_noise;
        
        con.implied_volatility = findImpliedVolatility(con.market_price, S0, K + i, r, T, isCallOption);
        
        chain.push_back(con);
    }
    return chain;
}

int main() {
    // Option parameters
    double S0 = 100.0;   // Initial stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double sigma = 0.2;  // Volatility
    double T = 1;        // Time to maturity (in years)

    // Calculate option prices
    auto callChain = blackScholesOptionPricing(S0, K, r, sigma, T, true);
    auto putChain = blackScholesOptionPricing(S0, K, r, sigma, T, false);

    // Output the results
    for (const auto& con : callChain) {
        std::cout << "Strike: " << con.strike << " European Call Option Price: " << con.premium
                  << ", Market Price: " << con.market_price << ", dte: " << con.dte
                  << ", delta: " << con.delta << ", gamma: " << con.gamma
                  << ", theta: " << con.theta << ", vega: " << con.vega
                  << ", rho: " << con.rho << ", implied volatility: " << con.implied_volatility
                  << ", intrinsic value: " << con.intrinsic_value << std::endl;
    }

    std::cout << "" << std::endl;

    for (const auto& con : putChain) {
        std::cout << "Strike: " << con.strike << " European Put Option Price: " << con.premium
                  << ", Market Price: " << con.market_price << ", dte: " << con.dte
                  << ", delta: " << con.delta << ", gamma: " << con.gamma
                  << ", theta: " << con.theta << ", vega: " << con.vega
                  << ", rho: " << con.rho << ", implied volatility: " << con.implied_volatility
                  << ", intrinsic value: " << con.intrinsic_value << std::endl;
    }

    return 0;
}